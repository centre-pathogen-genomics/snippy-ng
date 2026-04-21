from pathlib import Path
from typing import Optional, List
from pydantic import Field
from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.pipelines import PipelineBuilder, SnippyPipeline
from snippy_ng.stages.clean_reads import FastpCleanReads
from snippy_ng.stages.reporting import PrintVcfHistogram
from snippy_ng.stages.stats import SeqKitReadStatsBasic, VcfStats
from snippy_ng.stages.alignment import BWAMEMShortReadAligner, Minimap2ShortReadAligner
from snippy_ng.stages.filtering import SamtoolsFilter
from snippy_ng.stages.vcf import VcfFilterShort, AddDeletionsToVCF, VcfPassFilter
from snippy_ng.stages.calling import FreebayesCaller
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
from snippy_ng.stages.compression import CramCompressor, VcfCompressor
from snippy_ng.stages.masks import ApplyMask, DepthBedsFromBam, ApplyDepthMaskToFasta
from snippy_ng.stages.copy import FinaliseFasta
from snippy_ng.pipelines.common import load_or_prepare_reference
from snippy_ng.utils.gather import guess_sample_id


class ShortPipelineBuilder(PipelineBuilder):
    """Builder for short-read SNP calling pipeline."""
    reference: Path = Field(..., description="Reference genome file path")
    reads: List[Path] = Field(..., description="Short read files (FASTQ, R1 and optionally R2)")
    prefix: str = Field(default="snippy", description="Output file prefix")
    bam: Optional[Path] = Field(default=None, description="Pre-aligned BAM/CRAM file")
    clean_reads: bool = Field(default=False, description="Clean reads with fastp")
    downsample: Optional[float] = Field(default=None, description="Target coverage for downsampling")
    aligner: str = Field(default="minimap2", description="Aligner to use (minimap2 or bwamem)")
    aligner_opts: str = Field(default="", description="Additional aligner options")
    caller_opts: str = Field(default="", description="Additional caller options")
    mask: Optional[str] = Field(default=None, description="BED file with regions to mask")
    depth_mask: int = Field(default=10, description="Mask regions in the output fasta with Ns if the read depth is below this threshold")
    min_qual: float = Field(default=100, description="Mark variants below this QUAL threshold as LowQual in the output VCF")
    min_mapping_quality: int = Field(default=30, description="Minimum mapping quality for FreeBayes calls and depth masks")
    sample_name: Optional[str] = Field(default=None, description="Optional sample name override for output tables")
    add_deletions_to_vcf: bool = Field(default=True, description="Add zero-depth regions to VCF as symbolic deletion blocks")


    def build(self) -> SnippyPipeline:
        """Build and return the short-read pipeline."""
        stages = []
        globals = {'prefix': self.prefix}
        stats_tsv = None
        sample_name = self.sample_name
        
        # Setup reference (load existing or prepare new)
        setup = load_or_prepare_reference(
            reference_path=self.reference,
            output_directory=Path("reference"),
        )
        reference_file = setup.output.reference
        features_file = setup.output.gff
        reference_index = setup.output.reference_index
        ref_metadata = ReferenceMetadata(setup.output.metadata)
        stages.append(setup)
        
        # Track current reads through potential cleaning and downsampling
        current_reads = self.reads.copy() if self.reads else []
        
        if self.downsample and current_reads:
            from snippy_ng.stages.downsample_reads import RasusaDownsampleReadsByCoverage
            
            # We need the genome length at run time (once we know the reference)
            downsample_stage = RasusaDownsampleReadsByCoverage(
                ref_metadata=ref_metadata,
                coverage=self.downsample,
                reads=current_reads,
                **globals
            )
            # Update reads to use downsampled reads
            current_reads = [downsample_stage.output.downsampled_r1]
            if downsample_stage.output.downsampled_r2:
                current_reads.append(downsample_stage.output.downsampled_r2)
            stages.append(downsample_stage)
        
        # Clean reads (optional)
        if self.clean_reads and current_reads:
            clean_reads_stage = FastpCleanReads(
                reads=current_reads,
                **globals
            )
            # Update reads to use cleaned reads
            current_reads = [clean_reads_stage.output.cleaned_r1]
            if clean_reads_stage.output.cleaned_r2:
                current_reads.append(clean_reads_stage.output.cleaned_r2)
            stages.append(clean_reads_stage)

        if self.bam:
            if sample_name is None:
                sample_name = guess_sample_id(Path(self.bam).name)
            aligned_reads = Path(self.bam).absolute()  
        else:
            # SeqKit read statistics
            stats_stage = SeqKitReadStatsBasic(
                reads=current_reads,
                sample_name=sample_name,
                **globals
            )
            stages.append(stats_stage)
            stats_tsv = stats_stage.output.stats_tsv
            # Aligner
            if self.aligner == "bwamem":
                aligner_stage = BWAMEMShortReadAligner(
                    reads=current_reads,
                    reference=reference_file,
                    reference_index=reference_index,
                    aligner_opts=self.aligner_opts,
                    **globals
                )
            else:
                # Minimap2
                aligner_stage = Minimap2ShortReadAligner(
                    reads=current_reads,
                    reference=reference_file,
                    aligner_opts=self.aligner_opts,
                    **globals
                )
            
            if current_reads and sample_name is None:
                sample_name = guess_sample_id(Path(current_reads[0]).name)
            aligned_reads = aligner_stage.output.bam
            stages.append(aligner_stage)
        
        # Filter alignment
        align_filter = SamtoolsFilter(
            bam=aligned_reads,
            reference=reference_file,
            **globals
        )
        aligned_reads = align_filter.output.bam
        stages.append(align_filter)
        
        # SNP calling
        caller = FreebayesCaller(
            bam=aligned_reads,
            bam_index=align_filter.output.bai,
            reference=reference_file,
            reference_index=reference_index,
            fbopt=self.caller_opts,
            min_mapping_quality=self.min_mapping_quality,
            **globals
        )
        stages.append(caller)
        
        # Filter VCF
        variant_filter = VcfFilterShort(
            vcf=caller.output.vcf,
            reference=reference_file,
            min_qual=self.min_qual,
            min_depth=self.depth_mask,
            **globals
        )
        stages.append(variant_filter)
        variants_file = variant_filter.output.vcf

        depth_beds = DepthBedsFromBam(
            bam=aligned_reads,
            min_depth=self.depth_mask,
            min_base_quality=13,
            min_mapping_quality=self.min_mapping_quality,
            **globals
        )
        stages.append(depth_beds)

        if self.add_deletions_to_vcf:
            # Add zero-depth regions to VCF as symbolic deletion blocks
            add_deletions = AddDeletionsToVCF(
                zero_depth_bed=depth_beds.output.zero_depth_bed,
                vcf=variants_file,
                reference=reference_file,
                **globals
            )
            stages.append(add_deletions)
            variants_file = add_deletions.output.vcf
        
        # Consequences calling
        consequences = BcftoolsConsequencesCaller(
            variants=variants_file,
            features=features_file,
            reference=reference_file,
            **globals
        )
        stages.append(consequences)

        vcf_stats = VcfStats(
            vcf=consequences.output.annotated_vcf,
            sample_name=sample_name,
            **globals
        )
        stages.append(vcf_stats)

        # Compress VCF
        gzip_vcf = VcfCompressor(
            input=consequences.output.annotated_vcf,
            **globals
        )
        stages.append(gzip_vcf)
        
        # Filter to PASS-only variants
        pass_filter = VcfPassFilter(
            vcf=consequences.output.annotated_vcf,
            **globals
        )
        stages.append(pass_filter)

        # Pseudo-alignment
        pseudo = BcftoolsPseudoAlignment(
            ref_metadata=ref_metadata,
            vcf=pass_filter.output.vcf,
            reference=reference_file,
            **globals
        )
        stages.append(pseudo)
        
        # Track the current reference/fasta through the masking stages
        current_fasta = pseudo.output.fasta

        # Apply minimum-depth masking after consensus so the reference bases still
        # match VCF REF alleles while bcftools consensus is running.
        if self.depth_mask > 0:
            depth_mask = ApplyDepthMaskToFasta(
                fasta=current_fasta,
                mask_bed=depth_beds.output.min_depth_bed,
                min_depth=self.depth_mask,
                **globals
            )
            stages.append(depth_mask)
            current_fasta = depth_mask.output.masked_fasta

        
        
        # Apply user mask if provided
        if self.mask:
            user_mask = ApplyMask(
                fasta=current_fasta,
                mask_bed=Path(self.mask),
                **globals
            )
            stages.append(user_mask)
            current_fasta = user_mask.output.masked_fasta

        # Copy final masked consensus to standard output location
        copy_final = FinaliseFasta(
            input=current_fasta,
            output_path=f"{self.prefix}.pseudo.fna",
            **globals
        )
        stages.append(copy_final)

        # Compress BAM to CRAM with embedded reference
        cram_compressor = CramCompressor(
            input=aligned_reads,
            reference=reference_file,
        )
        stages.append(cram_compressor)

        # Print VCF histogram to terminal
        vcf_histogram = PrintVcfHistogram(
            vcf_path=variants_file,
            **globals
        )
        stages.append(vcf_histogram)

        keep_files = [
            copy_final.output.fasta, 
            gzip_vcf.output.gz,
            pass_filter.output.vcf, 
            cram_compressor.output.cram,
            vcf_stats.output.summary_tsv,
            vcf_stats.output.breakdown_tsv,
        ]
        if stats_tsv is not None:
            keep_files.append(stats_tsv)
        return SnippyPipeline(stages=stages, outputs_to_keep=keep_files)
