from pathlib import Path
from typing import Optional
from pydantic import Field
from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.pipelines import PipelineBuilder, SnippyPipeline
from snippy_ng.stages.reporting import PrintVcfHistogram
from snippy_ng.stages.stats import SeqKitReadStatsBasic, VcfStats
from snippy_ng.stages.alignment import Minimap2LongReadAligner
from snippy_ng.stages.filtering import SamtoolsFilter
from snippy_ng.stages.vcf import VcfFilterLong, AddDeletionstoVCF
from snippy_ng.stages.compression import CramCompressor, VcfCompressor
from snippy_ng.stages.clean_reads import SeqkitCleanLongReads
from snippy_ng.stages.calling import FreebayesCallerLong, Clair3Caller
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
from snippy_ng.stages.masks import DepthMask, ApplyMask, ZeroDepthBedFromBam
from snippy_ng.stages.copy import FinaliseFasta
from snippy_ng.pipelines.common import load_or_prepare_reference
from snippy_ng.utils.gather import guess_sample_id


class LongPipelineBuilder(PipelineBuilder):
    """Builder for long-read SNP calling pipeline."""
    reference: Path = Field(..., description="Reference genome file path")
    reads: Optional[Path] = Field(default=None, description="Long reads file (FASTQ)")
    bam: Optional[Path] = Field(default=None, description="Pre-aligned BAM/CRAM file")
    prefix: str = Field(default="snippy", description="Output file prefix")
    clean_reads: bool = Field(default=False, description="Clean reads with fastp")
    downsample: Optional[float] = Field(default=None, description="Target coverage for downsampling")
    aligner: str = Field(default="minimap2", description="Aligner to use")
    aligner_opts: str = Field(default="", description="Additional aligner options")
    minimap_preset: str = Field(default="map-ont", description="Minimap2 preset")
    caller: str = Field(default="clair3", description="Variant caller (clair3 or freebayes)")
    caller_opts: str = Field(default="", description="Additional caller options")
    clair3_model: Optional[Path] = Field(default=None, description="Clair3 model path")
    clair3_fast_mode: bool = Field(default=False, description="Use Clair3 fast mode")
    min_read_len: int = Field(default=1000, description="Minimum read length")
    min_read_qual: float = Field(default=10, description="Minimum read quality")
    mask: Optional[str] = Field(default=None, description="BED file with regions to mask")
    depth_mask: int = Field(default=10, description="Mask regions in the output fasta with Ns if the read depth is below this threshold")
    min_qual: float = Field(default=2, description="Mark variants below this QUAL threshold as LowQual in the output VCF")
    sample_name: Optional[str] = Field(default=None, description="Optional sample name override for output tables")

    def build(self) -> SnippyPipeline:
        """Build and return the long-read pipeline."""
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
        # Always keep as strings for Pydantic validation
        current_reads = [str(self.reads)] if self.reads else []
        
        if self.downsample and current_reads:
            from snippy_ng.stages.downsample_reads import RasusaDownsampleReadsByCoverage
            
            # We need the genome length at run time (once we know the reference)
            downsample_stage = RasusaDownsampleReadsByCoverage(
                ref_metadata=ref_metadata,
                coverage=self.downsample,
                reads=current_reads,
                **globals
            )
            # Update reads to use downsampled reads (convert to strings for Pydantic)
            current_reads = [str(downsample_stage.output.downsampled_r1)]
            if downsample_stage.output.downsampled_r2:
                current_reads.append(str(downsample_stage.output.downsampled_r2))
            stages.append(downsample_stage)
        
        # Clean reads
        if self.clean_reads and current_reads:
            clean_reads_stage = SeqkitCleanLongReads(
                reads=str(current_reads[0]),
                min_length=self.min_read_len,
                min_qscore=self.min_read_qual,
                **globals
            )
            # Update reads to use cleaned reads (convert to strings for Pydantic)
            current_reads = [str(clean_reads_stage.output.cleaned_reads)]
            stages.append(clean_reads_stage)

        # Aligner
        if self.bam:
            if sample_name is None:
                sample_name = guess_sample_id(Path(self.bam).name)
            aligned_reads = self.bam
        else:
            # SeqKit read statistics
            stats_stage = SeqKitReadStatsBasic(
                reads=current_reads,
                sample_name=sample_name,
                **globals
            )
            stages.append(stats_stage)
            stats_tsv = getattr(stats_stage.output, "stats_tsv", None)
            # Minimap2
            if self.aligner == "minimap2":
                aligner_stage = Minimap2LongReadAligner(
                    reads=current_reads,
                    reference=reference_file,
                    aligner_opts=self.aligner_opts,
                    minimap_preset=self.minimap_preset,
                    **globals
                )
            else:
                raise ValueError(f"Unsupported aligner '{self.aligner}'")
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
        if self.caller == "clair3":
            assert self.clair3_model is not None, "Clair3 model must be provided when using Clair3 caller."
            assert Path(self.clair3_model).is_absolute(), f"Clair3 model path '{self.clair3_model}' must be an absolute path."
            platform = 'ont'
            if 'hifi' in str(self.clair3_model).lower():
                platform = 'hifi'
            caller_stage = Clair3Caller(
                bam=aligned_reads,
                reference=reference_file,
                reference_index=reference_index,
                clair3_model=self.clair3_model,
                fast_mode=self.clair3_fast_mode,
                platform=platform,
                **globals
            )
        else:
            caller_stage = FreebayesCallerLong(
                bam=aligned_reads,
                bam_index=align_filter.output.bai,
                reference=reference_file,
                reference_index=reference_index,
                fbopt=self.caller_opts,
                **globals
            )
        stages.append(caller_stage)
        
        # Filter VCF
        variant_filter = VcfFilterLong(
            vcf=caller_stage.output.vcf,
            reference=reference_file,
            reference_index=reference_index,
            min_qual=self.min_qual,
            min_depth=self.depth_mask,
            **globals
        )
        stages.append(variant_filter)
        variants_file = variant_filter.output.vcf

        zero_depth_bed = ZeroDepthBedFromBam(
            bam=aligned_reads,
            **globals
        )
        stages.append(zero_depth_bed)

        # Add zero-depth regions to VCF as symbolic deletion blocks
        add_deletions = AddDeletionstoVCF(
            zero_depth_bed=zero_depth_bed.output.zero_depth_bed,
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
        gzip = VcfCompressor(
            input=consequences.output.annotated_vcf,
            **globals
        )
        stages.append(gzip)
        
        # Pseudo-alignment
        pseudo = BcftoolsPseudoAlignment(
            ref_metadata=ref_metadata,
            vcf_gz=gzip.output.gz,
            reference=reference_file,
            **globals
        )
        stages.append(pseudo)
        
        # Track the current reference/fasta through the masking stages
        current_fasta = pseudo.output.fasta
        
        # Apply minimum-depth masking
        if self.depth_mask > 0:
            depth_mask = DepthMask(
                bam=aligned_reads,
                fasta=current_fasta,
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
            consequences.output.annotated_vcf, 
            cram_compressor.output.cram,
            vcf_stats.output.summary_tsv,
            vcf_stats.output.breakdown_tsv,
        ]
        if stats_tsv is not None:
            keep_files.append(stats_tsv)
        return SnippyPipeline(stages=stages, outputs_to_keep=keep_files)
