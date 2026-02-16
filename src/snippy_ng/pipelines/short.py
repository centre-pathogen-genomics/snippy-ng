from pathlib import Path
from typing import Optional, List
from pydantic import Field
from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.pipelines import PipelineBuilder, SnippyPipeline
from snippy_ng.stages.clean_reads import FastpCleanReads
from snippy_ng.stages.reporting import PrintVcfHistogram
from snippy_ng.stages.stats import SeqKitReadStatsBasic
from snippy_ng.stages.alignment import BWAMEMShortReadAligner, Minimap2ShortReadAligner
from snippy_ng.stages.filtering import SamtoolsFilter, VcfFilterShort
from snippy_ng.stages.calling import FreebayesCaller
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
from snippy_ng.stages.compression import BgzipCompressor
from snippy_ng.stages.masks import DepthMask, ApplyMask, HetMask
from snippy_ng.stages.copy import CopyFasta
from snippy_ng.pipelines.common import load_or_prepare_reference


class ShortPipelineBuilder(PipelineBuilder):
    """Builder for short-read SNP calling pipeline."""
    reference: Path = Field(..., description="Reference genome file path")
    reads: List[Path] = Field(..., description="Short read files (FASTQ, R1 and optionally R2)")
    prefix: str = Field(default="snps", description="Output file prefix")
    bam: Optional[Path] = Field(default=None, description="Pre-aligned BAM/CRAM file")
    clean_reads: bool = Field(default=False, description="Clean reads with fastp")
    downsample: Optional[float] = Field(default=None, description="Target coverage for downsampling")
    aligner: str = Field(default="minimap2", description="Aligner to use (minimap2 or bwamem)")
    aligner_opts: str = Field(default="", description="Additional aligner options")
    caller_opts: str = Field(default="", description="Additional caller options")
    mask: Optional[str] = Field(default=None, description="BED file with regions to mask")
    min_depth: int = Field(default=10, description="Minimum variant depth")
    min_qual: float = Field(default=100, description="Minimum variant quality")
    tmpdir: Optional[Path] = Field(default=None, description="Temporary directory")
    cpus: int = Field(default=1, description="Number of CPUs to use")
    ram: int = Field(default=8, description="RAM in GB")

    def build(self) -> SnippyPipeline:
        """Build and return the short-read pipeline."""
        stages = []
        globals = {'prefix': self.prefix, 'cpus': self.cpus, 'ram': self.ram, 'tmpdir': self.tmpdir}
        
        # Setup reference (load existing or prepare new)
        setup = load_or_prepare_reference(
            reference_path=self.reference
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
            aligned_reads = Path(self.bam).resolve()  
        else:
            # SeqKit read statistics
            stats_stage = SeqKitReadStatsBasic(
                reads=current_reads,
                **globals
            )
            stages.append(stats_stage)
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
            
            aligned_reads = aligner_stage.output.cram
            stages.append(aligner_stage)
        
        # Filter alignment
        align_filter = SamtoolsFilter(
            bam=aligned_reads,
            **globals
        )
        aligned_reads = align_filter.output.cram
        stages.append(align_filter)
        
        # SNP calling
        caller = FreebayesCaller(
            bam=aligned_reads,
            reference=reference_file,
            reference_index=reference_index,
            fbopt=self.caller_opts,
            mincov=self.min_depth,
            **globals
        )
        stages.append(caller)
        
        # Filter VCF
        variant_filter = VcfFilterShort(
            vcf=caller.output.vcf,
            reference=reference_file,
            min_depth=self.min_depth,
            min_qual=self.min_qual,
            **globals
        )
        stages.append(variant_filter)
        variants_file = variant_filter.output.vcf
        
        # Consequences calling
        consequences = BcftoolsConsequencesCaller(
            variants=variants_file,
            features=features_file,
            reference=reference_file,
            **globals
        )
        stages.append(consequences)
        
        # Compress VCF
        gzip = BgzipCompressor(
            input=consequences.output.annotated_vcf,
            suffix="gz",
            **globals
        )
        stages.append(gzip)
        
        # Pseudo-alignment
        pseudo = BcftoolsPseudoAlignment(
            ref_metadata=ref_metadata,
            vcf_gz=gzip.output.compressed,
            reference=reference_file,
            **globals
        )
        stages.append(pseudo)
        
        # Track the current reference/fasta through the masking stages
        current_fasta = pseudo.output.fasta
        
        # Apply depth masking
        depth_mask = DepthMask(
            bam=aligned_reads,
            fasta=current_fasta,
            min_depth=self.min_depth,
            **globals
        )
        stages.append(depth_mask)
        current_fasta = depth_mask.output.masked_fasta

        # Apply heterozygous and low quality sites masking
        het_mask = HetMask(
            vcf=caller.output.vcf,  # Use raw VCF for complete site information
            fasta=current_fasta,
            min_qual=self.min_qual,
            **globals
        )
        stages.append(het_mask)
        current_fasta = het_mask.output.masked_fasta
        
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
        copy_final = CopyFasta(
            input=current_fasta,
            output_path=f"{self.prefix}.pseudo.fna",
            **globals
        )
        stages.append(copy_final)

        # Print VCF histogram to terminal
        vcf_histogram = PrintVcfHistogram(
            vcf_path=variants_file,
            **globals
        )
        stages.append(vcf_histogram)

        return SnippyPipeline(stages=stages)
