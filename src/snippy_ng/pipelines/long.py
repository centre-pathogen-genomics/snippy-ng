from pathlib import Path
from typing import Optional
from pydantic import Field
from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.pipelines import PipelineBuilder, SnippyPipeline
from snippy_ng.stages.reporting import PrintVcfHistogram
from snippy_ng.stages.stats import SeqKitReadStatsBasic
from snippy_ng.stages.alignment import Minimap2LongReadAligner
from snippy_ng.stages.filtering import SamtoolsFilter, VcfFilterLong
from snippy_ng.stages.clean_reads import SeqkitCleanLongReads
from snippy_ng.stages.calling import FreebayesCallerLong, Clair3Caller
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
from snippy_ng.stages.compression import BgzipCompressor
from snippy_ng.stages.masks import DepthMask, ApplyMask, HetMask
from snippy_ng.stages.copy import CopyFasta
from snippy_ng.pipelines.common import load_or_prepare_reference


class LongPipelineBuilder(PipelineBuilder):
    """Builder for long-read SNP calling pipeline."""
    reference: Path = Field(..., description="Reference genome file path")
    reads: Optional[Path] = Field(default=None, description="Long reads file (FASTQ)")
    bam: Optional[Path] = Field(default=None, description="Pre-aligned BAM/CRAM file")
    prefix: str = Field(default="snps", description="Output file prefix")
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
    min_qual: float = Field(default=100, description="Minimum variant quality")
    min_depth: int = Field(default=10, description="Minimum variant depth")
    mask: Optional[str] = Field(default=None, description="BED file with regions to mask")

    def build(self) -> SnippyPipeline:
        """Build and return the long-read pipeline."""
        stages = []
        globals = {'prefix': self.prefix}
        
        # Setup reference (load existing or prepare new)
        setup = load_or_prepare_reference(
            reference_path=self.reference,
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
        if current_reads and current_reads:
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
            aligned_reads = self.bam
        else:
            # SeqKit read statistics
            stats_stage = SeqKitReadStatsBasic(
                reads=current_reads,
                **globals
            )
            stages.append(stats_stage)
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
                reference=reference_file,
                reference_index=reference_index,
                fbopt=self.caller_opts,
                mincov=2,
                **globals
            )
        stages.append(caller_stage)
        
        # Filter VCF
        variant_filter = VcfFilterLong(
            vcf=caller_stage.output.vcf,
            reference=reference_file,
            reference_index=reference_index,
            min_qual=self.min_qual,
            min_depth=self.min_depth,
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
            vcf=caller_stage.output.vcf,  # Use raw VCF for complete site information
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
