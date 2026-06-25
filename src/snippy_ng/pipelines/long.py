from pathlib import Path
from typing import Optional
from pydantic import Field
from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.pipelines import PipelineBuilder, SnippyPipeline
from snippy_ng.stages.reporting import PrintVcfHistogram, SampleReport
from snippy_ng.stages.stats import SeqKitReadStatsBasic, VcfStats
from snippy_ng.stages.alignment import Minimap2LongReadAligner
from snippy_ng.stages.filtering import BamReferenceValidator, SamtoolsFilter
from snippy_ng.stages.vcf import VcfFilterLong, AddDeletionsToVCF, VcfPassFilter, CollapseDiploidGenotypes, VcfToTab
from snippy_ng.stages.compression import CramCompressor, VcfCompressor
from snippy_ng.stages.clean_reads import SeqkitCleanLongReads
from snippy_ng.stages.calling import FreebayesCallerLong, Clair3Caller, LongbowClair3ModelSelector
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
from snippy_ng.stages.masks import DepthBedsFromBam, ApplyDepthMaskToFasta, ApplyMask, MaskMixedSites
from snippy_ng.stages.copy import CopyFile, FinaliseFasta
from snippy_ng.pipelines.common import download_assembly, load_or_prepare_reference
from snippy_ng.utils.gather import strip_bio_suffixes


class LongPipelineBuilder(PipelineBuilder):
    """Builder for long-read SNP calling pipeline."""
    reference: Optional[Path] = Field(default=None, description="Reference genome file path")
    reference_accession: Optional[str] = Field(default=None, description="Reference assembly accession to download")
    reads: Optional[Path] = Field(default=None, description="Long reads file (FASTQ)")
    bam: Optional[Path] = Field(default=None, description="Pre-aligned BAM/CRAM file")
    prefix: str = Field(default="snippy", description="Output file prefix")
    clean_reads: bool = Field(default=False, description="Clean reads with SeqKit before alignment")
    min_read_len: int = Field(default=1000, description="Minimum read length")
    min_read_qual: float = Field(default=10, description="Minimum read quality")
    downsample: Optional[float] = Field(default=None, description="Target coverage for downsampling")
    aligner: str = Field(default="minimap2", description="Aligner to use")
    aligner_opts: str = Field(default="", description="Additional aligner options")
    minimap_preset: str = Field(default="map-ont", description="Minimap2 preset")
    caller: str = Field(default="clair3", description="Variant caller (clair3 or freebayes)")
    caller_opts: str = Field(default="", description="Additional caller options")
    clair3_model: Optional[Path] = Field(default=None, description="Clair3 model path")
    mask: Optional[Path] = Field(default=None, description="BED file with regions to mask")
    depth_mask: int = Field(default=10, description="Mask regions in the output fasta with Ns if the read depth is below this threshold")
    min_qual: Optional[float] = Field(default=None, description="Mark variants below this QUAL threshold as LowQual in the output VCF")
    min_mapping_quality: int = Field(default=10, description="Minimum mapping quality for FreeBayes calls and depth masks")
    sample_name: Optional[str] = Field(default=None, description="Optional sample name override for output tables")
    add_deletions_to_vcf: bool = Field(default=True, description="Add zero-depth regions to VCF as symbolic deletion blocks")
    haploid: bool = Field(default=True, description="Collapse diploid genotypes to haploid genotypes after consequence calling")
    report: bool = Field(default=True, description="Create per-sample HTML report")
    report_scope: str = Field(default="all", description="Variant scope for the sample report: pass or all")
    report_window_size: int = Field(default=100, description="Base pairs of context around each report variant")


    def build(self) -> SnippyPipeline:
        """Build and return the long-read pipeline."""
        stages = []
        globals = {'prefix': self.prefix}
        stats_tsv = None
        sample_name = self.sample_name
        reference_input = self.reference

        if self.reference_accession:
            reference_input = download_assembly(
                self.reference_accession,
                stages,
                output_directory=Path("reference"),
            )
        if reference_input is None:
            raise ValueError("Reference genome path or reference accession must be provided.")
        
        # Setup reference (load existing or prepare new)
        setup = load_or_prepare_reference(
            reference_path=reference_input,
            output_directory=Path("reference"),
        )
        reference_file = setup.output.reference
        features_file = setup.output.gff
        reference_index = setup.output.reference_index
        ref_metadata = ReferenceMetadata(setup.output.metadata)
        stages.append(setup)
        
        # Track current reads through potential cleaning and downsampling
        # Always keep as strings for Pydantic validation
        current_reads: list[Path] = [self.reads] if self.reads else []
        
        # Clean reads
        if self.clean_reads and current_reads:
            clean_reads_stage = SeqkitCleanLongReads(
                reads=str(current_reads[0]),
                min_length=self.min_read_len,
                min_qscore=self.min_read_qual,
                **globals
            )
            # Update reads to use cleaned reads (convert to strings for Pydantic)
            current_reads = [clean_reads_stage.output.cleaned_reads]
            stages.append(clean_reads_stage)

        if self.downsample and current_reads:
            from snippy_ng.stages.downsample import RasusaDownsampleReadsByCoverage
            
            # We need the genome length at run time (once we know the reference)
            downsample_stage = RasusaDownsampleReadsByCoverage(
                ref_metadata=ref_metadata,
                coverage=self.downsample,
                reads=current_reads,
                **globals
            )
            # Update reads to use downsampled reads (convert to strings for Pydantic)
            current_reads = [downsample_stage.output.downsampled_r1]
            if downsample_stage.output.downsampled_r2:
                current_reads.append(downsample_stage.output.downsampled_r2)
            stages.append(downsample_stage)

        # Aligner
        if self.bam:
            if sample_name is None:
                sample_name = strip_bio_suffixes(Path(self.bam).name)
            aligned_reads = self.bam
            stages.append(BamReferenceValidator(
                bam=aligned_reads,
                reference=reference_file,
                reference_index=reference_index,
                **globals,
            ))
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
                sample_name = strip_bio_suffixes(Path(current_reads[0]).name)
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
            clair3_model = self.clair3_model
            if clair3_model is None:
                if not current_reads:
                    raise ValueError("Clair3 model can not be auto-detected without reads. Provide --clair3-model when using BAM/CRAM input.")
                longbow_stage = LongbowClair3ModelSelector(
                    reads=Path(current_reads[0]),
                    **globals
                )
                stages.append(longbow_stage)
                clair3_model = longbow_stage.output.clair3_model
            assert clair3_model is not None, "Clair3 model must be provided or resolved when using Clair3 caller."
            platform = 'ont'
            if 'hifi' in str(clair3_model).lower():
                platform = 'hifi'
            caller_stage = Clair3Caller(
                bam=aligned_reads,
                reference=reference_file,
                reference_index=reference_index,
                clair3_model=clair3_model,
                additional_options=self.caller_opts,
                platform=platform,
                **globals
            )
        else:
            caller_stage = FreebayesCallerLong(
                bam=aligned_reads,
                bam_index=align_filter.output.bai,
                reference=reference_file,
                reference_index=reference_index,
                min_mapping_quality=self.min_mapping_quality,
                additional_options=self.caller_opts,
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

        depth_beds = DepthBedsFromBam(
            bam=aligned_reads,
            min_depth=self.depth_mask,
            min_mapping_quality=self.min_mapping_quality if self.caller == "freebayes" else 0,
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

        variants_file = consequences.output.annotated_vcf
        if self.haploid:
            collapse_genotypes = CollapseDiploidGenotypes(
                vcf=variants_file,
                **globals,
            )
            stages.append(collapse_genotypes)
            variants_file = collapse_genotypes.output.vcf

        final_vcf = CopyFile(
            input=variants_file,
            output_path=Path(f"{self.prefix}.all.vcf"),
        )
        stages.append(final_vcf)
        variants_file = final_vcf.output.copied_file

        vcf_stats = VcfStats(
            vcf=variants_file,
            sample_name=sample_name,
            **globals
        )
        stages.append(vcf_stats)

        # Compress VCF
        gzip_vcf = VcfCompressor(
            input=variants_file,
            **globals
        )
        stages.append(gzip_vcf)
        
        # Filter to PASS-only variants
        pass_filter = VcfPassFilter(
            vcf=variants_file,
            **globals
        )
        stages.append(pass_filter)

        variants_tab = VcfToTab(
            vcf=pass_filter.output.vcf,
            **globals
        )
        stages.append(variants_tab)
        
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

        # Mask sites flagged with the MixedSite VCF filter using the full-call VCF
        mixed_sites = MaskMixedSites(
            fasta=current_fasta,
            vcf=variants_file,
            **globals,
        )
        stages.append(mixed_sites)
        current_fasta = mixed_sites.output.masked_fasta

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
                mask_bed=self.mask,
                **globals
            )
            stages.append(user_mask)
            current_fasta = user_mask.output.masked_fasta

        # Copy final masked consensus to standard output location
        copy_final = FinaliseFasta(
            input=current_fasta,
            output_path=f"{self.prefix}.fna",
            **globals
        )
        stages.append(copy_final)

        # Compress BAM to CRAM with embedded reference
        cram_compressor = CramCompressor(
            input=aligned_reads,
            reference=reference_file,
        )
        stages.append(cram_compressor)

        sample_report_stage = None
        if self.report:
            sample_report_stage = SampleReport(
                vcf=variants_file,
                alignment=aligned_reads,
                reference=reference_file,
                reference_index=reference_index,
                title=f"{sample_name or self.prefix} Sample Report",
                sample_name=sample_name,
                variant_scope=self.report_scope,
                window_size=self.report_window_size,
                **globals,
            )
            stages.append(sample_report_stage)

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
            variants_tab.output.tab,
            cram_compressor.output.cram,
            vcf_stats.output.summary_tsv,
            vcf_stats.output.breakdown_tsv,
        ]
        keep_files.extend(setup.output.paths)
        if sample_report_stage is not None:
            keep_files.append(sample_report_stage.output.rendered)
        if stats_tsv is not None:
            keep_files.append(stats_tsv)
        return SnippyPipeline(stages=stages, outputs_to_keep=keep_files)
