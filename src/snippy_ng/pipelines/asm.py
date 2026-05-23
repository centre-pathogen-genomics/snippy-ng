from pathlib import Path
from typing import Literal, Optional
from pydantic import Field
from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.pipelines import PipelineBuilder, SnippyPipeline
from snippy_ng.stages.vcf import AddDeletionsToVCF, AssemblyVariantContextFilter, VcfFilterAsm, VcfPassFilter, CollapseDiploidGenotypes
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
from snippy_ng.stages.compression import VcfCompressor
from snippy_ng.stages.masks import ApplyMask, MaskMixedSites
from snippy_ng.stages.copy import CopyFile, FinaliseFasta
from snippy_ng.pipelines.common import load_or_prepare_reference
from snippy_ng.stages.alignment import AssemblyAligner, AssemblyNucmerAligner
from snippy_ng.stages.calling import PAFCaller, ShowSnpsCaller
from snippy_ng.stages.reporting import PrintVcfHistogram, SampleReport
from snippy_ng.stages.stats import VcfStats
from snippy_ng.utils.gather import strip_bio_suffixes


class AsmPipelineBuilder(PipelineBuilder):
    """Builder for assembly-based SNP calling pipeline."""
    reference: Path = Field(..., description="Reference genome (FASTA or GenBank) or prepared reference directory")
    assembly: Path = Field(..., description="Assembly file path")
    prefix: str = Field(default="snippy", description="Output file prefix")
    caller: Literal["paftools", "nucmer"] = Field(default="nucmer", description="Caller to use for assembly-based SNP calling")
    caller_opts: str = Field(default="", description="Additional caller options")
    mask: Optional[str] = Field(default=None, description="BED file with regions to mask")
    sample_name: Optional[str] = Field(default=None, description="Optional sample name override for output tables")
    add_deletions_to_vcf: bool = Field(default=True, description="Add zero-depth regions to VCF as symbolic deletion blocks")
    minimap_preset: Literal["asm5", "asm10", "asm20"] = Field(default="asm20", description="Minimap2 preset for assembly alignment")
    min_qual: int = Field(default=60, description="Minimum QUAL score for variants to retain in VCF")
    haploid: bool = Field(default=True, description="Collapse diploid genotypes to haploid genotypes after consequence calling")
    report: bool = Field(default=True, description="Create a per-sample HTML report")
    report_scope: str = Field(default="all", description="Variant scope for the sample report: pass or all")

    def build(self) -> SnippyPipeline:
        """Build and return the assembly pipeline."""
        stages = []
        sample_name = self.sample_name or strip_bio_suffixes(Path(self.assembly).name)
        
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
        
        if self.caller == "nucmer":
            aligner = AssemblyNucmerAligner(
                reference=reference_file,
                assembly=Path(self.assembly),
                prefix=self.prefix,
            )
            stages.append(aligner)
            caller = ShowSnpsCaller(
                delta=aligner.output.delta,
                ref_dict=setup.output.reference_dict,
                assembly=Path(self.assembly),
                reference=reference_file,
                reference_index=reference_index,
                additional_options=self.caller_opts,
                prefix=self.prefix,
            )
        else:
            aligner = AssemblyAligner(
                reference=reference_file,
                assembly=Path(self.assembly),
                minimap_preset=self.minimap_preset,
                prefix=self.prefix,
            )
            stages.append(aligner)
            caller = PAFCaller(
                paf=aligner.output.paf,
                ref_dict=setup.output.reference_dict,
                reference=reference_file,
                reference_index=reference_index,
                additional_options=self.caller_opts,
                prefix=self.prefix,
            )
        stages.append(caller)

        # Filter VCF
        variant_filter = VcfFilterAsm(
            vcf=caller.output.vcf,
            reference=reference_file,
            # hard code for asm-based calling
            min_qual=self.min_qual,
            prefix=self.prefix
        )
        stages.append(variant_filter)
        variants_file = variant_filter.output.vcf

        if self.add_deletions_to_vcf:
            # Add zero-depth regions to VCF as symbolic deletion blocks
            add_deletions = AddDeletionsToVCF(
                zero_depth_bed=caller.output.missing_bed,
                vcf=variants_file,
                reference=reference_file,
                prefix=self.prefix
            )
            stages.append(add_deletions)
            variants_file = add_deletions.output.vcf

        context_filter = AssemblyVariantContextFilter(
            vcf=variants_file,
            prefix=self.prefix,
        )
        stages.append(context_filter)

        # Consequences calling
        consequences = BcftoolsConsequencesCaller(
            variants=context_filter.output.vcf,
            features=features_file,
            reference=reference_file,
            prefix=self.prefix
        )
        stages.append(consequences)

        variants_file = consequences.output.annotated_vcf
        if self.haploid:
            collapse_genotypes = CollapseDiploidGenotypes(
                vcf=variants_file,
                prefix=self.prefix,
            )
            stages.append(collapse_genotypes)
            variants_file = collapse_genotypes.output.vcf

        final_vcf = CopyFile(
            input=variants_file,
            output_path=f"{self.prefix}.all.vcf",
        )
        stages.append(final_vcf)
        variants_file = final_vcf.output.copied_file

        vcf_stats = VcfStats(
            vcf=variants_file,
            sample_name=sample_name,
            prefix=self.prefix
        )
        stages.append(vcf_stats)
        
        # Filter to PASS-only variants
        pass_filter = VcfPassFilter(
            vcf=variants_file,
            prefix=self.prefix
        )
        stages.append(pass_filter)

        # Compress VCF
        gzip_vcf = VcfCompressor(
            input=variants_file,
            prefix=self.prefix,
        )
        stages.append(gzip_vcf)
        
        # Pseudo-alignment
        pseudo = BcftoolsPseudoAlignment(
            ref_metadata=ref_metadata,
            vcf=pass_filter.output.vcf,
            reference=reference_file,
            prefix=self.prefix
        )
        stages.append(pseudo)
        
        # Track the current reference/fasta through the masking stages
        current_fasta = pseudo.output.fasta

        # Mask sites flagged with the MixedSite VCF filter using the full-call VCF
        # I think this is very unlikely in assembly-based calling but there is a
        # theoretical possibility of contigs being collapsed in the assembly and showing 
        # up as heterozygous sites in the VCF, so we will mask them just in case.
        mixed_sites = MaskMixedSites(
            fasta=current_fasta,
            vcf=variants_file,
            prefix=self.prefix
        )
        stages.append(mixed_sites)
        current_fasta = mixed_sites.output.masked_fasta

        # Apply user mask if provided
        if self.mask:
            user_mask = ApplyMask(
                fasta=current_fasta,
                mask_bed=Path(self.mask),
                prefix=self.prefix
            )
            stages.append(user_mask)
            current_fasta = user_mask.output.masked_fasta

        # Copy final masked consensus to standard output location
        copy_final = FinaliseFasta(
            input=current_fasta,
            output_path=f"{self.prefix}.fna",
            prefix=self.prefix
        )
        stages.append(copy_final)

        # Print VCF histogram to terminal
        vcf_histogram = PrintVcfHistogram(
            vcf_path=variants_file,
            prefix=self.prefix
        )
        stages.append(vcf_histogram)

        sample_report_stage = None
        if self.report:
            sample_report_stage = SampleReport(
                vcf=variants_file,
                title="Snippy-NG Sample Report",
                sample_name=sample_name,
                variant_scope=self.report_scope,
                window_size=100,
                prefix=self.prefix,
            )
            stages.append(sample_report_stage)
        
        keep_files = [
            copy_final.output.fasta, 
            gzip_vcf.output.gz,
            pass_filter.output.vcf,
            vcf_stats.output.summary_tsv,
            vcf_stats.output.breakdown_tsv,
        ]
        if sample_report_stage is not None:
            keep_files.append(sample_report_stage.output.rendered)
        return SnippyPipeline(stages=stages, outputs_to_keep=keep_files)
