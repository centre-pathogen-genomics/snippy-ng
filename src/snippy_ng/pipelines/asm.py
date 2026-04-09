from pathlib import Path
from typing import Optional
from pydantic import Field
from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.pipelines import PipelineBuilder, SnippyPipeline
from snippy_ng.stages.vcf import AddDeletionstoVCF, VcfFilterAsm, VcfPassFilter
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
from snippy_ng.stages.compression import VcfCompressor
from snippy_ng.stages.masks import ApplyMask, QualMask
from snippy_ng.stages.copy import FinaliseFasta
from snippy_ng.pipelines.common import load_or_prepare_reference
from snippy_ng.stages.alignment import AssemblyAligner
from snippy_ng.stages.calling import PAFCaller
from snippy_ng.stages.reporting import PrintVcfHistogram
from snippy_ng.stages.stats import VcfStats
from snippy_ng.utils.gather import guess_sample_id


class AsmPipelineBuilder(PipelineBuilder):
    """Builder for assembly-based SNP calling pipeline."""
    reference: Path = Field(..., description="Reference genome file path")
    assembly: Path = Field(..., description="Assembly file path")
    prefix: str = Field(default="snippy", description="Output file prefix")
    mask: Optional[str] = Field(default=None, description="BED file with regions to mask")
    sample_name: Optional[str] = Field(default=None, description="Optional sample name override for output tables")

    def build(self) -> SnippyPipeline:
        """Build and return the assembly pipeline."""
        stages = []
        sample_name = self.sample_name or guess_sample_id(Path(self.assembly).name)
        
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
        
        # Aligner 
        aligner = AssemblyAligner(
            reference=reference_file,
            assembly=Path(self.assembly),
            prefix=self.prefix
        )
        stages.append(aligner)
        
        # Call variants
        caller = PAFCaller(
            paf=aligner.output.paf,
            ref_dict=setup.output.reference_dict,
            reference=reference_file,
            reference_index=reference_index,
            prefix=self.prefix
        )
        stages.append(caller)
        
        # Filter VCF
        variant_filter = VcfFilterAsm(
            vcf=caller.output.vcf,
            reference=reference_file,
            # hard code for asm-based calling
            min_qual=60,
            prefix=self.prefix
        )
        stages.append(variant_filter)
        variants_file = variant_filter.output.vcf

        # Add zero-depth regions to VCF as symbolic deletion blocks
        add_deletions = AddDeletionstoVCF(
            zero_depth_bed=caller.output.missing_bed,
            vcf=variants_file,
            reference=reference_file,
            prefix=self.prefix
        )
        stages.append(add_deletions)
        variants_file = add_deletions.output.vcf

        # Consequences calling
        consequences = BcftoolsConsequencesCaller(
            variants=variants_file,
            features=features_file,
            reference=reference_file,
            prefix=self.prefix
        )
        stages.append(consequences)

        vcf_stats = VcfStats(
            vcf=consequences.output.annotated_vcf,
            sample_name=sample_name,
            prefix=self.prefix
        )
        stages.append(vcf_stats)
        
        # Filter to PASS-only variants for consensus generation
        pass_filter = VcfPassFilter(
            vcf=consequences.output.annotated_vcf,
            prefix=self.prefix
        )
        stages.append(pass_filter)

        # Compress VCF
        gzip_vcf = VcfCompressor(
            input=consequences.output.annotated_vcf,
            prefix=self.prefix
        )
        stages.append(gzip_vcf)
        
        # Pseudo-alignment
        pseudo = BcftoolsPseudoAlignment(
            ref_metadata=ref_metadata,
            vcf_gz=gzip_vcf.output.gz,
            reference=reference_file,
            prefix=self.prefix
        )
        stages.append(pseudo)
        
        # Track the current reference/fasta through the masking stages
        current_fasta = pseudo.output.fasta

        # Apply heterozygous and low quality sites masking
        het_mask = QualMask(
            vcf=caller.output.vcf,  # Use raw VCF for complete site information
            fasta=current_fasta,
            prefix=self.prefix
        )
        stages.append(het_mask)
        current_fasta = het_mask.output.masked_fasta
        
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
            output_path=f"{self.prefix}.pseudo.fna",
            prefix=self.prefix
        )
        stages.append(copy_final)

        # Print VCF histogram to terminal
        vcf_histogram = PrintVcfHistogram(
            vcf_path=variants_file,
            prefix=self.prefix
        )
        stages.append(vcf_histogram)
        
        keep_files = [
            copy_final.output.fasta, 
            gzip_vcf.output.gz,
            pass_filter.output.vcf,
            vcf_stats.output.summary_tsv,
            vcf_stats.output.breakdown_tsv,
        ]
        return SnippyPipeline(stages=stages, outputs_to_keep=keep_files)
