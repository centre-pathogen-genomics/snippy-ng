from pathlib import Path
from typing import Optional
from snippy_ng.stages.filtering import VcfFilter
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
from snippy_ng.stages.compression import BgzipCompressor
from snippy_ng.stages.masks import ApplyMask, HetMask
from snippy_ng.stages.copy import CopyFasta
from snippy_ng.pipelines.common import load_or_prepare_reference
from snippy_ng.stages.alignment import AssemblyAligner
from snippy_ng.stages.calling import PAFCaller


def create_asm_pipeline_stages(
    reference: str,
    assembly: str,
    prefix: str = "snps",
    mask: Optional[str] = None,
    tmpdir: Path = Path("/tmp"),
    cpus: int = 1,
    ram: int = 8,
) -> list:
    stages = []
    globals = {'prefix': prefix, 'cpus': cpus, 'ram': ram, 'tmpdir': tmpdir}
    
    # Setup reference (load existing or prepare new)
    setup = load_or_prepare_reference(
        reference_path=reference,
        reference_prefix=prefix
    )
    reference_file = setup.output.reference
    features_file = setup.output.gff
    reference_index = setup.output.reference_index
    stages.append(setup)
    
    # Aligner 
    aligner = AssemblyAligner(
        reference=reference_file,
        assembly=Path(assembly),
        **globals
    )
    stages.append(aligner)
    
    # Call variants
    caller = PAFCaller(
        paf=aligner.output.paf,
        ref_dict=setup.output.reference_dict,
        reference=reference_file,
        reference_index=reference_index,
        **globals
    )
    stages.append(caller)
    
    # Filter VCF
    variant_filter = VcfFilter(
        vcf=caller.output.vcf,
        reference=reference_file,
        # hard code for asm-based calling
        min_depth=1,
        min_qual=60,
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
        vcf_gz=gzip.output.compressed,
        reference=reference_file,
        **globals
    )
    stages.append(pseudo)
    
    # Track the current reference/fasta through the masking stages
    current_fasta = pseudo.output.fasta
    
    # Apply depth masking
    missing_mask = ApplyMask(
        fasta=current_fasta,
        mask_bed=caller.output.missing_bed,
        mask_char="-",
        **globals
    )
    stages.append(missing_mask)
    current_fasta = missing_mask.output.masked_fasta

    # Apply heterozygous and low quality sites masking
    het_mask = HetMask(
        vcf=caller.output.vcf,  # Use raw VCF for complete site information
        fasta=current_fasta,
        **globals
    )
    stages.append(het_mask)
    current_fasta = het_mask.output.masked_fasta
    
    # Apply user mask if provided
    if mask:
        user_mask = ApplyMask(
            fasta=current_fasta,
            mask_bed=Path(mask),
            **globals
        )
        stages.append(user_mask)
        current_fasta = user_mask.output.masked_fasta

    # Copy final masked consensus to standard output location
    copy_final = CopyFasta(
        input=current_fasta,
        output_path=f"{prefix}.pseudo.fna",
        **globals
    )
    stages.append(copy_final)
    
    return stages