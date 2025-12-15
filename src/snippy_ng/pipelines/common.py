"""Utilities for reference handling in CLI commands."""

from pathlib import Path
from typing import Optional
from snippy_ng.stages.setup import LoadReference, PrepareReference
from snippy_ng.stages.filtering import VcfFilter
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
from snippy_ng.stages.compression import BgzipCompressor
from snippy_ng.stages.masks import ApplyMask, HetMask
from snippy_ng.stages.copy import CopyFasta
from snippy_ng.seq_utils import guess_format
from snippy_ng.cli.utils import error


def load_or_prepare_reference(reference_path, reference_prefix="ref") -> PrepareReference | LoadReference:
    """
    Load an existing reference directory or prepare a new reference from a FASTA/GenBank file.
    
    Args:
        reference_path: Path to reference file or directory
        reference_prefix: Prefix for output reference files
        
    Returns:
        An instance of LoadReference or PrepareReference stage.
        
    Raises:
        SystemExit: If reference format cannot be determined
    """
    if Path(reference_path).is_dir():
        setup = LoadReference(
            reference_dir=reference_path,
            reference_prefix=reference_prefix,
        )
    else:
        reference_format = guess_format(reference_path)
        if not reference_format:
            error(f"Could not determine format of reference file '{reference_path}'")

        # Determine reference directory - use outdir/reference if outdir provided, otherwise just "reference"
        reference_dir = Path("reference")

        setup = PrepareReference(
            input=reference_path,
            ref_fmt=reference_format,
            reference_prefix=reference_prefix,
            reference_dir=reference_dir,
        )
    
    return setup


def filter_annotate_and_generate_consensus(
    stages: list,
    vcf: Path,
    reference: Path,
    reference_index: Path,
    features: Path,
    header: Optional[str],
    min_depth: int,
    min_qual: float,
    globals: dict,
) -> Path:
    """
    Filter VCF, annotate with consequences, compress, and generate pseudo-alignment.
    
    This is common logic shared between asm and short pipelines for processing
    variant calls into a consensus sequence.
    
    Args:
        stages: List to append stages to
        vcf: Raw VCF file from variant caller
        reference: Reference file
        reference_index: Reference index file
        features: GFF features file
        header: Optional header for the output FASTA
        min_depth: Minimum depth for filtering
        min_qual: Minimum quality for filtering
        globals: Dictionary containing prefix, cpus, ram, outdir, tmpdir
        
    Returns:
        Path to the generated pseudo-alignment FASTA file
    """
    # Filter VCF
    variant_filter = VcfFilter(
        vcf=vcf,
        reference=reference,
        min_depth=min_depth,
        min_qual=min_qual,
        **globals
    )
    stages.append(variant_filter)
    variants_file = variant_filter.output.vcf
    
    # Consequences calling
    consequences = BcftoolsConsequencesCaller(
        variants=variants_file,
        features=features,
        reference=reference,
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
        reference=reference,
        reference_index=reference_index,
        header=header,
        **globals
    )
    stages.append(pseudo)
    
    return pseudo.output.fasta


def apply_consensus_masking(
    stages: list,
    fasta: Path,
    vcf: Path,
    reference: Path,
    mask: Optional[str],
    min_qual: float,
    globals: dict,
) -> list:
    """
    Apply heterozygous/low-quality site masking, optional user mask, and copy final consensus.
    
    This is common logic shared between asm and short pipelines after their respective
    depth/missing masking stages.
    
    Args:
        stages: List to append stages to
        fasta: Current FASTA file (after any depth masking)
        vcf: Raw VCF file for heterozygous site detection
        reference: Reference file
        mask: Optional user mask BED file path
        min_qual: Minimum quality threshold for masking
        globals: Dictionary containing prefix, cpus, ram, outdir, tmpdir
        
    Returns:
        Updated stages list
    """
    current_fasta = fasta
    
    # Apply heterozygous and low quality sites masking
    het_mask = HetMask(
        vcf=vcf,  # Use raw VCF for complete site information
        fasta=current_fasta,
        reference=reference,
        min_qual=min_qual,
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
    prefix = globals['prefix']
    copy_final = CopyFasta(
        input=current_fasta,
        output_path=f"{prefix}.pseudo.fna",
        **globals
    )
    stages.append(copy_final)
    
    return stages
