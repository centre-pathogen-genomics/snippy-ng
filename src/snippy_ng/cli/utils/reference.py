"""Utilities for reference handling in CLI commands."""

from pathlib import Path
from snippy_ng.stages.setup import LoadReference, PrepareReference
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
