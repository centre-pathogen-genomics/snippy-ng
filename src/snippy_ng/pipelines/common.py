from __future__ import annotations

from pathlib import Path
from typing import Optional
from snippy_ng.stages.setup import LoadReferenceFromMetadataFile, PrepareReference
from snippy_ng.stages.vcf import VcfFilterShort
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
from snippy_ng.stages.consensus import BcftoolsPseudoAlignment
from snippy_ng.stages.compression import BgzipCompressor
from snippy_ng.stages.masks import ApplyMask, QualMask
from snippy_ng.stages.copy import CopyFasta
from snippy_ng.utils.seq import guess_reference_format
from snippy_ng.exceptions import InvalidReferenceError


def load_or_prepare_reference(
        reference_path,
        output_directory: Path = Path("reference"),
        ) -> PrepareReference | LoadReferenceFromMetadataFile:
    """
    Load an existing reference directory or prepare a new reference from a FASTA/GenBank file.
    
    Args:
        reference_path: Path to reference file, directory, or metadata.json file.
        
    Returns:
        An instance of LoadReferenceFromMetadataFile or PrepareReference stage.
        
    Raises:
        SystemExit: If reference format cannot be determined
    """
    if Path(reference_path).is_dir():
        # check for metadata.json in directory
        metadata = Path(reference_path) / "metadata.json"
        if not metadata.exists():
            raise InvalidReferenceError(f"No metadata.json found in reference directory '{reference_path}'. Ensure you are providing a valid reference directory.")
        setup = LoadReferenceFromMetadataFile(
            metadata=metadata
        )
    elif Path(reference_path).suffix.lower() == ".json":
        setup = LoadReferenceFromMetadataFile(
            metadata=Path(reference_path)
        )
    else:
        setup = prepare_reference(reference_path, output_directory)
    
    return setup


def prepare_reference(reference_path, output_directory) -> PrepareReference:
    """
    Prepare a new reference from a FASTA/GenBank file.
    
    Args:
        reference_path: Path to reference file.
    Returns:
        An instance of PrepareReference stage.
    reference_format = guess_format(reference_path)
    """
    reference_format = guess_reference_format(reference_path)
    if not reference_format:
        raise InvalidReferenceError(f"Could not determine reference format for '{reference_path}'. Supported formats are FASTA and GenBank.")

    setup = PrepareReference(
        input=reference_path,
        ref_fmt=reference_format,
        directory=output_directory,
    )
    
    return setup

