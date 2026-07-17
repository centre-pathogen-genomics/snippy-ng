from __future__ import annotations

from pathlib import Path
import re

from snippy_ng.stages.setup import LoadReferenceFromMetadataFile, PrepareReference
from snippy_ng.utils.seq import guess_reference_format
from snippy_ng.exceptions import InvalidReferenceError

NCBI_ASSEMBLY_ACCESSION_RE = re.compile(r"^(GC[AF]_\d{9})(?:\.\d+)?$")
ATB_ASSEMBLY_ACCESSION_RE = re.compile(r"^(?:SAM(E|D|N)[A-Z]?)[0-9]+$")
SRA_ACCESSION_RE = re.compile(r"^(SRR|ERR|DRR)\d{6,}$")


def is_ncbi_assembly_accession(reference) -> bool:
    return NCBI_ASSEMBLY_ACCESSION_RE.fullmatch(str(reference)) is not None


def is_atb_assembly_accession(reference) -> bool:
    return ATB_ASSEMBLY_ACCESSION_RE.fullmatch(str(reference)) is not None


def is_reference_accession(reference) -> bool:
    return is_ncbi_assembly_accession(reference) or is_atb_assembly_accession(reference)


def is_sra_accession(read_accession) -> bool:
    """Check if a string is a valid SRA accession (SRR, ERR, or DRR format)."""
    return SRA_ACCESSION_RE.fullmatch(str(read_accession)) is not None


def download_reads(read_accession, stages: list, output_directory: Path | None = None) -> Path:
    """Download reads from SRA accession using sracha-rs."""
    from snippy_ng.stages.download import DownloadSraReads

    if is_sra_accession(read_accession):
        download_stage = DownloadSraReads(
            accession=str(read_accession),
            output_directory=output_directory,
        )
        stages.append(download_stage)
        return download_stage.output.reads

    raise InvalidReferenceError(
        f"Unsupported SRA accession '{read_accession}'. Supported formats are SRR, ERR, and DRR accessions."
    )


def get_download_stage_outputs(stages: list) -> list:
    """Return output paths from any download stages in the pipeline."""
    from snippy_ng.stages.download import DownloadAtbAssemblyReference, DownloadNcbiAssemblyFasta, DownloadSraReads
    outputs = []
    for stage in stages:
        if isinstance(stage, (DownloadNcbiAssemblyFasta, DownloadAtbAssemblyReference, DownloadSraReads)):
            outputs.extend(stage.output.paths)
    return outputs


def download_assembly(reference_accession, stages: list, output_directory: Path | None = None) -> Path:
    from snippy_ng.stages.download import DownloadAtbAssemblyReference, DownloadNcbiAssemblyFasta

    if is_ncbi_assembly_accession(reference_accession):
        download_reference = DownloadNcbiAssemblyFasta(
            accession=str(reference_accession),
            output_directory=output_directory,
            genbank=True,
        )
        stages.append(download_reference)
        return download_reference.output.fasta

    if is_atb_assembly_accession(reference_accession):
        download_reference = DownloadAtbAssemblyReference(
            accession=str(reference_accession),
            output_directory=output_directory,
        )
        stages.append(download_reference)
        return download_reference.output.fasta

    raise InvalidReferenceError(
        f"Unsupported assembly accession '{reference_accession}'. Supported accessions are NCBI GCF/GCA and AllTheBacteria SAMN/SAMEA/SAMD IDs."
    )


def download_assembly_fasta(assembly_accession, stages: list, output_directory: Path | None = None) -> Path:
    from snippy_ng.stages.download import DownloadAtbAssemblyReference, DownloadNcbiAssemblyFasta

    if is_ncbi_assembly_accession(assembly_accession):
        download_stage = DownloadNcbiAssemblyFasta(
            accession=str(assembly_accession),
            output_directory=output_directory,
        )
        stages.append(download_stage)
        return download_stage.output.fasta

    if is_atb_assembly_accession(assembly_accession):
        download_stage = DownloadAtbAssemblyReference(
            accession=str(assembly_accession),
            output_directory=output_directory,
        )
        stages.append(download_stage)
        return download_stage.output.fasta

    raise InvalidReferenceError(
        f"Unsupported assembly accession '{assembly_accession}'. Supported accessions are NCBI GCF/GCA and AllTheBacteria SAMN/SAMEA/SAMD IDs."
    )


def load_or_prepare_reference(
        reference_path,
        output_directory: Path | None = None,
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
        setup = prepare_reference(reference_path, output_directory=output_directory)
    
    return setup


def prepare_reference(reference_path, output_directory: Path | None = None) -> PrepareReference:
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
        output_directory=output_directory,
    )
    
    return setup
