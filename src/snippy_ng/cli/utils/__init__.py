import click
from pathlib import Path
import re

class AbsolutePath(click.Path):
    """Click path type that always returns absolute `Path` values without resolving symlinks."""

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("path_type", Path)
        super().__init__(*args, **kwargs)

    def convert(self, value, param, ctx):
        p = super().convert(value, param, ctx)
        if p is None:
            return None
        return Path(p).absolute()

def error(msg):
    click.echo(f"Error: {msg}", err=True)
    raise click.Abort()


def absolute_path(path):
    """Convert a path or sequence of paths to absolute Path objects."""
    if path is None:
        return None

    # click passes tuple/list for nargs=-1 arguments
    if isinstance(path, tuple):
        return tuple(absolute_path(p) for p in path)
    if isinstance(path, list):
        return [absolute_path(p) for p in path]

    if not path:
        return path
    return Path(path).absolute()

def absolute_path_callback(ctx, param, value):
    if ctx.resilient_parsing:
        return value
    return absolute_path(value)


def resolve_cli_input(option_value, arg_value, *, option_name: str, arg_name: str):
    """Resolve a value that may be provided by either an option or a positional argument."""
    if option_value is not None and arg_value is not None and option_value != arg_value:
        raise click.UsageError(
            f"Conflicting values provided for {option_name} and {arg_name}. Please use only one or provide the same value."
        )
    return option_value if option_value is not None else arg_value


def _path_or_accession_callback(kind: str, accession_hint: str | None = None):
    def callback(ctx, param, value):
        if ctx.resilient_parsing or value is None:
            return value

        path = absolute_path(value)
        if path.exists():
            return path

        from snippy_ng.pipelines.common import is_reference_accession

        if is_reference_accession(value):
            return value

        hint = f" {accession_hint}" if accession_hint else ""
        raise click.BadParameter(
            f"{kind} file '{path}' does not exist.{hint}",
            ctx=ctx,
            param=param,
        )

    return callback


reference_or_accession_callback = _path_or_accession_callback(
    "Reference",
    accession_hint="If you meant to provide an assembly accession, use NCBI GCF/GCA or AllTheBacteria SAMN/SAMEA format (e.g. GCF_000174395.2).",
)
assembly_or_accession_callback = _path_or_accession_callback(
    "Assembly",
    accession_hint="If you meant to provide an assembly accession, use NCBI GCF/GCA or AllTheBacteria SAMN/SAMEA format (e.g. GCF_000174395.2).",
)


SRA_ACCESSION_RE = re.compile(r"^(SRR|ERR|DRR)\d{6,}$")

def is_sra_accession(read_accession) -> bool:
    """Check if a string is a valid SRA accession (SRR, ERR, or DRR format)."""
    return SRA_ACCESSION_RE.fullmatch(str(read_accession)) is not None


def reads_or_accession_callback(ctx, param, value):
    """Callback for read files or SRA accessions."""
    if ctx.resilient_parsing or value is None:
        return value

    # If it's a tuple (multiple files), convert each element individually
    if isinstance(value, tuple):
        return tuple(reads_or_accession_callback(ctx, param, v) for v in value)

    # Try as path first
    path = absolute_path(value)
    if path.exists():
        return path

    if is_sra_accession(value):
        return value

    raise click.BadParameter(
        f"Reads file '{path}' does not exist. If you meant to provide an SRA accession, use SRR/ERR/DRR format (e.g. SRR1234567).",
        ctx=ctx,
        param=param,
    )
