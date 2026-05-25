import click
from pathlib import Path


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


def reference_or_accession_callback(ctx, param, value):
    if ctx.resilient_parsing or value is None:
        return value

    reference = absolute_path(value)
    if reference.exists():
        return reference

    from snippy_ng.pipelines.common import is_reference_accession

    if is_reference_accession(value):
        return value
    
    raise click.BadParameter(
        f"Reference path does not exist: {reference}",
        ctx=ctx,
        param=param,
    )
