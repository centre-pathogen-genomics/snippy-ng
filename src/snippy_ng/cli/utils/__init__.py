import click
from pathlib import Path

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