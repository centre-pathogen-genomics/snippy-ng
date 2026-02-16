import click 
from pathlib import Path
from typing import Optional, Union

def error(msg):
    click.echo(f"Error: {msg}", err=True)
    raise click.Abort()

def absolute_path(path: Optional[Union[str, Path]]) -> Optional[Path]:
    if not path:
        return path
    return Path(path).absolute()

def absolute_path_callback(ctx, param, value):
    if ctx.resilient_parsing:
        return
    return absolute_path(value)