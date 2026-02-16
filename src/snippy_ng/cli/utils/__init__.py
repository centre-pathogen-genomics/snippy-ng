import click 
from pathlib import Path

def error(msg):
    click.echo(f"Error: {msg}", err=True)
    raise click.Abort()

class AbsolutePath(type(Path())):
    """
    pathlib path type that always resolves to an absolute path.

    Used as ``click.Path(path_type=AbsolutePath)``.
    """

    def __init__(self, *args, **kwargs):
        try:
            resolved = Path(*args, **kwargs).expanduser().absolute()
        except Exception as e:
            value = args[0] if args else ""
            raise click.BadParameter(f"Invalid path '{value}': {e}")

        super().__init__(str(resolved))