import click
import tempfile
from pathlib import Path


class GlobalOption(click.Option):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_global = True

GLOBAL_DEFS = [
     {
        "param_decls": ("--outdir", "-o"),
        "attrs": {
            "type": click.Path(writable=True, readable=True, file_okay=False, dir_okay=True, path_type=Path),
            "required": True,
            "default": Path("out"),
            "help": "Where to put everything",
        },
    },
    {
        "param_decls": ("--tmpdir", "-t"),
        "attrs": {
            "type": click.Path(writable=True, readable=True, file_okay=False, dir_okay=True, path_type=Path),
            "default": tempfile.gettempdir(),
            "help": "Temporary directory for fast storage",
        },
    },
    {
        "param_decls": ("--cpus", "-c"),
        "attrs": {
            "type": int,
            "default": 1,
            "help": "Max cores to use",
        },
    },
    {
        "param_decls": ("--ram", "-r"),
        "attrs": {
            "type": int,
            "default": 8,
            "help": "Try and keep RAM under this many GB",
        },
    },
    {
        "param_decls": ("--force", "-f"),
        "attrs": {
            "is_flag": True,
            "default": False,
            "help": "Overwrite existing output directory",
        },
    },
    {
        "param_decls": ("--quiet", "-q"),
        "attrs": {
            "is_flag": True,
            "default": False,
            "help": "Capture tool output",
        },
    },
    {
        "param_decls": ("--skip-check", "-k"),
        "attrs": {
            "is_flag": True,
            "default": False,
            "help": "Skip dependency checks",
        },
    },
    {
        "param_decls": ("--check", "-d"),
        "attrs": {
            "is_flag": True,
            "default": False,
            "help": "Check dependencies are installed then exit",
        },
    },
    {
        "param_decls": ("--continue",),
        "attrs": {
            "is_flag": True,
            "default": False,
            "help": "Continue from the last run, skipping completed stages",
        },
    },
    {
        "param_decls": ("--keep-incomplete",),
        "attrs": {
            "is_flag": True,
            "default": False,
            "help": "Keep outputs from incomplete stages if an error occurs",
        },
    },
]


def snippy_global_options(f):
    """
    Decorator that prepends each GLOBAL_DEFS entry as
    @click.option(..., cls=GlobalOption, **attrs).
    """
    for entry in reversed(GLOBAL_DEFS):
        param_decls = entry["param_decls"]
        attrs = entry["attrs"]
        option_decorator = click.option(*param_decls, cls=GlobalOption, **attrs)
        f = option_decorator(f)
    return f


class CommandWithGlobals(click.Command):
    def format_options(self, ctx, formatter):
        global_opts = []
        other_opts = []
        for param in self.params:
            if getattr(param, 'is_global', False):
                global_opts.append(param)
            else:
                other_opts.append(param)

        if global_opts:
            with formatter.section('Globals'):
                rows = [p.get_help_record(ctx) for p in global_opts if p.get_help_record(ctx)]
                formatter.write_dl(rows)

        if other_opts:
            with formatter.section('Options'):
                rows = [p.get_help_record(ctx) for p in other_opts if p.get_help_record(ctx)]
                formatter.write_dl(rows)
