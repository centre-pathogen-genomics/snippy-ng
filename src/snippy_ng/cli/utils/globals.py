import click
import tempfile
from pathlib import Path
import os

from snippy_ng.cli.utils import absolute_path


class GlobalOption(click.Option):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_global = True

def debug_callback(ctx, param, value):
    if value is None or ctx.resilient_parsing:
        return
    if not value:
        # don't update env
        return value
    os.environ["SNIPPY_NG_DEBUG"] = "1"
    return value

def create_outdir_callback(ctx, param, value):
    if ctx.resilient_parsing:
        return
    if value is None:
        return value
    value = absolute_path(value)
    if value.exists() and not ctx.params.get("force", False):
        raise click.UsageError(f"Output folder '{value}' already exists! Use --force to overwrite.")
    if not value.exists():
        value.mkdir(parents=True, exist_ok=True)
    return value

def not_implemented_callback(ctx, param, value):
    if ctx.resilient_parsing:
        return
    if value:
        raise NotImplementedError(f"The option '{param.name}' is not yet implemented.")
    return value

def cap_cpus_callback(ctx, param, value):
    if ctx.resilient_parsing:
        return
    available = os.cpu_count()
    return min(available, value) if available else value

GLOBAL_DEFS = [
    {
        "param_decls": ("--outdir", "-o"),
        "attrs": {
            "type": click.Path(writable=True, readable=True, file_okay=False, dir_okay=True),
            "default": Path("out"),
            "help": "Where to put everything",
            "callback": create_outdir_callback,
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
        "param_decls": ("--prefix", "-p"),
        "attrs": {
            "type": click.STRING,
            "default": "snps",
            "help": "Prefix for output files",
        },
    },
    {
        "param_decls": ("--cpus", "-c"),
        "attrs": {
            "type": int,
            "default": 1,
            "help": "Max cores to use",
            "callback": cap_cpus_callback,
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
            "is_eager": True,
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
        "param_decls": ("--debug",),
        "attrs": {
            "is_flag": True,
            "default": False,
            "help": "Print debug output",
            "callback": debug_callback, 
        }, 
    },
    {
        "param_decls": ("--skip-check",),
        "attrs": {
            "is_flag": True,
            "default": False,
            "help": "Skip dependency checks",
        },
    },
    {
        "param_decls": ("--check",),
        "attrs": {
            "is_flag": True,
            "default": False,
            "help": "Check dependencies are installed then exit",
        },
    },
    {
        "param_decls": ("--create-missing",),
        "attrs": {
            "is_flag": True,
            "default": False,
            "help": "Continue from last run by creating missing outputs",
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
    {
        "param_decls": ("--no-cleanup",),
        "attrs": {
            "is_flag": True,
            "default": False,
            "help": "Do not delete temporary files after run",
            "callback": not_implemented_callback
        },
    }
]


def add_snippy_global_options(exclude: list = None):
    """
    Decorator that prepends each GLOBAL_DEFS entry as
    @click.option(..., cls=GlobalOption, **attrs).
    """
    if exclude is None:
        exclude = []
    def wraps(f):
        # Click stores params after decoration, so inspect __click_params__ first
        existing = {
            param.name
            for param in getattr(f, "__click_params__", [])
            if isinstance(param, click.Option)
        }

        for entry in reversed(GLOBAL_DEFS):
            param_decls = entry["param_decls"]
            attrs = entry["attrs"]

            # Click derives the internal name like this:
            option = click.Option(param_decls)
            option_name = option.name

            if option_name in existing or option_name in exclude:
                continue

            f = click.option(*param_decls, cls=GlobalOption, **attrs)(f)

        return f
    return wraps


class CommandWithGlobals(click.Command):
    def format_options(self, ctx, formatter):
        global_opts = []
        other_opts = []
        for param in self.params:
            # Only Options belong in format_options()
            if not isinstance(param, click.Option):
                continue
            if getattr(param, 'is_global', False):
                global_opts.append(param)
            else:
                other_opts.append(param)

        if global_opts:
            order = {opts['param_decls'][0]: i for i, opts in enumerate(GLOBAL_DEFS)}
            global_opts.sort(key=lambda p: order.get(p.opts[0], float('inf')))
            with formatter.section('Globals'):
                rows = [p.get_help_record(ctx) for p in global_opts if p.get_help_record(ctx)]
                formatter.write_dl(rows)

        if other_opts:
            with formatter.section('Options'):
                rows = [p.get_help_record(ctx) for p in other_opts if p.get_help_record(ctx)]
                formatter.write_dl(rows)
