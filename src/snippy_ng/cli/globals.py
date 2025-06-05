import click
import tempfile


class GlobalOption(click.Option):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_global = True

GLOBAL_DEFS = [
    {
        "param_decls": ("--threads",),
        "attrs": {
            "type": int,
            "default": 1,
            "help": "Max cores to use",
        },
    },
    {
        "param_decls": ("--ram",),
        "attrs": {
            "type": int,
            "default": None,
            "help": "Try and keep RAM under this many GB",
        },
    },
    {
        "param_decls": ("--outdir", "-o"),
        "attrs": {
            "type": click.Path(),
            "default": '.',
            "help": "Where to put everything",
        },
    },
    {
        "param_decls": ("--tmpdir",),
        "attrs": {
            "type": click.Path(),
            "default": tempfile.gettempdir(),
            "help": "Temporary directory for fast storage (e.g., local SSD)",
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
            "help": "Minimal screen output",
        },
    }
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
