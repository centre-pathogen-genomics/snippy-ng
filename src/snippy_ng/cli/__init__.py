import click

from snippy_ng.__about__ import __version__, EXE
from snippy_ng.cli.run import run

def show_citation(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(f"Please cite {EXE} in your research: ...")
    ctx.exit()

def version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(f"{EXE} version {__version__}")
    ctx.exit()


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(version=__version__, prog_name=EXE)
@click.option("--citation", is_flag=True, callback=show_citation, expose_value=False, help="Print citation for referencing Snippy-NG.")
def snippy_ng():
    """
    Snippy-NG: The Next Generation of Variant Calling.
    """

########################
# Register Subcommands #
########################
snippy_ng.add_command(run)