import click
import os

from snippy_ng.__about__ import __version__, EXE, GITHUB_URL
from snippy_ng.cli.og_cli import og
from snippy_ng.cli.short_cli import short
from snippy_ng.cli.asm_cli import asm
from snippy_ng.cli.utils.bug_catcher import BugCatchingGroup


def show_citation(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(f"Please cite '{EXE}' in your research: …")
    ctx.exit()


def bug_report(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    import webbrowser
    url = f"{GITHUB_URL}/issues/new?template=bug_report.md&labels=bug&type=bug"
    click.echo(f"Please report bugs at: {url}")
    webbrowser.open(url, new=2)
    ctx.exit()


def version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(__version__)
    ctx.exit()


@click.group(
    cls=BugCatchingGroup,
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 120},
)
@click.option("--version", "-v", is_flag=True, callback=version, expose_value=False, help="Show version and exit.")
@click.option("--citation", "-c", is_flag=True, callback=show_citation, expose_value=False,
              help="Print citation for referencing Snippy-NG.")
@click.option("--bug", "-b", is_flag=True, callback=bug_report, expose_value=False,
              help="Report a bug or issue with Snippy-NG.")
def snippy_ng():
    """
    Snippy-NG: The Next Generation of Variant Calling.
    """
    pass

########################
# Register Subcommands #
########################
if os.getenv("SNIPPY_NG_OG_ENABLED", "0") == "1":
    snippy_ng.add_command(og)
snippy_ng.add_command(short)
snippy_ng.add_command(asm)