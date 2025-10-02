import click
from snippy_ng.__about__ import __version__, EXE, GITHUB_URL 

class BugCatchingGroup(click.Group):
    """
    A Click group that catches exceptions and provides a bug-report URL
    plus a pre-filled bug report template.
    """

    def main(self, *args, **kwargs):
        try:
            # Pass standalone_mode=False so that Click will not swallow exceptions,
            # and will let them bubble up here.
            return super().main(*args, **kwargs)

        except Exception as e:
            import sys
            import traceback
            import platform
            from urllib.parse import quote
            from snippy_ng.logging import logger
            
            traceback.print_exc(file=sys.stderr)

            url = (
                f"{GITHUB_URL}/issues/new"
                f"?template=bug_report.md"
                f"&labels=cli,bug&type=bug"
                f"&title={quote(str(e))}"
            )
            
            try:
                cmd_string = " ".join([EXE] + sys.argv[1:])
            except Exception:
                cmd_string = "<unable to retrieve sys.argv>"

            snippy_version = __version__ if __version__ else "<unknown>"

            try:
                os_info = platform.platform()
            except Exception:
                os_info = "<unable to retrieve platform info>"

            try:
                full_tb = traceback.format_exc()
                if not full_tb.strip():
                    # format_exc() returns an empty string if called a second time immediately.
                    # In that rare case, we can re-capture from the exception object:
                    full_tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
            except Exception:
                full_tb = "<unable to capture traceback>"

            logger.echo(err=True)
            logger.horizontal_rule("Bug report template (copy/paste the sections below into GitHub)", color='red')

            logger.echo("\n**Describe the bug**", err=True)
            logger.echo("> A clear and concise description of what the bug is.\n", err=True)

            logger.echo("You encountered an exception. Here is the exception message:\n", err=True)
            logger.echo(f"```\n{type(e).__name__}: {e}\n```\n", err=True)

            logger.echo("The command you ran was:\n", err=True)
            logger.echo(f"```\n{cmd_string}\n```\n", err=True)

            logger.echo("**Environment**", err=True)
            logger.echo("> Please provide your OS and Snippy-NG version. If you are running in conda, "
                       "you can add 'conda list' output here as well.\n", err=True)
            logger.echo(f"- OS: `{os_info}`", err=True)
            logger.echo(f"- Snippy-NG version: `{snippy_version}`\n", err=True)
            
            logger.echo("**Additional context**", err=True)
            logger.echo("> Add any other context about the problem here. (e.g., input files, "
                       "exact command arguments, steps to reproduce, etc.)\n", err=True)

            logger.echo("**Backtrace**", err=True)
            logger.echo("> If possible, please include the complete error log below.\n", err=True)
            logger.echo("```", err=True)
            logger.echo(full_tb.rstrip("\n"), err=True)
            logger.echo("```\n", err=True)

            logger.horizontal_rule()
            logger.echo(
                "\nOh no! You broke Snippy-NG... Congrats! Please use the following URL to report this bug:",
                err=True,
            )
            logger.echo(f"\n{url}\n", err=True)
            logger.echo("Above is a pre-filled bug report template. "
                       "Please copy/paste it into the GitHub issue form.\n", err=True)
            
            sys.exit(1)
