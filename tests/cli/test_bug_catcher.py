import platform

from urllib.parse import quote

import click
from click.testing import CliRunner


from snippy_ng.cli.utils.bug_catcher import BugCatchingGroup
from snippy_ng.__about__ import GITHUB_URL


def test_exception_triggers_bug_report(capsys):
    """
    Ensure that when a subcommand raises an exception, BugCatchingGroup.main()
    prints the bug-report template to stderr and exits with code 1.
    """

    # Define a dummy CLI group that uses BugCatchingGroup as its base
    @click.group(cls=BugCatchingGroup)
    def cli():
        pass

    # Add a subcommand that simply raises an exception
    @cli.command()
    @click.argument("x", type=int)
    def explode(x):
        raise ValueError("test error")

    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(cli, ["explode", "42"])

    assert result.exit_code == 1

    stderr = result.stderr

    assert "ValueError: test error" in stderr

    assert GITHUB_URL in stderr
