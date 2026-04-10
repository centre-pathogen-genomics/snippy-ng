import click
from click.testing import CliRunner
from pathlib import Path


from snippy_ng.cli.utils.bug_catcher import BugCatchingGroup
from snippy_ng.__about__ import GITHUB_URL
from snippy_ng.logging import logger


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


def test_exception_bug_report_is_written_to_log_file(tmp_path):
    @click.group(cls=BugCatchingGroup)
    @click.option("--outdir", type=click.Path(path_type=Path), default=tmp_path)
    def cli(outdir):
        pass

    @cli.command()
    def explode():
        logger.set_log_path(tmp_path / "LOG.txt")
        logger.set_log_path(None)
        raise AssertionError("BOOM!")

    runner = CliRunner(mix_stderr=False)
    result = runner.invoke(cli, ["--outdir", str(tmp_path), "explode"])

    assert result.exit_code == 1

    log_contents = (tmp_path / "LOG.txt").read_text()
    assert "AssertionError: BOOM!" in log_contents
    assert "Traceback (most recent call last):" in log_contents
