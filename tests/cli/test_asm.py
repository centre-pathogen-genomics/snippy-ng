import pytest
from click.testing import CliRunner

from snippy_ng.cli import snippy_ng
import snippy_ng.pipelines.pipeline_runner as _pl


@pytest.fixture(autouse=True)
def stub_everything(stub_pipeline, stub_reference_format, stub_common_stages, stub_asm_stages):
    """
    Combine all required stubs for assembly CLI tests.
    This fixture runs automatically for every test.
    """
    pass


##############################################################################
#                               2.  TEST DATA                                 #
##############################################################################

@pytest.mark.parametrize(
    "case_name, extra, expect_exit, expect_run",
    [
        (
            "assembly_ok",
            lambda p: [
                "--reference", p["ref"],
                "--assembly",  p["asm"],
                "--outdir",    p["out"],
                "--skip-check",
            ],
            0,
            True,
        ),
        (
            "check_only",
            lambda p: [
                "--reference", p["ref"],
                "--assembly",  p["asm"],
                "--outdir",    p["out"],
                "--check",
                "--skip-check",
            ],
            0,
            False,
        ),
        (
            "outdir_exists",
            lambda p: [
                "--reference", p["ref"],
                "--assembly",  p["asm"],
                "--outdir",    p["out"],
                "--skip-check",
            ],
            2,
            False,
        ),
        (
            "bad_reference",
            lambda p: [
                "--reference", p["ref"],
                "--assembly",  p["asm"],
                "--outdir",    p["out"],
                "--skip-check",
            ],
            1,
            False,
        ),
    ],
)
def test_asm_cli(monkeypatch, tmp_path, case_name, extra, expect_exit, expect_run):
    """
    Parameterised test for the `asm` command.
    """

    # --------------- Arrange --------------------------------------------------
    paths = {
        "ref": tmp_path / "ref.fa",
        "asm": tmp_path / "assembly.fa",
        "out": tmp_path / "output",
    }
    for f in ["ref", "asm"]:
        paths[f].write_text(">dummy\nA")

    if case_name == "outdir_exists":
        paths["out"].mkdir()

    if case_name == "bad_reference":
        monkeypatch.setattr("snippy_ng.pipelines.common.guess_reference_format", lambda _: None)

    args = ["asm"] + extra(paths)
    runner = CliRunner()

    # --------------- Act ------------------------------------------------------
    result = runner.invoke(snippy_ng, args)

    # --------------- Assert ---------------------------------------------------
    assert result.exit_code == expect_exit, result.output

    # Did we create / run a pipeline?
    last_pipeline = _pl.Snippy.last        # may be None if creation failed early

    if expect_run:
        assert last_pipeline and last_pipeline.ran is True
    else:
        if last_pipeline:
            assert last_pipeline.ran is False
