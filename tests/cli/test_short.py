import pytest
from tests.cli.helpers import apply_cli_case_overrides, assert_cli_result, get_bad_reference_target, run_cli_command, write_dummy_files


@pytest.fixture(autouse=True)
def stub_everything(stub_pipeline, stub_reference_format, stub_common_stages, stub_short_stages):
    """
    Combine all required stubs for short-read CLI tests.
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
            "reads_ok",
            lambda p: [
                "--reference", p["ref"],
                "--R1",        p["r1"],
                "--R2",        p["r2"],
                "--outdir",    p["out"],
                "--skip-check",
            ],
            0,
            True,
        ),
        (
            "bam_ok",
            lambda p: [
                "--reference", p["ref"],
                "--bam",       p["bam"],
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
                "--R1",        p["r1"],
                "--R2",        p["r2"],
                "--outdir",    p["out"],
                "--check",
            ],
            0,
            False,
        ),
        (
            "outdir_exists",
            lambda p: [
                "--reference", p["ref"],
                "--R1",        p["r1"],
                "--R2",        p["r2"],
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
                "--R1",        p["r1"],
                "--R2",        p["r2"],
                "--outdir",    p["out"],
                "--skip-check",
            ],
            1,
            False,
        ),
    ],
)
def test_short_cli(monkeypatch, tmp_path, case_name, extra, expect_exit, expect_run):
    """
    Parameterised test for the `run` command.
    """

    # --------------- Arrange --------------------------------------------------
    paths = {
        "ref": tmp_path / "ref.fa",
        "r1":  tmp_path / "reads_1.fq",
        "r2":  tmp_path / "reads_2.fq",
        "bam": tmp_path / "reads.bam",
        "out": tmp_path / "output",
    }
    write_dummy_files(paths, ["ref", "r1", "r2", "bam"])
    apply_cli_case_overrides(
        monkeypatch,
        case_name,
        paths["out"],
        bad_reference_target=get_bad_reference_target("short"),
    )

    args = ["short"] + extra(paths)
    result = run_cli_command(args)
    assert_cli_result(result, expect_exit, expect_run)
