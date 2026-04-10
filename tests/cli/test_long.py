import pytest
from tests.cli.helpers import apply_cli_case_overrides, assert_cli_result, get_bad_reference_target, run_cli_command, write_dummy_files


@pytest.fixture(autouse=True)
def stub_everything(stub_pipeline, stub_reference_format, stub_common_stages, stub_long_stages):
    """
    Combine all required stubs for long-read CLI tests.
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
                "--reads",     p["reads"],
                "--outdir",    p["out"],
                "--clair3-model", p["model"],
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
                "--clair3-model", p["model"],
                "--skip-check",
            ],
            0,
            True,
        ),
        (
            "check_only",
            lambda p: [
                "--reference", p["ref"],
                "--reads",     p["reads"],
                "--outdir",    p["out"],
                "--clair3-model", p["model"],
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
                "--reads",     p["reads"],
                "--outdir",    p["out"],
                "--clair3-model", p["model"],
                "--skip-check",
            ],
            2,
            False,
        ),
        (
            "bad_reference",
            lambda p: [
                "--reference", p["ref"],
                "--reads",     p["reads"],
                "--outdir",    p["out"],
                "--clair3-model", p["model"],
                "--skip-check",
            ],
            1,
            False,
        ),
        (
            "reads_no_model",
            lambda p: [
                "--reference", p["ref"],
                "--reads",     p["reads"],
                "--outdir",    p["out"],
                "--caller",    "clair3",
                "--skip-check",
            ],
            0,
            True,
        ),
        (
            "freebayes_caller",
            lambda p: [
                "--reference", p["ref"],
                "--reads",     p["reads"],
                "--outdir",    p["out"],
                "--caller",    "freebayes",
                "--skip-check",
            ],
            0,
            True,
        ),
        (
            "bam_no_model",
            lambda p: [
                "--reference", p["ref"],
                "--bam",       p["bam"],
                "--outdir",    p["out"],
                "--caller",    "clair3",
                "--skip-check",
            ],
            2,
            False,
        ),

    ],
)
def test_long_cli(monkeypatch, tmp_path, case_name, extra, expect_exit, expect_run):
    """
    Parameterised test for the `long` command.
    """

    # --------------- Arrange --------------------------------------------------
    paths = {
        "ref": tmp_path / "ref.fa",
        "reads": tmp_path / "long_reads.fq",
        "bam": tmp_path / "reads.bam",
        "out": tmp_path / "output",
        "model": tmp_path / "clair3_model",
    }
    write_dummy_files(paths, ["ref", "reads", "bam"])
    apply_cli_case_overrides(
        monkeypatch,
        case_name,
        paths["out"],
        bad_reference_target=get_bad_reference_target("long"),
    )

    args = ["long"] + extra(paths)
    result = run_cli_command(args)
    assert_cli_result(result, expect_exit, expect_run)
