import types

import pytest
from click.testing import CliRunner

from snippy_ng.cli import snippy_ng         # the click *group*
import snippy_ng.pipeline as _pl           # <-- real module we patch


##############################################################################
#                           1.  DUMMY IMPLEMENTATIONS                         #
##############################################################################

@pytest.fixture(autouse=True)
def stub_everything(monkeypatch, tmp_path):
    """
    Replace the heavy pipeline and stage classes with tiny stand‑ins.
    The fixture runs automatically for every test.
    """

    # ---------- Dummy Pipeline ------------------------------------------------
    class DummyPipeline:
        """Light-weight stand-in that records what happened."""
        last = None                        # will hold most‑recent instance

        def __init__(self, *_, **__):
            self.validated = False
            self.ran       = False
            DummyPipeline.last = self      # remember myself

        def welcome(self):                 pass
        def validate_dependencies(self):   self.validated = True
        def set_working_directory(self, *_): pass
        def run(self, quiet=False, continue_last_run=False, keep_incomplete=False):        self.ran = True
        def cleanup(self):                 pass
        def goodbye(self):                 pass
        def error(self, *_):               pass

    # patch the real module, not the click group
    monkeypatch.setattr(_pl, "Pipeline", DummyPipeline)

    # ---------- Dummy Stage objects ------------------------------------------
    def _stage_factory(output):
        class _Stage:
            def __init__(self, *_, **__):
                self.output = types.SimpleNamespace(**{out_key: out_val for out_key, out_val in output.items()})
        return _Stage

    monkeypatch.setattr(
        "snippy_ng.stages.setup.PrepareReference",
        _stage_factory({"reference": tmp_path / "ref.fa", "gff": tmp_path / "ref.gff", "reference_index": tmp_path / "ref.fa.fai"}),
    )
    monkeypatch.setattr(
        "snippy_ng.stages.alignment.BWAMEMReadsAligner",
        _stage_factory({"bam": tmp_path / "align.bam"}),
    )
    monkeypatch.setattr(
        "snippy_ng.stages.alignment.PreAlignedReads",
        _stage_factory({"bam": tmp_path / "align.bam"}),
    )
    monkeypatch.setattr(
        "snippy_ng.stages.calling.FreebayesCaller",
        _stage_factory({"vcf": tmp_path / "calls.vcf"}),
    )

    # ---------- Always recognise the reference format ------------------------
    monkeypatch.setattr("snippy_ng.seq_utils.guess_format", lambda _: "fasta")


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
                "--skip-check",
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
            1,
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
def test_run_cli(monkeypatch, tmp_path, case_name, extra, expect_exit, expect_run):
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
    for f in ["ref", "r1", "r2", "bam"]:
        paths[f].write_text(">dummy\nA")

    if case_name == "outdir_exists":
        paths["out"].mkdir()

    if case_name == "bad_reference":
        monkeypatch.setattr("snippy_ng.cli.utils.reference.guess_format", lambda _: None)

    args = ["short"] + extra(paths)
    runner = CliRunner()

    # --------------- Act ------------------------------------------------------
    result = runner.invoke(snippy_ng, args)

    # --------------- Assert ---------------------------------------------------
    assert result.exit_code == expect_exit, result.output

    # Did we create / run a pipeline?
    last_pipeline = _pl.Pipeline.last        # may be None if creation failed early

    if expect_run:
        assert last_pipeline and last_pipeline.ran is True
    else:
        if last_pipeline:
            assert last_pipeline.ran is False
