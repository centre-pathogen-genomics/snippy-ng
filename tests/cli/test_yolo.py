from pathlib import Path
from types import SimpleNamespace

from click.testing import CliRunner

from snippy_ng.cli import snippy_ng


def _stub_common_yolo_dependencies(monkeypatch, tmp_path, samples):
    (tmp_path / "reference.fasta").write_text(">ref\nACGT\n")

    monkeypatch.setattr(
        "snippy_ng.utils.gather.gather_samples_config",
        lambda **_: {"samples": samples},
    )

    ref_dir = tmp_path / "prepared_reference"
    ref_dir.mkdir(parents=True, exist_ok=True)
    ref_file = ref_dir / "genomic.fa"
    ref_file.write_text(">ref\nACGT\n")

    def fake_load_or_prepare_reference(reference_path, output_directory):
        return SimpleNamespace(output=SimpleNamespace(reference=ref_file))

    monkeypatch.setattr(
        "snippy_ng.pipelines.common.load_or_prepare_reference",
        fake_load_or_prepare_reference,
    )

    class DummySnippyPipeline:
        def __init__(self, stages=None, outputs_to_keep=None):
            self.stages = stages or []

        def run(self, _ctx):
            return 0

    monkeypatch.setattr("snippy_ng.pipelines.SnippyPipeline", DummySnippyPipeline)

    class DummyCorePipeline:
        def __init__(self):
            self.stages = [
                SimpleNamespace(
                    output=SimpleNamespace(
                        soft_core=Path("core.095.aln"),
                        constant_sites=Path("core.095.fconst"),
                    )
                )
            ]

        def run(self, ctx):
            outdir = Path(ctx.outdir)
            outdir.mkdir(parents=True, exist_ok=True)
            (outdir / "core.095.aln").write_text(">s1\nACGT\n")
            (outdir / "core.095.fconst").write_text("1,2,3,4\n")
            return 0

    class DummyCorePipelineBuilder:
        def __init__(self, **kwargs):
            pass

        def build(self):
            return DummyCorePipeline()

    monkeypatch.setattr(
        "snippy_ng.pipelines.core.CorePipelineBuilder",
        DummyCorePipelineBuilder,
    )


def test_yolo_uses_soft_core_output_for_tree(monkeypatch, tmp_path):
    # YOLO requires >=3 samples to proceed to tree stage
    samples = {
        "s1": {"type": "short"},
        "s2": {"type": "short"},
        "s3": {"type": "short"},
    }
    _stub_common_yolo_dependencies(monkeypatch, tmp_path, samples)

    monkeypatch.setattr("snippy_ng.pipelines.multi.run_multi_pipeline", lambda **_: (list(samples.keys()), []))

    captured = {}

    class DummyTreePipelineBuilder:
        def __init__(self, aln, fconst, fast_mode):
            captured["aln"] = Path(aln)
            captured["fconst"] = fconst
            captured["fast_mode"] = fast_mode

        def build(self):
            return SimpleNamespace(run=lambda _ctx: 0, stages=[SimpleNamespace(output=SimpleNamespace(tree=Path("tree.newick")))])

    monkeypatch.setattr(
        "snippy_ng.pipelines.tree.TreePipelineBuilder",
        DummyTreePipelineBuilder,
    )

    # Act
    outdir = tmp_path / "out"
    runner = CliRunner()
    result = runner.invoke(
        snippy_ng,
        ["yolo", str(tmp_path), "-o", str(outdir), "--skip-check"],
    )

    # Assert
    assert result.exit_code == 0, result.output
    assert captured["aln"] == outdir / "core" / "core.095.aln"
    assert captured["fconst"] == "1,2,3,4"
    assert captured["fast_mode"] is False


def test_yolo_skips_tree_when_less_than_three_samples(monkeypatch, tmp_path):
    samples = {
        "s1": {"type": "short"},
        "s2": {"type": "short"},
    }
    _stub_common_yolo_dependencies(monkeypatch, tmp_path, samples)
    monkeypatch.setattr("snippy_ng.pipelines.multi.run_multi_pipeline", lambda **_: (list(samples.keys()), []))

    class ShouldNotBeCalledTreePipelineBuilder:
        def __init__(self, *args, **kwargs):
            raise AssertionError("TreePipelineBuilder should not be called when sample count < 3")

    monkeypatch.setattr(
        "snippy_ng.pipelines.tree.TreePipelineBuilder",
        ShouldNotBeCalledTreePipelineBuilder,
    )

    outdir = tmp_path / "out"
    runner = CliRunner()
    result = runner.invoke(
        snippy_ng,
        ["yolo", str(tmp_path), "-o", str(outdir), "--skip-check"],
    )

    assert result.exit_code == 0, result.output


def test_yolo_sets_long_samples_to_freebayes_and_cpus_per_sample(monkeypatch, tmp_path):
    samples = {
        "short_1": {"type": "short"},
        "long_1": {"type": "long"},
        "short_2": {"type": "short"},
    }
    _stub_common_yolo_dependencies(monkeypatch, tmp_path, samples)

    captured_multi = {}

    def fake_run_multi_pipeline(**kwargs):
        captured_multi.update(kwargs)
        return list(samples.keys()), []

    monkeypatch.setattr("snippy_ng.pipelines.multi.run_multi_pipeline", fake_run_multi_pipeline)

    class DummyTreePipelineBuilder:
        def __init__(self, aln, fconst, fast_mode):
            pass

        def build(self):
            return SimpleNamespace(run=lambda _ctx: 0, stages=[SimpleNamespace(output=SimpleNamespace(tree=Path("tree.newick")))])

    monkeypatch.setattr(
        "snippy_ng.pipelines.tree.TreePipelineBuilder",
        DummyTreePipelineBuilder,
    )

    outdir = tmp_path / "out"
    runner = CliRunner()
    result = runner.invoke(
        snippy_ng,
        ["yolo", str(tmp_path), "-o", str(outdir), "-c", "8", "--skip-check"],
    )

    assert result.exit_code == 0, result.output
    assert captured_multi["samples"]["long_1"]["caller"] == "freebayes"
    assert captured_multi["cpus_per_sample"] == 4


def test_yolo_errors_when_no_reference_found(tmp_path):
    runner = CliRunner()
    result = runner.invoke(
        snippy_ng,
        ["yolo", str(tmp_path), "-o", str(tmp_path / "out"), "--skip-check"],
    )

    assert result.exit_code != 0
    assert "No reference file found" in result.output
