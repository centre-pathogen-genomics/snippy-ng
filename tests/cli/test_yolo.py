from pathlib import Path
from types import SimpleNamespace

from click.testing import CliRunner

from snippy_ng.cli import snippy_ng
from snippy_ng.stages.alignment_filter import FilterAlignmentByAlignedPercentage
from snippy_ng.stages.core import SoftCoreFilter
from snippy_ng.stages.trees import ScaleTreeToSNPs
from tests.cli.helpers import make_prepared_reference, stub_load_or_prepare_reference


def _stub_common_yolo_dependencies(monkeypatch, tmp_path, samples):
    (tmp_path / "reference.fasta").write_text(">ref\nACGT\n")

    monkeypatch.setattr(
        "snippy_ng.utils.gather.gather",
        lambda **_: {"samples": samples},
    )

    _, ref_file = make_prepared_reference(tmp_path)
    stub_load_or_prepare_reference(monkeypatch, ref_file)
    captured = {}

    class DummySnippyPipeline:
        def __init__(self, stages=None, outputs_to_keep=None):
            self.stages = stages or []

        def get_stage(self, stage_class):
            for stage in reversed(self.stages):
                if isinstance(stage, stage_class):
                    return stage
            return None

        def run(self, _ctx):
            captured.setdefault("run_contexts", []).append(_ctx.model_copy(deep=True))
            return 0

    monkeypatch.setattr("snippy_ng.pipelines.SnippyPipeline", DummySnippyPipeline)

    class DummyCorePipeline:
        def __init__(self):
            self.stages = [
                FilterAlignmentByAlignedPercentage(
                    aln=Path("core.full.aln"),
                    alignment_stats=Path("core.aligned.tsv"),
                    prefix="core",
                ),
                SoftCoreFilter(
                    aln=Path("core.full.aln"),
                    core_threshold=0.95,
                    prefix="core",
                )
            ]

        def get_stage(self, stage_class):
            for stage in reversed(self.stages):
                if isinstance(stage, stage_class):
                    return stage
            return None

        def run(self, ctx):
            outdir = Path(ctx.outdir)
            outdir.mkdir(parents=True, exist_ok=True)
            (outdir / "core.filtered.aln").write_text(">s1\nACGT\n")
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
    return captured


def _run_yolo(tmp_path, *extra_args):
    outdir = tmp_path / "out"
    result = CliRunner().invoke(
        snippy_ng,
        ["yolo", str(tmp_path), "-o", str(outdir), *extra_args, "--skip-check"],
    )
    return outdir, result


def test_yolo_uses_sample_filtered_full_alignment_for_clonalframe(monkeypatch, tmp_path):
    # YOLO requires >=3 samples to proceed to tree stage
    samples = {
        "s1": {"type": "short"},
        "s2": {"type": "short"},
        "s3": {"type": "short"},
    }
    captured_pipeline = _stub_common_yolo_dependencies(monkeypatch, tmp_path, samples)

    monkeypatch.setattr("snippy_ng.pipelines.multi.run_multi_pipeline", lambda **_: (list(samples.keys()), []))

    captured = {}

    class DummyTreePipelineBuilder:
        def __init__(self, aln, fast_mode, clonalframe):
            captured["aln"] = Path(aln)
            captured["fast_mode"] = fast_mode
            captured["clonalframe"] = clonalframe

        def build(self):
            stage = ScaleTreeToSNPs(tree=Path("tree.treefile"), aln=Path("core.filtered.aln"), prefix="tree")
            return SimpleNamespace(
                run=lambda _ctx: 0,
                stages=[stage],
                get_stage=lambda stage_class: stage if isinstance(stage, stage_class) else None,
            )

    monkeypatch.setattr(
        "snippy_ng.pipelines.tree.TreePipelineBuilder",
        DummyTreePipelineBuilder,
    )

    # Act
    outdir, result = _run_yolo(tmp_path)

    # Assert
    assert result.exit_code == 0, result.output
    assert captured["aln"] == outdir / "core" / "core.filtered.aln"
    assert captured["fast_mode"] is True
    assert captured["clonalframe"] is True
    assert captured_pipeline["run_contexts"][0].log_path == (outdir / "reference" / "LOG.txt").absolute()
    assert captured_pipeline["run_contexts"][0].outdir == outdir / "reference"


def test_yolo_skips_tree_when_less_than_three_samples(monkeypatch, tmp_path):
    samples = {
        "s1": {"type": "short"},
        "s2": {"type": "short"},
    }
    captured_pipeline = _stub_common_yolo_dependencies(monkeypatch, tmp_path, samples)
    monkeypatch.setattr("snippy_ng.pipelines.multi.run_multi_pipeline", lambda **_: (list(samples.keys()), []))

    class ShouldNotBeCalledTreePipelineBuilder:
        def __init__(self, *args, **kwargs):
            raise AssertionError("TreePipelineBuilder should not be called when sample count < 3")

    monkeypatch.setattr(
        "snippy_ng.pipelines.tree.TreePipelineBuilder",
        ShouldNotBeCalledTreePipelineBuilder,
    )

    _, result = _run_yolo(tmp_path)

    assert result.exit_code == 0, result.output
    assert captured_pipeline["run_contexts"][0].log_path == ((tmp_path / "out") / "reference" / "LOG.txt").absolute()


def test_yolo_accepts_reference_accession(monkeypatch, tmp_path):
    samples = {
        "s1": {"type": "short"},
        "s2": {"type": "short"},
    }
    _stub_common_yolo_dependencies(monkeypatch, tmp_path, samples)

    captured_gather = {}

    def fake_gather(**kwargs):
        captured_gather.update(kwargs)
        return {"samples": samples}

    monkeypatch.setattr("snippy_ng.utils.gather.gather", fake_gather)

    captured_download = {}

    def fake_download_assembly(reference_accession, stages, output_directory):
        captured_download["reference_accession"] = reference_accession
        captured_download["output_directory"] = output_directory
        downloaded_reference = Path(output_directory) / f"{reference_accession}.fa"
        stages.append(SimpleNamespace(output=SimpleNamespace(fasta=downloaded_reference)))
        return downloaded_reference

    monkeypatch.setattr("snippy_ng.pipelines.common.download_assembly", fake_download_assembly)
    monkeypatch.setattr("snippy_ng.pipelines.multi.run_multi_pipeline", lambda **_: (list(samples.keys()), []))

    class ShouldNotBeCalledTreePipelineBuilder:
        def __init__(self, *args, **kwargs):
            raise AssertionError("TreePipelineBuilder should not be called when sample count < 3")

    monkeypatch.setattr(
        "snippy_ng.pipelines.tree.TreePipelineBuilder",
        ShouldNotBeCalledTreePipelineBuilder,
    )

    outdir, result = _run_yolo(tmp_path, "--reference", "SAMN123456")

    assert result.exit_code == 0, result.output
    assert captured_gather["reference"] is None
    assert captured_download == {
        "reference_accession": "SAMN123456",
        "output_directory": outdir / "reference",
    }


def test_yolo_preserves_long_sample_caller_and_delegates_auto_cpu_allocation(monkeypatch, tmp_path):
    samples = {
        "short_1": {"type": "short"},
        "long_1": {"type": "long", "caller": "clair3", "clair3_model": None},
        "short_2": {"type": "short"},
    }
    _stub_common_yolo_dependencies(monkeypatch, tmp_path, samples)

    captured_multi = {}

    def fake_run_multi_pipeline(**kwargs):
        captured_multi.update(kwargs)
        return list(samples.keys()), []

    monkeypatch.setattr("snippy_ng.pipelines.multi.run_multi_pipeline", fake_run_multi_pipeline)

    class DummyTreePipelineBuilder:
        def __init__(self, aln, fast_mode, clonalframe):
            pass

        def build(self):
            stage = ScaleTreeToSNPs(tree=Path("tree.treefile"), aln=Path("core.filtered.aln"), prefix="tree")
            return SimpleNamespace(
                run=lambda _ctx: 0,
                stages=[stage],
                get_stage=lambda stage_class: stage if isinstance(stage, stage_class) else None,
            )

    monkeypatch.setattr(
        "snippy_ng.pipelines.tree.TreePipelineBuilder",
        DummyTreePipelineBuilder,
    )

    _, result = _run_yolo(tmp_path, "-c", "8")

    assert result.exit_code == 0, result.output
    assert captured_multi["samples"]["long_1"]["caller"] == "clair3"
    assert captured_multi["cpus_per_sample"] is None
    assert captured_multi["run_ctx"].cpus == 8
    assert captured_multi["snippy_reference_dir"] == tmp_path / "prepared_reference"


def test_yolo_errors_when_no_reference_found(tmp_path):
    runner = CliRunner()
    result = runner.invoke(
        snippy_ng,
        ["yolo", str(tmp_path), "-o", str(tmp_path / "out"), "--skip-check"],
    )

    assert result.exit_code != 0
    assert "No reference file found" in result.output
