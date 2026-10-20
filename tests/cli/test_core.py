from click.testing import CliRunner

from snippy_ng.cli import snippy_ng


def test_core_cli_accepts_ref_alias_and_passes_reference(monkeypatch, tmp_path):
    ref_file = tmp_path / "ref.fa"
    ref_file.write_text(">ref\nATCG\n")

    snippy_dir = tmp_path / "sample1"
    snippy_dir.mkdir()

    captured = {}

    class DummyCorePipeline:
        def run(self, _ctx):
            return 0

    class DummyCorePipelineBuilder:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def build(self):
            return DummyCorePipeline()

    monkeypatch.setattr(
        "snippy_ng.pipelines.core.CorePipelineBuilder", DummyCorePipelineBuilder
    )

    result = CliRunner().invoke(
        snippy_ng,
        [
            "core",
            "--ref",
            str(ref_file),
            str(snippy_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["reference"] == ref_file
    assert captured["snippy_dirs"] == (snippy_dir,)
    assert captured["mask"] is None


def test_core_cli_passes_mask(monkeypatch, tmp_path):
    ref_file = tmp_path / "ref.fa"
    ref_file.write_text(">ref\nATCG\n")
    mask_file = tmp_path / "mask.bed"
    mask_file.write_text("ref\t0\t2\n")
    snippy_dir = tmp_path / "sample1"
    snippy_dir.mkdir()

    captured = {}

    class DummyCorePipelineBuilder:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def build(self):
            class DummyCorePipeline:
                def run(self, _ctx):
                    return 0

            return DummyCorePipeline()

    monkeypatch.setattr(
        "snippy_ng.pipelines.core.CorePipelineBuilder", DummyCorePipelineBuilder
    )

    result = CliRunner().invoke(
        snippy_ng,
        ["core", "--ref", str(ref_file), "--mask", str(mask_file), str(snippy_dir)],
    )

    assert result.exit_code == 0, result.output
    assert captured["mask"] == mask_file


def _capture_core_builder(monkeypatch):
    captured = {}

    class DummyCorePipelineBuilder:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def build(self):
            class DummyCorePipeline:
                def run(self, _ctx):
                    return 0

            return DummyCorePipeline()

    monkeypatch.setattr(
        "snippy_ng.pipelines.core.CorePipelineBuilder", DummyCorePipelineBuilder
    )
    return captured


def test_core_cli_keeps_alignment_outliers_by_default(monkeypatch, tmp_path):
    ref_file = tmp_path / "ref.fa"
    ref_file.write_text(">ref\nATCG\n")
    snippy_dir = tmp_path / "sample1"
    snippy_dir.mkdir()
    captured = _capture_core_builder(monkeypatch)

    result = CliRunner().invoke(
        snippy_ng, ["core", "--ref", str(ref_file), str(snippy_dir)]
    )

    assert result.exit_code == 0, result.output
    assert captured["remove_outliers"] is False
    assert captured["outlier_threshold"] == 0.20


def test_core_cli_passes_alignment_outlier_options(monkeypatch, tmp_path):
    ref_file = tmp_path / "ref.fa"
    ref_file.write_text(">ref\nATCG\n")
    snippy_dir = tmp_path / "sample1"
    snippy_dir.mkdir()
    captured = _capture_core_builder(monkeypatch)

    result = CliRunner().invoke(
        snippy_ng,
        [
            "core",
            "--ref",
            str(ref_file),
            "--remove-alignment-outliers",
            "--alignment-outlier-threshold",
            "0.5",
            str(snippy_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["remove_outliers"] is True
    assert captured["outlier_threshold"] == 0.5
