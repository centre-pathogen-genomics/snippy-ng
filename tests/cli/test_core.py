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

    monkeypatch.setattr("snippy_ng.pipelines.core.CorePipelineBuilder", DummyCorePipelineBuilder)

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