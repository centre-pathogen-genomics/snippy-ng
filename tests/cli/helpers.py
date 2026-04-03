import json
from pathlib import Path
from typing import Iterable
from types import SimpleNamespace

from click.testing import CliRunner

from snippy_ng.cli import snippy_ng
import snippy_ng.pipelines as _pl


BAD_REFERENCE_TARGETS = {
    "short": "snippy_ng.pipelines.short.load_or_prepare_reference",
    "long": "snippy_ng.pipelines.long.load_or_prepare_reference",
    "asm": "snippy_ng.pipelines.asm.load_or_prepare_reference",
    "multi": "snippy_ng.pipelines.common.load_or_prepare_reference",
}


def write_dummy_files(paths: dict[str, Path], keys: Iterable[str], contents: str = ">dummy\nA") -> None:
    for key in keys:
        paths[key].write_text(contents)


def apply_cli_case_overrides(monkeypatch, case_name: str, outdir: Path, bad_reference_target: str | None = None) -> None:
    if case_name == "outdir_exists":
        outdir.mkdir()

    if case_name == "bad_reference":
        monkeypatch.setattr("snippy_ng.pipelines.common.guess_reference_format", lambda _: None)
        if bad_reference_target is not None:
            monkeypatch.setattr(
                bad_reference_target,
                lambda *args, **kwargs: (_ for _ in ()).throw(ValueError("Could not determine reference format")),
            )


def get_bad_reference_target(command_name: str) -> str:
    return BAD_REFERENCE_TARGETS[command_name]


def run_cli_command(args: list[str]):
    return CliRunner().invoke(snippy_ng, args)


def assert_cli_result(result, expect_exit: int, expect_run: bool) -> None:
    assert result.exit_code == expect_exit, result.output

    last_pipeline = _pl.SnippyPipeline.last
    if expect_run:
        assert last_pipeline and last_pipeline.ran is True
    elif last_pipeline:
        assert last_pipeline.ran is False


def make_prepared_reference(tmp_path: Path, dirname: str = "prepared_reference") -> tuple[Path, Path]:
    ref_dir = tmp_path / dirname
    ref_dir.mkdir(parents=True, exist_ok=True)
    ref_file = ref_dir / "genomic.fa"
    ref_file.write_text(">ref\nACGT\n")
    (ref_dir / "genomic.gff").write_text("##gff-version 3\n")
    (ref_dir / "genomic.fa.fai").write_text("ref\t4\t5\t4\t5\n")
    (ref_dir / "genomic.dict").write_text("@SQ\tSN:ref\tLN:4\n")
    (ref_dir / "metadata.json").write_text(json.dumps({
        "reference": "genomic.fa",
        "format": "fasta",
        "num_sequences": 1,
        "total_length": 4,
        "num_features": 0,
        "prefix": "genomic",
        "datetime": "2026-01-01T00:00:00",
        "version": "test",
    }))
    return ref_dir, ref_file


def stub_load_or_prepare_reference(monkeypatch, ref_file: Path, target: str = "snippy_ng.pipelines.common.load_or_prepare_reference") -> None:
    monkeypatch.setattr(
        target,
        lambda *args, **kwargs: SimpleNamespace(output=SimpleNamespace(reference=ref_file)),
    )
