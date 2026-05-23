from pathlib import Path

import pytest

from snippy_ng.context import Context
from snippy_ng.exceptions import StageExecutionError
from snippy_ng.stages.filtering import BamReferenceValidator


def _write_fai(path: Path, records: dict[str, int]) -> None:
    offset = 0
    lines = []
    for name, length in records.items():
        lines.append(f"{name}\t{length}\t{offset}\t80\t81\n")
        offset += length + 1
    path.write_text("".join(lines), encoding="utf-8")


def _make_stage(
    tmp_path: Path,
    reference_index: Path,
    header: str,
    mapped_count: str = "1\n",
) -> BamReferenceValidator:
    stage = BamReferenceValidator(
        bam=tmp_path / "sample.bam",
        reference=tmp_path / "ref.fa",
        reference_index=reference_index,
        prefix=str(tmp_path / "snippy"),
    )
    stage.reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    stage.output.header.write_text(header, encoding="utf-8")
    stage.output.mapped_count.write_text(mapped_count, encoding="utf-8")
    return stage


def test_bam_reference_validator_commands_use_samtools(tmp_path):
    reference = tmp_path / "ref.fa"
    reference_index = tmp_path / "ref.fa.fai"
    stage = BamReferenceValidator(
        bam=tmp_path / "sample.bam",
        reference=reference,
        reference_index=reference_index,
        prefix=str(tmp_path / "snippy"),
    )

    commands = stage.create_commands(Context())

    assert commands[0].command == ["samtools", "view", "-H", str(stage.bam)]
    assert commands[0].output_file == stage.output.header
    assert commands[1].command == [
        "samtools",
        "view",
        "-c",
        "-F",
        "4",
        "-T",
        str(reference),
        str(stage.bam),
    ]
    assert commands[1].output_file == stage.output.mapped_count
    assert commands[2].func == stage.validate_bam_reference


def test_bam_reference_validator_accepts_matching_aligned_bam(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    reference_index = tmp_path / "ref.fa.fai"
    _write_fai(reference_index, {"chr1": 4, "plasmid": 2})
    stage = _make_stage(
        tmp_path,
        reference_index,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\n@SQ\tSN:plasmid\tLN:2\n",
    )

    stage.validate_bam_reference()


def test_bam_reference_validator_rejects_unaligned_bam_without_sq(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    reference_index = tmp_path / "ref.fa.fai"
    _write_fai(reference_index, {"chr1": 4})
    stage = _make_stage(tmp_path, reference_index, "@HD\tVN:1.6\n")

    with pytest.raises(StageExecutionError, match="no @SQ reference records"):
        stage.validate_bam_reference()


def test_bam_reference_validator_rejects_bam_with_no_mapped_reads(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    reference_index = tmp_path / "ref.fa.fai"
    _write_fai(reference_index, {"chr1": 4})
    stage = _make_stage(tmp_path, reference_index, "@SQ\tSN:chr1\tLN:4\n", "0\n")

    with pytest.raises(StageExecutionError, match="no mapped reads"):
        stage.validate_bam_reference()


def test_bam_reference_validator_rejects_missing_reference_contig(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n>plasmid\nAA\n", encoding="utf-8")
    reference_index = tmp_path / "ref.fa.fai"
    _write_fai(reference_index, {"chr1": 4, "plasmid": 2})
    stage = _make_stage(tmp_path, reference_index, "@SQ\tSN:chr1\tLN:4\n")

    with pytest.raises(
        StageExecutionError,
        match="reference contigs absent from BAM header: plasmid",
    ):
        stage.validate_bam_reference()


def test_bam_reference_validator_rejects_extra_bam_contig(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    reference_index = tmp_path / "ref.fa.fai"
    _write_fai(reference_index, {"chr1": 4})
    stage = _make_stage(
        tmp_path,
        reference_index,
        "@SQ\tSN:chr1\tLN:4\n@SQ\tSN:plasmid00002\tLN:2\n",
    )

    with pytest.raises(
        StageExecutionError,
        match="BAM contigs absent from reference: plasmid00002",
    ):
        stage.validate_bam_reference()


def test_bam_reference_validator_rejects_length_mismatch(tmp_path):
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    reference_index = tmp_path / "ref.fa.fai"
    _write_fai(reference_index, {"chr1": 4})
    stage = _make_stage(tmp_path, reference_index, "@SQ\tSN:chr1\tLN:5\n")

    with pytest.raises(StageExecutionError, match=r"chr1 \(BAM=5, reference=4\)"):
        stage.validate_bam_reference()
