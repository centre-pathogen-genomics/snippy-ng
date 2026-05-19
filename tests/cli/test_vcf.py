import subprocess
from pathlib import Path

from click.testing import CliRunner

from snippy_ng.cli import snippy_ng


def test_vcf_strand_bias_cli_annotates_output(monkeypatch, tmp_path):
    vcf = tmp_path / "input.vcf"
    bam = tmp_path / "sample.bam"
    reference = tmp_path / "reference.fa"
    output = tmp_path / "with.strand-bias.vcf"

    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tG\t60\tPASS\t.\n",
        encoding="utf-8",
    )
    bam.write_text("bam", encoding="utf-8")
    reference.write_text(
        ">chr1\nA\n",
        encoding="utf-8",
    )

    captured = {}

    def fake_run(command, check, stdout, stderr, text):
        captured["command"] = command
        assert check is True
        assert stdout == subprocess.PIPE
        assert stderr == subprocess.PIPE
        assert text is True
        bed_path = command[7]
        assert command[:7] == ["samtools", "mpileup", "-aa", "-B", "-f", str(reference.resolve()), "-l"]
        assert command[8] == str(bam.resolve())
        assert Path(bed_path).exists()
        return subprocess.CompletedProcess(
            command,
            0,
            stdout="chr1\t10\tA\t6\t..,,Gg\t~~~~~~\n",
            stderr="",
        )

    monkeypatch.setattr("snippy_ng.utils.strand_bias.subprocess.run", fake_run)

    result = CliRunner().invoke(
        snippy_ng,
        [
            "utils",
            "vcf",
            "strand-bias",
            "--alignment",
            str(bam),
            "--reference",
            str(reference),
            str(vcf),
            "-o",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["command"][:5] == ["samtools", "mpileup", "-aa", "-B", "-f"]
    output_text = output.read_text(encoding="utf-8")
    assert "SB_REF_FWD=2" in output_text
    assert "SB_ALT_REV=1" in output_text


def test_vcf_strand_bias_cli_accepts_bam_alias(monkeypatch, tmp_path):
    vcf = tmp_path / "input.vcf"
    bam = tmp_path / "sample.bam"
    reference = tmp_path / "reference.fa"
    output = tmp_path / "with.strand-bias.vcf"

    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tG\t60\tPASS\t.\n",
        encoding="utf-8",
    )
    bam.write_text("bam", encoding="utf-8")
    reference.write_text(
        ">chr1\nA\n",
        encoding="utf-8",
    )

    def fake_run(command, check, stdout, stderr, text):
        return subprocess.CompletedProcess(
            command,
            0,
            stdout="chr1\t10\tA\t6\t..,,Gg\t~~~~~~\n",
            stderr="",
        )

    monkeypatch.setattr("snippy_ng.utils.strand_bias.subprocess.run", fake_run)

    result = CliRunner().invoke(
        snippy_ng,
        [
            "utils",
            "vcf",
            "strand-bias",
            "--bam",
            str(bam),
            "--ref",
            str(reference),
            str(vcf),
            "-o",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    assert output.exists()


def test_vcf_strand_bias_cli_writes_to_stdout_by_default(monkeypatch, tmp_path):
    vcf = tmp_path / "input.vcf"
    bam = tmp_path / "sample.bam"
    reference = tmp_path / "reference.fa"

    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tG\t60\tPASS\t.\n",
        encoding="utf-8",
    )
    bam.write_text("bam", encoding="utf-8")
    reference.write_text(
        ">chr1\nA\n",
        encoding="utf-8",
    )

    def fake_run(command, check, stdout, stderr, text):
        return subprocess.CompletedProcess(
            command,
            0,
            stdout="chr1\t10\tA\t6\t..,,Gg\t~~~~~~\n",
            stderr="",
        )

    monkeypatch.setattr("snippy_ng.utils.strand_bias.subprocess.run", fake_run)

    result = CliRunner().invoke(
        snippy_ng,
        [
            "utils",
            "vcf",
            "strand-bias",
            "--alignment",
            str(bam),
            "--reference",
            str(reference),
            str(vcf),
        ],
    )

    assert result.exit_code == 0, result.output
    assert result.output.startswith("##fileformat=VCFv4.2\n")
    assert "SB_REF_FWD=2" in result.output
    assert "SB_ALT_REV=1" in result.output


def test_vcf_strand_bias_cli_requires_reference(tmp_path):
    vcf = tmp_path / "input.vcf"
    bam = tmp_path / "sample.bam"

    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tG\t60\tPASS\t.\n",
        encoding="utf-8",
    )
    bam.write_text("bam", encoding="utf-8")

    result = CliRunner().invoke(
        snippy_ng,
        [
            "utils",
            "vcf",
            "strand-bias",
            "--alignment",
            str(bam),
            str(vcf),
        ],
    )

    assert result.exit_code != 0
    assert "--reference" in result.output or "--ref" in result.output
