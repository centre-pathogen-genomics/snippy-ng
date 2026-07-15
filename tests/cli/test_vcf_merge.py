from pathlib import Path
from subprocess import CompletedProcess

from click.testing import CliRunner

from snippy_ng.cli.vcf.merge_cli import merge


def _write_vcf(path: Path, records: list[str]) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "\n".join(records) + "\n",
        encoding="utf-8",
    )


def test_merge_cli_runs_bcftools_and_writes_upset_plot(tmp_path, monkeypatch):
    first, second = tmp_path / "first.vcf", tmp_path / "second.vcf"
    output, plot = tmp_path / "merged.vcf", tmp_path / "overlap.svg"
    _write_vcf(first, ["chr1\t10\t.\tA\tT\t60\tPASS\t.", "chr1\t20\t.\tC\tG\t60\tPASS\t."])
    _write_vcf(second, ["chr1\t10\t.\tA\tT\t60\tPASS\t.", "chr1\t30\t.\tG\tA\t60\tPASS\t."])
    commands = []

    def fake_run(command, **_kwargs):
        commands.append(command)
        output.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
        return CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("snippy_ng.utils.vcf_merge.subprocess.run", fake_run)
    result = CliRunner().invoke(merge, [str(first), str(second), "-o", str(output), "--upset-plot", str(plot)])

    assert result.exit_code == 0, result.output
    assert commands == [["bcftools", "merge", "--no-index", "--force-samples", "--output-type", "v", "--output", str(output), str(first), str(second)]]
    svg = plot.read_text(encoding="utf-8")
    assert "VCF call-set overlap" in svg
    assert ">1</text>" in svg


def test_merge_cli_accepts_vcf_file_list(tmp_path, monkeypatch):
    first, second = tmp_path / "first.vcf", tmp_path / "second.vcf"
    output, file_list = tmp_path / "merged.vcf", tmp_path / "vcfs.txt"
    _write_vcf(first, [])
    _write_vcf(second, [])
    file_list.write_text(f"# inputs\n{first}\n{second}\n", encoding="utf-8")
    monkeypatch.setattr("snippy_ng.utils.vcf_merge.subprocess.run", lambda command, **_kwargs: CompletedProcess(command, 0, "", ""))

    result = CliRunner().invoke(merge, ["--file-list", str(file_list), "-o", str(output)])

    assert result.exit_code == 0, result.output
