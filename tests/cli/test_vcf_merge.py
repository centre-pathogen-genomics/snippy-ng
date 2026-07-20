from pathlib import Path
from subprocess import CompletedProcess

from click.testing import CliRunner

from snippy_ng.cli.vcf.merge_cli import merge
from snippy_ng.cli.vcf.plot_cli import plot
from snippy_ng.utils.vcf_merge import _unique_names


def _write_vcf(path: Path, records: list[str]) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "\n".join(records) + "\n",
        encoding="utf-8",
    )


def test_merge_cli_runs_bcftools(tmp_path, monkeypatch):
    first, second = tmp_path / "first.vcf", tmp_path / "second.vcf"
    output = tmp_path / "merged.vcf"
    _write_vcf(first, ["chr1\t10\t.\tA\tT\t60\tPASS\t.", "chr1\t20\t.\tC\tG\t60\tPASS\t."])
    _write_vcf(second, ["chr1\t10\t.\tA\tT\t60\tPASS\t.", "chr1\t30\t.\tG\tA\t60\tPASS\t."])
    commands = []

    def fake_run(command, **_kwargs):
        commands.append(command)
        output.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
        return CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("snippy_ng.utils.vcf_merge.subprocess.run", fake_run)
    result = CliRunner().invoke(merge, [str(first), str(second), "-o", str(output)])

    assert result.exit_code == 0, result.output
    assert commands == [["bcftools", "merge", "--no-index", "--force-samples", "--output-type", "v", "--output", str(output), str(first), str(second)]]


def test_merge_cli_accepts_vcf_file_list(tmp_path, monkeypatch):
    first, second = tmp_path / "first.vcf", tmp_path / "second.vcf"
    output, file_list = tmp_path / "merged.vcf", tmp_path / "vcfs.txt"
    _write_vcf(first, [])
    _write_vcf(second, [])
    file_list.write_text(f"# inputs\n{first}\n{second}\n", encoding="utf-8")
    monkeypatch.setattr("snippy_ng.utils.vcf_merge.subprocess.run", lambda command, **_kwargs: CompletedProcess(command, 0, "", ""))

    result = CliRunner().invoke(merge, ["--file-list", str(file_list), "-o", str(output)])

    assert result.exit_code == 0, result.output


def test_merge_cli_omits_output_path_and_streams_to_stdout(tmp_path, monkeypatch):
    first, second = tmp_path / "first.vcf", tmp_path / "second.vcf"
    _write_vcf(first, [])
    _write_vcf(second, [])
    commands = []

    def fake_run(command, **_kwargs):
        commands.append(command)
        return CompletedProcess(command, 0, "##fileformat=VCFv4.2\n", "")

    monkeypatch.setattr("snippy_ng.utils.vcf_merge.subprocess.run", fake_run)
    result = CliRunner().invoke(merge, [str(first), str(second)])

    assert result.exit_code == 0, result.output
    assert commands == [["bcftools", "merge", "--no-index", "--force-samples", "--output-type", "v", str(first), str(second)]]


def test_plot_upset_cli_builds_svg_from_vcf_paths(tmp_path):
    first, second = tmp_path / "first.vcf", tmp_path / "second.vcf"
    output = tmp_path / "overlap.svg"
    _write_vcf(first, ["chr1\t10\t.\tA\tT\t60\tPASS\t.", "chr1\t20\t.\tC\tG\t60\tPASS\t."])
    _write_vcf(second, ["chr1\t10\t.\tA\tT\t60\tPASS\t.", "chr1\t30\t.\tG\tA\t60\tPASS\t."])

    result = CliRunner().invoke(plot, ["upset", str(first), str(second), "-o", str(output)])

    assert result.exit_code == 0, result.output
    svg = output.read_text(encoding="utf-8")
    assert "VCF call-set overlap" in svg
    assert ">first</text>" in svg
    assert ">second</text>" in svg


def test_plot_upset_cli_accepts_merged_vcf_on_stdin(tmp_path):
    output = tmp_path / "stdin-overlap.svg"
    merged_vcf = "\n".join(
        [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcaller_a\tcaller_b",
            "chr1\t10\t.\tA\tT\t60\tPASS\t.\tGT\t1/1\t0/0",
            "chr1\t20\t.\tC\tG\t60\tPASS\t.\tGT\t1/1\t1/1",
            "chr1\t30\t.\tG\tA\t60\tPASS\t.\tGT\t0/0\t1/1",
        ]
    ) + "\n"

    result = CliRunner().invoke(plot, ["upset", "-o", str(output)], input=merged_vcf)

    assert result.exit_code == 0, result.output
    svg = output.read_text(encoding="utf-8")
    assert "VCF call-set overlap" in svg
    assert ">caller_a</text>" in svg
    assert ">caller_b</text>" in svg


def test_unique_names_prefers_sample_name_from_vcf_header(tmp_path):
    first = tmp_path / "alpha.vcf"
    second = tmp_path / "beta.vcf"
    first.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\n",
        encoding="utf-8",
    )
    second.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_B\n",
        encoding="utf-8",
    )

    assert _unique_names([first, second]) == ["SAMPLE_A", "SAMPLE_B"]


def test_unique_names_falls_back_to_filename_without_sample_columns(tmp_path):
    first = tmp_path / "first.vcf"
    second = tmp_path / "first.vcf.gz"
    first.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        encoding="utf-8",
    )
    second.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        encoding="utf-8",
    )

    assert _unique_names([first, second]) == ["first", "first (2)"]
