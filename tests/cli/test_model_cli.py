from pathlib import Path

from click.testing import CliRunner

from snippy_ng.cli import snippy_ng


def test_model_download_cli_uses_default_model_root(monkeypatch, tmp_path):
    captured = {}

    def fake_download(model, output_directory, force=False):
        captured["model"] = model
        captured["outdir"] = output_directory
        captured["force"] = force
        return Path(output_directory) / "model-dir"

    monkeypatch.setattr(
        "snippy_ng.utils.clair3_models.download_clair3_model",
        fake_download,
    )
    monkeypatch.setenv("CLAIR3_MODELS", str(tmp_path))

    result = CliRunner().invoke(
        snippy_ng,
        ["utils", "model", "download", "dna_r10.4.1_e8.2_400bps_hac@v4.2.0"],
    )

    assert result.exit_code == 0, result.output
    assert captured["model"] == "dna_r10.4.1_e8.2_400bps_hac@v4.2.0"
    assert captured["outdir"] == tmp_path
    assert captured["force"] is False
    assert "model-dir" in result.output


def test_model_download_cli_outdir_and_force(monkeypatch, tmp_path):
    captured = {}

    def fake_download(model, output_directory, force=False):
        captured["outdir"] = output_directory
        captured["force"] = force
        return Path(output_directory) / model

    monkeypatch.setattr(
        "snippy_ng.utils.clair3_models.download_clair3_model",
        fake_download,
    )
    outdir = tmp_path / "models"

    result = CliRunner().invoke(
        snippy_ng,
        [
            "utils",
            "model",
            "download",
            "r1041_e82_400bps_sup_v520",
            "--outdir",
            str(outdir),
            "--force",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["outdir"] == outdir
    assert captured["force"] is True


def test_model_check_uses_header_model(monkeypatch, tmp_path):
    models = tmp_path / "models"
    expected = models / "r1041_e82_400bps_hac_v520"
    expected.mkdir(parents=True)
    reads = tmp_path / "reads.fastq"
    reads.write_text(
        "@read1 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.2.0\nACGT\n+\nIIII\n"
    )
    monkeypatch.setenv("CLAIR3_MODELS", str(models))

    result = CliRunner().invoke(
        snippy_ng,
        ["utils", "model", "check", str(reads), "-o", str(tmp_path)],
    )

    assert result.exit_code == 0, result.output
    assert str(expected.resolve()) in result.output
    assert (tmp_path / "snippy.clair3_model.txt").read_text().strip() == str(expected.resolve())


def test_model_check_falls_back_to_longbow(monkeypatch, tmp_path):
    from snippy_ng.stages.calling import PostDownloadClair3ModelSelector

    captured = {}

    def fake_select(reads, prediction_json, model_file, cpus):
        captured["reads"] = reads
        captured["cpus"] = cpus
        model_file.write_text("/models/longbow_model\n")

    monkeypatch.setattr(
        PostDownloadClair3ModelSelector,
        "select_from_header_or_predict",
        staticmethod(fake_select),
    )
    reads = tmp_path / "reads.fastq"
    reads.write_text("@read1\nACGT\n+\nIIII\n")

    result = CliRunner().invoke(
        snippy_ng,
        ["utils", "model", "check", str(reads), "-o", str(tmp_path), "--cpus", "3"],
    )

    assert result.exit_code == 0, result.output
    assert captured["reads"] == reads
    assert captured["cpus"] == 3
    assert "/models/longbow_model" in result.output


def test_model_check_rejects_sra_accession(tmp_path):
    result = CliRunner().invoke(
        snippy_ng,
        ["utils", "model", "check", "SRR1234567", "-o", str(tmp_path)],
    )

    assert result.exit_code == 1
    assert "not an SRA accession" in result.output
