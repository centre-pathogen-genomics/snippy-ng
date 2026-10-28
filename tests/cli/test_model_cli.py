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
