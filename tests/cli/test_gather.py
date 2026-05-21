import csv
import json
from io import StringIO

from click.testing import CliRunner

from snippy_ng.cli import snippy_ng


def test_gather_default_applies_to_csv_rows(monkeypatch, tmp_path):
    monkeypatch.setattr("snippy_ng.logging.logger.info", lambda *_args, **_kwargs: None)

    def fake_gather_samples_config(**kwargs):
        defaults = kwargs.get("defaults") or {}
        samples = {
            "sample_a": {"type": "short", "reads1": "a_R1.fastq.gz"},
            "sample_b": {"type": "long", "reads": "b.fastq.gz"},
        }
        for sample_data in samples.values():
            for key, value in defaults.items():
                if sample_data.get(key) is None:
                    sample_data[key] = value
        return {"samples": samples}

    monkeypatch.setattr("snippy_ng.utils.gather.gather_samples_config", fake_gather_samples_config)

    result = CliRunner().invoke(
        snippy_ng,
        ["utils", "gather", str(tmp_path), "--default", "platform", "illumina"],
    )

    assert result.exit_code == 0, result.output

    rows = list(csv.DictReader(StringIO(result.output)))
    assert len(rows) == 2
    assert rows[0]["platform"] == "illumina"
    assert rows[1]["platform"] == "illumina"


def test_gather_default_keeps_existing_value_and_supports_multiple(monkeypatch, tmp_path):
    monkeypatch.setattr("snippy_ng.logging.logger.info", lambda *_args, **_kwargs: None)

    def fake_gather_samples_config(**kwargs):
        defaults = kwargs.get("defaults") or {}
        samples = {
            "sample_a": {"type": "short", "platform": "nanopore"},
            "sample_b": {"type": "long"},
        }
        for sample_data in samples.values():
            for key, value in defaults.items():
                if sample_data.get(key) is None:
                    sample_data[key] = value
        return {"samples": samples}

    monkeypatch.setattr("snippy_ng.utils.gather.gather_samples_config", fake_gather_samples_config)

    result = CliRunner().invoke(
        snippy_ng,
        [
            "utils",
            "gather",
            str(tmp_path),
            "--json",
            "--default",
            "platform",
            "illumina",
            "--default",
            "country",
            "uk",
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    samples = payload["samples"]

    assert samples["sample_a"]["platform"] == "nanopore"
    assert samples["sample_b"]["platform"] == "illumina"
    assert samples["sample_a"]["country"] == "uk"
    assert samples["sample_b"]["country"] == "uk"


def test_gather_default_replaces_none_values(monkeypatch, tmp_path):
    monkeypatch.setattr("snippy_ng.logging.logger.info", lambda *_args, **_kwargs: None)

    def fake_gather_samples_config(**kwargs):
        defaults = kwargs.get("defaults") or {}
        samples = {
            "sample_a": {"type": "long", "clair3_model": None},
            "sample_b": {"type": "long", "clair3_model": "custom_model"},
        }
        for sample_data in samples.values():
            for key, value in defaults.items():
                if sample_data.get(key) is None:
                    sample_data[key] = value
        return {"samples": samples}

    monkeypatch.setattr("snippy_ng.utils.gather.gather_samples_config", fake_gather_samples_config)

    result = CliRunner().invoke(
        snippy_ng,
        [
            "utils",
            "gather",
            str(tmp_path),
            "--json",
            "--default",
            "clair3_model",
            "clair3_models/r1041_e82_400bps_sup_v520",
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    samples = payload["samples"]

    assert samples["sample_a"]["clair3_model"] == "clair3_models/r1041_e82_400bps_sup_v520"
    assert samples["sample_b"]["clair3_model"] == "custom_model"


def test_gather_accepts_reference_directory(monkeypatch, tmp_path):
    monkeypatch.setattr("snippy_ng.logging.logger.info", lambda *_args, **_kwargs: None)
    reference_dir = tmp_path / "reference"
    reference_dir.mkdir()

    captured = {}

    def fake_gather_samples_config(**kwargs):
        captured.update(kwargs)
        return {"reference": str(reference_dir), "samples": {}}

    monkeypatch.setattr("snippy_ng.utils.gather.gather_samples_config", fake_gather_samples_config)

    result = CliRunner().invoke(
        snippy_ng,
        ["utils", "gather", str(tmp_path), "--ref", str(reference_dir), "--json"],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["reference"] == str(reference_dir)
    assert captured["reference"] == reference_dir
