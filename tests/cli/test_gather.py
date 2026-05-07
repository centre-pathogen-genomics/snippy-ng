import csv
import json
from io import StringIO

from click.testing import CliRunner

from snippy_ng.cli import snippy_ng


def test_gather_default_applies_to_csv_rows(monkeypatch, tmp_path):
    monkeypatch.setattr("snippy_ng.logging.logger.info", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        "snippy_ng.utils.gather.gather_samples_config",
        lambda **_: {
            "samples": {
                "sample_a": {"type": "short", "reads1": "a_R1.fastq.gz"},
                "sample_b": {"type": "long", "reads": "b.fastq.gz"},
            }
        },
    )

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
    monkeypatch.setattr(
        "snippy_ng.utils.gather.gather_samples_config",
        lambda **_: {
            "samples": {
                "sample_a": {"type": "short", "platform": "nanopore"},
                "sample_b": {"type": "long"},
            }
        },
    )

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
    monkeypatch.setattr(
        "snippy_ng.utils.gather.gather_samples_config",
        lambda **_: {
            "samples": {
                "sample_a": {"type": "long", "clair3_model": None},
                "sample_b": {"type": "long", "clair3_model": "custom_model"},
            }
        },
    )

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
