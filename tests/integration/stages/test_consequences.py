from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from snippy_ng.context import Context
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller


pytestmark = pytest.mark.integration_sim


def _parse_info_field(info: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for field in info.split(";"):
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
            values[key] = value
        else:
            values[field] = ""
    return values


def _annotated_records(path: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            records.append(
                {
                    "chrom": cols[0],
                    "pos": cols[1],
                    "ref": cols[3],
                    "alt": cols[4],
                    "filter": cols[6],
                    "info": cols[7],
                    "format": cols[8] if len(cols) > 8 else "",
                    "sample": cols[9] if len(cols) > 9 else "",
                }
            )
    return records


def test_bcftools_csq_compound_reverse_strand_bug_reproducer_is_annotated(tmp_path: Path):
    if shutil.which("bcftools") is None:
        pytest.skip("bcftools required for this consequence regression test")

    project_root = Path(__file__).resolve().parents[3]
    fixture_dir = project_root / "bcftools-csq-compound-bug"
    stage = BcftoolsConsequencesCaller(
        prefix=str(tmp_path / "compound-csq"),
        reference=fixture_dir / "mini.fa",
        variants=fixture_dir / "mini.vcf",
        features=fixture_dir / "mini.gff",
    )

    stage.run(Context(cpus=1, quiet=True))
    stage.run_tests()

    records = _annotated_records(stage.output.annotated_vcf)

    assert len(records) == 2
    assert [(record["chrom"], record["pos"], record["ref"], record["alt"]) for record in records] == [
        ("mini", "2214", "T", "G"),
        ("mini", "2216", "A", "AT"),
    ]
    assert all(record["filter"] == "PASS" for record in records)
    assert all(record["format"] == "GT:BCSQ" for record in records)
    assert all(record["sample"] == "1/1:3" for record in records)

    annotations = [_parse_info_field(record["info"])["BCSQ"] for record in records]
    assert annotations == [
        "missense|hypF|BFV63_RS17260_transcript_1|protein_coding|-|2S>2R|2214T>G",
        "start_lost|hypF|BFV63_RS17260_transcript_1|protein_coding|-",
    ]