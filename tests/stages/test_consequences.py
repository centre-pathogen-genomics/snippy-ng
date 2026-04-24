from pathlib import Path

from snippy_ng.context import Context
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller


def test_consequences_uses_local_csq_by_default(tmp_path, monkeypatch):
    monkeypatch.delenv("SNIPPY_NG_LOCAL_BCFTOOLS_CSQ", raising=False)
    reference = tmp_path / "reference.fa"
    variants = tmp_path / "variants.vcf"
    features = tmp_path / "reference.gff"
    reference.write_text(">ref\nA\n")
    variants.write_text("##fileformat=VCFv4.2\n")
    features.write_text("ref\tsnippy-ng\tCDS\t1\t1\t.\t+\t0\tID=cds1\n")

    stage = BcftoolsConsequencesCaller(
        prefix=str(tmp_path / "snippy"),
        reference=reference,
        variants=variants,
        features=features,
    )

    commands = stage.create_commands(Context(cpus=4))

    assert commands[0].command[:5] == ["bcftools", "csq", "--local-csq", "--threads", "4"]


def test_consequences_can_disable_local_csq_via_env_flag(tmp_path, monkeypatch):
    monkeypatch.setenv("SNIPPY_NG_LOCAL_BCFTOOLS_CSQ", "0")
    reference = tmp_path / "reference.fa"
    variants = tmp_path / "variants.vcf"
    features = tmp_path / "reference.gff"
    reference.write_text(">ref\nA\n")
    variants.write_text("##fileformat=VCFv4.2\n")
    features.write_text("ref\tsnippy-ng\tCDS\t1\t1\t.\t+\t0\tID=cds1\n")

    stage = BcftoolsConsequencesCaller(
        prefix=str(tmp_path / "snippy"),
        reference=reference,
        variants=variants,
        features=features,
    )

    commands = stage.create_commands(Context(cpus=2))

    assert commands[0].command[:4] == ["bcftools", "csq", "--threads", "2"]
    assert "--local-csq" not in commands[0].command