from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from snippy_ng.__about__ import __version__
from snippy_ng.stages.setup import PrepareReference


def _make_genbank(path: Path) -> None:
    record = SeqRecord(
        Seq("A" * 100),
        id="mini",
        name="mini",
        description="Mini organism",
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["date"] = "01-JAN-2026"
    record.features.append(
        SeqFeature(FeatureLocation(0, 100), type="source", qualifiers={"organism": "Mini organism"})
    )
    record.features.append(
        SeqFeature(FeatureLocation(9, 60), type="gene", qualifiers={"locus_tag": "mini_0001"})
    )
    record.features.append(
        SeqFeature(FeatureLocation(9, 60), type="CDS", qualifiers={"locus_tag": "mini_0001"})
    )
    # standalone feature (no gene/transcript hierarchy) -> separate GFF writer branch
    record.features.append(
        SeqFeature(FeatureLocation(70, 85), type="tRNA", qualifiers={"locus_tag": "mini_trna1"})
    )
    SeqIO.write([record], path, "genbank")


def test_prepare_reference_gff_source_is_snippy_ng_version(tmp_path: Path):
    gbk = tmp_path / "mini.gbk"
    _make_genbank(gbk)

    stage = PrepareReference(input=gbk, ref_fmt="genbank", output_directory=tmp_path)
    out_fasta = tmp_path / "mini.fa"
    out_gff = tmp_path / "mini.gff"
    stage.process_reference(gbk, "genbank", out_fasta, out_gff)

    expected_source = f"snippy-ng:{__version__}"
    rows = [
        line.split("\t")
        for line in out_gff.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    # gene, transcript, CDS, and the standalone tRNA's gene row
    assert len(rows) >= 4
    for row in rows:
        assert len(row) == 9, row
        assert row[1] == expected_source, row

    types = [row[2] for row in rows]
    assert types.count("gene") >= 2  # regular gene + standalone tRNA gene
    assert "transcript" in types
    assert "CDS" in types
    # metadata must be written inside the output directory
    assert (tmp_path / "metadata.json").exists()