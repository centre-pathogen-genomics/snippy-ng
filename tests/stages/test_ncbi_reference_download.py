import io
import gzip
import json
import zipfile
from pathlib import Path

from snippy_ng.pipelines.common import (
    is_atb_assembly_accession,
    is_ncbi_assembly_accession,
    is_reference_accession,
)
from snippy_ng.pipelines.asm import AsmPipelineBuilder
from snippy_ng.stages.setup import DownloadAtbAssemblyReference, DownloadNcbiGenbankReference, PrepareReference


def _ncbi_package() -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr(
            "ncbi_dataset/data/GCF_000000001.1/GCF_000000001.1_genomic.gbff",
            "LOCUS       chr1                       4 bp    DNA     linear   BCT 01-JAN-2000\n"
            "FEATURES             Location/Qualifiers\n"
            "ORIGIN\n"
            "        1 acgt\n"
            "//\n",
        )
        archive.writestr(
            "ncbi_dataset/data/assembly_data_report.jsonl",
            json.dumps({"assembly_info": {"assembly_name": "ASM1"}}) + "\n",
        )
    return buffer.getvalue()


def test_detects_ncbi_refseq_and_genbank_assembly_accessions():
    assert is_ncbi_assembly_accession("GCF_000000001.1")
    assert is_ncbi_assembly_accession("GCA_000000001")
    assert not is_ncbi_assembly_accession("NC_000001.11")
    assert not is_ncbi_assembly_accession("reference.fasta")


def test_detects_allthebacteria_assembly_accessions():
    assert is_atb_assembly_accession("SAMN123456")
    assert is_atb_assembly_accession("SAMEA123456")
    assert is_atb_assembly_accession("SAMD123456")
    assert not is_atb_assembly_accession("GCF_000000001.1")
    assert not is_atb_assembly_accession("SAMX123456")


def test_detects_supported_reference_accessions():
    assert is_reference_accession("GCF_000000001.1")
    assert is_reference_accession("SAMN123456")
    assert not is_reference_accession("reference.fasta")


def test_ncbi_reference_download_extracts_genbank_file(monkeypatch, tmp_path):
    class Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            self.close()

    def fake_urlopen(request):
        assert request.full_url.startswith("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000000001.1/download?")
        assert "include_annotation_type=GENOME_GBFF" in request.full_url
        return Response(_ncbi_package())

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)

    stage = DownloadNcbiGenbankReference(
        accession="GCF_000000001.1",
        output_directory=tmp_path / "reference",
    )

    stage.download_genbank_reference(
        stage.accession,
        stage.output.genbank,
    )

    assert stage.output.genbank == tmp_path / "reference" / "GCF_000000001.1.gbff"
    assert "LOCUS       chr1" in stage.output.genbank.read_text(encoding="utf-8")


def test_assembly_pipeline_downloads_ncbi_genbank_before_prepare_reference(tmp_path):
    assembly = tmp_path / "assembly.fa"
    assembly.write_text(">contig1\nACGT\n", encoding="utf-8")

    pipeline = AsmPipelineBuilder(
        reference_accession="GCF_000000001.1",
        assembly=assembly,
        prefix="sample",
    ).build()

    assert isinstance(pipeline.stages[0], DownloadNcbiGenbankReference)
    assert pipeline.stages[0].output.genbank == Path("reference/GCF_000000001.1.gbff")
    assert isinstance(pipeline.stages[1], PrepareReference)
    assert pipeline.stages[1].input == pipeline.stages[0].output.genbank
    assert pipeline.stages[1].ref_fmt == "genbank"


def test_atb_reference_download_extracts_fasta_file(monkeypatch, tmp_path):
    class Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            self.close()

    def fake_urlopen(request):
        assert request.full_url == "https://allthebacteria-assemblies.s3.eu-west-2.amazonaws.com/SAMN123456.fa.gz"
        return Response(gzip.compress(b">contig1\nACGT\n"))

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)

    stage = DownloadAtbAssemblyReference(
        accession="SAMN123456",
        output_directory=tmp_path / "reference",
    )

    stage.download_assembly_reference(
        stage.accession,
        stage.output.fasta,
    )

    assert stage.output.fasta == tmp_path / "reference" / "SAMN123456.fa"
    assert stage.output.fasta.read_text(encoding="utf-8") == ">contig1\nACGT\n"


def test_assembly_pipeline_downloads_atb_fasta_before_prepare_reference(tmp_path):
    assembly = tmp_path / "assembly.fa"
    assembly.write_text(">contig1\nACGT\n", encoding="utf-8")

    pipeline = AsmPipelineBuilder(
        reference_accession="SAMN123456",
        assembly=assembly,
        prefix="sample",
    ).build()

    assert isinstance(pipeline.stages[0], DownloadAtbAssemblyReference)
    assert pipeline.stages[0].output.fasta == Path("reference/SAMN123456.fa")
    assert isinstance(pipeline.stages[1], PrepareReference)
    assert pipeline.stages[1].input == pipeline.stages[0].output.fasta
    assert pipeline.stages[1].ref_fmt == "fasta"
