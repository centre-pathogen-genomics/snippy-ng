import io
import gzip
import json
import zipfile
from pathlib import Path

from snippy_ng.pipelines.common import download_assembly_fasta, is_atb_assembly_accession, is_ncbi_assembly_accession, is_reference_accession
from snippy_ng.pipelines.asm import AsmPipelineBuilder
from snippy_ng.stages.download import DownloadAtbAssemblyReference, DownloadNcbiAssemblyFasta
from snippy_ng.stages.setup import PrepareReference


def _ncbi_fasta_package() -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr(
            "ncbi_dataset/data/GCF_000000001.1/GCF_000000001.1_genomic.fna",
            ">contig1\nACGT\n",
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


def test_assembly_pipeline_downloads_ncbi_fasta_before_prepare_reference(tmp_path):
    assembly = tmp_path / "assembly.fa"
    assembly.write_text(">contig1\nACGT\n", encoding="utf-8")

    pipeline = AsmPipelineBuilder(
        reference_accession="GCF_000000001.1",
        assembly=assembly,
        prefix="sample",
    ).build()

    assert isinstance(pipeline.stages[0], DownloadNcbiAssemblyFasta)
    assert pipeline.stages[0].output.fasta == Path("reference/GCF_000000001.1.fa")
    assert isinstance(pipeline.stages[1], PrepareReference)
    assert pipeline.stages[1].input == pipeline.stages[0].output.fasta
    assert pipeline.stages[1].ref_fmt == "fasta"


def test_ncbi_assembly_download_extracts_fasta_file(monkeypatch, tmp_path):
    class Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            self.close()

    def fake_urlopen(request):
        assert request.full_url.startswith("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000000001.1/download?")
        assert "include_annotation_type=GENOME_FASTA" in request.full_url
        return Response(_ncbi_fasta_package())

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)

    stage = DownloadNcbiAssemblyFasta(
        accession="GCF_000000001.1",
        output_directory=tmp_path / "data",
    )

    stage.download_fasta(
        stage.accession,
        stage.output.fasta,
    )

    assert stage.output.fasta == tmp_path / "data" / "GCF_000000001.1.fa"
    assert stage.output.fasta.read_text(encoding="utf-8") == ">contig1\nACGT\n"


def test_download_assembly_fasta_uses_data_output_directory(tmp_path):
    stages = []

    fasta = download_assembly_fasta(
        "GCF_000000001.1",
        stages,
        output_directory=Path("data"),
    )

    assert isinstance(stages[0], DownloadNcbiAssemblyFasta)
    assert fasta == Path("data/GCF_000000001.1.fa")


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


def test_assembly_pipeline_downloads_assembly_accession_to_data_before_alignment(tmp_path):
    reference = tmp_path / "reference.fa"
    reference.write_text(">ref\nACGT\n", encoding="utf-8")

    pipeline = AsmPipelineBuilder(
        reference=reference,
        assembly_accession="GCF_000000001.1",
        prefix="sample",
    ).build()

    assert isinstance(pipeline.stages[0], DownloadNcbiAssemblyFasta)
    assert pipeline.stages[0].output.fasta == Path("data/GCF_000000001.1.fa")
    assert isinstance(pipeline.stages[1], PrepareReference)
    assert pipeline.stages[2].assembly == Path("data/GCF_000000001.1.fa")
