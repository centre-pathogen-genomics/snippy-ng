"""Stages for downloading data from remote sources."""
from pathlib import Path
import gzip
import json
import os
import shutil
import tempfile
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pydantic import Field
from snippy_ng.stages import BaseStage, BaseOutput
from snippy_ng.__about__ import __version__
from snippy_ng.dependencies import sracha
from snippy_ng.exceptions import StageExecutionError


class ReadDownloadOutput(BaseOutput):
    reads: Path = Field(..., description="Downloaded reads file")


class NcbiAssemblyFastaOutput(BaseOutput):
    fasta: Path = Field(..., description="Downloaded NCBI assembly FASTA")


class DownloadNcbiAssemblyFasta(BaseStage):
    _skip_if_outputs_exist = True
    accession: str = Field(..., description="NCBI RefSeq/GenBank assembly accession")
    output_directory: Path | None = Field(default=None, description="Optional directory for downloaded assembly FASTA outputs")
    genbank: bool = Field(default=False, description="Download in GenBank format instead of FASTA")

    @property
    def reference_prefix(self) -> str:
        return self.accession

    @property
    def output(self) -> NcbiAssemblyFastaOutput:
        output_directory = self.output_directory or Path(".")
        ext = ".gbff" if self.genbank else ".fa"
        return NcbiAssemblyFastaOutput(
            fasta=output_directory / f"{self.reference_prefix}{ext}",
        )

    def create_commands(self, ctx):
        description = f"Download NCBI assembly: {self.accession}"
        return [
            self.python_cmd(
                func=self.download_fasta,
                args=(self.accession, self.output.fasta, self.genbank),
                description=description,
            ),
        ]

    def download_fasta(
        self,
        accession: str,
        output_fasta_path: Path,
        genbank: bool = False,
    ) -> None:
        output_fasta_path.parent.mkdir(parents=True, exist_ok=True)

        url = self._download_url(accession, genbank=genbank)
        request = urllib.request.Request(url, headers={"User-Agent": f"snippy-ng/{__version__}"})

        with tempfile.TemporaryDirectory(prefix=f"snippy-ng-{accession}-") as tmpdir:
            zip_path = Path(tmpdir) / f"{accession}.zip"
            try:
                with urllib.request.urlopen(request) as response, open(zip_path, "wb") as handle:
                    shutil.copyfileobj(response, handle)
            except urllib.error.HTTPError as exc:
                raise StageExecutionError(
                    f"NCBI Datasets API could not download assembly FASTA '{accession}' "
                    f"(HTTP {exc.code}). Check the accession and try again."
                ) from exc
            except urllib.error.URLError as exc:
                raise StageExecutionError(
                    f"Could not reach NCBI Datasets API while downloading assembly FASTA '{accession}': {exc.reason}"
                ) from exc

            try:
                with zipfile.ZipFile(zip_path) as archive:
                    if genbank:
                        member = DownloadNcbiAssemblyFasta._find_archive_member(archive, "_genomic.gbff")
                        if member is None:
                            member = DownloadNcbiAssemblyFasta._find_archive_member(archive, ".gbff")
                        if member is None:
                            raise StageExecutionError(f"NCBI assembly package for '{accession}' did not contain a GenBank file (wrong id?).")
                    else:
                        member = DownloadNcbiAssemblyFasta._find_archive_member(archive, "_genomic.fna")
                        if member is None:
                            member = DownloadNcbiAssemblyFasta._find_archive_member(archive, ".fna")
                        if member is None:
                            raise StageExecutionError(f"NCBI assembly package for '{accession}' did not contain a FASTA file (wrong id?).")
                    with archive.open(member) as source, open(output_fasta_path, "wb") as target:
                        shutil.copyfileobj(source, target)
                    report = DownloadNcbiAssemblyFasta._load_assembly_report(archive)
            except zipfile.BadZipFile as exc:
                raise StageExecutionError(f"NCBI returned an invalid zip package for assembly '{accession}' (wrong id?).") from exc

        assembly_name = report.get("assembly_info", {}).get("assembly_name") if report else None
        suffix = f" ({assembly_name})" if assembly_name else ""
        print(f"Downloaded NCBI assembly {accession}{suffix} to {output_fasta_path}")

    @staticmethod
    def _download_url(accession: str, genbank: bool = False) -> str:
        annotation_type = "GENOME_GBFF" if genbank else "GENOME_FASTA"
        query = urllib.parse.urlencode(
            [
                ("include_annotation_type", annotation_type),
            ]
        )
        api_key = os.environ.get("NCBI_API_KEY")
        if api_key:
            query += "&" + urllib.parse.urlencode({"api_key": api_key})
        return f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{urllib.parse.quote(accession)}/download?{query}"

    @staticmethod
    def _find_archive_member(archive: zipfile.ZipFile, suffix: str) -> str | None:
        for name in archive.namelist():
            if name.endswith(suffix) and not name.endswith("/"):
                return name
        return None

    @staticmethod
    def _load_assembly_report(archive: zipfile.ZipFile) -> dict:
        member = DownloadNcbiAssemblyFasta._find_archive_member(archive, "assembly_data_report.jsonl")
        if member is None:
            return {}
        with archive.open(member) as handle:
            first_line = handle.readline().decode("utf-8").strip()
        if not first_line:
            return {}
        try:
            return json.loads(first_line)
        except json.JSONDecodeError:
            return {}


class AtbAssemblyReferenceOutput(BaseOutput):
    fasta: Path = Field(..., description="Downloaded AllTheBacteria assembly FASTA reference")


class DownloadAtbAssemblyReference(BaseStage):
    _skip_if_outputs_exist = True
    accession: str = Field(..., description="AllTheBacteria sample accession")
    output_directory: Path | None = Field(default=None, description="Optional directory for downloaded reference outputs")

    @property
    def reference_prefix(self) -> str:
        return self.accession

    @property
    def output(self) -> AtbAssemblyReferenceOutput:
        output_directory = self.output_directory or Path(".")
        return AtbAssemblyReferenceOutput(
            fasta=output_directory / f"{self.reference_prefix}.fa",
        )

    def create_commands(self, ctx):
        return [
            self.python_cmd(
                func=self.download_assembly_reference,
                args=(self.accession, self.output.fasta),
                description=f"Download AllTheBacteria assembly: {self.accession}",
            ),
        ]

    def download_assembly_reference(
        self,
        accession: str,
        output_fasta_path: Path,
    ) -> None:
        output_fasta_path.parent.mkdir(parents=True, exist_ok=True)

        url = self._download_url(accession)
        request = urllib.request.Request(url, headers={"User-Agent": f"snippy-ng/{__version__}"})

        try:
            with urllib.request.urlopen(request) as response:
                compressed = response.read()
        except urllib.error.HTTPError as exc:
            raise StageExecutionError(
                f"AllTheBacteria S3 could not download assembly '{accession}' "
                f"(HTTP {exc.code}). Check the accession and try again."
            ) from exc
        except urllib.error.URLError as exc:
            raise StageExecutionError(
                f"Could not reach AllTheBacteria S3 while downloading assembly '{accession}': {exc.reason}"
            ) from exc

        try:
            fasta_bytes = gzip.decompress(compressed)
        except OSError as exc:
            raise StageExecutionError(f"AllTheBacteria assembly '{accession}' was not valid gzip-compressed FASTA.") from exc

        output_fasta_path.write_bytes(fasta_bytes)
        print(f"Downloaded AllTheBacteria assembly {accession} to {output_fasta_path}")

    @staticmethod
    def _download_url(accession: str) -> str:
        return f"https://allthebacteria-assemblies.s3.eu-west-2.amazonaws.com/{urllib.parse.quote(accession)}.fa.gz"


class DownloadSraReads(BaseStage):
    """Download reads from SRA using sracha-rs."""

    _dependencies = [sracha]
    _skip_if_outputs_exist = True

    accession: str = Field(..., description="SRA accession (SRR, ERR, or DRR format)")
    output_directory: Path | None = Field(default=None, description="Optional directory for downloaded reads outputs")

    @property
    def read_prefix(self) -> str:
        return self.accession

    @property
    def output(self) -> ReadDownloadOutput:
        output_directory = self.output_directory or Path(".")
        # SRA downloads produce FASTQ files; naming depends on single-end vs paired-end
        # We'll use the accession as the base name and sracha-rs will expand it
        return ReadDownloadOutput(
            reads=output_directory / f"{self.read_prefix}.fastq.gz",
        )

    def create_commands(self, ctx):
        output_directory = self.output_directory or Path(".")
        return [
            self.shell_cmd(
                ["sracha", "get", str(self.accession), "--output-dir", str(output_directory), "--prefer-ena"],
                description=f"Download SRA reads: {self.accession}",
            ),
        ]
