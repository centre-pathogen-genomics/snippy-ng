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
from typing import Annotated
from pydantic import Field
from snippy_ng.stages import BaseStage, BaseOutput
from snippy_ng.__about__ import __version__
from snippy_ng.dependencies import sracha
from snippy_ng.envvars import EnvVarField
from snippy_ng.exceptions import StageExecutionError


class BaseDownloadStage(BaseStage):
    """Base class for all download stages.

    Provides common fields and helpers shared across download stages:
    - ``_skip_if_outputs_exist`` is always ``True`` so completed downloads are not re-run.
    - ``output_directory`` declares where downloaded files are placed (defaults to the
      current working directory when ``None``).
    - ``resolved_output_directory`` returns the effective ``Path`` to use in ``output``
      properties and ``create_commands`` implementations.
    """

    _skip_if_outputs_exist = True
    output_directory: Annotated[Path, EnvVarField(
        Path("data"),
        env_var="DOWNLOAD_DATA_DIR",
        description="Directory for downloaded output files (default: data, override with SNIPPY_NG_DOWNLOAD_DATA_DIR)",
        parser=lambda v: Path(v),
    )]

    @property
    def resolved_output_directory(self) -> Path:
        """Return the effective output directory."""
        return self.output_directory


class ReadDownloadOutput(BaseOutput):
    r1: Path = Field(..., description="Downloaded reads R1 file (or single-end reads file)")
    r2: Path | None = Field(default=None, description="Downloaded reads R2 file (paired-end only)")


class NcbiAssemblyFastaOutput(BaseOutput):
    fasta: Path = Field(..., description="Downloaded NCBI assembly FASTA")


class DownloadNcbiAssemblyFasta(BaseDownloadStage):
    accession: str = Field(..., description="NCBI RefSeq/GenBank assembly accession")
    genbank: bool = Field(default=False, description="Download in GenBank format instead of FASTA")

    @property
    def reference_prefix(self) -> str:
        return self.accession

    @property
    def output(self) -> NcbiAssemblyFastaOutput:
        ext = ".gbff" if self.genbank else ".fa"
        return NcbiAssemblyFastaOutput(
            fasta=self.resolved_output_directory / f"{self.reference_prefix}{ext}",
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


class DownloadAtbAssemblyFasta(BaseDownloadStage):
    accession: str = Field(..., description="AllTheBacteria sample accession")

    @property
    def reference_prefix(self) -> str:
        return self.accession

    @property
    def output(self) -> AtbAssemblyReferenceOutput:
        return AtbAssemblyReferenceOutput(
            fasta=self.resolved_output_directory / f"{self.reference_prefix}.fa",
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


class DownloadSraReads(BaseDownloadStage):
    """Download reads from SRA using sracha-rs."""

    _dependencies = [sracha]

    accession: str = Field(..., description="SRA accession (SRR, ERR, or DRR format)")

    @property
    def read_prefix(self) -> str:
        return self.accession

    @property
    def output(self) -> ReadDownloadOutput:
        # sracha outputs paired-end reads as {acc}_1.fastq.gz / {acc}_2.fastq.gz
        # and single-end reads as {acc}.fastq.gz
        return ReadDownloadOutput(
            r1=self.resolved_output_directory / f"{self.read_prefix}_1.fastq.gz",
            r2=self.resolved_output_directory / f"{self.read_prefix}_2.fastq.gz",
        )

    def error_if_outputs_missing(self, include_temporary: bool = False):
        """Override to also accept single-end output ({acc}.fastq.gz)."""
        from snippy_ng.exceptions import MissingOutputError
        single_end = self.resolved_output_directory / f"{self.read_prefix}.fastq.gz"
        paired_r1 = self.resolved_output_directory / f"{self.read_prefix}_1.fastq.gz"
        if single_end.exists() or paired_r1.exists():
            return
        raise MissingOutputError(
            f"Expected output files are missing: r1 ({paired_r1}) or single-end reads ({single_end})"
        )

    def _normalize_sracha_output(self):
        """Rename single-end {acc}.fastq.gz → {acc}_1.fastq.gz for consistent output naming."""
        single_end = self.resolved_output_directory / f"{self.read_prefix}.fastq.gz"
        r1 = self.resolved_output_directory / f"{self.read_prefix}_1.fastq.gz"
        if single_end.exists() and not r1.exists():
            shutil.move(str(single_end), str(r1))

    def create_commands(self, ctx):
        return [
            self.shell_cmd(
                ["sracha", "get", str(self.accession), "--output-dir", str(self.resolved_output_directory), "--prefer-ena"],
                description=f"Download SRA reads: {self.accession}",
            ),
            self.python_cmd(
                self._normalize_sracha_output,
                description=f"Normalize sracha output paths for {self.accession}",
            ),
        ]
