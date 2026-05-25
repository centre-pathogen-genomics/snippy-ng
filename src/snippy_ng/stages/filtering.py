from pathlib import Path
from typing import List, Optional
from snippy_ng.stages import BaseStage, BaseOutput, ShellCommand, TempPath
from snippy_ng.dependencies import samtools
from snippy_ng.exceptions import StageExecutionError
from pydantic import Field, field_validator


class BamReferenceValidationOutput(BaseOutput):
    header: TempPath = Field(..., description="Temporary BAM/CRAM header")
    mapped_count: TempPath = Field(..., description="Temporary mapped-read count")


class BamReferenceValidator(BaseStage):
    """
    Validate that a user-provided BAM/CRAM is aligned and uses the selected reference.
    """

    bam: Path = Field(..., description="Input BAM/CRAM file to validate")
    reference_index: Path = Field(..., description="Reference FASTA index (.fai)")
    reference: Path = Field(..., description="Reference FASTA file")

    _dependencies = [samtools]

    @property
    def output(self) -> BamReferenceValidationOutput:
        return BamReferenceValidationOutput(
            header=f"{self.prefix}.bam_header.txt",
            mapped_count=f"{self.prefix}.mapped_count.txt",
        )

    def create_commands(self, ctx) -> List:
        return [
            self.shell_cmd(
                ["samtools", "view", "-H", str(self.bam)],
                description=f"Read BAM/CRAM header: {self.bam}",
                output_file=self.output.header,
            ),
            self.shell_cmd(
                ["samtools", "view", "-c", "-F", "4", "-T", str(self.reference), str(self.bam)],
                description=f"Count mapped reads in BAM/CRAM: {self.bam}",
                output_file=self.output.mapped_count,
            ),
            self.python_cmd(
                func=self.validate_bam_reference,
                description=f"Validate BAM/CRAM alignment and reference compatibility: {self.bam}",
            )
        ]

    def validate_bam_reference(self) -> None:
        reference_sequences = self._read_fai(self.reference_index)
        bam_sequences = self._read_bam_header_sequences(self.output.header)

        if not bam_sequences:
            raise StageExecutionError(
                f"Input alignment '{self.bam}' has no @SQ reference records. "
                "This looks like an unaligned BAM/CRAM; provide FASTQ reads or align "
                "the BAM to the selected reference first."
            )

        self._validate_reference_sequences_match(bam_sequences, reference_sequences)

        mapped_count = self._read_mapped_count(self.output.mapped_count)
        if mapped_count == 0:
            raise StageExecutionError(
                f"Input alignment '{self.bam}' has no mapped reads. "
                "Provide an aligned BAM/CRAM, or provide reads so Snippy-NG can align them."
            )

    @staticmethod
    def _read_fai(reference_index: Path) -> dict[str, int]:
        sequences: dict[str, int] = {}
        try:
            with open(reference_index, "r", encoding="utf-8") as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) >= 2:
                        sequences[fields[0]] = int(fields[1])
        except FileNotFoundError as exc:
            raise StageExecutionError(f"Reference index does not exist: {reference_index}") from exc
        except ValueError as exc:
            raise StageExecutionError(f"Reference index has invalid sequence lengths: {reference_index}") from exc

        if not sequences:
            raise StageExecutionError(f"Reference index has no sequences: {reference_index}")
        return sequences

    @staticmethod
    def _read_bam_header_sequences(header_path: Path) -> dict[str, int]:
        sequences: dict[str, int] = {}
        try:
            with open(header_path, "r", encoding="utf-8") as handle:
                for line in handle:
                    if not line.startswith("@SQ\t"):
                        continue
                    fields = {}
                    for item in line.rstrip("\n").split("\t")[1:]:
                        if ":" not in item:
                            continue
                        key, value = item.split(":", 1)
                        fields[key] = value
                    name = fields.get("SN")
                    length = fields.get("LN")
                    if name is None or length is None:
                        continue
                    try:
                        sequences[name] = int(length)
                    except ValueError as exc:
                        raise StageExecutionError(f"BAM/CRAM header has invalid length for contig '{name}'") from exc
        except FileNotFoundError as exc:
            raise StageExecutionError(f"BAM/CRAM header output does not exist: {header_path}") from exc
        return sequences

    @staticmethod
    def _read_mapped_count(mapped_count_path: Path) -> int:
        try:
            with open(mapped_count_path, "r", encoding="utf-8") as handle:
                count_text = handle.read().strip()
        except FileNotFoundError as exc:
            raise StageExecutionError(f"Mapped-read count output does not exist: {mapped_count_path}") from exc

        try:
            return int(count_text)
        except ValueError as exc:
            raise StageExecutionError(
                f"samtools returned an invalid mapped-read count in '{mapped_count_path}'"
            ) from exc

    @staticmethod
    def _validate_reference_sequences_match(
        bam_sequences: dict[str, int],
        reference_sequences: dict[str, int],
    ) -> None:
        bam_names = set(bam_sequences)
        reference_names = set(reference_sequences)
        missing_from_bam = sorted(reference_names - bam_names)
        missing_from_reference = sorted(bam_names - reference_names)
        length_mismatches = sorted(
            name
            for name in bam_names & reference_names
            if bam_sequences[name] != reference_sequences[name]
        )

        if not missing_from_bam and not missing_from_reference and not length_mismatches:
            return

        details = []
        if missing_from_bam:
            details.append(
                "reference contigs absent from BAM header: "
                f"{BamReferenceValidator._summarise_names(missing_from_bam)}"
            )
        if missing_from_reference:
            details.append(
                "BAM contigs absent from reference: "
                f"{BamReferenceValidator._summarise_names(missing_from_reference)}"
            )
        if length_mismatches:
            examples = [
                f"{name} (BAM={bam_sequences[name]}, reference={reference_sequences[name]})"
                for name in length_mismatches[:5]
            ]
            if len(length_mismatches) > 5:
                examples.append(f"... {len(length_mismatches) - 5} more")
            details.append(f"contig length mismatches: {', '.join(examples)}")

        raise StageExecutionError(
            "Input BAM/CRAM was not aligned to the selected reference. "
            + "; ".join(details)
            + ". Re-align the reads to this reference or rerun with the matching reference."
        )

    @staticmethod
    def _summarise_names(names: list[str]) -> str:
        shown = names[:5]
        suffix = f", ... {len(names) - 5} more" if len(names) > 5 else ""
        return ", ".join(shown) + suffix


class SamtoolsFilterOutput(BaseOutput):
    bam: Path = Field(..., description="Filtered alignment file in BAM format")
    bai: Path = Field(..., description="Index file for the filtered BAM")
    stats: Path = Field(..., description="Flagstat output file for the BAM file")


class SamtoolsFilter(BaseStage):
    """
    Filter BAM files using Samtools to remove unwanted alignments.
    """
    
    bam: Path = Field(..., description="Input BAM file to filter")
    reference: Path = Field(..., description="Reference FASTA file for BAM output")
    min_mapq: int = Field(20, description="Minimum mapping quality")
    exclude_flags: int = Field(1796, description="SAM flags to exclude (default: unmapped, secondary, qcfail, duplicate)")
    include_flags: Optional[int] = Field(None, description="SAM flags to include")
    regions: Optional[str] = Field(None, description="Regions to include (BED file or region string)")
    additional_filters: str = Field("", description="Additional samtools view options")
    
    _dependencies = [samtools]
    
    @property
    def output(self) -> SamtoolsFilterOutput:
        filtered_bam = f"{self.prefix}.filtered.bam"
        return SamtoolsFilterOutput(
            bam=filtered_bam,
            bai=f"{filtered_bam}.bai",
            stats=f"{filtered_bam}.flagstat.txt"
        )
    
    def build_filter_command(self, ctx) -> ShellCommand:
        """Constructs the samtools view command for filtering."""
        cmd_parts = [
            "samtools",
            "view",
            "-O",
            "bam",
            "--reference", str(self.reference),
            "-o",
            str(self.output.bam),
        ]
        
        # Add threading
        if ctx.cpus > 1:
            cmd_parts.extend(["--threads", str(ctx.cpus - 1)])
        
        # Add mapping quality filter
        if self.min_mapq > 0:
            cmd_parts.extend(["--min-MQ", str(self.min_mapq)])
        
        # Add flag filters
        if self.exclude_flags:
            cmd_parts.extend(["-F", str(self.exclude_flags)])
        
        if self.include_flags is not None:
            cmd_parts.extend(["-f", str(self.include_flags)])
        
        # Add regions if specified as BED file
        if self.regions and Path(self.regions).exists():
            cmd_parts.extend(["-L", str(self.regions)])
        
        # Add additional filters (split if it contains spaces)
        if self.additional_filters:
            import shlex
            cmd_parts.extend(shlex.split(self.additional_filters))
        
        # Add input file
        cmd_parts.append(str(self.bam))
        
        # Add region string if not a file
        if self.regions and not Path(self.regions).exists():
            cmd_parts.append(str(self.regions))
        
        filter_cmd = self.shell_cmd(
            command=cmd_parts,
            description=f"Filter BAM file with MAPQ>={self.min_mapq}, flags={self.exclude_flags}",
        )
        
        return filter_cmd
    
    def build_index_command(self):
        """Returns the samtools index command."""
        return self.shell_cmd([
            "samtools", "index", str(self.output.bam)
        ], description=f"Index filtered BAM file: {self.output.bam}")
    
    def build_flagstat_command(self):
        """Returns the samtools flagstat command."""
        return self.shell_cmd(
            ["samtools", "flagstat", str(self.output.bam)],
            description=f"Generate alignment statistics for {self.output.bam}",
            output_file=Path(self.output.stats),
        )

    def test_mapped_reads_not_zero(self):
        """Test that the BAM file contains mapped reads."""
        flagstat_path = self.output.stats
        if not flagstat_path.exists():
            raise FileNotFoundError(
                f"Expected flagstat output file {flagstat_path} was not created"
            )

        with open(flagstat_path) as f:
            for line in f:
                if "mapped (" in line:
                    mapped_reads = int(line.split()[0])
                    if mapped_reads == 0:
                        raise ValueError(
                            "No reads were mapped in BAM file after filtering. Did you use the correct reference?"
                        )
                    break

    def create_commands(self, ctx) -> List:
        """Constructs the filtering commands."""
        filter_cmd = self.build_filter_command(ctx)
        index_cmd = self.build_index_command()
        stats_cmd = self.build_flagstat_command()
        return [filter_cmd, index_cmd, stats_cmd]


class SamtoolsFilterByRegion(SamtoolsFilter):
    """
    Filter BAM file to include only alignments in specified regions.
    """
    
    regions: str = Field(..., description="Regions to include (BED file or region string)")
    
    @field_validator("regions")
    @classmethod
    def validate_regions(cls, v):
        if not v:
            raise ValueError("Regions must be specified")
        return v


class SamtoolsFilterByQuality(SamtoolsFilter):
    """
    Filter BAM file based on mapping quality and alignment flags.
    """
    
    min_mapq: int = Field(30, description="Minimum mapping quality (higher than default)")
    exclude_flags: int = Field(3844, description="SAM flags to exclude (default + supplementary)")


class SamtoolsFilterProperPairs(SamtoolsFilter):
    """
    Filter BAM file to include only properly paired reads.
    """
    
    include_flags: int = Field(2, description="Include only properly paired reads")
    exclude_flags: int = Field(1796, description="Exclude unmapped, secondary, qcfail, duplicate")
