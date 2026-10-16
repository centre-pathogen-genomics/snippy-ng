from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pydantic import Field

from snippy_ng.dependencies import biopython, core_snp_filter, distle
from snippy_ng.exceptions import (
    MSAValidationError,
    MissingInputError,
    StageExecutionError,
)
from snippy_ng.stages import BaseOutput, BaseStage


class CombineFastaFileOutput(BaseOutput):
    aln: Path = Field(
        ..., description="Combined multi-sample alignment in FASTA format"
    )


class CombineFastaFile(BaseStage):
    """
    Build a concatenated multi-sample alignment by collapsing multi-contig
    pseudo-alignments into a single supercontig per sample.
    """

    snippy_dirs: list[Path] = Field(
        ..., description="List of snippy output directories to combine"
    )
    reference: Path = Field(
        ..., description="Reference FASTA used to define contig order"
    )
    reference_id: str = Field(
        default="reference",
        description="ID to use for reference sequence in output MSA",
    )
    mask_bed: Optional[Path] = Field(
        default=None,
        description="BED file (reference coordinates) of regions to mask in every sequence of the alignment",
    )
    mask_char: str = Field(
        default="N",
        min_length=1,
        max_length=1,
        description="Character to use for masking",
    )

    _dependencies = [biopython]

    @property
    def output(self) -> CombineFastaFileOutput:
        return CombineFastaFileOutput(aln=Path(f"{self.prefix}.full.aln"))

    def create_commands(self, ctx):
        return [
            self.python_cmd(
                func=self.build_concatenated_alignment,
                args=(self.snippy_dirs, self.reference, self.output.aln),
                description="Combine FASTAs into a concatenated MSA",
            )
        ]

    @staticmethod
    def load_mask_intervals(
        mask_bed: Path,
        contig_order: list[str],
        contig_lengths: Dict[str, int],
    ) -> list[tuple[int, int]]:
        """
        Read a BED file and convert it to sorted, non-overlapping half-open
        intervals in the coordinates of the concatenated alignment.
        """
        offsets: Dict[str, int] = {}
        running = 0
        for contig in contig_order:
            offsets[contig] = running
            running += contig_lengths[contig]

        intervals: list[tuple[int, int]] = []
        with mask_bed.open("r") as handle:
            for line_number, line in enumerate(handle, start=1):
                line = line.strip()
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.split("\t") if "\t" in line else line.split()
                if len(fields) < 3:
                    raise MSAValidationError(
                        f"{mask_bed}:{line_number}: expected at least 3 BED columns, got {len(fields)}"
                    )
                contig = fields[0]
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    raise MSAValidationError(
                        f"{mask_bed}:{line_number}: BED start and end must be integers"
                    )
                if contig not in offsets:
                    raise MSAValidationError(
                        f"{mask_bed}:{line_number}: contig '{contig}' is not in the reference"
                    )
                if start < 0 or end < start or end > contig_lengths[contig]:
                    raise MSAValidationError(
                        f"{mask_bed}:{line_number}: interval {start}-{end} is invalid for contig "
                        f"'{contig}' (length {contig_lengths[contig]})"
                    )
                if end > start:
                    intervals.append((offsets[contig] + start, offsets[contig] + end))

        intervals.sort()
        merged: list[tuple[int, int]] = []
        for start, end in intervals:
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        return merged

    @staticmethod
    def apply_mask(seq: str, intervals: list[tuple[int, int]], mask_char: str) -> str:
        """Replace each half-open interval in seq with mask_char."""
        if not intervals:
            return seq
        parts = []
        pos = 0
        for start, end in intervals:
            parts.append(seq[pos:start])
            parts.append(mask_char * (end - start))
            pos = end
        parts.append(seq[pos:])
        return "".join(parts)

    def _find_pseudo_fasta(self, snippy_dir: Path) -> Path:
        pseudo_fna_files = sorted(snippy_dir.glob("*.fna"))
        if not pseudo_fna_files:
            raise MissingInputError(f"No *.fna file found in {snippy_dir}")
        return pseudo_fna_files[0]

    def build_concatenated_alignment(
        self,
        snippy_dirs: list[Path],
        reference: Path,
        msa_out: Path,
    ) -> None:
        # Read contig order once
        ref_records = list(SeqIO.parse(str(reference), "fasta"))
        contig_order = [rec.id for rec in ref_records]

        if not contig_order:
            raise MSAValidationError(
                f"No contigs found in reference FASTA: {reference}"
            )

        contig_set = set(contig_order)

        # Prepare output MSA file
        msa_out.unlink(missing_ok=True)

        # Add reference as first record in MSA
        ref_str = "".join(str(rec.seq) for rec in ref_records)
        expected_length = len(ref_str)

        mask_intervals: list[tuple[int, int]] = []
        if self.mask_bed:
            mask_intervals = self.load_mask_intervals(
                self.mask_bed,
                contig_order,
                {rec.id: len(rec.seq) for rec in ref_records},
            )
            ref_str = self.apply_mask(ref_str, mask_intervals, self.mask_char)
        ref_record = SeqRecord(
            Seq(ref_str),
            id=self.reference_id,
            name=self.reference_id,
            description="",
        )
        with msa_out.open("w") as msa_handle:
            SeqIO.write([ref_record], msa_handle, "fasta-2line")

        # Process one sample at a time
        for snippy_dir in snippy_dirs:
            sample = snippy_dir.name
            fasta_path = self._find_pseudo_fasta(snippy_dir)

            # Load only THIS sample's contigs
            sample_map: Dict[str, SeqRecord] = {}
            for rec in SeqIO.parse(str(fasta_path), "fasta"):
                sample_map[rec.id] = rec

            # Validate contigs
            missing = [c for c in contig_order if c not in sample_map]
            if missing:
                raise MSAValidationError(
                    f"{sample}: missing contigs in {fasta_path}: "
                    + ", ".join(missing[:10])
                    + (" ..." if len(missing) > 10 else "")
                )

            extra = [cid for cid in sample_map.keys() if cid not in contig_set]
            if extra:
                raise MSAValidationError(
                    f"{sample}: extra contigs not in reference in {fasta_path}: "
                    + ", ".join(extra[:10])
                    + (" ..." if len(extra) > 10 else "")
                )

            # Concatenate contigs in order
            super_seq_str = "".join(str(sample_map[c].seq) for c in contig_order)

            # Validate alignment length consistency
            L = len(super_seq_str)
            if L != expected_length:
                raise MSAValidationError(
                    f"Alignment length mismatch: {sample} has {L}, expected {expected_length}"
                )

            super_seq_str = self.apply_mask(
                super_seq_str, mask_intervals, self.mask_char
            )

            super_record = SeqRecord(
                Seq(super_seq_str),
                id=sample,
                name=sample,
                description="",
            )

            # Append to final MSA
            with msa_out.open("a") as msa_handle:
                SeqIO.write([super_record], msa_handle, "fasta-2line")


class SoftCoreFilterOutput(BaseOutput):
    soft_core: Path = Field(
        ..., description="Filtered MSA containing only soft core positions"
    )
    constant_sites: Path = Field(
        ...,
        description="File containing constant-site counts for phylogenetic model correction",
    )


class DistleDistanceMatrixOutput(BaseOutput):
    phylip: Path = Field(..., description="Pairwise distance matrix in PHYLIP format")


class DistleDistanceMatrix(BaseStage):
    """Convert an alignment into a PHYLIP distance matrix using distle."""

    aln: Path = Field(..., description="Input multiple sequence alignment")
    format: str = Field(
        default="tabular", description="Output format for pairwise snip distance matrix"
    )

    _dependencies = [distle]

    @property
    def output(self) -> DistleDistanceMatrixOutput:
        suffix = ".distance.phylip" if self.format == "phylip" else ".distance.tsv"
        return DistleDistanceMatrixOutput(
            phylip=self.aln.with_suffix(suffix),
        )

    def create_commands(self, ctx):
        return [
            self.shell_cmd(
                [
                    "distle",
                    "--threads",
                    str(ctx.cpus),
                    "-o",
                    self.format,
                    str(self.aln),
                    str(self.output.phylip),
                ],
                description="Calculate pairwise distances from MSA in PHYLIP format",
            )
        ]


class SoftCoreError(StageExecutionError):
    pass


class SoftCoreFilter(BaseStage):
    """
    Filter a multiple sequence alignment to retain only positions present
    in a specified fraction of samples (soft core).
    """

    aln: Path = Field(..., description="Input multiple sequence alignment")
    core_threshold: float = Field(
        ..., description="Fraction of samples a position must be present in"
    )

    _dependencies = [core_snp_filter]

    @property
    def output(self) -> SoftCoreFilterOutput:
        aln = Path(f"{self.prefix}.{round(self.core_threshold * 100):03d}.aln")
        return SoftCoreFilterOutput(
            soft_core=aln, constant_sites=aln.with_suffix(".fconst")
        )

    def test_soft_core_is_not_empty(self):
        first_record = next(SeqIO.parse(str(self.output.soft_core), "fasta"), None)
        if first_record is None or len(first_record.seq) == 0:
            raise SoftCoreError(
                f"Soft core MSA has no sites: {self.output.soft_core}. One or more samples may have too few aligned or called bases for the requested core threshold. Check each sample's alignment percentage and the alignment-filter decisions."
            )

    def create_commands(self, ctx):
        return [
            self.shell_cmd(
                [
                    "coresnpfilter",
                    "-C",
                    str(self.aln),
                ],
                description="Extract constant sites from MSA",
                output_file=self.output.constant_sites,
            ),
            self.shell_cmd(
                [
                    "coresnpfilter",
                    "--core",
                    str(self.core_threshold),
                    "--exclude_invariant",
                    str(self.aln),
                ],
                description="Filter MSA to soft core positions",
                output_file=self.output.soft_core,
            ),
        ]
