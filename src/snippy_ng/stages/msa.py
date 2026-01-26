from __future__ import annotations

from pathlib import Path
from typing import Dict

from pydantic import Field
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from snippy_ng.stages.base import BaseStage, BaseOutput
from snippy_ng.dependencies import biopython, core_snp_filter


class CombineFastaFileOutput(BaseOutput):
    aln: Path


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

    _dependencies = [biopython]

    @property
    def output(self) -> CombineFastaFileOutput:
        return CombineFastaFileOutput(aln=Path(f"{self.prefix}.full.aln"))

    @property
    def commands(self):
        return [
            self.python_cmd(
                func=self.build_concatenated_alignment,
                args=(self.snippy_dirs, self.reference, self.output.aln, self.prefix),
                description="Combine FASTAs into a concatenated MSA",
            )
        ]

    def _find_pseudo_fasta(self, snippy_dir: Path) -> Path:
        pseudo_fna_files = sorted(snippy_dir.glob("*.pseudo.fna"))
        if not pseudo_fna_files:
            raise FileNotFoundError(f"No *.pseudo.fna file found in {snippy_dir}")
        return pseudo_fna_files[0]

    def build_concatenated_alignment(
        self,
        snippy_dirs: list[Path],
        reference: Path,
        msa_out: Path,
        prefix: str,
    ) -> None:
        # -------------------------
        # 1) Read contig order once
        # -------------------------
        ref_records = list(SeqIO.parse(str(reference), "fasta"))
        contig_order = [rec.id for rec in ref_records]

        if not contig_order:
            raise ValueError(f"No contigs found in reference FASTA: {reference}")

        contig_set = set(contig_order)

        # -------------------------
        # 2) Prepare output MSA file
        # -------------------------
        msa_out.unlink(missing_ok=True)

        expected_length = None  # enforce alignment consistency

        # -------------------------
        # 3) Process one sample at a time
        # -------------------------
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
                raise ValueError(
                    f"{sample}: missing contigs in {fasta_path}: "
                    + ", ".join(missing[:10])
                    + (" ..." if len(missing) > 10 else "")
                )

            extra = [cid for cid in sample_map.keys() if cid not in contig_set]
            if extra:
                raise ValueError(
                    f"{sample}: extra contigs not in reference in {fasta_path}: "
                    + ", ".join(extra[:10])
                    + (" ..." if len(extra) > 10 else "")
                )

            # -------------------------
            # 4) Concatenate contigs in order
            # -------------------------
            super_seq_str = "".join(str(sample_map[c].seq) for c in contig_order)

            # Validate alignment length consistency
            L = len(super_seq_str)
            if expected_length is None:
                expected_length = L
            elif L != expected_length:
                raise ValueError(
                    f"Alignment length mismatch: {sample} has {L}, expected {expected_length}"
                )

            super_record = SeqRecord(
                Seq(super_seq_str),
                id=sample,
                name=sample,
                description="",
            )

            # -------------------------
            # 6) Append to final MSA (streaming)
            # -------------------------
            with msa_out.open("a") as msa_handle:
                SeqIO.write([super_record], msa_handle, "fasta-2line")


class SoftCoreFilterOutput(BaseOutput):
    aln: Path
    constant_sites: Path


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
        return SoftCoreFilterOutput(
            aln=Path(f"{self.prefix}.aln"),
            constant_sites=Path(f"{self.prefix}.aln.sites")
        )

    @property
    def commands(self):
        return [
            self.shell_pipeline(
                [
                    self.shell_cmd(
                        [
                            "coresnpfilter",
                            "-C",
                            str(self.aln),
                        ],
                        description="Extract constant sites from MSA",
                    )
                ],
                output_file=self.output.constant_sites,
                description="Write constant sites file",
            ), 
            self.shell_pipeline(
                [
                    self.shell_cmd(
                        [
                            "coresnpfilter",
                            "--core",
                            str(self.core_threshold),
                            "--exclude_invariant",
                            str(self.aln),
                        ],
                        description="Filter MSA to soft core positions",
                    )
                ],
                output_file=self.output.aln,
                description="Write filtered core alignment",
            )
        ]
