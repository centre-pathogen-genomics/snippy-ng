from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict

import numpy as np
from pydantic import Field
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from snippy_ng.exceptions import StageExecutionError, MissingInputError
from snippy_ng.logging import logger
from snippy_ng.stages import BaseStage, BaseOutput
from snippy_ng.dependencies import biopython, core_snp_filter, numpy


class MSAValidationError(StageExecutionError):
    pass


class CombineFastaFileOutput(BaseOutput):
    aln: Path = Field(..., description="Combined multi-sample alignment in FASTA format")


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
    reference_id: str = Field(default="reference", description="ID to use for reference sequence in output MSA")

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

    def _find_pseudo_fasta(self, snippy_dir: Path) -> Path:
        pseudo_fna_files = sorted(snippy_dir.glob("*.pseudo.fna"))
        if not pseudo_fna_files:
            raise MissingInputError(f"No *.pseudo.fna file found in {snippy_dir}")
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
    soft_core: Path = Field(..., description="Filtered MSA containing only soft core positions")
    constant_sites: Path = Field(..., description="File containing constant-site counts for phylogenetic model correction")


class AlignmentSampleFilterOutput(BaseOutput):
    filtered_aln: Path = Field(..., description="MSA with low-alignment samples removed before soft core filtering")


class FilterAlignmentByAlignedPercentage(BaseStage):
    """
    Remove samples whose aligned percentage is unlikely to belong to the main sample cluster.
    """

    aln: Path = Field(..., description="Input multiple sequence alignment")
    alignment_stats: Path = Field(..., description="TSV file of per-sequence aligned percentages")
    inclusion_threshold: float = Field(0.20, description="Posterior probability threshold for retaining membership in the main alignment cluster")
    identical_alignment_spread: float = Field(10.0, description="If the spread of aligned percentages is less than or equal to this value, keep all samples to avoid filtering when there are no clear outliers")

    _dependencies = [biopython, numpy]

    @property
    def output(self) -> AlignmentSampleFilterOutput:
        return AlignmentSampleFilterOutput(
            filtered_aln=Path(f"{self.prefix}.filtered.aln"),
        )

    def create_commands(self, ctx):
        return [
            self.python_cmd(
                func=self.filter_alignment,
                args=(self.aln, self.alignment_stats, self.output.filtered_aln, self.inclusion_threshold, self.identical_alignment_spread),
                description="Filter low-alignment samples from the MSA before soft core filtering",
            )
        ]

    @classmethod
    def _fit_gmm(
        cls,
        values: list[float],
        max_iter: int = 100,
        tol: float = 1e-6,
        n_components: int = 2,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Fit a one-dimensional Gaussian mixture model with EM.

        The input is treated as a single numeric feature, and the algorithm
        estimates overlapping normal distributions that could have produced
        those values. Each iteration alternates between:
        1. estimating each sample's responsibility for each component, and
        2. updating the component weights, means, and variances from those
           responsibilities.

        This helper is used to separate the main alignment-quality cluster from
        one or two tails without imposing a hard cutoff on aligned percentage.
        """
        x = np.asarray(values, dtype=float).reshape(-1, 1)
        n = x.shape[0]
        if n_components < 2:
            raise ValueError("Need at least two components")
        if n < n_components:
            raise ValueError(f"Need at least {n_components} values to fit a {n_components}-component mixture")
        if np.var(x[:, 0]) < 1e-12:
            raise ValueError("Data are almost identical; mixture model is not identifiable")

        sorted_values = np.sort(x[:, 0])
        quantile_positions = np.linspace(0, n - 1, num=n_components)
        means = np.array([sorted_values[int(round(pos))] for pos in quantile_positions], dtype=float)
        if len(np.unique(np.round(means, 6))) < n_components:
            unique_values = np.unique(sorted_values)
            if unique_values.size < n_components:
                raise ValueError(f"Need at least {n_components} distinct values to fit a {n_components}-component mixture")
            means = np.linspace(unique_values[0], unique_values[-1], num=n_components, dtype=float)
        overall_var = float(np.var(x[:, 0])) if n > 1 else 1.0
        variances = np.full(n_components, max(overall_var, 1.0), dtype=float)
        weights = np.full(n_components, 1.0 / n_components, dtype=float)
        responsibilities = np.full((n, n_components), 1.0 / n_components, dtype=float)
        previous_ll = None

        for _ in range(max_iter):
            safe_variances = np.maximum(variances, 1e-6)
            scales = np.sqrt(2.0 * np.pi * safe_variances)
            exponent = -((x - means) ** 2) / (2.0 * safe_variances)
            pdf = np.exp(exponent) / scales
            weighted_pdf = pdf * weights
            totals = np.maximum(weighted_pdf.sum(axis=1, keepdims=True), 1e-12)
            responsibilities = weighted_pdf / totals
            log_likelihood = float(np.log(totals).sum())

            if previous_ll is not None and abs(log_likelihood - previous_ll) < tol:
                break
            previous_ll = log_likelihood

            component_mass = np.maximum(responsibilities.sum(axis=0), 1e-6)
            weights = component_mass / n
            means = (responsibilities * x).sum(axis=0) / component_mass
            variances = np.maximum((responsibilities * ((x - means) ** 2)).sum(axis=0) / component_mass, 1e-6)

        return weights, means, variances, responsibilities

    @staticmethod
    def _stats_fieldnames() -> list[str]:
        return [
            "sequence",
            "aligned",
            "probability_component_0",
            "probability_component_1",
            "probability_component_2",
            "probability_main",
            "removed",
        ]

    @classmethod
    def _empty_stats_row(cls, sequence: str, aligned: float) -> dict[str, str]:
        return {
            "sequence": sequence,
            "aligned": f"{aligned:.2f}",
            "probability_component_0": "",
            "probability_component_1": "",
            "probability_component_2": "",
            "probability_main": "",
            "removed": "false",
        }

    @staticmethod
    def _rewrite_alignment_stats(
        alignment_stats: Path,
        rows: list[dict[str, str]],
    ) -> None:
        with alignment_stats.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=FilterAlignmentByAlignedPercentage._stats_fieldnames(),
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rows)

    @classmethod
    def filter_alignment(
        cls,
        aln: Path,
        alignment_stats: Path,
        filtered_aln: Path,
        inclusion_threshold: float = 0.50,
        identical_alignment_spread: float = 10.0,
    ) -> None:
        records = list(SeqIO.parse(str(aln), "fasta"))
        if not records:
            raise MSAValidationError(f"No sequences found in alignment: {aln}")

        aligned_by_sequence: dict[str, float] = {}
        with alignment_stats.open("r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                aligned_by_sequence[row["sequence"]] = float(row["aligned"])

        reference_id = records[0].id
        sample_records = records[1:]
        stats_rows = [cls._empty_stats_row(reference_id, aligned_by_sequence[reference_id])]

        if len(sample_records) < 3:
            kept_records = records
        else:
            sample_values = []
            for record in sample_records:
                if record.id not in aligned_by_sequence:
                    raise MSAValidationError(f"Missing aligned percentage for sample '{record.id}' in {alignment_stats}")
                sample_values.append(aligned_by_sequence[record.id])

            if max(sample_values) - min(sample_values) <= identical_alignment_spread:
                kept_records = records
            else:
                n_components = 3 if len(sample_records) >= 5 else 2
                weights, _, _, responsibilities = cls._fit_gmm(sample_values, n_components=n_components)
                main_component = int(np.argmax(weights))

                keep_ids = {reference_id}
                removed_ids: list[str] = []

                for record, aligned, probs in zip(sample_records, sample_values, responsibilities):
                    probability_main = float(probs[main_component])
                    removed = probability_main < inclusion_threshold
                    row = cls._empty_stats_row(record.id, aligned)
                    for component_idx, probability in enumerate(probs):
                        row[f"probability_component_{component_idx}"] = f"{float(probability):.4f}"
                    row["probability_main"] = f"{probability_main:.4f}"
                    row["removed"] = str(removed).lower()
                    stats_rows.append(
                        row
                    )
                    if not removed:
                        keep_ids.add(record.id)
                    else:
                        removed_ids.append(record.id)

                kept_records = [record for record in records if record.id in keep_ids]
                if removed_ids:
                    logger.warning(
                        "Filtered out "
                        f"{len(removed_ids)} sample(s) from alignment {aln} "
                        f"using inclusion threshold {inclusion_threshold:.2f}: "
                        + ", ".join(removed_ids)
                    )

        if len(stats_rows) == 1:
            for record in sample_records:
                stats_rows.append(cls._empty_stats_row(record.id, aligned_by_sequence[record.id]))

        with filtered_aln.open("w") as handle:
            SeqIO.write(kept_records, handle, "fasta-2line")
        cls._rewrite_alignment_stats(alignment_stats, stats_rows)


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
        aln=Path(f"{self.prefix}.{round(self.core_threshold*100):03d}.aln")
        return SoftCoreFilterOutput(
            soft_core=aln,
            constant_sites=aln.with_suffix(".fconst")
        )

    def test_soft_core_is_not_empty(self):
        first_record = next(SeqIO.parse(str(self.output.soft_core), "fasta"), None)
        if first_record is None or len(first_record.seq) == 0:
            logger.warning(
                f"Soft core MSA has no sites: {self.output.soft_core}. You likely have samples with no variant sites. Check the % alignment for each sample to the reference."
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
