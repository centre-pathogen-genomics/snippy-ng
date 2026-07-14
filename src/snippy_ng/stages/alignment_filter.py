from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from Bio import SeqIO
from pydantic import Field
from sklearn.metrics import adjusted_mutual_info_score
from sklearn.mixture import GaussianMixture

from snippy_ng.dependencies import biopython, numpy, scikit_learn
from snippy_ng.exceptions import MSAValidationError
from snippy_ng.logging import logger
from snippy_ng.stages import BaseOutput, BaseStage


class AlignmentSampleFilterOutput(BaseOutput):
    filtered_aln: Path = Field(..., description="MSA with low-alignment samples removed before soft core filtering")
    filter_stats: Path = Field(..., description="Per-sample alignment filtering decisions and mixture-model diagnostics")


class AlignmentClusterTechnicalCheckOutput(BaseOutput):
    summary_tsv: Path = Field(..., description="Association between alignment mixture components and sample pipeline types")


class FilterAlignmentByAlignedPercentage(BaseStage):
    """Remove clear low-alignment outliers using a mixture model or small-cohort MAD rule."""

    aln: Path = Field(..., description="Input multiple sequence alignment")
    alignment_stats: Path = Field(..., description="TSV file of per-sequence aligned percentages")
    inclusion_threshold: float = Field(
        0.20,
        ge=0.0,
        le=1.0,
        description="Posterior probability threshold for membership in the retained alignment clusters",
    )
    minimum_cluster_separation: float = Field(
        10.0,
        ge=0.0,
        le=100.0,
        description="Minimum percentage-point separation between fitted component means required to filter samples",
    )
    minimum_mixture_samples: int = Field(
        10,
        ge=2,
        description="Minimum number of samples required to fit mixture models; smaller cohorts use a one-sided MAD rule",
    )
    minimum_samples_per_component: int = Field(5, ge=2, description="Minimum cohort size per candidate mixture component")
    bic_improvement: float = Field(6.0, ge=0.0, description="Minimum BIC improvement required to select a more complex mixture model")
    mad_threshold: float = Field(
        3.5,
        gt=0.0,
        description="One-sided robust z-score threshold for low-alignment outliers in cohorts too small for mixture modelling",
    )

    _dependencies = [biopython, numpy, scikit_learn]

    @property
    def output(self) -> AlignmentSampleFilterOutput:
        return AlignmentSampleFilterOutput(
            filtered_aln=Path(f"{self.prefix}.filtered.aln"),
            filter_stats=Path(f"{self.prefix}.alignment-filter.tsv"),
        )

    def create_commands(self, ctx):
        return [
            self.python_cmd(
                func=self.filter_alignment,
                args=(
                    self.aln,
                    self.alignment_stats,
                    self.output.filtered_aln,
                    self.output.filter_stats,
                    self.inclusion_threshold,
                    self.minimum_cluster_separation,
                    self.minimum_mixture_samples,
                    self.minimum_samples_per_component,
                    self.bic_improvement,
                    self.mad_threshold,
                ),
                description="Filter low-alignment samples from the MSA before soft core filtering",
            )
        ]

    @classmethod
    def _fit_gmm(
        cls,
        values: list[float],
        n_components: int = 2,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
        """Fit a deterministic one-dimensional Gaussian mixture."""
        x = np.asarray(values, dtype=float).reshape(-1, 1)
        n = x.shape[0]
        if n_components < 1:
            raise ValueError("Need at least one component")
        if n < n_components:
            raise ValueError(f"Need at least {n_components} values to fit a {n_components}-component mixture")
        if not np.all(np.isfinite(x)):
            raise ValueError("Mixture model input contains non-finite values")

        model = GaussianMixture(
            n_components=n_components,
            covariance_type="full",
            reg_covar=1e-6,
            n_init=10,
            random_state=0,
        ).fit(x)
        if not model.converged_:
            raise ValueError("Mixture model did not converge")

        weights = model.weights_
        means = model.means_.reshape(-1)
        variances = model.covariances_.reshape(-1)
        responsibilities = model.predict_proba(x)

        if not np.all(np.isfinite(responsibilities)):
            raise ValueError("Mixture model produced non-finite responsibilities")
        if not np.all(np.isfinite(means)):
            raise ValueError("Mixture model produced non-finite means")
        if not np.all(np.isfinite(weights)):
            raise ValueError("Mixture model produced non-finite weights")
        if not np.all(np.isfinite(variances)):
            raise ValueError("Mixture model produced non-finite variances")
        if np.any(weights < 0.05):
            raise ValueError("Mixture model contains a near-empty component")

        return weights, means, variances, responsibilities, float(model.bic(x))

    @classmethod
    def _select_gmm(
        cls,
        values: list[float],
        *,
        max_components: int = 3,
        minimum_samples_per_component: int = 5,
        bic_improvement: float = 6.0,
    ) -> tuple[
        int,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        dict[int, float],
    ]:
        """Select a one- to three-component mixture using guarded BIC improvement."""
        allowed_components = min(max_components, len(values) // minimum_samples_per_component)
        allowed_components = max(1, allowed_components)
        fits: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]] = {}
        bics: dict[int, float] = {}

        for n_components in range(1, allowed_components + 1):
            try:
                fit = cls._fit_gmm(values, n_components=n_components)
                bics[n_components] = fit[4]
                fits[n_components] = fit
            except ValueError as exc:
                if n_components == 1:
                    raise
                logger.warning(
                    f"Skipping {n_components}-component alignment mixture: {exc}"
                )

        selected_components = 1
        for n_components in sorted(component for component in fits if component > 1):
            if bics[n_components] <= bics[selected_components] - bic_improvement:
                selected_components = n_components

        weights, means, variances, responsibilities, _bic = fits[selected_components]
        return (
            selected_components,
            weights,
            means,
            variances,
            responsibilities,
            bics,
        )

    @staticmethod
    def _lower_mad_outliers(
        values: list[float],
        *,
        mad_threshold: float = 3.5,
        minimum_separation: float = 10.0,
    ) -> tuple[np.ndarray, float, float, float]:
        """Identify extreme lower-tail values using a conservative MAD cutoff."""
        x = np.asarray(values, dtype=float)
        if x.size == 0:
            return np.zeros(0, dtype=bool), float("nan"), float("nan"), float("nan")

        median = float(np.median(x))
        scaled_mad = float(1.4826 * np.median(np.abs(x - median)))
        required_deviation = max(minimum_separation, mad_threshold * scaled_mad)
        cutoff = median - required_deviation
        return x < cutoff, median, scaled_mad, cutoff

    @staticmethod
    def _stats_fieldnames() -> list[str]:
        return [
            "sequence",
            "aligned",
            "assigned_component",
            "probability_component_0",
            "probability_component_1",
            "probability_component_2",
            "probability_retained",
            "probability_main",
            "removed",
            "filter_reason",
            "selected_components",
            "retained_components",
            "largest_component_gap",
            "component_mean_0",
            "component_mean_1",
            "component_mean_2",
            "component_weight_0",
            "component_weight_1",
            "component_weight_2",
            "bic_1",
            "bic_2",
            "bic_3",
            "mad_median",
            "mad_scaled",
            "mad_cutoff",
        ]

    @classmethod
    def _empty_stats_row(
        cls,
        sequence: str,
        aligned: float,
        filter_reason: str,
        model_stats: dict[str, str] | None = None,
    ) -> dict[str, str]:
        row = {
            "sequence": sequence,
            "aligned": f"{aligned:.2f}",
            "assigned_component": "",
            "probability_component_0": "",
            "probability_component_1": "",
            "probability_component_2": "",
            "probability_retained": "",
            "probability_main": "",
            "removed": "false",
            "filter_reason": filter_reason,
            "selected_components": "",
            "retained_components": "",
            "largest_component_gap": "",
            "component_mean_0": "",
            "component_mean_1": "",
            "component_mean_2": "",
            "component_weight_0": "",
            "component_weight_1": "",
            "component_weight_2": "",
            "bic_1": "",
            "bic_2": "",
            "bic_3": "",
            "mad_median": "",
            "mad_scaled": "",
            "mad_cutoff": "",
        }
        if model_stats is not None:
            row.update(model_stats)
        return row

    @staticmethod
    def _aligned_percentage(
        sequence: str,
        aligned_by_sequence: dict[str, float],
        alignment_stats: Path,
    ) -> float:
        aligned = aligned_by_sequence.get(sequence)
        if aligned is not None:
            return aligned
        raise ValueError(
            f"Missing aligned percentage for sequence {sequence!r} in {alignment_stats}"
        )

    @staticmethod
    def _rewrite_alignment_stats(
        filter_stats: Path,
        rows: list[dict[str, str]],
    ) -> None:
        with filter_stats.open("w", newline="") as handle:
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
        filter_stats: Path | None = None,
        inclusion_threshold: float = 0.20,
        minimum_cluster_separation: float = 10.0,
        minimum_mixture_samples: int = 10,
        minimum_samples_per_component: int = 5,
        bic_improvement: float = 6.0,
        mad_threshold: float = 3.5,
    ) -> None:
        if filter_stats is None:
            filter_stats = filtered_aln.with_suffix(".stats.tsv")
        records = list(SeqIO.parse(str(aln), "fasta"))
        if not records:
            raise MSAValidationError(f"No sequences found in alignment: {aln}")
        record_ids = [record.id for record in records]
        if len(record_ids) != len(set(record_ids)):
            raise MSAValidationError(f"Duplicate sequence IDs found in alignment: {aln}")

        aligned_by_sequence: dict[str, float] = {}
        with alignment_stats.open("r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"Alignment statistics file has no header: {alignment_stats}")
            required_columns = {"sequence", "aligned"}
            missing = required_columns - set(reader.fieldnames)
            if missing:
                raise ValueError(
                    f"Alignment statistics file is missing columns: {sorted(missing)}"
                )
            for row in reader:
                sequence = row["sequence"]
                if sequence in aligned_by_sequence:
                    raise ValueError(f"Duplicate sequence {sequence!r} in {alignment_stats}")
                aligned = float(row["aligned"])
                if not np.isfinite(aligned) or not 0.0 <= aligned <= 100.0:
                    raise ValueError(
                        f"Invalid aligned percentage for {sequence!r}: {aligned}"
                    )
                aligned_by_sequence[sequence] = aligned

        reference_id = records[0].id
        sample_records = records[1:]
        reference_aligned = cls._aligned_percentage(
            reference_id,
            aligned_by_sequence,
            alignment_stats,
        )
        sample_values = [
            cls._aligned_percentage(record.id, aligned_by_sequence, alignment_stats)
            for record in sample_records
        ]
        stats_rows: list[dict[str, str]] = []
        model_stats = {
            "selected_components": "",
            "retained_components": "",
            "largest_component_gap": "",
            "component_mean_0": "",
            "component_mean_1": "",
            "component_mean_2": "",
            "component_weight_0": "",
            "component_weight_1": "",
            "component_weight_2": "",
            "bic_1": "",
            "bic_2": "",
            "bic_3": "",
            "mad_median": "",
            "mad_scaled": "",
            "mad_cutoff": "",
        }
        model_responsibilities: np.ndarray | None = None
        model_retained_components: np.ndarray | None = None
        skip_reason: str | None = None
        if not sample_records:
            kept_records = records
            skip_reason = "too_few_samples"
        elif len(sample_records) < minimum_mixture_samples:
            removed_by_mad, mad_median, mad_scaled, mad_cutoff = cls._lower_mad_outliers(
                sample_values,
                mad_threshold=mad_threshold,
                minimum_separation=minimum_cluster_separation,
            )
            model_stats["mad_median"] = f"{mad_median:.4f}"
            model_stats["mad_scaled"] = f"{mad_scaled:.4f}"
            model_stats["mad_cutoff"] = f"{mad_cutoff:.4f}"
            keep_ids = {reference_id}
            removed_ids: list[str] = []
            for record, aligned, removed in zip(sample_records, sample_values, removed_by_mad):
                row = cls._empty_stats_row(
                    record.id,
                    aligned,
                    "removed_by_mad" if removed else "retained_by_mad",
                    model_stats,
                )
                row["removed"] = str(bool(removed)).lower()
                stats_rows.append(row)
                if removed:
                    removed_ids.append(record.id)
                else:
                    keep_ids.add(record.id)

            kept_records = [record for record in records if record.id in keep_ids]
            if removed_ids:
                logger.warning(
                    "Filtered out "
                    f"{len(removed_ids)} sample(s) from alignment {aln} "
                    "using the small-cohort lower-tail MAD rule: "
                    + ", ".join(removed_ids)
                )
        elif max(sample_values) - min(sample_values) <= minimum_cluster_separation:
            kept_records = records
            skip_reason = "low_spread"
        else:
            (
                selected_components,
                weights,
                means,
                _variances,
                responsibilities,
                bics,
            ) = cls._select_gmm(
                sample_values,
                max_components=3,
                minimum_samples_per_component=minimum_samples_per_component,
                bic_improvement=bic_improvement,
            )
            model_responsibilities = responsibilities
            model_stats["selected_components"] = str(selected_components)
            for component_idx, (mean, weight) in enumerate(zip(means, weights)):
                model_stats[f"component_mean_{component_idx}"] = f"{float(mean):.4f}"
                model_stats[f"component_weight_{component_idx}"] = f"{float(weight):.4f}"
            for component_count, bic in bics.items():
                model_stats[f"bic_{component_count}"] = f"{bic:.4f}"

            if selected_components == 1:
                kept_records = records
                skip_reason = "single_component"
                model_retained_components = np.array([0])
                model_stats["retained_components"] = "0"
            else:
                component_order = np.argsort(means)
                ordered_gaps = np.diff(means[component_order])
                split_idx = int(np.argmax(ordered_gaps))
                largest_gap = float(ordered_gaps[split_idx])
                model_stats["largest_component_gap"] = f"{largest_gap:.4f}"
                retained_components = component_order[split_idx + 1:]
                model_stats["retained_components"] = ",".join(
                    str(int(component)) for component in retained_components
                )

                if largest_gap <= minimum_cluster_separation:
                    kept_records = records
                    skip_reason = "low_component_separation"
                    model_retained_components = component_order
                    model_stats["retained_components"] = ",".join(
                        str(int(component)) for component in component_order
                    )
                else:
                    model_retained_components = retained_components
                    keep_ids = {reference_id}
                    removed_ids: list[str] = []

                    for record, aligned, probs in zip(sample_records, sample_values, responsibilities):
                        probability_retained = float(probs[retained_components].sum())
                        removed = probability_retained < inclusion_threshold
                        filter_reason = "removed_by_gmm" if removed else "retained_by_gmm"
                        row = cls._empty_stats_row(
                            record.id,
                            aligned,
                            filter_reason,
                            model_stats,
                        )
                        for component_idx, probability in enumerate(probs):
                            row[f"probability_component_{component_idx}"] = f"{float(probability):.4f}"
                        row["assigned_component"] = str(int(np.argmax(probs)))
                        row["probability_retained"] = f"{probability_retained:.4f}"
                        row["probability_main"] = row["probability_retained"]
                        row["removed"] = str(removed).lower()
                        stats_rows.append(row)
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

        if skip_reason is not None:
            for sample_idx, record in enumerate(sample_records):
                row = cls._empty_stats_row(
                    record.id,
                    aligned_by_sequence[record.id],
                    skip_reason,
                    model_stats,
                )
                if model_responsibilities is not None and model_retained_components is not None:
                    probs = model_responsibilities[sample_idx]
                    for component_idx, probability in enumerate(probs):
                        row[f"probability_component_{component_idx}"] = f"{float(probability):.4f}"
                    row["assigned_component"] = str(int(np.argmax(probs)))
                    probability_retained = float(probs[model_retained_components].sum())
                    row["probability_retained"] = f"{probability_retained:.4f}"
                    row["probability_main"] = row["probability_retained"]
                stats_rows.append(row)

        stats_rows.insert(
            0,
            cls._empty_stats_row(
                reference_id,
                reference_aligned,
                "reference",
                model_stats,
            ),
        )

        with filtered_aln.open("w") as handle:
            SeqIO.write(kept_records, handle, "fasta-2line")
        cls._rewrite_alignment_stats(filter_stats, stats_rows)


class CheckAlignmentClustersByPipelineType(BaseStage):
    """Check whether alignment mixture components track technical input type."""

    filter_stats: Path = Field(..., description="Alignment filter statistics containing component assignments")
    qc_files: list[Path] = Field(..., description="Per-sample QC TSV files containing pipeline_type")
    warning_threshold: float = Field(0.8, ge=0.0, le=1.0, description="Adjusted mutual information threshold for warning")

    _dependencies = [scikit_learn]

    @property
    def output(self) -> AlignmentClusterTechnicalCheckOutput:
        return AlignmentClusterTechnicalCheckOutput(
            summary_tsv=Path(f"{self.prefix}.alignment-cluster-technical-check.tsv")
        )

    def create_commands(self, ctx):
        return [
            self.python_cmd(
                func=self.check_association,
                args=(self.filter_stats, self.qc_files, self.output.summary_tsv, self.warning_threshold),
                description="Check alignment mixture components for association with sample pipeline type",
            )
        ]

    @staticmethod
    def check_association(
        filter_stats: Path,
        qc_files: list[Path],
        output_tsv: Path,
        warning_threshold: float = 0.8,
    ) -> None:
        pipeline_types: dict[str, str] = {}
        for qc_file in qc_files:
            with qc_file.open("r", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                required = {"sample", "pipeline_type"}
                if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                    raise ValueError(f"QC TSV is missing sample or pipeline_type columns: {qc_file}")
                for row in reader:
                    sample = row["sample"].strip()
                    pipeline_type = row["pipeline_type"].strip()
                    if sample and pipeline_type:
                        if sample in pipeline_types:
                            raise ValueError(f"Duplicate sample {sample!r} across QC files")
                        pipeline_types[sample] = pipeline_type

        matched_types: list[str] = []
        assigned_types: list[str] = []
        assigned_components: list[str] = []
        with filter_stats.open("r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"sequence", "assigned_component"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValueError(f"Alignment filter statistics are missing required columns: {filter_stats}")
            for row in reader:
                pipeline_type = pipeline_types.get(row["sequence"])
                component = row["assigned_component"].strip()
                if pipeline_type:
                    matched_types.append(pipeline_type)
                if pipeline_type and component:
                    assigned_types.append(pipeline_type)
                    assigned_components.append(component)

        type_count = len(set(matched_types))
        component_count = len(set(assigned_components))
        score: float | None = None
        reason = ""
        if not matched_types:
            reason = "no_matched_samples"
        elif not assigned_components:
            reason = "no_component_assignments"
        elif len(assigned_components) < 2:
            reason = "too_few_assigned_samples"
        elif type_count < 2:
            reason = "single_pipeline_type"
        elif component_count < 2:
            reason = "single_component"
        else:
            score = float(adjusted_mutual_info_score(assigned_types, assigned_components))

        strong_association = score is not None and score >= warning_threshold
        if strong_association:
            logger.warning(
                "Alignment mixture components are strongly associated with sample pipeline type "
                f"(adjusted mutual information={score:.2f}); inspect technical input effects "
                "before interpreting clusters as biological"
            )

        fieldnames = [
            "matched_samples",
            "component_assigned_samples",
            "pipeline_type_count",
            "component_count",
            "adjusted_mutual_information",
            "strong_association",
            "reason",
        ]
        with output_tsv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerow({
                "matched_samples": len(matched_types),
                "component_assigned_samples": len(assigned_components),
                "pipeline_type_count": type_count,
                "component_count": component_count,
                "adjusted_mutual_information": "" if score is None else f"{score:.4f}",
                "strong_association": str(strong_association).lower(),
                "reason": reason,
            })
