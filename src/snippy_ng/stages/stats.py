import csv
import gzip
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

from Bio import AlignIO
from snippy_ng.stages import BaseStage, ShellCommand, BaseOutput, TempPath
from snippy_ng.stages import PythonCommand
from snippy_ng.dependencies import seqkit
from snippy_ng.dependencies import biopython
from snippy_ng.dependencies import samtools
from snippy_ng.context import Context 
from snippy_ng.utils.gather import strip_bio_suffixes
from pydantic import Field, field_validator


def compute_aligned_percentage(seq: str) -> float:
    seq = seq.upper()

    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    gaps = seq.count("-")

    nongap_length = len(seq) - gaps
    denominator = gaps + nongap_length

    if denominator == 0:
        return 0.0

    return (a + c + g + t) * 100.0 / denominator


class AlignmentStatsOutput(BaseOutput):
    aligned_tsv: Path = Field(..., description="TSV file containing per-sequence aligned percentages from a multiple sequence alignment")


class VcfStatsOutput(BaseOutput):
    summary_tsv: Path = Field(..., description="One-row TSV summary of the final VCF for MultiQC ingestion")
    breakdown_tsv: Path = Field(..., description="Long-format TSV with VCF breakdown counts for MultiQC ingestion")


class AlignmentQcStatsOutput(BaseOutput):
    flagstat_txt: Path = Field(..., description="Raw samtools flagstat output")
    stats_txt: Path = Field(..., description="Raw samtools stats output")
    coverage_tsv: Path = Field(..., description="Raw samtools coverage output")
    summary_tsv: Path = Field(..., description="One-row TSV with normalized alignment QC metrics")


class FastaCompositionStatsOutput(BaseOutput):
    summary_tsv: Path = Field(..., description="One-row TSV with final FASTA composition metrics")


class SampleQcSummaryOutput(BaseOutput):
    qc_tsv: Path = Field(..., description="One-row combined sample QC TSV")


SAMPLE_QC_COLUMNS = [
    "sample",
    "pipeline_type",
    "read_files",
    "read_num_seqs",
    "read_sum_len",
    "read_min_len",
    "read_avg_len",
    "read_max_len",
    "alignment_total",
    "alignment_primary",
    "alignment_secondary",
    "alignment_supplementary",
    "alignment_duplicates",
    "alignment_mapped",
    "alignment_mapped_percent",
    "alignment_properly_paired",
    "alignment_properly_paired_percent",
    "alignment_reads_mapped",
    "alignment_reads_unmapped",
    "alignment_bases_mapped",
    "alignment_error_rate",
    "alignment_average_length",
    "alignment_average_quality",
    "alignment_insert_size_average",
    "alignment_mean_depth",
    "alignment_mean_breadth",
    "alignment_contigs",
    "vcf_total",
    "vcf_pass",
    "vcf_filtered",
    "vcf_lowqual",
    "vcf_lowdepth",
    "vcf_snp",
    "vcf_mnp",
    "vcf_indel",
    "vcf_complex",
    "vcf_symbolic",
    "vcf_contigs_with_variants",
    "vcf_bcsq_annotated",
    "vcf_transitions",
    "vcf_transversions",
    "vcf_ts_tv_ratio",
    "vcf_qual_mean",
    "vcf_depth_mean",
    "vcf_depth_min",
    "vcf_depth_max",
    "vcf_allele_fraction_mean",
    "final_length",
    "final_acgt",
    "final_A",
    "final_C",
    "final_G",
    "final_T",
    "final_N",
    "final_n",
    "final_gap",
    "final_other",
    "final_acgt_fraction",
    "final_N_fraction",
    "final_n_fraction",
    "final_gap_fraction",
    "core_aligned_percent",
]


class AlignmentAlignedPercentage(BaseStage):
    """Compute per-sequence aligned percentage from a multiple sequence alignment."""

    alignment: Path = Field(..., description="Input multiple sequence alignment file")
    alignment_format: str = Field("fasta", description="Alignment format understood by Biopython AlignIO")

    _dependencies = [biopython]

    @property
    def output(self) -> AlignmentStatsOutput:
        return AlignmentStatsOutput(
            aligned_tsv=Path(f"{self.prefix}.aligned.tsv")
        )

    def create_commands(self, ctx) -> List[PythonCommand]:
        return [
            self.python_cmd(
                func=self.write_aligned_percentages,
                args=[self.alignment, self.output.aligned_tsv, self.alignment_format],
                description="Compute per-sequence aligned percentages from MSA",
            )
        ]

    @staticmethod
    def write_aligned_percentages(
        alignment_path: Path,
        output_path: Path,
        alignment_format: str = "fasta",
    ) -> None:
        alignment = AlignIO.read(str(alignment_path), alignment_format)

        with output_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sequence", "aligned"])
            for record in alignment:
                writer.writerow(
                    [record.id, f"{compute_aligned_percentage(str(record.seq)):.2f}"]
                )


class SamtoolsAlignmentQcStats(BaseStage):
    """Generate normalized alignment QC metrics from a BAM/CRAM with samtools."""

    bam: Path = Field(..., description="Input BAM/CRAM alignment")
    reference: Optional[Path] = Field(default=None, description="Reference FASTA for CRAM input")
    sample_name: Optional[str] = Field(default=None, description="Optional sample name override for output tables")

    _dependencies = [samtools]

    @property
    def output(self) -> AlignmentQcStatsOutput:
        return AlignmentQcStatsOutput(
            flagstat_txt=Path(f"{self.prefix}.alignment.flagstat.txt"),
            stats_txt=Path(f"{self.prefix}.alignment.stats.txt"),
            coverage_tsv=Path(f"{self.prefix}.alignment.coverage.tsv"),
            summary_tsv=Path(f"{self.prefix}.alignment.tsv"),
        )

    def create_commands(self, ctx) -> List[ShellCommand | PythonCommand]:
        stats_cmd = ["samtools", "stats"]
        coverage_cmd = ["samtools", "coverage"]
        if self.reference is not None:
            stats_cmd.extend(["--reference", str(self.reference)])
        stats_cmd.append(str(self.bam))
        coverage_cmd.append(str(self.bam))
        return [
            self.shell_cmd(
                ["samtools", "flagstat", str(self.bam)],
                description=f"Generate flagstat QC metrics for {self.bam}",
                output_file=self.output.flagstat_txt,
            ),
            self.shell_cmd(
                stats_cmd,
                description=f"Generate samtools stats QC metrics for {self.bam}",
                output_file=self.output.stats_txt,
            ),
            self.shell_cmd(
                coverage_cmd,
                description=f"Generate samtools coverage QC metrics for {self.bam}",
                output_file=self.output.coverage_tsv,
            ),
            self.python_cmd(
                func=self.write_alignment_qc,
                args=[
                    self.output.flagstat_txt,
                    self.output.stats_txt,
                    self.output.coverage_tsv,
                    self.output.summary_tsv,
                    self.sample_name or self.prefix,
                ],
                description="Normalize alignment QC metrics",
            ),
        ]

    @staticmethod
    def _as_number(value: str) -> int | float | str:
        cleaned = value.strip().replace(",", "")
        if cleaned in {"", "."}:
            return ""
        try:
            number = float(cleaned)
        except ValueError:
            return value.strip()
        return int(number) if number.is_integer() else number

    @staticmethod
    def _parse_flagstat(path: Path) -> dict[str, int | float | str]:
        metrics: dict[str, int | float | str] = {}
        if not path.exists():
            return metrics

        for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            line = raw_line.strip()
            if not line:
                continue
            count_text = line.split(" + ", 1)[0].strip()
            try:
                count = int(count_text)
            except ValueError:
                continue
            parts = line.split(maxsplit=3)
            label = parts[3] if len(parts) > 3 else ""
            percent = ""
            if "(" in line and "%" in line:
                percent_text = line.split("(", 1)[1].split("%", 1)[0].strip()
                try:
                    percent = float(percent_text)
                except ValueError:
                    percent = ""

            if label.startswith("in total"):
                metrics["alignment_total"] = count
            elif label == "primary":
                metrics["alignment_primary"] = count
            elif label == "secondary":
                metrics["alignment_secondary"] = count
            elif label == "supplementary":
                metrics["alignment_supplementary"] = count
            elif label == "duplicates":
                metrics["alignment_duplicates"] = count
            elif label.startswith("mapped ("):
                metrics["alignment_mapped"] = count
                metrics["alignment_mapped_percent"] = percent
            elif label.startswith("properly paired"):
                metrics["alignment_properly_paired"] = count
                metrics["alignment_properly_paired_percent"] = percent
        return metrics

    @staticmethod
    def _parse_samtools_stats(path: Path) -> dict[str, int | float | str]:
        wanted = {
            "reads mapped": "alignment_reads_mapped",
            "reads unmapped": "alignment_reads_unmapped",
            "bases mapped": "alignment_bases_mapped",
            "error rate": "alignment_error_rate",
            "average length": "alignment_average_length",
            "average quality": "alignment_average_quality",
            "insert size average": "alignment_insert_size_average",
        }
        metrics: dict[str, int | float | str] = {}
        if not path.exists():
            return metrics

        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith("SN\t"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                label = parts[1].rstrip(":")
                key = wanted.get(label)
                if key is not None:
                    metrics[key] = SamtoolsAlignmentQcStats._as_number(parts[2])
        return metrics

    @staticmethod
    def _parse_samtools_coverage(path: Path) -> dict[str, int | float | str]:
        if not path.exists():
            return {}

        rows: list[dict[str, str]] = []
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                rows.append(row)
        if not rows:
            return {}

        total_length = 0.0
        weighted_meandepth = 0.0
        weighted_coverage = 0.0
        for row in rows:
            length = float(str(row.get("endpos", 0)).replace(",", "") or 0) - float(str(row.get("startpos", 1)).replace(",", "") or 1) + 1
            if length <= 0:
                continue
            total_length += length
            weighted_meandepth += float(str(row.get("meandepth", 0)).replace(",", "") or 0) * length
            weighted_coverage += float(str(row.get("coverage", 0)).replace(",", "") or 0) * length

        if total_length == 0:
            return {"alignment_contigs": len(rows)}
        return {
            "alignment_mean_depth": round(weighted_meandepth / total_length, 4),
            "alignment_mean_breadth": round(weighted_coverage / total_length, 4),
            "alignment_contigs": len(rows),
        }

    @staticmethod
    def write_alignment_qc(
        flagstat_txt: Path,
        stats_txt: Path,
        coverage_tsv: Path,
        output_tsv: Path,
        sample: str,
    ) -> None:
        row: dict[str, int | float | str] = {
            "sample": sample,
            **SamtoolsAlignmentQcStats._parse_flagstat(flagstat_txt),
            **SamtoolsAlignmentQcStats._parse_samtools_stats(stats_txt),
            **SamtoolsAlignmentQcStats._parse_samtools_coverage(coverage_tsv),
        }
        fieldnames = ["sample"] + [
            column
            for column in SAMPLE_QC_COLUMNS
            if column.startswith("alignment_")
        ]
        with output_tsv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerow({field: row.get(field, "") for field in fieldnames})


class FastaCompositionStats(BaseStage):
    """Summarize final pseudo-alignment FASTA composition for QC."""

    fasta: Path = Field(..., description="Final per-sample FASTA")
    sample_name: Optional[str] = Field(default=None, description="Optional sample name override for output tables")

    @property
    def output(self) -> FastaCompositionStatsOutput:
        return FastaCompositionStatsOutput(
            summary_tsv=Path(f"{self.prefix}.fasta.tsv")
        )

    def create_commands(self, ctx) -> List[PythonCommand]:
        return [
            self.python_cmd(
                func=self.write_fasta_composition,
                args=[
                    self.fasta,
                    self.output.summary_tsv,
                    self.sample_name or self.prefix,
                ],
                description="Summarize final FASTA composition",
            )
        ]

    @staticmethod
    def count_fasta_composition(fasta: Path) -> dict[str, int | float]:
        counts: Counter = Counter()
        total = 0
        with fasta.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    continue
                seq = line.strip()
                total += len(seq)
                counts.update(seq)

        upper_counts = {
            "A": counts.get("A", 0) + counts.get("a", 0),
            "C": counts.get("C", 0) + counts.get("c", 0),
            "G": counts.get("G", 0) + counts.get("g", 0),
            "T": counts.get("T", 0) + counts.get("t", 0),
        }
        n_upper = counts.get("N", 0)
        n_lower = counts.get("n", 0)
        gaps = counts.get("-", 0)
        acgt = sum(upper_counts.values())
        known = acgt + n_upper + n_lower + gaps
        denominator = total or 1
        return {
            "final_length": total,
            "final_acgt": acgt,
            "final_A": upper_counts["A"],
            "final_C": upper_counts["C"],
            "final_G": upper_counts["G"],
            "final_T": upper_counts["T"],
            "final_N": n_upper,
            "final_n": n_lower,
            "final_gap": gaps,
            "final_other": total - known,
            "final_acgt_fraction": round(acgt / denominator, 6),
            "final_N_fraction": round(n_upper / denominator, 6),
            "final_n_fraction": round(n_lower / denominator, 6),
            "final_gap_fraction": round(gaps / denominator, 6),
        }

    @staticmethod
    def write_fasta_composition(fasta: Path, output_tsv: Path, sample: str) -> None:
        row = {"sample": sample, **FastaCompositionStats.count_fasta_composition(fasta)}
        fieldnames = ["sample"] + [
            column
            for column in SAMPLE_QC_COLUMNS
            if column.startswith("final_")
        ]
        with output_tsv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerow({field: row.get(field, "") for field in fieldnames})


class VcfStats(BaseStage):
    """Summarize the final VCF into TSV outputs that MultiQC can ingest."""

    vcf: Path = Field(..., description="Final annotated VCF file to summarize")
    sample_name: Optional[str] = Field(default=None, description="Optional sample name override for output tables")

    @property
    def output(self) -> VcfStatsOutput:
        return VcfStatsOutput(
            summary_tsv=Path(f"{self.vcf.name}.summary.tsv"),
            breakdown_tsv=Path(f"{self.vcf.name}.breakdown.tsv"),
        )

    def create_commands(self, ctx) -> List[PythonCommand]:
        return [
            self.python_cmd(
                func=self.write_vcf_stats,
                args=[
                    self.vcf,
                    self.output.summary_tsv,
                    self.output.breakdown_tsv,
                    self.sample_name,
                ],
                description="Summarize final VCF into MultiQC-friendly TSV outputs",
            )
        ]

    @staticmethod
    def _open_text(path: Path):
        if path.suffix == ".gz":
            return gzip.open(path, "rt", encoding="utf-8", errors="replace")
        return path.open("r", encoding="utf-8", errors="replace")

    @staticmethod
    def _parse_info(info: str) -> Dict[str, str]:
        parsed: Dict[str, str] = {}
        if not info or info == ".":
            return parsed
        for field in info.split(";"):
            if not field:
                continue
            if "=" in field:
                key, value = field.split("=", 1)
                parsed[key] = value
            else:
                parsed[field] = ""
        return parsed

    @staticmethod
    def _parse_sample_map(fmt: str, sample: str) -> Dict[str, str]:
        if not fmt or not sample or fmt == "." or sample == ".":
            return {}
        return dict(zip(fmt.split(":"), sample.split(":")))

    @staticmethod
    def _as_float(value: Optional[str]) -> Optional[float]:
        if value in (None, "", "."):
            return None
        try:
            return float(value)
        except ValueError:
            return None

    @staticmethod
    def _extract_depth(sample_map: Dict[str, str], info_map: Dict[str, str]) -> Optional[float]:
        for key in ("DP",):
            value = VcfStats._as_float(sample_map.get(key))
            if value is not None:
                return value
        return VcfStats._as_float(info_map.get("DP"))

    @staticmethod
    def _extract_allele_fraction(sample_map: Dict[str, str], info_map: Dict[str, str], depth: Optional[float]) -> Optional[float]:
        af = VcfStats._as_float(sample_map.get("AF"))
        if af is not None:
            return af
        af = VcfStats._as_float(info_map.get("AF"))
        if af is not None:
            return af

        ao = sample_map.get("AO")
        if ao and depth not in (None, 0):
            try:
                alt_depths = [float(x) for x in ao.split(",") if x not in ("", ".")]
            except ValueError:
                alt_depths = []
            if alt_depths:
                return sum(alt_depths) / depth
        return None

    @staticmethod
    def _infer_type(ref: str, alt: str) -> str:
        if alt.startswith("<") and alt.endswith(">"):
            if alt == "<DEL>":
                return "DEL"
            return "SYMBOLIC"
        if len(ref) == len(alt):
            if len(ref) == 1:
                return "SNP"
            return "MNP"
        if len(ref) < len(alt):
            return "INS"
        if len(ref) > len(alt):
            return "DEL"
        return "COMPLEX"

    @staticmethod
    def _transition_or_transversion(ref: str, alt: str) -> Optional[str]:
        if len(ref) != 1 or len(alt) != 1:
            return None
        pair = f"{ref}>{alt}"
        if pair in {"A>G", "G>A", "C>T", "T>C"}:
            return "transition"
        if ref in "ACGT" and alt in "ACGT":
            return "transversion"
        return None

    @staticmethod
    def _write_breakdown_rows(
        writer: csv.writer,
        sample: str,
        section: str,
        counts: Counter,
    ) -> None:
        for key in sorted(counts):
            writer.writerow([sample, section, key, counts[key]])

    @staticmethod
    def write_vcf_stats(
        vcf_path: Path,
        summary_tsv: Path,
        breakdown_tsv: Path,
        sample_name: Optional[str] = None,
    ) -> None:
        sample = sample_name or vcf_path.stem
        filter_counts: Counter = Counter()
        type_counts: Counter = Counter()
        substitution_counts: Counter = Counter()
        consequence_counts: Counter = Counter()

        total_records = 0
        pass_records = 0
        filtered_records = 0
        bcsq_annotated_records = 0
        contigs_with_variants = set()

        qual_sum = 0.0
        qual_count = 0
        depth_sum = 0.0
        depth_count = 0
        min_depth = None
        max_depth = None
        af_sum = 0.0
        af_count = 0
        transitions = 0
        transversions = 0

        with VcfStats._open_text(vcf_path) as handle:
            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) > 9 and not sample_name:
                        sample = fields[9]
                    continue
                if not line.strip():
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue

                chrom, pos, _id, ref, alt_field, qual_field, filter_field, info_field = fields[:8]
                fmt = fields[8] if len(fields) > 8 else "."
                sample_field = fields[9] if len(fields) > 9 else "."

                total_records += 1
                contigs_with_variants.add(chrom)

                filter_labels = ["PASS"] if filter_field in (".", "PASS") else filter_field.split(";")
                if filter_labels == ["PASS"]:
                    pass_records += 1
                else:
                    filtered_records += 1
                for label in filter_labels:
                    filter_counts[label] += 1

                info_map = VcfStats._parse_info(info_field)
                sample_map = VcfStats._parse_sample_map(fmt, sample_field)

                qual = VcfStats._as_float(qual_field)
                if qual is not None:
                    qual_sum += qual
                    qual_count += 1

                depth = VcfStats._extract_depth(sample_map, info_map)
                if depth is not None:
                    depth_sum += depth
                    depth_count += 1
                    min_depth = depth if min_depth is None else min(min_depth, depth)
                    max_depth = depth if max_depth is None else max(max_depth, depth)

                allele_fraction = VcfStats._extract_allele_fraction(sample_map, info_map, depth)
                if allele_fraction is not None:
                    af_sum += allele_fraction
                    af_count += 1

                types = info_map.get("TYPE")
                type_list = [t for t in types.split(",") if t] if types else []

                alts = [a for a in alt_field.split(",") if a]
                if not type_list:
                    type_list = [VcfStats._infer_type(ref, alt) for alt in alts]
                for variant_type in type_list:
                    type_counts[variant_type] += 1

                if "BCSQ" in info_map:
                    bcsq_annotated_records += 1
                    for annotation in info_map["BCSQ"].split(","):
                        if annotation.startswith("@"):
                            continue
                        term_field = annotation.split("|", 1)[0]
                        for term in term_field.split("&"):
                            if term:
                                consequence_counts[term] += 1

                if len(alts) == 1:
                    change = VcfStats._transition_or_transversion(ref, alts[0])
                    if change is not None:
                        substitution = f"{ref}>{alts[0]}"
                        substitution_counts[substitution] += 1
                        if change == "transition":
                            transitions += 1
                        else:
                            transversions += 1

        summary_row = {
            "sample": sample,
            "total": total_records,
            "pass": pass_records,
            "filtered": filtered_records,
            "lowqual": filter_counts.get("LowQual", 0),
            "lowdepth": filter_counts.get("LowDepth", 0),
            "snp": type_counts.get("SNP", 0),
            "mnp": type_counts.get("MNP", 0),
            "indel": type_counts.get("INDEL", 0),
            "complex": type_counts.get("COMPLEX", 0),
            "symbolic": type_counts.get("SYMBOLIC", 0),
            "contigs_with_variants": len(contigs_with_variants),
            "bcsq_annotated": bcsq_annotated_records,
            "transitions": transitions,
            "transversions": transversions,
            "ts_tv_ratio": round(transitions / transversions, 4) if transversions else "",
            "qual_mean": round(qual_sum / qual_count, 4) if qual_count else "",
            "depth_mean": round(depth_sum / depth_count, 4) if depth_count else "",
            "depth_min": round(min_depth, 4) if min_depth is not None else "",
            "depth_max": round(max_depth, 4) if max_depth is not None else "",
            "allele_fraction_mean": round(af_sum / af_count, 4) if af_count else "",
        }

        with summary_tsv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summary_row.keys()), delimiter="\t")
            writer.writeheader()
            writer.writerow(summary_row)

        with breakdown_tsv.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "section", "label", "count"])
            VcfStats._write_breakdown_rows(writer, sample, "filter", filter_counts)
            VcfStats._write_breakdown_rows(writer, sample, "type", type_counts)
            VcfStats._write_breakdown_rows(writer, sample, "substitution", substitution_counts)
            VcfStats._write_breakdown_rows(writer, sample, "consequence", consequence_counts)


class SampleQcSummary(BaseStage):
    """Merge per-sample read, alignment, VCF, and final FASTA QC into one TSV."""

    sample_name: str = Field(..., description="Sample name for the QC row")
    pipeline_type: str = Field(..., description="Pipeline type: short, long, or asm")
    reads_tsv: Optional[Path] = Field(default=None, description="Optional SeqKit reads TSV")
    alignment_tsv: Optional[Path] = Field(default=None, description="Optional alignment QC TSV")
    vcf_summary_tsv: Optional[Path] = Field(default=None, description="Optional VCF summary TSV")
    fasta_tsv: Optional[Path] = Field(default=None, description="Optional final FASTA composition TSV")

    @property
    def output(self) -> SampleQcSummaryOutput:
        return SampleQcSummaryOutput(
            qc_tsv=Path(f"{self.prefix}.qc.tsv")
        )

    def create_commands(self, ctx) -> List[PythonCommand]:
        return [
            self.python_cmd(
                func=self.write_sample_qc_summary,
                args=[
                    self.output.qc_tsv,
                    self.sample_name,
                    self.pipeline_type,
                    self.reads_tsv,
                    self.alignment_tsv,
                    self.vcf_summary_tsv,
                    self.fasta_tsv,
                ],
                description="Merge sample QC metrics into one TSV",
            )
        ]

    @staticmethod
    def _read_first_row(path: Optional[Path]) -> dict[str, str]:
        if path is None or not path.exists():
            return {}
        with path.open("r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return next(reader, {}) or {}

    @staticmethod
    def _parse_numeric(value: str) -> Optional[float]:
        cleaned = str(value).strip().replace(",", "")
        if cleaned in {"", "."}:
            return None
        try:
            return float(cleaned)
        except ValueError:
            return None

    @staticmethod
    def _merge_read_stats(reads_tsv: Optional[Path]) -> dict[str, str | int | float]:
        if reads_tsv is None or not reads_tsv.exists():
            return {}

        rows: list[dict[str, str]] = []
        with reads_tsv.open("r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
        if not rows:
            return {}

        files = [row.get("file", "") for row in rows if row.get("file")]
        num_seqs_values = [SampleQcSummary._parse_numeric(row.get("num_seqs", "")) for row in rows]
        sum_len_values = [SampleQcSummary._parse_numeric(row.get("sum_len", "")) for row in rows]
        min_len_values = [SampleQcSummary._parse_numeric(row.get("min_len", "")) for row in rows]
        avg_len_values = [SampleQcSummary._parse_numeric(row.get("avg_len", "")) for row in rows]
        max_len_values = [SampleQcSummary._parse_numeric(row.get("max_len", "")) for row in rows]

        num_seqs = sum(value for value in num_seqs_values if value is not None)
        sum_len = sum(value for value in sum_len_values if value is not None)
        weighted_avg_denominator = sum(
            value
            for value in num_seqs_values
            if value is not None
        )
        weighted_avg_numerator = sum(
            (avg or 0) * (count or 0)
            for avg, count in zip(avg_len_values, num_seqs_values)
            if avg is not None and count is not None
        )
        min_lengths = [value for value in min_len_values if value is not None]
        max_lengths = [value for value in max_len_values if value is not None]

        return {
            "read_files": ",".join(files),
            "read_num_seqs": int(num_seqs) if num_seqs else "",
            "read_sum_len": int(sum_len) if sum_len else "",
            "read_min_len": int(min(min_lengths)) if min_lengths else "",
            "read_avg_len": round(weighted_avg_numerator / weighted_avg_denominator, 4) if weighted_avg_denominator else "",
            "read_max_len": int(max(max_lengths)) if max_lengths else "",
        }

    @staticmethod
    def write_sample_qc_summary(
        output_tsv: Path,
        sample: str,
        pipeline_type: str,
        reads_tsv: Optional[Path],
        alignment_tsv: Optional[Path],
        vcf_summary_tsv: Optional[Path],
        fasta_tsv: Optional[Path],
    ) -> None:
        row: dict[str, str | int | float] = {
            column: ""
            for column in SAMPLE_QC_COLUMNS
        }
        row["sample"] = sample
        row["pipeline_type"] = pipeline_type
        row.update(SampleQcSummary._merge_read_stats(reads_tsv))

        alignment_row = SampleQcSummary._read_first_row(alignment_tsv)
        row.update({key: value for key, value in alignment_row.items() if key in SAMPLE_QC_COLUMNS and key != "sample"})

        vcf_row = SampleQcSummary._read_first_row(vcf_summary_tsv)
        for key, value in vcf_row.items():
            target_key = f"vcf_{key}"
            if target_key in row:
                row[target_key] = value

        fasta_row = SampleQcSummary._read_first_row(fasta_tsv)
        row.update({key: value for key, value in fasta_row.items() if key in SAMPLE_QC_COLUMNS and key != "sample"})

        fieldnames = [
            column
            for column in SAMPLE_QC_COLUMNS
            if column in {"sample", "pipeline_type"} or row.get(column) not in ("", None)
        ]

        with output_tsv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerow({field: row.get(field, "") for field in fieldnames})


class SeqKitReadStatsOutput(BaseOutput):
    stats_tsv: Path = Field(..., description="Tab Separated Values (TSV) file containing read statistics")
    raw_tsv: TempPath = Field(..., description="Temporary raw TSV file produced by seqkit stats")


class SeqKitReadStats(BaseStage):
    """Generate read statistics using seqkit stats command.
    
    This stage analyzes FASTQ/FASTA files and produces a TSV file with 
    comprehensive statistics including sequence count, length metrics,
    quality scores, and GC content.
    
    Args:
        reads: List of input read files (FASTQ or FASTA).
        prefix: Output file prefix for the generated stats.tsv file.
        all_stats: Whether to output all statistics including quartiles and N50.
        tabular: Whether to output in machine-friendly tabular format.
        basename_only: Whether to only output basename of files in results.
        skip_errors: Whether to skip files with errors and show warnings.
        fastq_encoding: FASTQ quality encoding (sanger, solexa, illumina-1.3+, etc.).
        gap_letters: Gap letters to be counted in sequences.
        additional_options: Additional seqkit stats command-line options.
        cpus: Number of CPU cores to use for processing.
        tmpdir: Temporary directory for intermediate files.
        
    Returns:
        SeqKitReadStatsOutput: Output containing path to the generated stats.tsv file.
        
    Example:
        >>> from pathlib import Path
        >>> stage = SeqKitReadStats(
        ...     reads=["sample_R1.fastq", "sample_R2.fastq"],
        ...     prefix="sample_stats",
        ...     tmpdir=Path("/tmp"),
        ...     cpus=4
        ... )
        >>> print(stage.output.stats_tsv)
        sample_stats.reads.tsv
    """
    
    reads: List[Path] = Field(..., description="List of input read files (FASTQ or FASTA)")
    sample_name: Optional[str] = Field(default=None, description="Optional sample name override for output tables")
    all_stats: bool = Field(True, description="Output all statistics including quartiles and N50")
    tabular: bool = Field(True, description="Output in machine-friendly tabular format")
    basename_only: bool = Field(False, description="Only output basename of files")
    skip_errors: bool = Field(True, description="Skip files with errors and show warnings")
    fastq_encoding: str = Field("sanger", description="FASTQ quality encoding (sanger, solexa, illumina-1.3+, etc.)")
    gap_letters: str = Field("- .", description="Gap letters to be counted")
    additional_options: str = Field("", description="Additional seqkit stats options")
    
    _dependencies = [seqkit]
    
    @field_validator("reads")
    @classmethod
    def validate_reads(cls, v):
        """Validate that read files are provided and exist.
        
        Args:
            v: List of read file paths to validate.
            
        Returns:
            List[str]: Validated list of read file paths.
            
        Raises:
            ValueError: If no read files provided or if any file doesn't exist.
        """
        if not v or len(v) == 0:
            raise ValueError("At least one read file must be provided")
        return v
    
    @field_validator("fastq_encoding")
    @classmethod
    def validate_fastq_encoding(cls, v):
        """Validate FASTQ quality encoding format.
        
        Args:
            v: FASTQ encoding string to validate.
            
        Returns:
            str: Validated FASTQ encoding string.
            
        Raises:
            ValueError: If encoding is not in the list of supported formats.
        """
        valid_encodings = [
            "sanger", "solexa", "illumina-1.3+", "illumina-1.5+", "illumina-1.8+"
        ]
        if v not in valid_encodings:
            raise ValueError(f"Invalid FASTQ encoding. Must be one of: {', '.join(valid_encodings)}")
        return v
    
    @property
    def output(self) -> SeqKitReadStatsOutput:
        """Get the output specification for this stage.
        
        Returns:
            SeqKitReadStatsOutput: Object containing paths to output files.
        """
        return SeqKitReadStatsOutput(
            stats_tsv=Path(f"{self.prefix}.reads.tsv"),
            raw_tsv=Path(f"{self.prefix}.reads.raw.tsv"),
        )
    
    def build_seqkit_stats_command(self, ctx: Context) -> ShellCommand:
        """Constructs the seqkit stats command.
        
        Builds the complete seqkit stats command with all specified options
        including threading, output format, encoding, and file handling options.
        
        Returns:
            ShellCommand: Complete seqkit stats command ready for execution.
        """
        shell_cmd = self.shell_cmd(
            command=["seqkit", "stats"],
            description=f"Generate read statistics for {len(self.reads)} files using seqkit stats",
            output_file=Path(self.output.raw_tsv)
        )
        
        # Threading
        if ctx.cpus > 1:
            shell_cmd.command.extend(["-j", str(ctx.cpus)])
        
        # Output format options
        if self.tabular:
            shell_cmd.command.append("-T")
        
        if self.all_stats:
            shell_cmd.command.append("-a")
        
        if self.basename_only:
            shell_cmd.command.append("-b")
        
        if self.skip_errors:
            shell_cmd.command.append("-e")
        
        # FASTQ encoding
        if self.fastq_encoding != "sanger":
            shell_cmd.command.extend(["-E", self.fastq_encoding])
        
        # Gap letters
        if self.gap_letters != "- .":
            shell_cmd.command.extend(["-G", self.gap_letters])
        
        # Additional options (split if it contains spaces)
        if self.additional_options:
            import shlex
            shell_cmd.command.extend(shlex.split(self.additional_options))
        
        # Input files
        shell_cmd.command.extend(self.reads)
        
        # Create shell command with output file
        return shell_cmd
    
    @staticmethod
    def add_sample_column(input_tsv: Path, output_tsv: Path, sample_name: str) -> None:
        with input_tsv.open("r", newline="") as in_handle, output_tsv.open("w", newline="") as out_handle:
            reader = csv.reader(in_handle, delimiter="\t")
            writer = csv.writer(out_handle, delimiter="\t")

            header = next(reader, None)
            if header is None:
                writer.writerow(["sample"])
                return

            writer.writerow(["sample", *header])
            for row in reader:
                writer.writerow([sample_name, *row])

    def create_commands(self, ctx) -> List[ShellCommand | PythonCommand]:
        """Get the list of commands to execute for this stage.
        
        Returns:
            List[ShellCommand]: List containing the seqkit stats command.
        """
        sample_name = self.sample_name or strip_bio_suffixes(Path(self.reads[0]).name)
        return [
            self.build_seqkit_stats_command(ctx),
            self.python_cmd(
                func=self.add_sample_column,
                args=[self.output.raw_tsv, self.output.stats_tsv, sample_name],
                description=f"Add sample column to read statistics for {sample_name}",
            ),
        ]


class SeqKitReadStatsBasic(SeqKitReadStats):
    """Basic read statistics without extra information like quartiles and N50.
    
    This variant of SeqKitReadStats outputs only essential statistics
    without additional metrics like quartiles, N50, or other detailed measures.
    It's useful for quick quality checks when comprehensive statistics aren't needed.
    
    Args:
        all_stats: Set to False to output only basic statistics.
        
    Note:
        Inherits all other parameters from SeqKitReadStats.
        
    Example:
        >>> stage = SeqKitReadStatsBasic(
        ...     reads=["sample.fastq"],
        ...     prefix="basic_stats",
        ...     tmpdir=Path("/tmp")
        ... )
        >>> # Will generate command without -a flag for basic stats only
    """
    
    all_stats: bool = Field(False, description="Output basic statistics only")


class SeqKitReadStatsDetailed(SeqKitReadStats):
    """Detailed read statistics with all available metrics.
    
    This variant provides comprehensive statistics including all standard metrics
    plus additional N-statistics (like N75, N90, N95) that can be customized
    based on analysis requirements.
    
    Args:
        all_stats: Set to True to output all available statistics.
        additional_n_stats: List of integers (0-100) specifying additional 
            N-statistics to compute (e.g., [75, 90, 95] for N75, N90, N95).
            
    Note:
        Inherits all other parameters from SeqKitReadStats.
        
    Example:
        >>> stage = SeqKitReadStatsDetailed(
        ...     reads=["assembly.fasta"],
        ...     prefix="detailed_stats",
        ...     tmpdir=Path("/tmp"),
        ...     additional_n_stats=[75, 90, 95]
        ... )
        >>> # Will include N75, N90, N95 in addition to standard statistics
    """
    
    all_stats: bool = Field(True, description="Output all available statistics")
    additional_n_stats: List[int] = Field(
        default_factory=list, 
        description="Additional N-statistics to compute (e.g., [90, 95] for N90 and N95)"
    )
    
    @field_validator("additional_n_stats")
    @classmethod
    def validate_n_stats(cls, v):
        """Validate N-statistic values are within valid range.
        
        Args:
            v: List of N-statistic values to validate.
            
        Returns:
            List[int]: Validated list of N-statistic values.
            
        Raises:
            ValueError: If any N-statistic value is not between 0 and 100.
        """
        for stat in v:
            if not (0 <= stat <= 100):
                raise ValueError(f"N-statistic values must be between 0 and 100, got: {stat}")
        return v
    
    def build_seqkit_stats_command(self, ctx) -> ShellCommand:
        """Constructs the seqkit stats command with additional N-statistics.
        
        Builds the complete seqkit stats command including any additional
        N-statistics specified in the additional_n_stats parameter.
        
        Returns:
            ShellCommand: Complete seqkit stats command with N-statistics options.
        """
        shell_cmd = self.shell_cmd(
            command=["seqkit", "stats"],
            description=f"Generate detailed read statistics for {len(self.reads)} files using seqkit stats",
            output_file=Path(self.output.raw_tsv)
        ) 
        
        # Threading
        if ctx.cpus > 1:
            shell_cmd.command.extend(["-j", str(ctx.cpus)])
        
        # Output format options
        if self.tabular:
            shell_cmd.command.append("-T")
        
        if self.all_stats:
            shell_cmd.command.append("-a")
        
        # Additional N-statistics
        if self.additional_n_stats:
            n_stats_str = ",".join(map(str, self.additional_n_stats))
            shell_cmd.command.extend(["-N", n_stats_str])
        
        if self.basename_only:
            shell_cmd.command.append("-b")
        
        if self.skip_errors:
            shell_cmd.command.append("-e")
        
        # FASTQ encoding
        if self.fastq_encoding != "sanger":
            shell_cmd.command.extend(["-E", self.fastq_encoding])
        
        # Gap letters
        if self.gap_letters != "- .":
            shell_cmd.command.extend(["-G", self.gap_letters])
        
        # Additional options (split if it contains spaces)
        if self.additional_options:
            import shlex
            shell_cmd.command.extend(shlex.split(self.additional_options))
        
        # Input files
        shell_cmd.command.extend(self.reads)
        
        return shell_cmd
