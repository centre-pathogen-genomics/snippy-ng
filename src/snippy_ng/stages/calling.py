import json
import os
from bisect import bisect_right
from collections import defaultdict

# Concrete Alignment Strategies
from pathlib import Path
from typing import List, Annotated, Optional

from snippy_ng.stages import BaseStage, BaseOutput, TempPath
from snippy_ng.dependencies import freebayes, bcftools, bedtools, paftools, clair3, longbow, nucmer
from snippy_ng.exceptions import StageExecutionError
from snippy_ng.logging import logger
from snippy_ng.envvars import EnvVarField

from pydantic import Field, AfterValidator


MIN_SHORT_CHUNK_SIZE = 1000
MIN_LONG_CHUNK_SIZE = 10000


def estimate_reference_bases(reference: Path, reference_index: Path) -> int:
    """Estimate reference size in bases, preferring the FASTA index if available."""
    if reference_index.exists():
        total = 0
        with open(reference_index, "r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    continue
                try:
                    total += int(fields[1])
                except ValueError:
                    continue
        if total > 0:
            return total
    return max(1, reference.stat().st_size)


def get_calling_chunk_size(reference: Path, reference_index: Path, cpus: int, min_chunk_size: int) -> tuple[int, int]:
    """Choose a chunk size that oversamples slightly relative to the available CPUs."""
    refsize = estimate_reference_bases(reference, reference_index)
    num_chunks = 1 + 2 * (max(1, cpus) - 1)
    chunk_size = max(min_chunk_size, int(refsize / num_chunks))
    return num_chunks, chunk_size


def get_short_chunk_size(reference: Path, reference_index: Path, cpus: int) -> tuple[int, int]:
    """Determine short read chunk size based on reference size and available CPUs, with a minimum threshold."""
    return get_calling_chunk_size(reference, reference_index, cpus, MIN_SHORT_CHUNK_SIZE)

def get_long_chunk_size(reference: Path, reference_index: Path, cpus: int) -> tuple[int, int]:
    """Determine long read chunk size based on reference size and available CPUs, with a minimum threshold."""
    return get_calling_chunk_size(reference, reference_index, cpus, MIN_LONG_CHUNK_SIZE)

def no_spaces(v: str) -> str:
    """Ensure that a string contains no spaces."""
    if " " in v:
        raise ValueError(
            "Prefix must not contain spaces, please use underscores or hyphens instead."
        )
    return v


# Define the base Pydantic model for alignment parameters
class Caller(BaseStage):
    reference: Path = Field(
        ...,
        description="Reference file",
    )
    reference_index: Path = Field(..., description="Reference index file")
    prefix: Annotated[str, AfterValidator(no_spaces)] = Field(
        ..., description="Output file prefix"
    )
    additional_options: str = Field("", description="Additional options for the caller")

    def test_check_if_vcf_has_variants(self):
        """Test that the output VCF file is not empty."""
        vcf_path = self.output.vcf
        with open(vcf_path, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    return  # Found a non-header line, VCF is not empty
        # give warning instead of error as some callers may produce empty VCFs if no variants are found
        logger.warning(f"Output VCF file {vcf_path} has no variants (only header lines). Please check if this is expected based on your data and parameters.")

class BaseCallerOutput(BaseOutput):
    vcf: Path = Field(..., description="VCF file containing raw variant calls")

class FreebayesCallerOutput(BaseCallerOutput):
    vcf: Path = Field(..., description="VCF file containing raw variant calls")
    regions: TempPath = Field(..., description="BED file with regions used for parallel calling")


class FreebayesCaller(Caller):
    """
    Call variants using Freebayes.
    """

    bam: Path = Field(..., description="Input BAM file")
    bam_index: Path = Field(..., description="Index file for the input BAM")
    ploidy: int = Field(2, description="Ploidy for variant calling")
    min_mapping_quality: int = Field(30, description="Minimum mapping quality for FreeBayes to count reads")
    exclude_insertions: bool = Field(
        True,
        description="Exclude insertions from variant calls so the pseudo-alignment remains the same length as the reference",
    )

    _dependencies = [freebayes, bcftools]

    @property
    def output(self) -> FreebayesCallerOutput:
        return FreebayesCallerOutput(
            vcf=self.prefix + ".raw.vcf",
            regions=self.prefix + ".regions.txt",
        )
    
    def create_commands(self, ctx) -> List:
        """Constructs the Freebayes variant calling and postprocessing commands."""
        num_chunks, chunk_size = get_short_chunk_size(self.reference, self.reference_index, ctx.cpus)
        logger.info(
            f"Freebayes will process {num_chunks} chunks of {chunk_size} bp, {ctx.cpus} chunks at a time."
        )

        # 1) Regions for parallel FreeBayes
        generate_regions_cmd = self.shell_cmd(
            ["fasta_generate_regions.py", str(self.reference_index), str(chunk_size)],
            description="Generate genomic regions for parallel variant calling",
        )
        generate_regions_pipeline = self.shell_pipe(
            commands=[generate_regions_cmd],
            description="Generate regions file for parallel processing",
            output_file=Path(self.output.regions),
        )

        # 2) FreeBayes parallel call
        freebayes_cmd_parts = [
            "freebayes-parallel",
            str(self.output.regions),
            str(ctx.cpus),
            "--ploidy", str(self.ploidy),
            "--genotype-qualities",
            "--min-alternate-count", "2",
            "--min-repeat-entropy", "1.0",
            "--min-base-quality", "13",
            "--min-mapping-quality", str(self.min_mapping_quality),
            "--strict-vcf",
        ]
        if self.additional_options:
            import shlex
            freebayes_cmd_parts.extend(shlex.split(self.additional_options))
        
        freebayes_cmd_parts.extend(["-f", str(self.reference), str(self.bam)])

        freebayes_cmd = self.shell_cmd(
            freebayes_cmd_parts,
            description="Call variants with FreeBayes in parallel",
            output_file=Path(self.output.vcf),
        )
        return [generate_regions_pipeline, freebayes_cmd]


class FreebayesCallerLong(FreebayesCaller):
    """
    Call variants using Freebayes for long-read data.
    """
    # def test_freebayes_ran_successfully(self):
    min_mapping_quality: int = Field(30, description="Minimum mapping quality for FreeBayes to count reads")
        

    def create_commands(self, ctx) -> List:
        """Constructs the Freebayes variant calling and postprocessing commands."""
        num_chunks, chunk_size = get_long_chunk_size(self.reference, self.reference_index, ctx.cpus)
        logger.info(
            f"Freebayes will process {num_chunks} chunks of {chunk_size} bp, {ctx.cpus} chunks at a time."
        )

        # Regions for parallel FreeBayes
        generate_regions_cmd = self.shell_cmd(
            ["fasta_generate_regions.py", str(self.reference_index), str(chunk_size)],
            description="Generate genomic regions for parallel variant calling",
            output_file=Path(self.output.regions),
        )
        # FreeBayes parallel call
        freebayes_cmd_parts = [
            "freebayes-parallel",
            str(self.output.regions),
            str(ctx.cpus),
            "--ploidy", str(self.ploidy),
            "--genotype-qualities",
            "--haplotype-length", "-1",
            "--min-mapping-quality", str(self.min_mapping_quality),
            "--min-base-quality", "10",
        ]
        if self.additional_options:
            import shlex
            freebayes_cmd_parts.extend(shlex.split(self.additional_options))
        
        freebayes_cmd_parts.extend(["-f", str(self.reference), str(self.bam)])

        freebayes_cmd = self.shell_cmd(
            freebayes_cmd_parts,
            description="Call variants with FreeBayes in parallel",
            output_file=Path(self.output.vcf),
        )
        
        return [generate_regions_cmd, freebayes_cmd]

class PAFCallerOutput(BaseCallerOutput):
    vcf: Path = Field(..., description="VCF file with raw PAF-derived variant calls and annotations")
    tmp_vcf: TempPath = Field(..., description="Temporary intermediate VCF generated by paftools.js")
    missing_bed: Path = Field(..., description="BED file of unaligned (missing) reference regions")
    aln_bed: TempPath = Field(..., description="Temporary BED file of merged aligned reference intervals")


class PAFCaller(Caller):
    """
    Call variants from PAF alignments using paftools.js.
    """

    paf: Path = Field(..., description="Input PAF file")
    ref_dict: Path = Field(..., description="Reference FASTA dictionary file")
    min_mapping_quality: int = EnvVarField(30, "PAFTOOLS_MIN_MAPPING_QUAL", description="Mark PAF-derived variants below this mapping quality as LowQual")
    min_alignment_length_coverage: int = EnvVarField(10000, "PAFTOOLS_MIN_ALIGNMENT_LENGTH_COVERAGE", description="Minimum alignment length to compute coverage")
    min_alignment_length_variant_calling: int = EnvVarField(1000, "PAFFTOOLS_MIN_ALIGNMENT_LENGTH_VARIANT_CALLING", description="Minimum alignment length for variant calling")
    max_secondary_to_primary_ratio: float = EnvVarField(
        0.8,
        "PAFTOOLS_MAX_SECONDARY_TO_PRIMARY_RATIO",
        description="Mark variants as LowQual when minimap2 s2/s1 is at least this value. 0 disables this filter.",
    )
    max_gap_compressed_divergence: float = EnvVarField(
        0.05,
        "PAFTOOLS_MAX_GAP_COMPRESSED_DIVERGENCE",
        description="Mark variants as LowQual when minimap2 de is above this value. 0 disables this filter.",
    )
    max_sequence_divergence: float = EnvVarField(
        0.05,
        "PAFTOOLS_MAX_SEQUENCE_DIVERGENCE",
        description="Mark variants as LowQual when minimap2 dv is above this value. 0 disables this filter.",
    )
    max_repeat_query_fraction: float = EnvVarField(
        0.5,
        "PAFTOOLS_MAX_REPEAT_QUERY_FRACTION",
        description="Mark variants as LowQual when minimap2 rl/query-aligned-length is above this value. 0 disables this filter.",
    )
    min_chain_minimizers_per_kb: float = EnvVarField(
        0.0,
        "PAFTOOLS_MIN_CHAIN_MINIMIZERS_PER_KB",
        description="Mark variants as LowQual when minimap2 cm per aligned query kb is below this value. 0 disables this filter.",
    )

    _dependencies = [bedtools, bcftools, paftools]

    @property
    def output(self) -> PAFCallerOutput:
        return PAFCallerOutput(
            vcf=Path(f"{self.prefix}.raw.vcf"),
            aln_bed=Path(f"{self.prefix}.aln.bed"),
            missing_bed=Path(f"{self.prefix}.missing.bed"),
            tmp_vcf=Path(f"{self.prefix}.tmp.vcf"),
        )

    def create_commands(self, ctx) -> List:
        """Constructs the PAF processing and BED generation commands."""

        # 4) Convert PAF to merged aligned reference intervals (BED)
        # Keep primary or pseudo-primary hits: tp:A:P or tp:A:I
        paf_to_pipeline = self.shell_pipe(
            commands=[
                self.shell_cmd(
                    ["grep", "-E", "tp:A:[PI]", str(self.paf)],
                    description="Filter PAF for primary or pseudo-primary alignments",
                ),
                self.shell_cmd(
                    ["cut", "-f6,8,9"],
                    description="Extract relevant PAF fields (reference name, start, end)",
                ),
                self.shell_cmd(
                    ["bedtools", "sort", "-i", "-"],
                    description="Sort BED entries",
                ),
                self.shell_cmd(
                    ["bedtools", "merge", "-i", "-"],
                    description="Merge overlapping BED intervals",
                ),
            ],
            description="Convert PAF to merged aligned reference intervals (BED)",
            output_file=self.output.aln_bed,
        )

        # 5) Compute unaligned (missing) reference regions
        compute_missing_bed_cmd = self.shell_pipe(
            commands=[
                self.shell_cmd(
                    ["bedtools", "sort", "-g", str(self.ref_dict), "-i", str(self.output.aln_bed)],
                    description="Sort aligned reference intervals using reference dictionary order",
                ),
                self.shell_cmd(
                    ["bedtools", "complement", "-g", str(self.ref_dict), "-i", "-"],
                    description="Compute unaligned (missing) reference regions",
                ),
            ],
            description="Compute unaligned (missing) reference regions",
            output_file=self.output.missing_bed,
        )

        paftools_cmd_parts = [
            "paftools.js",
            "call",
            "-q",
            "0",
            "-L",
            str(self.min_alignment_length_coverage),
            "-l",
            str(self.min_alignment_length_variant_calling),
            "-s",
            self.prefix,
            "-f",
            str(self.reference),
            str(self.paf),
        ]
        if self.additional_options:
            import shlex
            paftools_cmd_parts.extend(shlex.split(self.additional_options))

        # variant calling
        paftools_cmd = self.shell_cmd(
                    paftools_cmd_parts,
                    description="Call variants from PAF using paftools.js",
                    output_file=self.output.tmp_vcf,
                )

        annotate_cmd = self.python_cmd(
            func=self.annotate_paf_vcf,
            args=[
                self.output.tmp_vcf,
                self.paf,
                self.output.vcf,
                self.min_mapping_quality,
                self.max_secondary_to_primary_ratio,
                self.max_gap_compressed_divergence,
                self.max_sequence_divergence,
                self.max_repeat_query_fraction,
                self.min_chain_minimizers_per_kb,
            ],
            description="Annotate PAF-derived variants with minimap2 mapping context",
        )

        return [
            paf_to_pipeline,
            compute_missing_bed_cmd,
            paftools_cmd,
            annotate_cmd,
        ]

    @staticmethod
    def _parse_paf_tags(fields: list[str]) -> dict[str, str]:
        tags: dict[str, str] = {}
        for field in fields[12:]:
            parts = field.split(":", 2)
            if len(parts) == 3:
                tags[parts[0]] = parts[2]
        return tags

    @staticmethod
    def _safe_int(value: str | None) -> int | None:
        if value is None:
            return None
        try:
            return int(value)
        except ValueError:
            return None

    @staticmethod
    def _safe_float(value: str | None) -> float | None:
        if value is None:
            return None
        try:
            return float(value)
        except ValueError:
            return None

    @classmethod
    def _read_paf_intervals(cls, paf: Path) -> dict[str, list[dict[str, object]]]:
        intervals: dict[str, list[dict[str, object]]] = defaultdict(list)
        with open(paf, "r", encoding="utf-8") as handle:
            for raw_line in handle:
                if not raw_line.strip():
                    continue
                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    continue
                tags = cls._parse_paf_tags(fields)
                query_start = int(fields[2])
                query_end = int(fields[3])
                target_start = int(fields[7])
                target_end = int(fields[8])
                mapq = int(fields[11])
                s1 = cls._safe_int(tags.get("s1"))
                intervals[fields[5]].append(
                    {
                        "query_name": fields[0],
                        "query_aligned_length": max(1, abs(query_end - query_start)),
                        "target_start": min(target_start, target_end),
                        "target_end": max(target_start, target_end),
                        "mapq": mapq,
                        "tp": tags.get("tp", ""),
                        "s1": s1,
                        "s2": cls._safe_int(tags.get("s2")),
                        "cm": cls._safe_int(tags.get("cm")),
                        "as": cls._safe_int(tags.get("AS")),
                        "dv": cls._safe_float(tags.get("dv")),
                        "de": cls._safe_float(tags.get("de")),
                        "rl": cls._safe_int(tags.get("rl")),
                    }
                )
        for chrom_intervals in intervals.values():
            chrom_intervals.sort(
                key=lambda item: (
                    int(item["target_start"]),
                    -int(item["target_end"]),
                    -int(item["mapq"]),
                    -(int(item["s1"]) if item["s1"] is not None else -1),
                )
            )
        return intervals

    @staticmethod
    def _best_paf_interval(intervals: list[dict[str, object]], pos0: int) -> dict[str, object] | None:
        candidates = [
            interval
            for interval in intervals
            if int(interval["target_start"]) <= pos0 < int(interval["target_end"])
        ]
        if not candidates:
            return None
        return max(
            candidates,
            key=lambda item: (
                int(item["mapq"]),
                int(item["s1"]) if item["s1"] is not None else -1,
                int(item["target_end"]) - int(item["target_start"]),
            ),
        )

    @staticmethod
    def _append_filter(filter_value: str, label: str) -> str:
        if filter_value in {"", ".", "PASS"}:
            return label
        labels = filter_value.split(";")
        if label not in labels:
            labels.append(label)
        return ";".join(labels)

    @staticmethod
    def _append_info(info_value: str, key: str, value: object) -> str:
        entry = f"{key}={value}"
        if not info_value or info_value == ".":
            return entry
        return f"{info_value};{entry}"

    @classmethod
    def _paf_lowqual_reasons(
        cls,
        interval: dict[str, object] | None,
        *,
        min_mapping_quality: int,
        max_secondary_to_primary_ratio: float,
        max_gap_compressed_divergence: float,
        max_sequence_divergence: float,
        max_repeat_query_fraction: float,
        min_chain_minimizers_per_kb: float,
    ) -> list[str]:
        if interval is None:
            return ["NO_PAF_CONTEXT"]

        reasons: list[str] = []
        mapq = int(interval["mapq"])
        if min_mapping_quality > 0 and mapq < min_mapping_quality:
            reasons.append("LOW_MAPQ")

        tp = str(interval["tp"])
        if tp and tp not in {"P", "I"}:
            reasons.append("NON_PRIMARY")

        s1 = interval["s1"]
        s2 = interval["s2"]
        if (
            max_secondary_to_primary_ratio > 0
            and s1 is not None
            and int(s1) > 0
            and s2 is not None
            and int(s2) / int(s1) >= max_secondary_to_primary_ratio
        ):
            reasons.append("REPETITIVE_CHAIN")

        de = interval["de"]
        if max_gap_compressed_divergence > 0 and de is not None and float(de) > max_gap_compressed_divergence:
            reasons.append("HIGH_DE")

        dv = interval["dv"]
        if max_sequence_divergence > 0 and dv is not None and float(dv) > max_sequence_divergence:
            reasons.append("HIGH_DV")

        aligned_query_length = int(interval["query_aligned_length"])
        rl = interval["rl"]
        if max_repeat_query_fraction > 0 and rl is not None and int(rl) / aligned_query_length > max_repeat_query_fraction:
            reasons.append("HIGH_RL")

        cm = interval["cm"]
        if min_chain_minimizers_per_kb > 0 and cm is not None:
            cm_per_kb = int(cm) / (aligned_query_length / 1000)
            if cm_per_kb < min_chain_minimizers_per_kb:
                reasons.append("LOW_CM_PER_KB")

        return reasons

    @classmethod
    def annotate_paf_vcf(
        cls,
        input_vcf: Path,
        paf: Path,
        output_vcf: Path,
        min_mapping_quality: int,
        max_secondary_to_primary_ratio: float,
        max_gap_compressed_divergence: float,
        max_sequence_divergence: float,
        max_repeat_query_fraction: float,
        min_chain_minimizers_per_kb: float,
    ) -> None:
        intervals_by_chrom = cls._read_paf_intervals(paf)
        with open(input_vcf, "r", encoding="utf-8") as in_handle, open(output_vcf, "w", encoding="utf-8") as out_handle:
            for raw_line in in_handle:
                if raw_line.startswith("##"):
                    out_handle.write(raw_line)
                    continue
                if raw_line.startswith("#CHROM"):
                    out_handle.write('##FILTER=<ID=LowQual,Description="Low-quality assembly variant call based on minimap2 PAF mapping context">\n')
                    out_handle.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
                    out_handle.write('##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count for each ALT">\n')
                    out_handle.write('##INFO=<ID=ASM_MAPQ,Number=1,Type=Integer,Description="minimap2 PAF mapping quality for the alignment covering this variant">\n')
                    out_handle.write('##INFO=<ID=ASM_TP,Number=1,Type=String,Description="minimap2 PAF tp tag for the alignment covering this variant">\n')
                    out_handle.write('##INFO=<ID=ASM_S1,Number=1,Type=Integer,Description="minimap2 PAF s1 chaining score">\n')
                    out_handle.write('##INFO=<ID=ASM_S2,Number=1,Type=Integer,Description="minimap2 PAF s2 secondary chaining score">\n')
                    out_handle.write('##INFO=<ID=ASM_S2_RATIO,Number=1,Type=Float,Description="minimap2 PAF s2/s1 chaining score ratio">\n')
                    out_handle.write('##INFO=<ID=ASM_CM,Number=1,Type=Integer,Description="minimap2 PAF cm minimizer count">\n')
                    out_handle.write('##INFO=<ID=ASM_CM_PER_KB,Number=1,Type=Float,Description="minimap2 PAF cm minimizer count per aligned query kb">\n')
                    out_handle.write('##INFO=<ID=ASM_AS,Number=1,Type=Integer,Description="minimap2 PAF AS alignment score">\n')
                    out_handle.write('##INFO=<ID=ASM_DV,Number=1,Type=Float,Description="minimap2 PAF dv approximate per-base sequence divergence">\n')
                    out_handle.write('##INFO=<ID=ASM_DE,Number=1,Type=Float,Description="minimap2 PAF de gap-compressed sequence divergence">\n')
                    out_handle.write('##INFO=<ID=ASM_RL,Number=1,Type=Integer,Description="minimap2 PAF rl length of query regions with repetitive seeds">\n')
                    out_handle.write('##INFO=<ID=ASM_RL_FRAC,Number=1,Type=Float,Description="minimap2 PAF rl divided by aligned query length">\n')
                    out_handle.write('##INFO=<ID=ASM_EDGE_DIST,Number=1,Type=Integer,Description="Distance in reference bases from this variant to the nearest PAF alignment edge">\n')
                    out_handle.write('##INFO=<ID=ASM_LOWQUAL_REASON,Number=.,Type=String,Description="Reasons this PAF-derived variant was marked LowQual">\n')
                    out_handle.write(raw_line)
                    continue
                if not raw_line.strip():
                    continue

                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    out_handle.write(raw_line)
                    continue

                chrom = fields[0]
                pos = int(fields[1])
                pos0 = pos - 1
                interval = cls._best_paf_interval(intervals_by_chrom.get(chrom, []), pos0)
                if interval is not None:
                    target_start = int(interval["target_start"])
                    target_end = int(interval["target_end"])
                    aligned_query_length = int(interval["query_aligned_length"])
                    fields[7] = cls._append_info(fields[7], "ASM_MAPQ", interval["mapq"])
                    if interval["tp"]:
                        fields[7] = cls._append_info(fields[7], "ASM_TP", interval["tp"])
                    for info_key, interval_key in (
                        ("ASM_S1", "s1"),
                        ("ASM_S2", "s2"),
                        ("ASM_CM", "cm"),
                        ("ASM_AS", "as"),
                        ("ASM_DV", "dv"),
                        ("ASM_DE", "de"),
                        ("ASM_RL", "rl"),
                    ):
                        value = interval[interval_key]
                        if value is not None:
                            fields[7] = cls._append_info(fields[7], info_key, value)
                    if interval["s1"] is not None and int(interval["s1"]) > 0 and interval["s2"] is not None:
                        fields[7] = cls._append_info(fields[7], "ASM_S2_RATIO", f"{int(interval['s2']) / int(interval['s1']):.6g}")
                    if interval["cm"] is not None:
                        fields[7] = cls._append_info(fields[7], "ASM_CM_PER_KB", f"{int(interval['cm']) / (aligned_query_length / 1000):.6g}")
                    if interval["rl"] is not None:
                        fields[7] = cls._append_info(fields[7], "ASM_RL_FRAC", f"{int(interval['rl']) / aligned_query_length:.6g}")
                    edge_distance = min(pos0 - target_start, target_end - pos0 - 1)
                    fields[7] = cls._append_info(fields[7], "ASM_EDGE_DIST", edge_distance)

                reasons = cls._paf_lowqual_reasons(
                    interval,
                    min_mapping_quality=min_mapping_quality,
                    max_secondary_to_primary_ratio=max_secondary_to_primary_ratio,
                    max_gap_compressed_divergence=max_gap_compressed_divergence,
                    max_sequence_divergence=max_sequence_divergence,
                    max_repeat_query_fraction=max_repeat_query_fraction,
                    min_chain_minimizers_per_kb=min_chain_minimizers_per_kb,
                )
                if reasons:
                    fields[6] = cls._append_filter(fields[6], "LowQual")
                    fields[7] = cls._append_info(fields[7], "ASM_LOWQUAL_REASON", ",".join(reasons))

                if len(fields) >= 10 and "DP" not in fields[8].split(":"):
                    fields[8] = f"{fields[8]}:DP:AO"
                    for sample_index in range(9, len(fields)):
                        fields[sample_index] = f"{fields[sample_index]}:1:1"

                out_handle.write("\t".join(fields) + "\n")


class ShowSnpsCallerOutput(BaseCallerOutput):
    vcf: Path = Field(..., description="VCF file with raw MUMmer-derived variant calls and annotations")
    tmp_vcf: TempPath = Field(..., description="Temporary intermediate VCF generated from the filtered delta")
    filtered_delta: TempPath = Field(..., description="Temporary filtered delta file used for variant calling")
    coords_tsv: TempPath = Field(..., description="Temporary show-coords table produced from the filtered delta")
    aln_bed: TempPath = Field(..., description="Temporary BED file of merged aligned reference intervals")
    missing_bed: Path = Field(..., description="BED file of unaligned (missing) reference regions")


class ShowSnpsCaller(Caller):
    """
    Call variants from nucmer alignments using a Python delta parser.
    """

    delta: Path = Field(..., description="Input delta file produced by nucmer")
    ref_dict: Path = Field(..., description="Reference FASTA dictionary file")
    assembly: Path = Field(..., description="Assembly FASTA used as the nucmer query")
    min_delta_identity: float = EnvVarField(
        96.05,
        "MUMMER_MIN_DELTA_IDENTITY",
        description="Minimum per-alignment delta identity percentage retained for MUMmer variant calling. 0 disables this filter.",
    )
    min_delta_block_length: int = EnvVarField(
        700,
        "MUMMER_MIN_DELTA_BLOCK_LENGTH",
        description="Minimum per-alignment delta block length retained for MUMmer variant calling. 0 disables this filter.",
    )
    max_delta_variants_per_kb: float = EnvVarField(
        37.5,
        "MUMMER_MAX_DELTA_VARIANTS_PER_KB",
        description="Maximum per-alignment variant density retained for MUMmer variant calling. 0 disables this filter.",
    )
    filter_ambiguous_delta_positions: bool = EnvVarField(
        True,
        "MUMMER_FILTER_AMBIGUOUS_DELTA_POSITIONS",
        description="Drop MUMmer variants at reference positions covered by multiple delta alignment blocks, similar to show-snps -C.",
    )

    _dependencies = [bedtools, bcftools, nucmer]

    @property
    def output(self) -> ShowSnpsCallerOutput:
        return ShowSnpsCallerOutput(
            vcf=Path(f"{self.prefix}.raw.vcf"),
            tmp_vcf=Path(f"{self.prefix}.mummer.tmp.vcf"),
            filtered_delta=Path(f"{self.prefix}.filtered.delta"),
            coords_tsv=Path(f"{self.prefix}.coords.tsv"),
            aln_bed=Path(f"{self.prefix}.aln.bed"),
            missing_bed=Path(f"{self.prefix}.missing.bed"),
        )

    def create_commands(self, ctx) -> List:
        filter_delta_cmd = self.shell_cmd(
            ["delta-filter", "-1", str(self.delta)],
            description="Filter nucmer alignments to the best one-to-one mappings",
            output_file=self.output.filtered_delta,
        )
        show_coords_cmd_parts = [
            "show-coords", "-T", "-H", "-r", "-c", "-l", str(self.output.filtered_delta)
            ]
        
        if self.additional_options:
            import shlex
            show_coords_cmd_parts.extend(shlex.split(self.additional_options))

        show_coords_cmd = self.shell_cmd(
            show_coords_cmd_parts,
            description="Summarize filtered nucmer alignments with show-coords",
            output_file=self.output.coords_tsv,
        )
        coords_to_bed_cmd = self.python_cmd(
            func=self.coords_to_bed,
            args=[self.output.coords_tsv, self.output.aln_bed],
            description="Convert show-coords output to merged BED intervals",
        )
        compute_missing_bed_cmd = self.shell_pipe(
            commands=[
                self.shell_cmd(
                    ["bedtools", "sort", "-g", str(self.ref_dict), "-i", str(self.output.aln_bed)],
                    description="Sort aligned reference intervals using reference dictionary order",
                ),
                self.shell_cmd(
                    ["bedtools", "complement", "-g", str(self.ref_dict), "-i", "-"],
                    description="Compute unaligned (missing) reference regions",
                ),
            ],
            description="Compute unaligned (missing) reference regions",
            output_file=self.output.missing_bed,
        )
        delta_to_vcf_cmd = self.python_cmd(
            func=self.delta_to_vcf,
            args=[
                self.output.filtered_delta,
                self.reference,
                self.assembly,
                self.output.tmp_vcf,
                self.min_delta_identity,
                self.min_delta_block_length,
                self.max_delta_variants_per_kb,
                self.filter_ambiguous_delta_positions,
            ],
            description="Convert filtered nucmer delta into a raw VCF",
        )
        finalize_vcf_cmd = self.python_cmd(
            func=self.postprocess_delta_vcf,
            args=[self.output.tmp_vcf, self.output.vcf, self.prefix],
            description="Add sample genotype fields to the MUMmer-derived VCF",
        )

        return [
            filter_delta_cmd,
            show_coords_cmd,
            coords_to_bed_cmd,
            compute_missing_bed_cmd,
            delta_to_vcf_cmd,
            finalize_vcf_cmd,
        ]

    @staticmethod
    def coords_to_bed(coords_tsv: Path, output_bed: Path) -> None:
        intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
        with open(coords_tsv, "r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 13:
                    continue
                start = min(int(fields[0]), int(fields[1])) - 1
                end = max(int(fields[0]), int(fields[1]))
                ref_name = fields[-2]
                intervals[ref_name].append((start, end))

        with open(output_bed, "w", encoding="utf-8") as handle:
            for ref_name in sorted(intervals):
                merged: list[list[int]] = []
                for start, end in sorted(intervals[ref_name]):
                    if not merged or start > merged[-1][1]:
                        merged.append([start, end])
                    else:
                        merged[-1][1] = max(merged[-1][1], end)
                for start, end in merged:
                    handle.write(f"{ref_name}\t{start}\t{end}\n")

    @staticmethod
    def _read_fasta_sequences(reference: Path) -> dict[str, str]:
        sequences: dict[str, list[str]] = {}
        current_name: str | None = None
        with open(reference, "r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    current_name = line[1:].split()[0]
                    sequences[current_name] = []
                    continue
                if current_name is None:
                    raise StageExecutionError(f"Invalid FASTA format in {reference}: sequence data found before a header")
                sequences[current_name].append(line.upper())
        return {name: "".join(parts) for name, parts in sequences.items()}

    @classmethod
    def _write_raw_vcf_record(cls, handle, chrom: str, pos: int, ref: str, alt: str, info: str = ".") -> None:
        handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t60\t.\t{info}\n")

    @classmethod
    def _reverse_complement(cls, sequence: str) -> str:
        return sequence[::-1].translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))

    @staticmethod
    def _delta_block_identity(errors: int, block_length: int) -> float:
        if block_length <= 0:
            return 0.0
        return 100.0 * (1.0 - (errors / block_length))

    @staticmethod
    def _delta_info(block_length: int, errors: int, identity: float, variants_per_kb: float) -> str:
        return (
            f"MUMMER_BLEN={block_length};"
            f"MUMMER_ERR={errors};"
            f"MUMMER_IDENTITY={identity:.5f};"
            f"MUMMER_VAR_KB={variants_per_kb:.5f}"
        )

    @staticmethod
    def _ambiguous_delta_intervals(
        block_outputs: list[dict[str, int | str | list[dict[str, int | str]]]]
    ) -> dict[str, list[tuple[int, int]]]:
        events_by_ref: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for block in block_outputs:
            ref_name = str(block["refname"])
            events_by_ref[ref_name].append((int(block["ref_start"]), 1))
            events_by_ref[ref_name].append((int(block["ref_end"]) + 1, -1))

        ambiguous: dict[str, list[tuple[int, int]]] = {}
        for ref_name, events in events_by_ref.items():
            depth = 0
            previous_position: int | None = None
            intervals: list[tuple[int, int]] = []
            for position, delta in sorted(events):
                if previous_position is not None and depth > 1 and previous_position < position:
                    intervals.append((previous_position, position - 1))
                depth += delta
                previous_position = position
            if intervals:
                ambiguous[ref_name] = intervals
        return ambiguous

    @staticmethod
    def _position_in_intervals(position: int, intervals: list[tuple[int, int]]) -> bool:
        starts = [start for start, _end in intervals]
        index = bisect_right(starts, position) - 1
        return index >= 0 and position <= intervals[index][1]

    @staticmethod
    def _record_info(record: dict[str, int | str]) -> str:
        parts = [str(record["info"])]
        if "mummer_edge_dist" in record:
            parts.append(f"MUMMER_EDGE_DIST={int(record['mummer_edge_dist'])}")
        return ";".join(part for part in parts if part and part != ".") or "."

    @classmethod
    def delta_to_vcf(
        cls,
        delta: Path,
        reference: Path,
        assembly: Path,
        output_vcf: Path,
        min_delta_identity: float = 0.0,
        min_delta_block_length: int = 0,
        max_delta_variants_per_kb: float = 0.0,
        filter_ambiguous_delta_positions: bool = True,
    ) -> None:
        reference_sequences = cls._read_fasta_sequences(reference)
        query_sequences = cls._read_fasta_sequences(assembly)
        query_sequences_rc = {
            name: cls._reverse_complement(sequence)
            for name, sequence in query_sequences.items()
        }
        block_outputs: list[dict[str, int | str | list[dict[str, int | str]]]] = []

        with open(delta, "r", encoding="utf-8") as handle:
            first_line = handle.readline().strip()
            if not first_line:
                raise StageExecutionError(f"Delta file {delta} is empty")
            header_fields = first_line.split()
            if len(header_fields) < 2:
                raise StageExecutionError(f"Delta file {delta} is missing the reference/query FASTA header")

            delta_type = handle.readline().strip()
            if delta_type != "NUCMER":
                raise StageExecutionError(f"Delta file {delta} is not a valid nucmer delta file")

            ref_name = None
            qry_name = None
            current_line = handle.readline()
            while current_line:
                line = current_line.strip()
                current_line = handle.readline()
                if not line:
                    continue

                fields = line.split()
                if len(fields) == 4 and line.startswith(">"):
                    ref_name = fields[0][1:]
                    qry_name = fields[1]
                    continue

                if len(fields) != 7:
                    continue
                if ref_name is None or qry_name is None:
                    raise StageExecutionError(f"Alignment block in {delta} appeared before a contig header")

                ref_start, ref_end, qry_start, qry_end = map(int, fields[:4])
                delta_errors = int(fields[4])
                alignment_block_length = max(abs(ref_end - ref_start) + 1, abs(qry_end - qry_start) + 1)
                block_identity = cls._delta_block_identity(delta_errors, alignment_block_length)
                delta_begins: list[int] = []
                delta_lengths: list[int] = []
                last_indel = 0

                while current_line:
                    delta_line = current_line.strip()
                    current_line = handle.readline()
                    if not delta_line:
                        continue
                    if delta_line == "0":
                        break
                    value = int(delta_line)
                    if value > 1 or value < -1:
                        if value > 0:
                            delta_begins.append(value + last_indel)
                            delta_lengths.append(1)
                            last_indel += value
                        else:
                            delta_begins.append(-(-value + last_indel))
                            delta_lengths.append(1)
                            last_indel += -value
                    else:
                        if not delta_lengths:
                            raise StageExecutionError(f"Malformed delta encoding in {delta}: extension without a preceding indel block")
                        delta_lengths[-1] += 1
                        last_indel += 1

                ref_seq = reference_sequences.get(ref_name)
                if ref_seq is None:
                    raise StageExecutionError(f"Reference contig {ref_name!r} from delta file was not found in the FASTA")

                if qry_start < qry_end:
                    qry_seq = query_sequences.get(qry_name)
                    if qry_seq is None:
                        raise StageExecutionError(f"Query contig {qry_name!r} from delta file was not found in the FASTA")
                else:
                    qry_seq = query_sequences_rc.get(qry_name)
                    if qry_seq is None:
                        raise StageExecutionError(f"Query contig {qry_name!r} from delta file was not found in the FASTA")
                    qry_start = len(qry_seq) - qry_start + 1
                    qry_end = len(qry_seq) - qry_end - 1

                ref_pos = ref_start
                qry_pos = qry_start
                alignment_coord = 1
                delta_index = 0
                block_records: list[dict[str, int | str]] = []

                while ref_pos <= ref_end:
                    at_indel = (
                        delta_index < len(delta_begins)
                        and alignment_coord == abs(delta_begins[delta_index])
                    )
                    if at_indel:
                        indel_length = delta_lengths[delta_index]
                        if delta_begins[delta_index] > 0:
                            for _ in range(indel_length):
                                ref_pos += 1
                                alignment_coord += 1
                            event_ref_pos = ref_pos - indel_length - 1
                            block_records.append(
                                {
                                    "refname": ref_name,
                                    "refpos": event_ref_pos,
                                    "ref": ref_seq[ref_pos - indel_length - 2:ref_pos - 1],
                                    "qry": qry_seq[qry_pos - 2:qry_pos - 1],
                                    "kind": "indel",
                                }
                            )
                        else:
                            for _ in range(indel_length):
                                qry_pos += 1
                                alignment_coord += 1
                            block_records.append(
                                {
                                    "refname": ref_name,
                                    "refpos": ref_pos - 1,
                                    "ref": ref_seq[ref_pos - 2:ref_pos - 1],
                                    "qry": qry_seq[qry_pos - indel_length - 2:qry_pos - 1],
                                    "kind": "indel",
                                }
                            )
                        delta_index += 1
                        continue

                    ref_base = ref_seq[ref_pos - 1]
                    qry_base = qry_seq[qry_pos - 1]
                    if ref_base.upper() != qry_base.upper():
                        block_records.append(
                            {
                                "refname": ref_name,
                                "refpos": ref_pos,
                                "ref": ref_base,
                                "qry": qry_base,
                                "kind": "snp",
                            }
                        )
                    ref_pos += 1
                    qry_pos += 1
                    alignment_coord += 1

                variants_per_kb = 0.0 if alignment_block_length <= 0 else (len(block_records) / alignment_block_length) * 1000.0
                if min_delta_identity > 0 and block_identity < min_delta_identity:
                    continue
                if min_delta_block_length > 0 and alignment_block_length < min_delta_block_length:
                    continue
                if max_delta_variants_per_kb > 0 and variants_per_kb > max_delta_variants_per_kb:
                    continue

                info = cls._delta_info(alignment_block_length, delta_errors, block_identity, variants_per_kb)
                for record in block_records:
                    record["info"] = info
                    record["mummer_edge_dist"] = min(
                        abs(int(record["refpos"]) - min(ref_start, ref_end)),
                        abs(max(ref_start, ref_end) - int(record["refpos"])),
                    )
                block_outputs.append(
                    {
                        "refname": ref_name,
                        "ref_start": min(ref_start, ref_end),
                        "ref_end": max(ref_start, ref_end),
                        "records": block_records,
                    }
                )

        ambiguous_intervals = (
            cls._ambiguous_delta_intervals(block_outputs)
            if filter_ambiguous_delta_positions
            else {}
        )
        out_records: list[dict[str, int | str]] = []
        for block in block_outputs:
            ref_name = str(block["refname"])
            intervals = ambiguous_intervals.get(ref_name, [])
            records = block["records"]
            if not isinstance(records, list):
                continue
            for record in records:
                if intervals and cls._position_in_intervals(int(record["refpos"]), intervals):
                    continue
                out_records.append(record)

        with open(output_vcf, "w", encoding="utf-8") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write("##source=snippy-ng-delta-parser\n")
            for contig, sequence in reference_sequences.items():
                handle.write(f"##contig=<ID={contig},length={len(sequence)}>\n")
            handle.write('##INFO=<ID=MUMMER_BLEN,Number=1,Type=Integer,Description="Length of the nucmer delta alignment block that produced this variant">\n')
            handle.write('##INFO=<ID=MUMMER_ERR,Number=1,Type=Integer,Description="MUMmer delta error count for the alignment block that produced this variant">\n')
            handle.write('##INFO=<ID=MUMMER_IDENTITY,Number=1,Type=Float,Description="Identity percentage estimated from the MUMmer delta alignment block error count">\n')
            handle.write('##INFO=<ID=MUMMER_VAR_KB,Number=1,Type=Float,Description="Variant records per kilobase emitted from this MUMmer delta alignment block">\n')
            handle.write('##INFO=<ID=MUMMER_EDGE_DIST,Number=1,Type=Integer,Description="Distance in reference bases from this variant to the nearest edge of its MUMmer delta alignment block">\n')
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for record in sorted(out_records, key=lambda item: (str(item["refname"]), int(item["refpos"]))):
                cls._write_raw_vcf_record(
                    handle,
                    str(record["refname"]),
                    int(record["refpos"]),
                    str(record["ref"]),
                    str(record["qry"]),
                    cls._record_info(record),
                )

    @staticmethod
    def postprocess_delta_vcf(input_vcf: Path, output_vcf: Path, sample_name: str) -> None:

        format_headers = [
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
            '##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count for each ALT">',
        ]

        with open(input_vcf, "r", encoding="utf-8") as src, open(output_vcf, "w", encoding="utf-8") as dst:
            saw_chrom_header = False
            for raw_line in src:
                if raw_line.startswith("##"):
                    dst.write(raw_line)
                    continue
                if raw_line.startswith("#CHROM"):
                    for header in format_headers:
                        dst.write(f"{header}\n")
                    dst.write(raw_line.rstrip("\n") + f"\tFORMAT\t{sample_name}\n")
                    saw_chrom_header = True
                    continue

                if not raw_line.strip():
                    continue

                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue
                if len(fields) > 8:
                    fields = fields[:8]
                if fields[5] in {"", "."}:
                    fields[5] = "60"
                fields.extend(["GT:DP:AO", "1/1:1:1"])
                dst.write("\t".join(fields) + "\n")

            if not saw_chrom_header:
                raise StageExecutionError(
                    f"delta2vcf output {input_vcf} did not contain a VCF header"
                )

class Clair3CallerOutput(BaseCallerOutput):
    vcf: Path = Field(..., description="VCF file containing raw variant calls produced by Clair3")


class Clair3ModelSelectorError(StageExecutionError):
    """Raised when the Clair3 model cannot be resolved based on Longbow predictions."""
    pass

class LongbowClair3ModelOutput(BaseOutput):
    prediction_json: TempPath = Field(..., description="Temporary Longbow prediction report")
    clair3_model: Path = Field(..., description="Text file containing the resolved Clair3 model path for the current run")


class LongbowClair3ModelSelector(BaseStage):
    reads: Path = Field(..., description="Input FASTQ file used for Longbow prediction")

    _dependencies = [longbow]

    @property
    def output(self) -> LongbowClair3ModelOutput:
        return LongbowClair3ModelOutput(
            prediction_json=Path(f"{self.prefix}.longbow.json"),
            clair3_model=Path(f"{self.prefix}.clair3_model.txt"),
        )

    def create_commands(self, ctx) -> List:
        logger.warning("Running Longbow to predict ONT basecalling configuration for Clair3 model selection. This may take a few minutes... Specify --clair3-model to skip this step and use a specific model.")
        return [
            self.shell_cmd(
                [
                    "longbow",
                    "-i", str(self.reads),
                    "-o", str(self.output.prediction_json),
                    "-t", str(ctx.cpus),
                ],
                description="Predict ONT basecalling configuration with Longbow",
            ),
            self.python_cmd(
                self.resolve_clair3_model,
                args=[self.output.prediction_json, self.output.clair3_model],
                description="Resolve the Clair3 model from Longbow output",
            ),
        ]

    @staticmethod
    def _normalize_string(value) -> Optional[str]:
        if value is None:
            return None
        if isinstance(value, (int, float)):
            value = str(value)
        if not isinstance(value, str):
            return None
        normalized = value.strip().lower()
        return normalized or None

    @classmethod
    def _get_nested_value(cls, data, *paths) -> Optional[str]:
        for path in paths:
            current = data
            for key in path:
                if not isinstance(current, dict):
                    current = None
                    break
                current = current.get(key)
            normalized = cls._normalize_string(current)
            if normalized is not None:
                return normalized
        return None

    @classmethod
    def _parse_longbow_prediction(cls, prediction_json: Path) -> dict[str, Optional[str]]:
        with open(prediction_json, "r", encoding="utf-8") as handle:
            data = json.load(handle)

        if not isinstance(data, dict):
            raise ValueError(f"Unexpected Longbow JSON format in '{prediction_json}'.")

        return {
            "basecaller": cls._get_nested_value(
                data,
                ("basecaller",),
                ("Software",),
                ("basecalling_software",),
                ("prediction", "basecaller"),
                ("result", "basecaller"),
            ),
            "flowcell": cls._get_nested_value(
                data,
                ("flowcell",),
                ("Flowcell",),
                ("flowcell_version",),
                ("prediction", "flowcell"),
                ("prediction", "flowcell_version"),
                ("result", "flowcell"),
            ),
            "major_version": cls._get_nested_value(
                data,
                ("major_version",),
                ("Version",),
                ("basecaller_version",),
                ("basecaller_major_version",),
                ("prediction", "major_version"),
                ("prediction", "basecaller_version"),
                ("result", "major_version"),
            ),
            "mode": cls._get_nested_value(
                data,
                ("mode",),
                ("Mode",),
                ("basecalling_mode",),
                ("prediction", "mode"),
                ("prediction", "basecalling_mode"),
                ("result", "mode"),
            ),
            "dorado_model_version": cls._get_nested_value(
                data,
                ("dorado_model_version",),
                ("prediction", "dorado_model_version"),
                ("result", "dorado_model_version"),
            ),
        }

    @staticmethod
    def _longbow_summary(prediction: dict[str, Optional[str]]) -> str:
        return ", ".join(
            [
                f"flowcell={prediction.get('flowcell') or 'unknown'}",
                f"basecaller={prediction.get('basecaller') or 'unknown'}",
                f"major_version={prediction.get('major_version') or 'unknown'}",
                f"mode={prediction.get('mode') or 'unknown'}",
            ]
        )

    @staticmethod
    def _candidate_roots() -> list[Path]:
        launch_dir = Path(os.environ.get("PWD", str(Path.cwd()))).expanduser()
        roots: list[Path] = []
        for env_var in ("CLAIR3_MODELS", "CLAIR3_MODEL_ROOT"):
            value = os.environ.get(env_var)
            if value:
                root = Path(value).expanduser()
                if not root.is_absolute():
                    root = launch_dir / root
                roots.append(root)

        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix:
            roots.extend(
                [
                    Path(conda_prefix) / "bin" / "models",
                    Path(conda_prefix) / "models",
                ]
            )

        roots.extend(
            [
                Path("/opt/models"),
                Path("/opt/models/clair3_models"),
            ]
        )

        unique_roots: list[Path] = []
        seen: set[Path] = set()
        for root in roots:
            root = root.expanduser()
            if root in seen:
                continue
            seen.add(root)
            unique_roots.append(root)
        return unique_roots

    @staticmethod
    def _container_roots() -> list[Path]:
        roots: list[Path] = []
        for env_var in ("CLAIR3_MODEL_CONTAINER_ROOT", "CLAIR3_MODELS_CONTAINER_ROOT"):
            value = os.environ.get(env_var)
            if value:
                roots.append(Path(value))
        return roots

    @staticmethod
    def _unique_preserving_order(items: list[str]) -> list[str]:
        seen: set[str] = set()
        unique: list[str] = []
        for item in items:
            if item in seen:
                continue
            seen.add(item)
            unique.append(item)
        return unique

    @staticmethod
    def _r10_model_series(mode: str) -> list[str]:
        return [
            f"r1041_e82_400bps_{mode}_v520",
            f"r1041_e82_400bps_{mode}_v500",
            f"r1041_e82_400bps_{mode}_v430",
            f"r1041_e82_400bps_{mode}_v410",
        ]

    @staticmethod
    def _r10_guppy56_model_series(mode: str) -> list[str]:
        if mode == "sup":
            return [
                "r1041_e82_400bps_sup_g615",
                "r1041_e82_260bps_sup_g632",
            ]
        if mode == "hac":
            return [
                "r1041_e82_400bps_hac_g632",
                "r1041_e82_400bps_hac_g615",
            ]
        if mode == "fast":
            return [
                "r1041_e82_400bps_fast_g632",
                "r1041_e82_400bps_fast_g615",
                "r1041_e82_260bps_fast_g632",
            ]
        return []

    @classmethod
    def _candidate_model_names(cls, prediction: dict[str, Optional[str]]) -> list[str]:
        basecaller = prediction.get("basecaller") or ""
        flowcell = prediction.get("flowcell") or ""
        major_version = prediction.get("major_version") or ""
        mode = prediction.get("mode") or ""
        is_guppy56 = "guppy" in basecaller and (
            "5or6" in major_version or "guppy5" in major_version or "guppy6" in major_version
        )

        if "r10" in flowcell:
            candidates: list[str] = []
            if "sup" in mode:
                candidates.extend(cls._r10_model_series("sup")[:2])
                if is_guppy56:
                    candidates.extend(cls._r10_guppy56_model_series("sup"))
                candidates.extend(cls._r10_model_series("sup")[2:])
            elif "hac" in mode:
                candidates.extend(cls._r10_model_series("hac")[:2])
                if is_guppy56:
                    candidates.extend(cls._r10_guppy56_model_series("hac"))
                candidates.extend(cls._r10_model_series("hac")[2:])
            elif "fast" in mode:
                candidates.extend(cls._r10_guppy56_model_series("fast"))
            else:
                candidates.extend(
                    [
                        "r1041_e82_400bps_sup_v520",
                        "r1041_e82_400bps_sup_v500",
                        "r1041_e82_400bps_hac_v520",
                        "r1041_e82_400bps_hac_v500",
                        "r1041_e82_400bps_sup_v410",
                        "r1041_e82_400bps_hac_v410",
                    ]
                )
            return cls._unique_preserving_order(candidates)

        if "r9" in flowcell and "guppy" in basecaller:
            if "guppy2" in major_version:
                return ["ont_guppy2", "r941_prom_hac_g238"]
            if "3or4" in major_version or "guppy3" in major_version or "guppy4" in major_version:
                return ["ont", "r941_prom_hac_g360+g422", "r941_prom_hac_g360+g422_1235"]
            if "5or6" in major_version or "guppy5" in major_version or "guppy6" in major_version:
                if "sup" in mode:
                    return ["ont_guppy5", "r941_prom_sup_g5014", "ont"]
                return ["ont_guppy5", "ont", "r941_prom_sup_g5014"]

        if "r9" in flowcell and ("dorado" in basecaller or "dorado0" in major_version or major_version == "0"):
            return ["ont_guppy5", "ont", "r941_prom_sup_g5014"]

        if "guppy2" in major_version:
            return ["ont_guppy2", "r941_prom_hac_g238"]
        if "guppy5" in major_version or "guppy6" in major_version or "dorado0" in major_version:
            return ["ont_guppy5", "r941_prom_sup_g5014", "ont"]
        if "guppy3" in major_version or "guppy4" in major_version:
            return ["ont", "r941_prom_hac_g360+g422", "r941_prom_hac_g360+g422_1235"]

        return ["ont_guppy5", "ont", "ont_guppy2", "r941_prom_sup_g5014", "r941_prom_hac_g360+g422", "r941_prom_hac_g238"]

    @staticmethod
    def _find_model_directory(root: Path, model_name: str) -> Optional[Path]:
        direct = root / model_name
        if direct.is_dir():
            return direct
        nested = root / "clair3_models" / model_name
        if nested.is_dir():
            return nested
        return None

    @classmethod
    def _prefer_v500_for_container(cls, candidates: list[str]) -> list[str]:
        def sort_key(model_name: str) -> tuple[int, int]:
            if "_v500" in model_name:
                return (0, 0)
            if "_v520" in model_name:
                return (0, 1)
            return (1, 0)

        head = sorted(
            [name for name in candidates if "_v500" in name or "_v520" in name],
            key=sort_key,
        )
        tail = [name for name in candidates if "_v500" not in name and "_v520" not in name]
        return cls._unique_preserving_order(head + tail)

    @classmethod
    def _resolve_local_model(cls, prediction: dict[str, Optional[str]]) -> Path:
        roots = [root for root in cls._candidate_roots() if root.exists()]
        candidates = cls._candidate_model_names(prediction)
        for model_name in candidates:
            for root in roots:
                resolved = cls._find_model_directory(root, model_name)
                if resolved is not None:
                    return resolved.resolve()

        container_roots = cls._container_roots()
        if container_roots:
            # TODO: the containers currently ship with only the v500 models
            # so we reorder the candidates to prefer the v500 models if 
            # container roots are being used. Remove once v520 is standard in containers.
            container_candidates = cls._prefer_v500_for_container(candidates)
            return container_roots[0] / container_candidates[0]

        if not roots:
            raise Clair3ModelSelectorError(
                "Could not find any Clair3 model roots. Set CLAIR3_MODELS for local models or CLAIR3_MODEL_CONTAINER_ROOT for container-only model paths."
            )

        raise Clair3ModelSelectorError(
            "Could not resolve a Clair3 model from Longbow prediction "
            f"({cls._longbow_summary(prediction)}). Searched roots: {', '.join(str(root) for root in roots)}. "
            f"Tried models: {', '.join(candidates)}. Provide --clair3-model to override."
        )

    @classmethod
    def resolve_clair3_model(cls, prediction_json: Path, output_path: Path) -> None:
        prediction = cls._parse_longbow_prediction(prediction_json)
        model_dir = cls._resolve_local_model(prediction)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(f"{model_dir}\n", encoding="utf-8")

        logger.info(
            f"Resolved Clair3 model '{model_dir}' from Longbow prediction ({cls._longbow_summary(prediction)})."
        )


class Clair3Caller(Caller):
    """
    Call variants using Clair3.
    """
    bam: Path = Field(..., description="Input BAM file")
    clair3_model: Path = Field(..., description="Path to Clair3 model")
    platform: str = Field("ont", description="Sequencing platform (e.g., ont, hifi)")
    min_mapping_quality: int = Field(30, description="Minimum mapping quality for Clair3 to count reads")

    _dependencies = [clair3]

    @property
    def output(self) -> Clair3CallerOutput:
        return Clair3CallerOutput(
            vcf=Path(f"{self.prefix}.raw.vcf"),
        )

    def _resolve_model_path_arg(self) -> Path:
        model_path = Path(self.clair3_model)
        if model_path.exists() and model_path.is_file():
            resolved = model_path.read_text(encoding="utf-8").strip()
            if not resolved:
                raise Clair3ModelSelectorError(f"Resolved Clair3 model file '{model_path}' is empty.")
            return Path(resolved)
        return model_path

    def create_commands(self, ctx) -> List:
        """Constructs the Clair3 variant calling commands."""
        model_path = self._resolve_model_path_arg()
        _, chunk_size = get_long_chunk_size(self.reference, self.reference_index, ctx.cpus)
        
        clair3_cmd_parts = [
            "run_clair3.sh",
            f"--model_path={model_path.absolute()}",
            f"--bam_fn={str(self.bam.absolute())}",
            f"--ref_fn={str(self.reference.absolute())}",
            f"--threads={str(ctx.cpus)}",
            f"--output={Path(self.prefix + '_clair3_out').absolute()}",
            f"--platform={self.platform}",
            f"--chunk_size={chunk_size}",
            f"--min_mq={self.min_mapping_quality}",
            "--include_all_ctgs",
            "--no_phasing_for_fa",
            "--enable_long_indel",
            "--enable_variant_calling_at_sequence_head_and_tail",
        ]
        if self.additional_options:
            import shlex
            clair3_cmd_parts.extend(shlex.split(self.additional_options))

        clair3_cmd = self.shell_cmd(
            clair3_cmd_parts,
            description="Call variants with Clair3",
        )
        unzip_cmd = self.shell_cmd(
            [
                "gunzip",
                f"{self.prefix}_clair3_out/merge_output.vcf.gz",
            ],
            description="Unzip Clair3 VCF output",
        )
        move_vcf_cmd = self.shell_cmd(
            [
                "mv",
                f"{self.prefix}_clair3_out/merge_output.vcf",
                str(self.output.vcf),
            ],
            description="Move Clair3 VCF to final output location",
        )
        commands = [clair3_cmd, unzip_cmd, move_vcf_cmd]

        if not ctx.no_cleanup:
            commands.append(self.shell_cmd(
                    [
                        "rm",
                        "-rf",
                        f"{self.prefix}_clair3_out",
                    ],
                    description="Clean up Clair3 output directory",
                )
            )

        return commands
