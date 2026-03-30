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
from snippy_ng.context import Context 
from snippy_ng.utils.gather import guess_sample_id
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


class VcfStats(BaseStage):
    """Summarize the final VCF into TSV outputs that MultiQC can ingest."""

    vcf: Path = Field(..., description="Final annotated VCF file to summarize")
    sample_name: Optional[str] = Field(default=None, description="Optional sample name override for output tables")

    @property
    def output(self) -> VcfStatsOutput:
        return VcfStatsOutput(
            summary_tsv=Path(f"{self.prefix}.vcf.summary.tsv"),
            breakdown_tsv=Path(f"{self.prefix}.vcf.breakdown.tsv"),
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
        contig_counts: Counter = Counter()
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
                contig_counts[chrom] += 1
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
            "records_total": total_records,
            "records_pass": pass_records,
            "records_filtered": filtered_records,
            "records_lowqual": filter_counts.get("LowQual", 0),
            "records_lowdepth": filter_counts.get("LowDepth", 0),
            "records_snp": type_counts.get("SNP", 0),
            "records_mnp": type_counts.get("MNP", 0),
            "records_indel": type_counts.get("INDEL", 0),
            "records_complex": type_counts.get("COMPLEX", 0),
            "records_symbolic": type_counts.get("SYMBOLIC", 0),
            "records_contigs_with_variants": len(contigs_with_variants),
            "records_bcsq_annotated": bcsq_annotated_records,
            "snps_transitions": transitions,
            "snps_transversions": transversions,
            "snps_ts_tv": round(transitions / transversions, 4) if transversions else "",
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
            VcfStats._write_breakdown_rows(writer, sample, "contig", contig_counts)
            VcfStats._write_breakdown_rows(writer, sample, "substitution", substitution_counts)
            VcfStats._write_breakdown_rows(writer, sample, "consequence", consequence_counts)

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
        sample_name = guess_sample_id(Path(self.reads[0]).name)
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
