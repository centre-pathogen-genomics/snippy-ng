import gzip
import os
import re
import csv
import base64
import sys
from pathlib import Path
import datetime
import json
from typing import Callable, Dict, Iterable, List, Optional, Tuple, Union
from snippy_ng.stages import BaseStage, BaseOutput, TempPath
from snippy_ng.stages import PythonCommand
from snippy_ng.logging import logger
from snippy_ng.exceptions import PipelineExecutionError
from snippy_ng.dependencies import biopython, phylocanvas, phylojs, bedtools, samtools
from snippy_ng.utils.files import load_metadata_as_json_str
from snippy_ng.__about__ import __version__
from snippy_ng.constants import NEWICK_BRANCH_LENGTH_FORMAT
from pydantic import Field


class PrintVcfSummary(BaseStage):
    vcf_files: List[Path] = Field(..., description="Input VCF files to summarize")

    @property
    def output(self) -> None:
        return BaseOutput()

    def create_commands(self, ctx) -> List[PythonCommand]:
        return [
            self.python_cmd(
                func=self.summarize_vcfs,
                args=[self.vcf_files],
                description="Print summary of VCF files"
            )
        ]

    def summarize_vcfs(vcf_files: List[Path]):
        import gzip
        from collections import Counter

        def open_vcf(path):
            if str(path).endswith(".gz"):
                return gzip.open(path, "rt")
            return open(path, "r")

        for vcf_path in vcf_files:
            type_counts = Counter()
            total_records = 0

            with open_vcf(vcf_path) as f:
                for line in f:
                    if line.startswith("#"):
                        continue

                    total_records += 1
                    fields = line.rstrip("\n").split("\t")
                    info = fields[7]

                    # Extract TYPE field
                    info_fields = info.split(";")
                    type_field = next(
                        (x for x in info_fields if x.startswith("TYPE=")),
                        None
                    )

                    if type_field is None:
                        type_counts["MISSING"] += 1
                        continue

                    types = type_field.replace("TYPE=", "").split(",")
                    for t in types:
                        type_counts[t] += 1

            # Print like a groupby().count()
            print(f"\n{vcf_path}")
            print(f"TOTAL_RECORDS={total_records}")
            for t, c in type_counts.items():
                print(f"TYPE={t}\tCOUNT={c}")


class PrintVcfHistogram(BaseStage):
    vcf_path: Path = Field(..., description="Input VCF file to print a terminal genome histogram for")
    height: int = Field(5, description="Histogram height (rows)")
    contig_order: str = Field("header", description="Contig order: 'header' or 'alpha'")
    margin_cols: int = Field(1, description="Spaces between contig blocks")
    min_cols_per_contig: int = Field(1, description="Minimum columns allocated per contig")
    draw_separators: bool = Field(True, description="Draw contig separators in the histogram")

    @property
    def output(self) -> None:
        return BaseOutput()

    def create_commands(self, ctx) -> List[PythonCommand]:
        return [
            self.python_cmd(
                func=self.print_histograms,
                args=[
                    self.vcf_path,
                    self.height,
                    self.contig_order,
                    self.margin_cols,
                    self.min_cols_per_contig,
                    self.draw_separators,
                    not ctx.debug,
                ],
                description="Print terminal genome-wide histogram for VCF files",
            )
        ]

    @staticmethod
    def print_histograms(
        vcf_path: Path,
        height: int = 5,
        contig_order: str = "header",
        margin_cols: int = 1,
        min_cols_per_contig: int = 1,
        draw_separators: bool = True,
        pass_only: bool = True,
    ) -> None:

        # --- VCF parsing (contigs + records) ---
        CONTIG_LINE_RE = re.compile(r"^##contig=<(.+)>")
        ID_RE = re.compile(r"(?:^|,)ID=([^,>]+)")
        LEN_RE = re.compile(r"(?:^|,)length=([0-9]+)")

        def open_text_maybe_gzip(path: Union[str, Path]):
            path = Path(path)
            if path.suffix == ".gz":
                return gzip.open(path, "rt", encoding="utf-8", errors="replace")
            return open(path, "rt", encoding="utf-8", errors="replace")

        def parse_contigs_from_vcf_header(vcf_path: Union[str, Path]) -> List[Tuple[str, int]]:
            contigs: List[Tuple[str, int]] = []
            seen: set[str] = set()

            with open_text_maybe_gzip(vcf_path) as f:
                for line in f:
                    if line.startswith("#CHROM"):
                        break
                    m = CONTIG_LINE_RE.match(line)
                    if not m:
                        continue

                    inside = m.group(1)
                    id_m = ID_RE.search(inside)
                    len_m = LEN_RE.search(inside)
                    if not id_m or not len_m:
                        continue

                    cid = id_m.group(1)
                    clen = int(len_m.group(1))
                    if clen <= 0:
                        continue

                    if cid not in seen:
                        contigs.append((cid, clen))
                        seen.add(cid)

            if not contigs:
                raise ValueError("No valid ##contig header lines with ID and length found in VCF.")

            return contigs

        def iter_vcf_records(vcf_path: Union[str, Path], pass_only: bool = True, snp_only: bool = True) -> Iterable[Tuple[str, int]]:
            with open_text_maybe_gzip(vcf_path) as f:
                for line in f:
                    if not line or line[0] == "#":
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 7:
                        continue
                    if pass_only and fields[6] != "PASS":
                        continue
                    if snp_only and fields[7].find("TYPE=SNP") == -1:
                        continue
                    chrom = fields[0]
                    try:
                        pos = int(fields[1])
                    except ValueError:
                        continue
                    yield chrom, pos

        # --- Width allocation (terminal columns -> contigs) ---
        def allocate_columns_proportional(
            contigs: List[Tuple[str, int]],
            total_columns: int,
            min_cols_per_contig: int = 1,
        ) -> List[int]:
            if total_columns <= 0:
                raise ValueError("total_columns must be > 0")
            if min_cols_per_contig < 0:
                raise ValueError("min_cols_per_contig must be >= 0")

            n = len(contigs)
            if n == 0:
                return []

            if min_cols_per_contig * n > total_columns:
                min_cols_per_contig = 0

            lengths = [clen for _, clen in contigs]
            genome_len = sum(lengths)
            if genome_len <= 0:
                raise ValueError("Sum of contig lengths must be > 0")

            cols = [min_cols_per_contig] * n
            remaining = total_columns - sum(cols)
            if remaining <= 0:
                return cols

            raw = [length / genome_len * remaining for length in lengths]
            add = [int(x) for x in raw]
            cols = [c + a for c, a in zip(cols, add)]

            remainder = total_columns - sum(cols)
            fracs = [r - int(r) for r in raw]
            for i in sorted(range(n), key=lambda i: fracs[i], reverse=True):
                if remainder <= 0:
                    break
                cols[i] += 1
                remainder -= 1

            diff = total_columns - sum(cols)
            if diff != 0:
                step = 1 if diff > 0 else -1
                diff = abs(diff)
                for i in range(n):
                    if diff == 0:
                        break
                    j = i if step > 0 else (n - 1 - i)
                    if step < 0 and cols[j] == 0:
                        continue
                    cols[j] += step
                    diff -= 1

            return cols

        # --- Counting: variable bins per contig -> concatenated genome bins ---
        def genome_wide_binned_counts_variable_width(
            vcf_path: Union[str, Path],
            contigs: List[Tuple[str, int]],
            columns_per_contig: List[int],
            pass_only: bool = True,
        ) -> List[int]:
            if len(contigs) != len(columns_per_contig):
                raise ValueError("contigs and columns_per_contig must have the same length")

            lengths: Dict[str, int] = {cid: clen for cid, clen in contigs}

            bins: Dict[str, List[int]] = {}
            for (cid, _), w in zip(contigs, columns_per_contig):
                if w < 0:
                    raise ValueError("columns_per_contig contains a negative value")
                if w > 0:
                    bins[cid] = [0] * w

            for chrom, pos in iter_vcf_records(vcf_path, pass_only=pass_only):
                if chrom not in bins:
                    continue
                L = lengths.get(chrom)
                if L is None or L <= 0:
                    continue
                if pos < 1 or pos > L:
                    continue

                w = len(bins[chrom])
                idx = ((pos - 1) * w) // L
                if idx >= w:
                    idx = w - 1
                bins[chrom][idx] += 1

            genome_counts: List[int] = []
            for (cid, _), w in zip(contigs, columns_per_contig):
                if w > 0:
                    genome_counts.extend(bins[cid])

            return genome_counts

        # --- Rendering: vertical terminal histogram ---
        def print_vertical_histogram(
            counts: List[int],
            height: int = 12,
            char: str = "█",
            separators_after: Optional[List[int]] = None,
            top_padding: int = 1,
        ) -> None:
            if height <= 0:
                raise ValueError("height must be > 0")
            if not counts:
                raise ValueError("counts must be non-empty")

            max_count = max(counts) if max(counts) > 0 else 1
            scaled = []
            for c in counts:
                if c == 0:
                    scaled.append(0)
                    continue
                s = (c / max_count) * height
                # Ensure that non-zero counts get at least 1 row, to be visible.
                scaled.append(max(1, int(round(s))))

            sep_set = set(separators_after or [])

            def render_row(level: int) -> str:
                out: List[str] = []
                for i, h in enumerate(scaled):
                    out.append(char if h >= level else " ")
                    if i in sep_set:
                        out.append("│")
                return "".join(out)

            print("\n" * top_padding, end="")
            for level in range(height, 0, -1):
                print(render_row(level))

            axis: List[str] = []
            for i in range(len(counts)):
                axis.append("─")
                if i in sep_set:
                    axis.append("┼")
            print("".join(axis))

        def print_contig_labels(
            contigs: List[Tuple[str, int]],
            columns_per_contig: List[int],
        ) -> None:
            if len(contigs) != len(columns_per_contig):
                raise ValueError("contigs and columns_per_contig must have the same length")

            pieces: List[str] = []
            for idx, ((cid, _), w) in enumerate(zip(contigs, columns_per_contig)):
                if w <= 0:
                    continue
                block = [" "] * w
                label = cid[:w]
                start = max(0, (w - len(label)) // 2)
                for j, ch in enumerate(label):
                    block[start + j] = ch
                pieces.append("".join(block))
                if idx != len(contigs) - 1:
                    pieces.append(" ")
            print("".join(pieces))

        # --- End-to-end per VCF ---
        def terminal_genome_histogram_from_vcf(
            vcf_path: Union[str, Path],
            height: int = 5,
            contig_order: str = "header",
            margin_cols: int = 1,
            min_cols_per_contig: int = 1,
            draw_separators: bool = True,
            pass_only: bool = True,
        ) -> None:
            contigs = parse_contigs_from_vcf_header(vcf_path)
            if contig_order == "alpha":
                contigs = sorted(contigs, key=lambda x: x[0])
            elif contig_order != "header":
                raise ValueError("contig_order must be 'header' or 'alpha'")

            try:
                term_cols = os.get_terminal_size().columns
            except OSError:
                term_cols = 80

            separator_cols = (len(contigs) - 1) if draw_separators else 0
            label_spacing_cols = max(0, len(contigs) - 1)
            rendered_overhead_cols = max(separator_cols, label_spacing_cols)
            usable_cols = max(40, term_cols - margin_cols - rendered_overhead_cols)

            cols_per_contig = allocate_columns_proportional(
                contigs, usable_cols, min_cols_per_contig=min_cols_per_contig
            )

            counts = genome_wide_binned_counts_variable_width(
                vcf_path,
                contigs,
                cols_per_contig,
                pass_only=pass_only,
            )

            separators_after: Optional[List[int]] = None
            if draw_separators:
                separators_after = []
                acc = 0
                for w in cols_per_contig[:-1]:
                    if w > 0:
                        acc += w
                        separators_after.append(acc - 1)
                if separators_after and separators_after[-1] == len(counts) - 1:
                    separators_after.pop()
            if not sum(counts):
                filter_label = "PASS variants" if pass_only else "variants"
                logger.warning(f"No {filter_label} found in VCF {vcf_path}, skipping histogram.")
                return
            print_vertical_histogram(counts, height=height, separators_after=separators_after)
            print_contig_labels(contigs, cols_per_contig)


        terminal_genome_histogram_from_vcf(
            vcf_path=vcf_path,
            height=height,
            contig_order=contig_order,
            margin_cols=margin_cols,
            min_cols_per_contig=min_cols_per_contig,
            draw_separators=draw_separators,
            pass_only=pass_only,
        )

class FormatHTMLReportTemplateOutput(BaseOutput):
    rendered: Path = Field(..., description="Rendered HTML report file")

ContextValue = Optional[Union[str, int, float, Path, Callable]]
Context = Dict[str, ContextValue]


class SampleReportOutput(BaseOutput):
    rendered: Path = Field(..., description="Rendered per-sample HTML report file")
    variants_json: TempPath = Field(..., description="Temporary JSON table data extracted from the VCF")
    regions_bed: TempPath = Field(..., description="Temporary BED file of variant windows used to subset alignments")
    merged_regions_bed: TempPath = Field(..., description="Temporary merged BED file of variant windows")
    alignment_index: TempPath = Field(..., description="Temporary index for the input alignment")
    cram: TempPath = Field(..., description="Temporary clipped CRAM embedded into the sample report")
    crai: TempPath = Field(..., description="Temporary CRAI index for the windowed CRAM")


class SampleReport(BaseStage):
    """Create a standalone per-sample HTML report with variants and an IGV.js alignment view."""

    template_path: Path = Path(__file__).resolve().parent.parent / "templates" / "sample-report" / "snippy-sample-report.html"
    vcf: Path = Field(..., description="Input VCF used to build the variant table")
    alignment: Optional[Path] = Field(default=None, description="Optional input BAM/CRAM alignment to window around variants")
    reference: Optional[Path] = Field(default=None, description="Reference FASTA used by the alignment")
    reference_index: Optional[Path] = Field(default=None, description="Reference FASTA index (.fai)")
    title: str = Field(default="Snippy-NG Sample Report", description="Title for the sample report")
    sample_name: Optional[str] = Field(default=None, description="Optional sample name override")
    variant_scope: str = Field(default="pass", description="Variant scope to include: pass or all")
    window_size: int = Field(default=1000, description="Base pairs of context to embed around each variant")

    _dependencies = [bedtools, samtools]

    def model_post_init(self, __context) -> None:
        self._dependencies = [bedtools, samtools] if self.alignment else []

    @property
    def output(self) -> SampleReportOutput:
        return SampleReportOutput(
            rendered=Path(f"{self.prefix}.report.html"),
            variants_json=Path(f"{self.prefix}.report.variants.json"),
            regions_bed=Path(f"{self.prefix}.report.regions.bed"),
            merged_regions_bed=Path(f"{self.prefix}.report.regions.merged.bed"),
            alignment_index=Path(f"{self.prefix}.report.alignment.index"),
            cram=Path(f"{self.prefix}.report.cram"),
            crai=Path(f"{self.prefix}.report.cram.crai"),
        )

    def create_commands(self, ctx) -> List:
        if self.alignment and not self.reference:
            raise PipelineExecutionError("reference is required when alignment is provided for sample report")
        reference_index = self.reference_index or (Path(f"{self.reference}.fai") if self.reference else None)
        commands = [
            self.python_cmd(
                func=self.write_variant_assets,
                args=[
                    self.vcf,
                    self.output.variants_json,
                    self.output.regions_bed,
                    self.variant_scope,
                    self.window_size,
                    self.sample_name,
                ],
                description="Create sample-report variant table and alignment windows",
            ),
        ]
        has_igv_assets = self.alignment is not None and self.reference is not None and reference_index is not None
        if has_igv_assets:
            commands.extend(
                [
                    self.shell_pipe(
                        [
                            self.shell_cmd(
                                [
                                    "sort",
                                    "-k1,1",
                                    "-k2,2n",
                                    str(self.output.regions_bed),
                                ],
                                description="Sort sample-report alignment windows",
                            ),
                            self.shell_cmd(
                                [
                                    "bedtools",
                                    "merge",
                                    "-i",
                                    "-",
                                ],
                                description="Merge sample-report alignment windows",
                            ),
                        ],
                        description="Merge sample-report alignment windows",
                        output_file=self.output.merged_regions_bed,
                    ),
                    self.shell_cmd(
                        [
                            "samtools",
                            "index",
                            str(self.alignment),
                            str(self.output.alignment_index),
                        ],
                        description="Index input alignment for sample-report window extraction",
                    ),
                    self.shell_pipe(
                        [
                            self.shell_cmd(
                                [
                                    "samtools",
                                    "view",
                                    "--threads",
                                    str(ctx.cpus),
                                    "-h",
                                    "-M",
                                    "-X",
                                    "-T",
                                    str(self.reference),
                                    "-L",
                                    str(self.output.merged_regions_bed),
                                    str(self.alignment),
                                    str(self.output.alignment_index),
                                ],
                                description="Create SAM stream for sample-report clipping",
                            ),
                            self.shell_cmd(
                                [
                                    sys.executable,
                                    "-m",
                                    "snippy_ng",
                                    "utils",
                                    "aln",
                                    "samcrop",
                                    "--bed",
                                    str(self.output.merged_regions_bed),
                                ],
                                description="Hard-crop reads to sample-report windows",
                            ),
                            self.shell_cmd(
                                [
                                    "samtools",
                                    "sort",
                                    "--threads",
                                    str(ctx.cpus),
                                    "-T",
                                    f"{self.prefix}.sample-report.sort.tmp",
                                    "-O",
                                    "cram,level=9",
                                    "--reference",
                                    str(self.reference),
                                    "-o",
                                    str(self.output.cram),
                                    "-",
                                ],
                                description="Sort cropped sample-report CRAM",
                            ),
                        ],
                        description="Hard-crop and sort sample-report reads to clipped CRAM",
                    ),
                    self.shell_cmd(
                        ["samtools", "index", str(self.output.cram), str(self.output.crai)],
                        description="Index clipped sample-report CRAM",
                    ),
                ]
            )
        commands.append(
            self.python_cmd(
                func=self.render_sample_report,
                args=[
                    self.template_path,
                    self.output.rendered,
                    self.output.variants_json,
                    self.output.cram if has_igv_assets else None,
                    self.output.crai if has_igv_assets else None,
                    self.vcf,
                    self.reference if has_igv_assets else None,
                    reference_index if has_igv_assets else None,
                    self.title,
                    self.sample_name,
                ],
                description="Render sample HTML report",
            ),
        )
        return commands

    @staticmethod
    def _open_text(path: Path):
        if str(path).endswith(".gz"):
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
    def _parse_structured_header_value(value: str) -> Dict[str, str]:
        parsed: Dict[str, str] = {}
        current = []
        in_quotes = False
        for char in value:
            if char == '"':
                in_quotes = not in_quotes
            if char == "," and not in_quotes:
                field = "".join(current)
                if "=" in field:
                    key, field_value = field.split("=", 1)
                    parsed[key] = field_value.strip('"')
                current = []
                continue
            current.append(char)
        field = "".join(current)
        if "=" in field:
            key, field_value = field.split("=", 1)
            parsed[key] = field_value.strip('"')
        return parsed

    @staticmethod
    def _parse_format_header(line: str) -> Optional[Dict[str, str]]:
        match = re.match(r"^##FORMAT=<(.+)>", line.rstrip("\n"))
        if not match:
            return None
        header = SampleReport._parse_structured_header_value(match.group(1))
        field_id = header.get("ID")
        if not field_id:
            return None
        return {
            "id": field_id,
            "number": header.get("Number", "."),
            "type": header.get("Type", "String"),
            "description": header.get("Description", ""),
        }

    @staticmethod
    def _parse_info_header(line: str) -> Optional[Dict[str, str]]:
        match = re.match(r"^##INFO=<(.+)>", line.rstrip("\n"))
        if not match:
            return None
        header = SampleReport._parse_structured_header_value(match.group(1))
        field_id = header.get("ID")
        if not field_id:
            return None
        return {
            "id": field_id,
            "field": f"INFO:{field_id}",
            "number": header.get("Number", "."),
            "type": header.get("Type", "String"),
            "description": header.get("Description", ""),
        }

    @staticmethod
    def _cast_vcf_value(value: Optional[str], field_type: str) -> object:
        if field_type == "Flag":
            return value is not None
        if value in (None, "", "."):
            return ""
        if "," in value:
            return value
        if field_type == "Integer":
            try:
                return int(value)
            except ValueError:
                return value
        if field_type == "Float":
            try:
                return float(value)
            except ValueError:
                return value
        return value

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
    def _extract_consequence(info_map: Dict[str, str]) -> str:
        bcsq = info_map.get("BCSQ", "")
        if not bcsq:
            return ""
        consequences = []
        for annotation in bcsq.split(","):
            if not annotation or annotation.startswith("@"):
                continue
            term_field = annotation.split("|", 1)[0]
            consequences.extend(term for term in term_field.split("&") if term)
        return ",".join(sorted(set(consequences)))

    @classmethod
    def parse_vcf_records(
        cls,
        vcf_path: Path,
        variant_scope: str = "pass",
        sample_name: Optional[str] = None,
    ) -> Tuple[List[Dict[str, object]], List[Dict[str, str]], List[Dict[str, str]]]:
        if variant_scope not in {"pass", "all"}:
            raise PipelineExecutionError("variant_scope must be one of: pass, all")

        records: List[Dict[str, object]] = []
        format_fields: List[Dict[str, str]] = []
        info_fields: List[Dict[str, str]] = []
        format_field_headers: Dict[str, Dict[str, str]] = {}
        info_field_headers: Dict[str, Dict[str, str]] = {}
        seen_format_fields: set[str] = set()
        seen_info_fields: set[str] = set()
        sample = sample_name or vcf_path.stem

        with cls._open_text(vcf_path) as handle:
            for line in handle:
                if line.startswith("##"):
                    format_header = cls._parse_format_header(line)
                    if format_header:
                        format_field_headers[format_header["id"]] = format_header
                    info_header = cls._parse_info_header(line)
                    if info_header:
                        info_field_headers[info_header["id"]] = info_header
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

                chrom, pos_field, record_id, ref, alt_field, qual_field, filter_field, info_field = fields[:8]
                if variant_scope == "pass" and filter_field not in (".", "PASS"):
                    continue

                try:
                    pos = int(pos_field)
                except ValueError:
                    continue

                fmt = fields[8] if len(fields) > 8 else "."
                sample_field = fields[9] if len(fields) > 9 else "."
                info_map = cls._parse_info(info_field)
                sample_map = cls._parse_sample_map(fmt, sample_field)
                alts = [alt for alt in alt_field.split(",") if alt]
                types = [t for t in info_map.get("TYPE", "").split(",") if t]
                if not types:
                    types = [cls._infer_type(ref, alt) for alt in alts]

                record = {
                    "sample": sample,
                    "chrom": chrom,
                    "pos": pos,
                    "id": "" if record_id == "." else record_id,
                    "ref": ref,
                    "alt": alt_field,
                    "qual": "" if qual_field == "." else cls._cast_vcf_value(qual_field, "Float"),
                    "filter": "PASS" if filter_field == "." else filter_field,
                    "type": ",".join(types),
                    "consequence": cls._extract_consequence(info_map),
                }
                for key, value in info_map.items():
                    if key not in seen_info_fields:
                        info_fields.append(
                            info_field_headers.get(
                                key,
                                {
                                    "id": key,
                                    "field": f"INFO:{key}",
                                    "number": ".",
                                    "type": "String",
                                    "description": "",
                                },
                            )
                        )
                        seen_info_fields.add(key)
                    record[f"INFO:{key}"] = cls._cast_vcf_value(
                        value,
                        info_field_headers.get(key, {}).get("type", "String"),
                    )
                for key, value in sample_map.items():
                    if key not in seen_format_fields:
                        format_fields.append(
                            format_field_headers.get(
                                key,
                                {
                                    "id": key,
                                    "number": ".",
                                    "type": "String",
                                    "description": "",
                                },
                            )
                        )
                        seen_format_fields.add(key)
                    record[key] = cls._cast_vcf_value(
                        value,
                        format_field_headers.get(key, {}).get("type", "String"),
                    )
                records.append(record)

        return records, format_fields, info_fields

    @classmethod
    def write_variant_assets(
        cls,
        vcf_path: Path,
        variants_json: Path,
        regions_bed: Path,
        variant_scope: str = "pass",
        window_size: int = 1000,
        sample_name: Optional[str] = None,
    ) -> None:
        records, format_fields, info_fields = cls.parse_vcf_records(vcf_path, variant_scope, sample_name)
        variants_json.write_text(
            json.dumps({"variants": records, "sampleFields": format_fields, "infoFields": info_fields}),
            encoding="utf-8",
        )
        with regions_bed.open("w", encoding="utf-8") as handle:
            for record in records:
                chrom = str(record["chrom"])
                pos = int(record["pos"])
                ref_len = max(len(str(record["ref"])), 1)
                start = max(0, pos - window_size - 1)
                end = pos + ref_len + window_size
                handle.write(f"{chrom}\t{start}\t{end}\n")

    @staticmethod
    def _base64_file(path: Path) -> str:
        if not path.exists():
            raise PipelineExecutionError(f"Required sample-report file does not exist: {path}")
        return base64.b64encode(path.read_bytes()).decode("ascii")

    @classmethod
    def render_sample_report(
        cls,
        template_path: Path,
        output_html: Path,
        variants_json: Path,
        cram: Optional[Path],
        crai: Optional[Path],
        vcf: Path,
        reference: Optional[Path],
        reference_index: Optional[Path],
        title: str,
        sample_name: Optional[str] = None,
    ) -> None:
        has_igv = cram is not None and crai is not None and reference is not None and reference_index is not None
        if has_igv and not reference_index.exists():
            raise PipelineExecutionError(f"Reference index for sample report does not exist: {reference_index}")

        context = {
            "REPORT_NAME": title,
            "SAMPLE_NAME": sample_name or "",
            "DATETIME": datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
            "VERSION": __version__,
            "VARIANTS_JSON_B64": cls._base64_file(variants_json),
            "VCF_B64": cls._base64_file(vcf),
            "REFERENCE_FASTA_B64": cls._base64_file(reference) if has_igv else "",
            "REFERENCE_INDEX_B64": cls._base64_file(reference_index) if has_igv else "",
            "CRAM_B64": cls._base64_file(cram) if has_igv else "",
            "CRAI_B64": cls._base64_file(crai) if has_igv else "",
            "HAS_IGV": "true" if has_igv else "false",
            "REFERENCE_NAME": reference.name if reference else "",
            "VCF_NAME": vcf.name,
            "CRAM_NAME": cram.name if cram else "",
            "CRAI_NAME": crai.name if crai else "",
        }

        template_content = template_path.read_text(encoding="utf-8")
        template_vars = set(re.findall(r"{{\s*(\w+)\s*}}", template_content))
        missing_template_vars = [v for v in template_vars if v not in context]
        if missing_template_vars:
            raise ValueError(
                f"Template variable(s) '{', '.join(missing_template_vars)}' not found in context for template {template_path}"
            )
        for key, value in context.items():
            template_content = re.sub(r"{{\s*" + re.escape(key) + r"\s*}}", str(value), template_content)
        output_html.write_text(template_content, encoding="utf-8")

class FormatHTMLReportTemplate(BaseStage):
    template_path: Path = Field(..., description="Path to the template file")
    context: Context = Field(default_factory=dict, description="Context variables for rendering the template")

    @property
    def output(self) -> FormatHTMLReportTemplateOutput:
        return FormatHTMLReportTemplateOutput(rendered=f"{self.prefix}.html")

    def create_commands(self, ctx) -> List[PythonCommand]:
        return [
            self.python_cmd(
                func=self.render_template,
                description="Render a template with provided context"
            )
        ]
    
    def validate_context(self, context: Context) -> None:
        """
        Override this method in subclasses to implement custom validation logic for the context variables before rendering the template.
        """
        pass

    def default_context(self) -> Context:
        return {
            'DATETIME': datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
            'VERSION': __version__,
            'USER': os.getenv("USER", "Unknown"),
        }

    def render_template(self) -> None:
        template = self.template_path
        context = {**self.default_context(), **self.context}
        # eval any callable context values (e.g. lambda functions) to get their actual value before rendering
        for k, v in context.items():
            if callable(v):
                context[k] = v()
            elif isinstance(context[k], Path):
                context[k] = v.read_text().strip()
        # Custom validation logic
        self.validate_context(context)
        # load template
        with open(template, "r") as f:
            template_content = f.read()
        template_vars = set(re.findall(r"{{\s*(\w+)\s*}}", template_content))
        # check to see if any of the template variables are missing from the context
        missing_template_vars = [v for v in template_vars if v not in context]
        if missing_template_vars:
            # if the template has variables that are not provided in the context, that likely 
            # indicates a mistake in the template or missing context variables, so we should 
            # raise an error instead of silently rendering with missing values
            raise ValueError(f"Template variable(s) '{', '.join(missing_template_vars)}' not found in context for template {template}")

        # Finally, render the template with the context values
        for k, v in context.items():
            template_content = re.sub(r"{{\s*" + re.escape(k) + r"\s*}}", str(v).strip(), template_content)

        with open(self.output.rendered, "w") as f:
            f.write(template_content)

class TreeReport(FormatHTMLReportTemplate):
    template_path: Path = Path(__file__).resolve().parent.parent / "templates" / "report-tree" / "snippy-report-tree.html"
    iso3166_2: Path = Path(__file__).resolve().parent.parent / "templates" / "report-tree" / "iso3166-2-export.csv"
    ladderize: bool = Field(default=False, description="Ladderize the tree in the report")
    mid_point_root: bool = Field(default=False, description="Mid-point root the tree in the report")
    remove_invalid_rows: bool = Field(default=True, description="Remove metadata rows that do not match any tip in the tree, instead of raising an error")

    _dependencies = [biopython, phylocanvas, phylojs]

    def validate_context(self, context: Dict[str, Union[str, int, float, Callable]]) -> None:
        required_keys = {"NEWICK", "REPORT_NAME", "METADATA", "LOGS"}
        missing_keys = [k for k in required_keys if k not in context]
        if missing_keys:
            raise ValueError(f"Context key(s) '{', '.join(missing_keys)}' are required for TreeReport but not found in context")
        
        from Bio import Phylo
        from io import StringIO

        if isinstance(context["NEWICK"], Path):
            # If NEWICK is a Path, read the content and replace
            context["NEWICK"] = context["NEWICK"].read_text().strip()
        try:
            newick_str = context["NEWICK"]
            tree = Phylo.read(StringIO(newick_str), "newick")
            if self.mid_point_root:
                tree.root_at_midpoint()
            if self.ladderize:
                tree.ladderize()
            handle = StringIO()
            # Biopython defaults to fixed-point formatting for branch lengths, which can
            # round very small branches to zero when we rewrite the Newick string.
            Phylo.write(tree, handle, "newick", format_branch_length=NEWICK_BRANCH_LENGTH_FORMAT)
            context["NEWICK"] = handle.getvalue().strip()
        except Exception as e:
            raise PipelineExecutionError(f"Invalid NEWICK string provided in context for TreeReport: {e}")
        
        if context.get("METADATA") is not None:
            if isinstance(context["METADATA"], str):
                # assume it's a JSON string 
                context["METADATA_JSON"] = context["METADATA"]
            elif isinstance(context["METADATA"], Path):
                try:
                    context["METADATA_JSON"] = load_metadata_as_json_str(context["METADATA"])
                except Exception as e:
                    raise PipelineExecutionError(f"Error loading metadata file provided in context for TreeReport: {e}")
            else:
                raise PipelineExecutionError("METADATA context variable for TreeReport must be either a JSON string or a Path to a metadata file")
        if "METADATA_JSON" not in context:
            context["METADATA_JSON"] = None
        if context["METADATA_JSON"] is not None:
            try:
                metadata_json_str = context["METADATA_JSON"]
                metadata_json = json.loads(metadata_json_str)
            except Exception as e:
                raise PipelineExecutionError(f"Invalid METADATA_JSON string provided in context for TreeReport: {e}")
            # check all the metadata entries have an id_column that matches a tip in the tree
            tree_tips = {tip.name for tip in tree.get_terminals()}
            id_column = None
            metadata = []
            for entry in metadata_json:
                if "id" in entry:
                    id_column = "id"
                elif "sample_id" in entry:
                    id_column = "sample_id"
                elif "sample" in entry:
                    id_column = "sample"
                elif "name" in entry:
                    id_column = "name"
                else:
                    raise PipelineExecutionError(f"Metadata entry {entry} is missing required 'id' field for TreeReport context")
                if entry[id_column] not in tree_tips:
                    if self.remove_invalid_rows:
                        print(f"Metadata {id_column} '{entry[id_column]}' does not match any tip in the NEWICK tree for TreeReport context, skipping this metadata entry")
                        continue
                    raise PipelineExecutionError(f"Metadata {id_column} '{entry[id_column]}' does not match any tip in the NEWICK tree for TreeReport context")
                metadata.append(entry)
            context["METADATA_JSON"] = json.dumps(metadata)
        # add iso3166-2 country code mapping to context for use in the template
        iso3166_mapping = {}
        with open(self.iso3166_2, "r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                subdivision_code = row.get("subdivisionCode")
                lat_lng_str = row.get("latLng")
                if not subdivision_code or not lat_lng_str:
                    continue

                try:
                    lat_lng = json.loads(lat_lng_str)
                    if isinstance(lat_lng, list) and len(lat_lng) == 2:
                        iso3166_mapping[subdivision_code] = lat_lng
                except Exception:
                    continue
                
        context["ISO3166_MAPPING_JSON"] = json.dumps(iso3166_mapping)

        # convert None to null for html
        for k, v in context.items():
            if v is None:
                context[k] = "null"
