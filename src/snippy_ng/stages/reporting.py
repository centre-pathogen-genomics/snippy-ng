import gzip
import os
import re
import csv
from pathlib import Path
import datetime
import json
from typing import Callable, Dict, Iterable, List, Optional, Tuple, Union
from snippy_ng.stages import BaseStage, BaseOutput
from snippy_ng.stages import PythonCommand
from snippy_ng.logging import logger
from snippy_ng.exceptions import PipelineExecutionError
from snippy_ng.dependencies import biopython
from snippy_ng.__about__ import __version__
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

        def iter_vcf_records(vcf_path: Union[str, Path]) -> Iterable[Tuple[str, int]]:
            with open_text_maybe_gzip(vcf_path) as f:
                for line in f:
                    if not line or line[0] == "#":
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 2:
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

            for chrom, pos in iter_vcf_records(vcf_path):
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
        ) -> None:
            contigs = parse_contigs_from_vcf_header(vcf_path)
            if contig_order == "alpha":
                contigs = sorted(contigs, key=lambda x: x[0])
            elif contig_order != "header":
                raise ValueError("contig_order must be 'header' or 'alpha'")

            try:
                term_cols = max(os.get_terminal_size().columns, 40)
            except OSError:
                term_cols = 80

            usable_cols = max(10, term_cols - max(0, margin_cols * (len(contigs) - 1)))

            cols_per_contig = allocate_columns_proportional(
                contigs, usable_cols, min_cols_per_contig=min_cols_per_contig
            )

            counts = genome_wide_binned_counts_variable_width(vcf_path, contigs, cols_per_contig)

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
                logger.warning(f"No variants found in VCF {vcf_path}, skipping histogram.")
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
        )

class FormatHTMLReportTemplateOutput(BaseOutput):
    rendered: Path

ContextValue = Optional[Union[str, int, float, Path, ]]
Context = Dict[str, ContextValue]

class FormatHTMLReportTemplate(BaseStage):
    template_path: Path = Field(..., description="Path to the template file")
    context: Context = Field(default_factory=dict, description="Context variables for rendering the template")

    @property
    def output(self) -> FormatHTMLReportTemplateOutput:
        return FormatHTMLReportTemplateOutput(rendered="report.html")

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
            # check value is of type Path and read the file content
            if isinstance(v, Path):
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

class EpiReport(FormatHTMLReportTemplate):
    template_path: Path = Path(__file__).resolve().parent.parent / "templates" / "epi-report" / "snippy-epi-report.html"
    iso3166_2: Path = Path(__file__).resolve().parent.parent / "templates" / "epi-report" / "iso3166-2-export.csv"
    preprocess_tree: bool = Field(default=True, description="Whether to pre-process the NEWICK tree by rooting at midpoint and ladderizing before rendering the report")

    _dependencies = [biopython]

    def validate_context(self, context: Dict[str, Union[str, int, float, Callable]]) -> None:
        required_keys = {"NEWICK", "REPORT_NAME", "METADATA_JSON", "LOGS"}
        missing_keys = [k for k in required_keys if k not in context]
        if missing_keys:
            raise ValueError(f"Context key(s) '{', '.join(missing_keys)}' are required for EpiReport but not found in context")
        
        from Bio import Phylo
        from io import StringIO
        try:
            newick_str = context["NEWICK"]
            tree = Phylo.read(StringIO(newick_str), "newick")
            if self.preprocess_tree:
                tree.root_at_midpoint()
                tree.ladderize()
            handle = StringIO()
            Phylo.write(tree, handle, "newick")
            context["NEWICK"] = handle.getvalue().strip()
        except Exception as e:
            raise PipelineExecutionError(f"Invalid NEWICK string provided in context for EpiReport: {e}")
        
        if context["METADATA_JSON"] is not None:
            try:
                metadata_json_str = context["METADATA_JSON"]
                metadata = json.loads(metadata_json_str)
            except Exception as e:
                raise PipelineExecutionError(f"Invalid METADATA_JSON string provided in context for EpiReport: {e}")
            # check all the metadata entries have an "id" field that matches a tip in the tree
            tree_tips = {tip.name for tip in tree.get_terminals()}
            for entry in metadata:
                if "id" not in entry:
                    raise PipelineExecutionError(f"Metadata entry {entry} is missing required 'id' field for EpiReport context")
                if entry["id"] not in tree_tips:
                    raise PipelineExecutionError(f"Metadata entry id '{entry['id']}' does not match any tip in the NEWICK tree for EpiReport context")

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