import gzip
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Union
from snippy_ng.stages import BaseStage, BaseOutput
from snippy_ng.stages import PythonCommand
from snippy_ng.logging import logger
from pydantic import Field


class PrintVcfSummary(BaseStage):
    vcf_files: List[Path] = Field(..., description="Input VCF files to summarize")

    @property
    def output(self) -> None:
        return BaseOutput()

    @property
    def commands(self) -> List[PythonCommand]:
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

    @property
    def commands(self) -> List[PythonCommand]:
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