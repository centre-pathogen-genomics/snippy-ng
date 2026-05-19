from __future__ import annotations

import math
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


class StrandBiasError(ValueError):
    pass


@dataclass(frozen=True)
class StrandCounts:
    ref_forward: int
    ref_reverse: int
    alt_forward: int
    alt_reverse: int

    @property
    def ref_depth(self) -> int:
        return self.ref_forward + self.ref_reverse

    @property
    def alt_depth(self) -> int:
        return self.alt_forward + self.alt_reverse


@dataclass(frozen=True)
class StrandBiasResult(StrandCounts):
    pvalue: float
    phred: float

    @property
    def alt_forward_fraction(self) -> float | None:
        total_alt = self.alt_forward + self.alt_reverse
        if total_alt == 0:
            return None
        return self.alt_forward / total_alt


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    row1 = a + b
    row2 = c + d
    col1 = a + c
    total = row1 + row2

    if total == 0 or row1 == 0 or row2 == 0 or col1 == 0 or col1 == total:
        return 1.0

    lower = max(0, col1 - row2)
    upper = min(row1, col1)

    observed = math.exp(_log_comb(col1, a) + _log_comb(total - col1, row1 - a) - _log_comb(total, row1))
    pvalue = 0.0
    epsilon = max(observed * 1e-12, 1e-15)

    for x in range(lower, upper + 1):
        probability = math.exp(
            _log_comb(col1, x)
            + _log_comb(total - col1, row1 - x)
            - _log_comb(total, row1)
        )
        if probability <= observed + epsilon:
            pvalue += probability

    return min(1.0, pvalue)


def _phred_score(pvalue: float) -> float:
    bounded = max(min(pvalue, 1.0), 1e-300)
    return -10.0 * math.log10(bounded)


def parse_mpileup_bases(read_bases: str) -> list[tuple[str, str | None, str | None]]:
    observations: list[tuple[str, str | None, str | None]] = []
    index = 0

    while index < len(read_bases):
        char = read_bases[index]

        if char == "^":
            index += 2
            continue
        if char == "$":
            index += 1
            continue
        if char in "*#<>":
            observations.append((char, None, None))
            index += 1
            continue
        if char not in ".,ACGTNacgtn":
            index += 1
            continue

        strand = "forward" if char in ".ACGTN" else "reverse"
        index += 1
        indel: str | None = None

        if index < len(read_bases) and read_bases[index] in "+-":
            sign = read_bases[index]
            index += 1
            digits_start = index
            while index < len(read_bases) and read_bases[index].isdigit():
                index += 1
            if digits_start == index:
                raise StrandBiasError(f"Invalid mpileup indel encoding: {read_bases}")
            indel_length = int(read_bases[digits_start:index])
            indel_sequence = read_bases[index:index + indel_length]
            if len(indel_sequence) != indel_length:
                raise StrandBiasError(f"Invalid mpileup indel length in: {read_bases}")
            indel = f"{sign}{indel_sequence.upper()}"
            index += indel_length

        observations.append((char, strand, indel))

    return observations


def _is_snp(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and not alt.startswith("<")


def _is_supported_indel(ref: str, alt: str) -> bool:
    if alt.startswith("<"):
        return False
    return (len(ref) != len(alt)) and (ref.startswith(alt) or alt.startswith(ref))


def count_strands_for_variant(ref: str, alt: str, read_bases: str) -> StrandCounts | None:
    observations = parse_mpileup_bases(read_bases)

    if _is_snp(ref, alt):
        ref_forward = sum(1 for base, _, _ in observations if base == ".")
        ref_reverse = sum(1 for base, _, _ in observations if base == ",")
        alt_forward = sum(1 for base, _, _ in observations if base == alt.upper())
        alt_reverse = sum(1 for base, _, _ in observations if base == alt.lower())
        return StrandCounts(
            ref_forward=ref_forward,
            ref_reverse=ref_reverse,
            alt_forward=alt_forward,
            alt_reverse=alt_reverse,
        )

    if not _is_supported_indel(ref, alt):
        return None

    if alt.startswith(ref):
        target_indel = f"+{alt[len(ref):].upper()}"
    else:
        target_indel = f"-{ref[len(alt):].upper()}"

    ref_forward = 0
    ref_reverse = 0
    alt_forward = 0
    alt_reverse = 0

    for base, strand, indel in observations:
        if strand is None:
            continue
        if indel == target_indel:
            if strand == "forward":
                alt_forward += 1
            else:
                alt_reverse += 1
            continue
        if base == ".":
            ref_forward += 1
        elif base == ",":
            ref_reverse += 1

    return StrandCounts(
        ref_forward=ref_forward,
        ref_reverse=ref_reverse,
        alt_forward=alt_forward,
        alt_reverse=alt_reverse,
    )


def calculate_strand_bias(ref: str, alt: str, read_bases: str) -> StrandBiasResult | None:
    counts = count_strands_for_variant(ref=ref, alt=alt, read_bases=read_bases)
    if counts is None:
        return None
    pvalue = fisher_exact_two_sided(
        counts.ref_forward,
        counts.ref_reverse,
        counts.alt_forward,
        counts.alt_reverse,
    )
    return StrandBiasResult(
        ref_forward=counts.ref_forward,
        ref_reverse=counts.ref_reverse,
        alt_forward=counts.alt_forward,
        alt_reverse=counts.alt_reverse,
        pvalue=pvalue,
        phred=_phred_score(pvalue),
    )


def _append_filter(filter_value: str, label: str) -> str:
    if filter_value in {"", ".", "PASS"}:
        return label
    labels = filter_value.split(";")
    if label not in labels:
        labels.append(label)
    return ";".join(labels)


def _append_info(info: str, key: str, value: str) -> str:
    if info in {"", "."}:
        return f"{key}={value}"
    return f"{info};{key}={value}"


def _format_float(value: float) -> str:
    return f"{value:.6g}"


def _format_pvalue_and_phred(pvalue: float, phred: float) -> tuple[str, str]:
    formatted_pvalue = _format_float(pvalue)
    if formatted_pvalue == "1":
        return formatted_pvalue, "0"
    return formatted_pvalue, _format_float(phred)


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return "."
    return _format_float(value)


def samtools_mpileup(alignment: Path, reference: Path, positions_bed: Path) -> str:
    command = ["samtools", "mpileup", "-aa", "-B", "-f", str(reference), "-l", str(positions_bed), str(alignment)]
    result = subprocess.run(
        command,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.stdout


def parse_mpileup_output(lines: Iterable[str]) -> dict[tuple[str, int], str]:
    pileups: dict[tuple[str, int], str] = {}

    for line in lines:
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            raise StrandBiasError(f"Invalid samtools mpileup row: {line}")
        try:
            position = int(fields[1])
        except ValueError as exc:
            raise StrandBiasError(f"Invalid mpileup coordinate: {line}") from exc
        pileups[(fields[0], position)] = fields[4] if len(fields) > 4 else ""

    return pileups


def _collect_variant_positions(vcf_lines: Iterable[str]) -> list[tuple[str, int]]:
    positions: list[tuple[str, int]] = []
    seen: set[tuple[str, int]] = set()

    for line in vcf_lines:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            continue
        chrom = fields[0]
        pos = int(fields[1])
        key = (chrom, pos)
        if key not in seen:
            positions.append(key)
            seen.add(key)

    return positions


def annotate_vcf_strand_bias(
    vcf: Path,
    bam: Path,
    reference: Path,
    output: Path | None = None,
    filter_pvalue: float | None = None,
) -> None:
    vcf_lines = vcf.read_text(encoding="utf-8").splitlines(keepends=True)
    positions = _collect_variant_positions(vcf_lines)

    pileups: dict[tuple[str, int], str] = {}
    if positions:
        with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=True) as bed_handle:
            for chrom, pos in positions:
                bed_handle.write(f"{chrom}\t{pos - 1}\t{pos}\n")
            bed_handle.flush()
            pileups = parse_mpileup_output(
                samtools_mpileup(bam, reference, Path(bed_handle.name)).splitlines()
            )

    if output is None:
        import sys

        handle: TextIO = sys.stdout
        close_handle = False
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        handle = output.open("w", encoding="utf-8")
        close_handle = True

    try:
        for line in vcf_lines:
            if line.startswith("#CHROM"):
                handle.write('##INFO=<ID=SB_REF_FWD,Number=1,Type=Integer,Description="Reference-supporting forward-strand reads from samtools mpileup">\n')
                handle.write('##INFO=<ID=SB_REF_REV,Number=1,Type=Integer,Description="Reference-supporting reverse-strand reads from samtools mpileup">\n')
                handle.write('##INFO=<ID=SB_ALT_FWD,Number=1,Type=Integer,Description="Alternate-supporting forward-strand reads from samtools mpileup">\n')
                handle.write('##INFO=<ID=SB_ALT_REV,Number=1,Type=Integer,Description="Alternate-supporting reverse-strand reads from samtools mpileup">\n')
                handle.write('##INFO=<ID=SB_ALT_FWD_FRAC,Number=1,Type=Float,Description="Fraction of alternate-supporting reads on the forward strand: SB_ALT_FWD / (SB_ALT_FWD + SB_ALT_REV)">\n')
                handle.write('##INFO=<ID=SB_PVALUE,Number=1,Type=Float,Description="Two-sided Fisher exact strand-bias p-value">\n')
                handle.write('##INFO=<ID=SB_PHRED,Number=1,Type=Float,Description="Phred-scaled strand-bias score derived from SB_PVALUE">\n')
                handle.write('##INFO=<ID=SB_METHOD,Number=1,Type=String,Description="Strand-bias method used to compute SB_PVALUE">\n')
                if filter_pvalue is not None:
                    handle.write(f'##FILTER=<ID=StrandBias,Description="Strand-bias p-value below {filter_pvalue:g}">\n')
                handle.write(line)
                continue

            if line.startswith("#"):
                handle.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                handle.write(line)
                continue

            chrom, pos_text, _id, ref, alt = fields[:5]
            alt_alleles = alt.split(",")
            if len(alt_alleles) != 1:
                handle.write(line)
                continue

            result = calculate_strand_bias(
                ref=ref,
                alt=alt_alleles[0],
                read_bases=pileups.get((chrom, int(pos_text)), ""),
            )
            if result is None:
                handle.write(line)
                continue

            formatted_pvalue, formatted_phred = _format_pvalue_and_phred(result.pvalue, result.phred)
            formatted_alt_forward_fraction = _format_optional_float(result.alt_forward_fraction)

            fields[7] = _append_info(fields[7], "SB_REF_FWD", str(result.ref_forward))
            fields[7] = _append_info(fields[7], "SB_REF_REV", str(result.ref_reverse))
            fields[7] = _append_info(fields[7], "SB_ALT_FWD", str(result.alt_forward))
            fields[7] = _append_info(fields[7], "SB_ALT_REV", str(result.alt_reverse))
            fields[7] = _append_info(fields[7], "SB_ALT_FWD_FRAC", formatted_alt_forward_fraction)
            fields[7] = _append_info(fields[7], "SB_PVALUE", formatted_pvalue)
            fields[7] = _append_info(fields[7], "SB_PHRED", formatted_phred)
            fields[7] = _append_info(fields[7], "SB_METHOD", "FisherExact")

            if filter_pvalue is not None and result.pvalue < filter_pvalue:
                fields[6] = _append_filter(fields[6], "StrandBias")

            handle.write("\t".join(fields) + "\n")
    finally:
        if close_handle:
            handle.close()
