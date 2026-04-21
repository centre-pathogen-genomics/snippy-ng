from pathlib import Path
from typing import List, Optional
import re

from snippy_ng.stages import BaseStage, BaseOutput
from snippy_ng.dependencies import bcftools
from snippy_ng.logging import logger
from pydantic import Field

class VcfFilterOutput(BaseOutput):
    vcf: Path = Field(..., description="Filtered and normalized VCF file")


class VcfPassFilterOutput(BaseOutput):
    vcf: Path = Field(..., description="VCF containing only PASS variants")


class VcfPassFilter(BaseStage):
    vcf: Path = Field(..., description="Input VCF file to subset to PASS variants")
    no_insertions: bool = Field(True, description="Remove insertions from the output VCF")

    _dependencies = [bcftools]

    @property
    def output(self) -> VcfPassFilterOutput:
        return VcfPassFilterOutput(
            vcf=Path(f"{self.prefix}.pass.vcf")
        )

    def create_commands(self, ctx) -> List:
        include_expr = 'FMT/GT="1/1" || FMT/GT="0/1"'
        if self.no_insertions:
            # Exclude insertions while retaining deletions, including symbolic DEL blocks.
            include_expr = f'({include_expr}) && (strlen(ALT)<=strlen(REF) || ALT="<DEL>")'
        return [
            self.shell_cmd(
                ["bcftools", "view", "-f", "PASS", "-i", include_expr, str(self.vcf)],
                description="Filter VCF to PASS variants only",
                output_file=self.output.vcf,
            )
        ]


class VcfFilter(BaseStage):
    vcf: Path = Field(..., description="Input VCF file to filter")
    reference: Path = Field(..., description="Reference FASTA file")

    _dependencies = [bcftools]

    @property
    def output(self) -> VcfFilterOutput:
        filtered_vcf = f"{self.prefix}.filtered.vcf"
        return VcfFilterOutput(vcf=filtered_vcf)
    
    def test_check_if_vcf_has_variants(self):
        """Test that the output VCF file is not empty."""
        vcf_path = self.output.vcf
        with open(vcf_path, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    return  # Found a non-header line, VCF is not empty
        # give warning instead of error as some callers may produce empty VCFs if no variants are found
        logger.warning(f"Output VCF file {vcf_path} has no variants (only header lines). Please check if this is expected based on your data and parameters.")

class VcfFilterShort(VcfFilter):
    """
    Filter VCF files using Samtools to remove unwanted variants.
    """

    min_qual: float = Field(100.0, description="Mark variants below this QUAL threshold as LowQual")
    min_depth: int = Field(10, description="Mark variants below this depth threshold as LowDepth")

    def create_commands(self, ctx) -> List:
        """Constructs the samtools view command for filtering."""

        # Build the post-norm filter. We filter AFTER splitting/normalizing and after recomputing TYPE.
        base_filter = " && ".join([
            'ALT!="*"'
        ])
        ram_gb_param = ["-m", f"{ctx.ram}G"] if ctx.ram else []
        commands = [
                self.shell_cmd(
                    ["bcftools", "view", "-Ou", "-i", 'GT="alt"', str(self.vcf)],
                    description="Remove non-alt alleles"
                ),
                self.shell_cmd(
                    [
                        "bcftools",
                        "norm",
                        "-f",
                        str(self.reference),
                        "--check-ref", "e",
                        "--remove-duplicates",
                        "-Ou",
                    ],
                    description="Normalize and split multiallelic variants",
                ),
                self.shell_cmd(
                    ["bcftools", "+fill-tags", "-Ou", "-", "--", "-t", "TYPE"],
                    description="Recompute TYPE from REF/ALT",
                ),
                self.shell_cmd(
                    ["bcftools", "sort", "-Ou", *ram_gb_param],
                    description="Sort VCF"
                ),
                self.shell_cmd(
                    ["bcftools", "view", "-Ou", "--include", base_filter, "-"],
                    description="Filter variants after normalization and TYPE recomputation",
                ),
                self.shell_cmd(
                    ["bcftools", "filter", "-s", "LowQual", "-m", "+", "-e", f"QUAL<{self.min_qual}", "-"],
                    description=f"Mark variants with QUAL<{self.min_qual} as LowQual and others as PASS",
                ),
                self.shell_cmd(
                    ["bcftools", "filter", "-s", "LowDepth", "-m", "+", "-e", f"FMT/DP<{self.min_depth}", "-"],
                    description=f"Mark variants with DP<{self.min_depth} as LowDepth and preserve existing FILTER labels",
                ),
            ]
        bcftools_pipeline = self.shell_pipe(
            commands=commands,
            description="Normalize, recompute TYPE, filter, and annotate variants",
            output_file=Path(self.output.vcf),
        )
        return [bcftools_pipeline]

class VcfFilterAsm(VcfFilterShort):
    """
    Filter VCF files for assemblies using bcftools to remove unwanted variants.
    """
    min_qual: float = Field(60, description="Minimum QUAL score for assembly-based calling")
    min_depth: int = Field(1, description="Minimum depth for assembly-based calling")


class VcfFilterLong(VcfFilter):
    """
    Filter VCF files for long-read variant calling using bcftools to remove unwanted variants.
    
    This pipeline handles long-read specific filtering including:
    - Reheadering with all reference contigs
    - Making heterozygous calls homozygous for the allele with most depth
    - Filtering out non-alt alleles and missing alleles
    - Normalizing and left-aligning indels
    - Removing long indels and duplicates
    - Converting to haploid genotypes
    """
    reference_index: Path = Field(..., description="Reference FASTA index file (.fai)")
    min_qual: Optional[float] = Field(None, description="Mark variants below this QUAL threshold as LowQual")
    min_depth: int = Field(10, description="Mark variants below this depth threshold as LowDepth")
    max_indel: int = Field(10000, description="Maximum indel length to keep")
    
    def create_commands(self, ctx) -> List:
        """Constructs the bcftools pipeline for long-read variant filtering."""
        
        # Create temp files for contigs and header
        contigs_file = f"{self.prefix}.contigs.txt"
        header_file = f"{self.prefix}.header.txt"
        
        # Generate contig lines from faidx
        create_contigs_cmd = self.shell_cmd(
            ["awk", '{print "##contig=<ID="$1",length="$2">"}', str(self.reference_index)],
            description="Generate contig lines from reference index",
            output_file=Path(contigs_file)
        )
        
        # Create new header with all contigs
        # This combines: bcftools view -h | grep -v "^##contig=" | sed -e "3r $contigs"
        create_header_pipeline = self.shell_pipe(
            commands=[
                self.shell_cmd(
                    ["bcftools", "view", "-h", str(self.vcf)],
                    description="Extract VCF header"
                ),
                self.shell_cmd(
                    ["grep", "-v", "^##contig="],
                    description="Drop existing contig headers"
                ),
                self.shell_cmd(
                    ["sed", "-e", f"3r {contigs_file}"],
                    description="Insert all reference contigs into header"
                ),
            ],
            description="Create VCF header with all contigs",
            output_file=Path(header_file),
        )
        
        # Build the main filtering pipeline
        pipeline_commands = [
            self.shell_cmd(
                ["bcftools", "reheader", "-h", header_file, str(self.vcf)],
                description="Replace VCF header with new header containing all contigs"
            ),
        ]
        
        # Continue with the filtering pipeline
        pipeline_commands.extend([
            self.shell_cmd(
                ["bcftools", "view", "-Ou", "-i", 'GT="alt"'],
                description="Remove non-alt alleles"
            ),
            self.shell_cmd(
                [
                    "bcftools",
                    "norm",
                    "-f",
                    str(self.reference),
                    "--check-ref", "e",
                    "--remove-duplicates",
                    "-Ou",
                ],
                description="Normalize and split multiallelic variants",
            ),
            self.shell_cmd(
                ["bcftools", "filter", "-Ou", "-e", f'abs(ILEN)>{self.max_indel} || ALT="*"'],
                description=f"Remove indels longer than {self.max_indel}bp or sites with unobserved alleles"
            ),
            self.shell_cmd(
                    ["bcftools", "+fill-tags", "-Ou", "-", "--", "-t", "TYPE"],
                    description="Recompute TYPE from REF/ALT",
            ),
            self.shell_cmd(
                ["bcftools", "sort", "-Ou", "-m", str(ctx.ram)],
                description="Sort VCF"
            ),
            self.shell_cmd(
                ["bcftools", "view", "-Ou", "-i", 'GT="alt"'],
                description="Remove non-alt alleles and output final VCF"
            ),
            self.shell_cmd(
                ["bcftools", "filter", "-s", "LowDepth", "-m", "+", "-e", f"FMT/DP<{self.min_depth}", "-"],
                description=f"Mark variants with DP<{self.min_depth} as LowDepth and preserve existing FILTER labels",
            ),
        ])

        if self.min_qual is not None:
            pipeline_commands.append( 
                self.shell_cmd(
                    ["bcftools", "filter", "-s", "LowQual", "-m", "+", "-e", f"QUAL<{self.min_qual}", "-"],
                    description=f"Mark variants with QUAL<{self.min_qual} as LowQual and others as PASS",
                )
            )

        main_pipeline = self.shell_pipe(
            commands=pipeline_commands,
            description="Long-read variant filtering pipeline",
            output_file=Path(self.output.vcf)
        )
        
        # Step 4: Cleanup temp files
        cleanup_cmd = self.shell_cmd(
            ["rm", "-f", contigs_file, header_file],
            description="Remove temporary files"
        )
        
        return [create_contigs_cmd, create_header_pipeline, main_pipeline, cleanup_cmd]

    




class AddDeletionsToVCFOutput(BaseOutput):
    """Output from the zero-depth deletion annotation stage."""
    vcf: Path = Field(..., description="VCF with zero-depth regions added as <DEL> blocks")


class AddDeletionsToVCF(BaseStage):
    """
    Add zero-depth regions as symbolic deletion blocks to a VCF.

    This stage detects contiguous depth==0 blocks from the BAM and appends
    `<DEL>` records to the input VCF with valid ALT/INFO headers.
    """
    zero_depth_bed: Path = Field(..., description="BED file with zero-depth regions")
    vcf: Path = Field(..., description="Input VCF file")
    reference: Path = Field(..., description="Reference FASTA file")

    @property
    def output(self) -> AddDeletionsToVCFOutput:
        return AddDeletionsToVCFOutput(
            vcf=Path(f"{self.prefix}.with_del.vcf")
        )

    def create_commands(self, ctx) -> List:
        """Generate zero-depth BED and merge symbolic deletion blocks into VCF."""
        return [
            self.python_cmd(
                func=self._merge_zero_depth_deletions_into_vcf,
                args=[self.vcf, self.zero_depth_bed, self.output.vcf, self.reference],
                description="Add zero-depth <DEL> blocks to VCF",
            ),
        ]

    @staticmethod
    def _extract_info_int(info_field: str, key: str):
        match = re.search(rf'(?:^|;){re.escape(key)}=([0-9]+)(?:;|$)', info_field)
        if not match:
            return None
        try:
            return int(match.group(1))
        except ValueError:
            return None

    @staticmethod
    def _existing_deletion_keys(chrom: str, pos: int, ref: str, alt: str, info_field: str) -> set[tuple[str, int, int]]:
        keys: set[tuple[str, int, int]] = set()

        if alt == "<DEL>" or "SVTYPE=DEL" in info_field:
            end_pos = AddDeletionsToVCF._extract_info_int(info_field, "END")
            if end_pos is not None:
                keys.add((chrom, pos, end_pos))

        if not ref or ref == ".":
            return keys

        for alt_allele in alt.split(","):
            if not alt_allele or alt_allele in {".", "*"}:
                continue
            if alt_allele.startswith("<") and alt_allele.endswith(">"):
                continue
            if "[" in alt_allele or "]" in alt_allele:
                continue

            alt_seq = "" if alt_allele == "-" else alt_allele
            if len(ref) <= len(alt_seq):
                continue
            if not ref.startswith(alt_seq):
                continue

            end_pos = pos + (len(ref) - len(alt_seq))
            keys.add((chrom, pos, end_pos))

        return keys

    @staticmethod
    def _variant_interval(pos: int, ref: str, info_field: str, alt: str) -> tuple[int, int]:
        start = pos
        end = pos

        if ref and ref != ".":
            end = pos + max(len(ref), 1) - 1

        if alt == "<DEL>" or "SVTYPE=DEL" in info_field:
            info_end = AddDeletionsToVCF._extract_info_int(info_field, "END")
            if info_end is not None:
                end = max(end, info_end)

        return (start, end)

    @staticmethod
    def _intervals_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
        return start_a <= end_b and start_b <= end_a

    @staticmethod
    def _overlaps_any_interval(
        chrom: str,
        start: int,
        end: int,
        occupied_intervals: dict[str, list[tuple[int, int]]],
    ) -> bool:
        for existing_start, existing_end in occupied_intervals.get(chrom, []):
            if AddDeletionsToVCF._intervals_overlap(start, end, existing_start, existing_end):
                return True
        return False

    @staticmethod
    def _load_reference_sequences(reference_fasta: Path) -> dict[str, str]:
        sequences: dict[str, str] = {}
        current_name: str | None = None
        current_chunks: list[str] = []

        with open(reference_fasta, "r") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_name is not None:
                        sequences[current_name] = "".join(current_chunks).upper()
                    current_name = line[1:].split()[0]
                    current_chunks = []
                    continue
                if current_name is not None:
                    current_chunks.append(line)

        if current_name is not None:
            sequences[current_name] = "".join(current_chunks).upper()

        return sequences

    @staticmethod
    def _merge_zero_depth_deletions_into_vcf(input_vcf: Path, zero_depth_bed: Path, output_vcf: Path, reference_fasta: Path) -> None:
        meta_lines: List[str] = []
        chrom_line = ""
        variant_rows: List[tuple[str, int, int, str]] = []
        contig_order: list[str] = []
        contig_rank: dict[str, int] = {}
        has_alt_del = False
        has_info_end = False
        has_info_svtype = False
        has_info_svlen = False
        has_info_type = False
        has_info_zd = False
        sample_count = 0
        reference_sequences = AddDeletionsToVCF._load_reference_sequences(reference_fasta)

        existing_del_keys: set[tuple[str, int, int]] = set()
        occupied_intervals: dict[str, list[tuple[int, int]]] = {}

        with open(input_vcf, "r") as handle:
            record_index = 0
            for raw_line in handle:
                if raw_line.startswith("##"):
                    meta_lines.append(raw_line)
                    if raw_line.startswith("##contig=<ID="):
                        contig_id = raw_line.split("##contig=<ID=", 1)[1].split(",", 1)[0].split(">", 1)[0]
                        if contig_id not in contig_rank:
                            contig_rank[contig_id] = len(contig_order)
                            contig_order.append(contig_id)
                    if raw_line.startswith("##ALT=<ID=DEL"):
                        has_alt_del = True
                    if raw_line.startswith("##INFO=<ID=END"):
                        has_info_end = True
                    if raw_line.startswith("##INFO=<ID=SVTYPE"):
                        has_info_svtype = True
                    if raw_line.startswith("##INFO=<ID=SVLEN"):
                        has_info_svlen = True
                    if raw_line.startswith("##INFO=<ID=TYPE"):
                        has_info_type = True
                    if raw_line.startswith("##INFO=<ID=ZERODEPTH"):
                        has_info_zd = True
                    continue

                if raw_line.startswith("#CHROM"):
                    chrom_line = raw_line
                    chrom_fields = raw_line.rstrip("\n").split("\t")
                    if len(chrom_fields) > 9:
                        sample_count = len(chrom_fields) - 9
                    continue

                stripped = raw_line.rstrip("\n")
                if not stripped:
                    continue

                cols = stripped.split("\t")
                if len(cols) < 8:
                    continue
                chrom = cols[0]
                try:
                    pos = int(cols[1])
                except ValueError:
                    continue

                ref = cols[3]
                alt = cols[4]
                info_field = cols[7]
                existing_del_keys.update(
                    AddDeletionsToVCF._existing_deletion_keys(chrom, pos, ref, alt, info_field)
                )
                interval_start, interval_end = AddDeletionsToVCF._variant_interval(pos, ref, info_field, alt)
                occupied_intervals.setdefault(chrom, []).append((interval_start, interval_end))

                variant_rows.append((chrom, pos, record_index, stripped))
                record_index += 1

        if not chrom_line:
            raise ValueError(f"Input VCF has no #CHROM header line: {input_vcf}")

        if not has_alt_del:
            meta_lines.append('##ALT=<ID=DEL,Description="Deletion">\n')
        if not has_info_end:
            meta_lines.append('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
        if not has_info_svtype:
            meta_lines.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        if not has_info_svlen:
            meta_lines.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">\n')
        if not has_info_type:
            meta_lines.append('##INFO=<ID=TYPE,Number=A,Type=String,Description="Allele type">\n')
        if not has_info_zd:
            meta_lines.append('##INFO=<ID=ZERODEPTH,Number=0,Type=Flag,Description="Constructed from a zero-depth alignment interval with no reads spanning the gap">\n')

        added_index = len(variant_rows)
        with open(zero_depth_bed, "r") as bed_handle:
            del_index = 1
            for bed_line in bed_handle:
                fields = bed_line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                chrom = fields[0]
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError:
                    continue

                if end <= start:
                    continue

                if start == 0:
                    # For contig-start deletions, use right-anchored explicit alleles:
                    # POS=1, REF=<deleted_bases + right_anchor>, ALT=<right_anchor>.
                    ref_seq = reference_sequences.get(chrom)
                    if ref_seq is None:
                        raise ValueError(f"Contig '{chrom}' not found in reference FASTA: {reference_fasta}")
                    if end >= len(ref_seq):
                        # A complete contig deletion has no valid left or right anchor in VCF.
                        # Emit a non-standard explicit deletion over the whole contig so downstream
                        # consensus generation can still remove the sequence instead of aborting.
                        pos = 1
                        ref_allele = ref_seq
                        alt_allele = "-"
                        end_pos = len(ref_seq)
                        svlen = -len(ref_seq)
                        logger.warning(
                            f"Encoding whole-contig deletion for '{chrom}' without an anchor: "
                            f"POS={pos}, REF=<len {len(ref_allele)}>, ALT={alt_allele}, END={end_pos}, SVLEN={svlen}. "
                            f"This is not strict VCF, but preserves deletion semantics for consensus generation."
                        )
                    else:
                        pos = 1
                        deleted_seq = ref_seq[:end]
                        anchor_base = ref_seq[end]
                        ref_allele = f"{deleted_seq}{anchor_base}"
                        # TODO: This is a hack to get bcftools consensus to properly apply the deletion
                        alt_allele = "-"
                        end_pos = end
                        svlen = -end
                        log_ref_allele = ref_allele[:100]
                        log_ref_suffix = "..." if len(ref_allele) > 100 else ""
                        logger.warning(
                            f"Encoding contig-start deletion for '{chrom}' with right anchor: "
                            f"POS={pos}, REF={log_ref_allele}{log_ref_suffix}, ALT={alt_allele}, END={end_pos}, SVLEN={svlen}. "
                            f"Please ensure your downstream tools can handle this encoding."
                        )
                else:
                    # BED is 0-based half-open [start, end).
                    # Symbolic DEL uses left-anchor POS with deleted region (POS, END].
                    pos = start
                    end_pos = end
                    svlen = -(end_pos - pos)
                    ref_allele = "N"
                    alt_allele = "<DEL>"

                key = (chrom, pos, end_pos)
                if key in existing_del_keys:
                    continue
                if AddDeletionsToVCF._overlaps_any_interval(chrom, pos, end_pos, occupied_intervals):
                    logger.debug(
                        f"Skipping synthetic zero-depth DEL at {chrom}:{pos}-{end_pos} because it overlaps an existing VCF variant"
                    )
                    continue
                existing_del_keys.add(key)
                occupied_intervals.setdefault(chrom, []).append((pos, end_pos))

                info = (
                    "TYPE=INDEL;"
                    f"ZERODEPTH;SVTYPE=DEL;END={end_pos};SVLEN={svlen}"
                )
                base_cols = [chrom, str(pos), f"DEL_{del_index}", ref_allele, alt_allele, ".", "PASS", info]

                if sample_count > 0:
                    base_cols.extend(["GT", *(["1/1"] * sample_count)])

                line = "\t".join(base_cols)
                variant_rows.append((chrom, pos, added_index, line))
                added_index += 1
                del_index += 1

        unknown_rank = len(contig_rank)

        def sort_key(row: tuple[str, int, int, str]):
            chrom, pos, idx, _ = row
            rank = contig_rank.get(chrom, unknown_rank)
            chrom_tiebreak = "" if chrom in contig_rank else chrom
            return (rank, chrom_tiebreak, pos, idx)

        variant_rows.sort(key=sort_key)

        with open(output_vcf, "w") as out:
            out.writelines(meta_lines)
            out.write(chrom_line)
            for _, _, _, variant_line in variant_rows:
                out.write(f"{variant_line}\n")
