from pathlib import Path
from typing import List
import re

from snippy_ng.stages import BaseStage, BaseOutput, TempPath
from snippy_ng.dependencies import bcftools, bedtools
from snippy_ng.logging import logger
from pydantic import Field

class VcfFilterOutput(BaseOutput):
    vcf: Path = Field(..., description="Filtered and normalized VCF file")

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

    # Keep only the tags you want; everything else is dropped.
    _keep_vcf_tags = ",".join(
        [f"^INFO/{tag}" for tag in ["TYPE", "DP", "RO", "AO", "AB"]]
        + [f"^FORMAT/{tag}" for tag in ["GT", "GQ", "DP", "RO", "AO", "QR", "QA", "GL"]]
    )

    def create_commands(self, ctx) -> List:
        """Constructs the samtools view command for filtering."""

        # Build the post-norm filter. We filter AFTER splitting/normalizing and after recomputing TYPE.
        base_filter = " && ".join([
            'ALT!="*"', 'FMT/GT="1/1"'
        ])
        commands = [
                self.shell_cmd(
                    ["cat", str(self.vcf)],
                    description="Read input VCF file",
                ),
                self.shell_cmd(
                    [
                        "bcftools",
                        "norm",
                        "-f",
                        str(self.reference),
                        "-m",
                        "-both",
                        "-Ob",
                    ],
                    description="Normalize and split multiallelic variants",
                ),
                self.shell_cmd(
                    ["bcftools", "+fill-tags", "-Ob", "-", "--", "-t", "TYPE"],
                    description="Recompute TYPE from REF/ALT",
                ),
                self.shell_cmd(
                    ["bcftools", "view", "--include", base_filter, "-"],
                    description="Filter variants after normalization and TYPE recomputation",
                ),
                self.shell_cmd(
                    ["bcftools", "+setGT", "-", "--", "-t", "a", "-n", "c:M"],
                    description="Make genotypes haploid (e.g., 1/1 -> 1)"
                ),
            ]
        if self._keep_vcf_tags:
            commands.append(
                self.shell_cmd(
                    ["bcftools", "annotate", "--remove", self._keep_vcf_tags, "-"],
                    description="Remove unnecessary VCF annotations",
                ),
            )
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
    min_qual: int = Field(60, description="Minimum QUAL score for assembly-based calling")
    # Keep only the tags you want; everything else is dropped.
    _keep_vcf_tags = None # keep all tags for assembly-based calling


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
        
        # Keep only the tags you want; everything else is dropped.
        keep_vcf_tags = ",".join(
            [f"^INFO/{tag}" for tag in ["TYPE", "DP", "RO", "AO", "AB"]]
            + [f"^FORMAT/{tag}" for tag in ["GT", "DP", "RO", "AO", "QR", "QA", "GL"]]
        ) 
        
        # Continue with the filtering pipeline
        pipeline_commands.extend([
            self.shell_cmd(
                ["bcftools", "view", "-i", 'GT="alt"'],
                description="Remove non-alt alleles"
            ),
            self.shell_cmd(
                ["bcftools", "norm", "-f", str(self.reference), "-a", "-c", "e", "-m", "-"],
                description="Normalize and left-align indels"
            ),
            self.shell_cmd(
                ["bcftools", "norm", "-aD"],
                description="Remove duplicates after normalization"
            ),
            self.shell_cmd(
                ["bcftools", "view", "--apply-filters", "PASS"],
                description="Keep only PASS variants"
            ),
            self.shell_cmd(
                ["bcftools", "filter", "-e", f'abs(ILEN)>{self.max_indel} || ALT="*"'],
                description=f"Remove indels longer than {self.max_indel}bp or sites with unobserved alleles"
            ),
            self.shell_cmd(
                ["bcftools", "+setGT", "-", "--", "-t", "a", "-n", "c:M"],
                description="Make genotypes haploid (e.g., 1/1 -> 1)"
            ),
            self.shell_cmd(
                    ["bcftools", "+fill-tags", "-", "--", "-t", "TYPE"],
                    description="Recompute TYPE from REF/ALT",
            ),
            self.shell_cmd(
                    ["bcftools", "annotate", "--remove", keep_vcf_tags, "-"],
                    description="Remove unnecessary VCF annotations",
            ),
            self.shell_cmd(
                ["bcftools", "sort"],
                description="Sort VCF"
            ),
            self.shell_cmd(
                ["bcftools", "view", "-i", 'GT="A"'],
                description="Remove non-alt alleles and output final VCF"
            ),
        ])
        
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

    




class AddDeletionstoVCFOutput(BaseOutput):
    """Output from the zero-depth deletion annotation stage."""
    vcf: Path = Field(..., description="VCF with zero-depth regions added as <DEL> blocks")
    zero_depth_bed: TempPath = Field(..., description="BED file with zero-depth regions")


class AddDeletionstoVCF(BaseStage):
    """
    Add zero-depth regions as symbolic deletion blocks to a VCF.

    This stage detects contiguous depth==0 blocks from the BAM and appends
    `<DEL>` records to the input VCF with valid ALT/INFO headers.
    """
    bam: Path = Field(..., description="Input BAM file")
    vcf: Path = Field(..., description="Input VCF file")
    reference: Path = Field(..., description="Reference FASTA file")

    _dependencies = [
        bedtools
    ]

    @property
    def output(self) -> AddDeletionstoVCFOutput:
        return AddDeletionstoVCFOutput(
            vcf=Path(f"{self.prefix}.with_del.vcf"),
            zero_depth_bed=Path(f"{self.prefix}.zerodepth.bed")
        )

    def create_commands(self, ctx) -> List:
        """Generate zero-depth BED and merge symbolic deletion blocks into VCF."""
        return [
            *self._generate_zero_depth_mask_commands(),
            self.python_cmd(
                func=self._merge_zero_depth_deletions_into_vcf,
                args=[self.vcf, self.output.zero_depth_bed, self.output.vcf, self.reference],
                description="Add zero-depth <DEL> blocks to VCF",
            ),
        ]

    def _generate_zero_depth_mask_commands(self) -> List:
        """Generate commands to create a zero-depth BED file."""
        genomecov_cmd = self.shell_cmd(
            ["bedtools", "genomecov", "-ibam", str(self.bam), "-bga"],
            description="Generate genome coverage in BED format"
        )
        awk_cmd = self.shell_cmd(
            ["awk", '$4==0 {print $1"\\t"$2"\\t"$3}'],
            description="Filter for regions with depth == 0"
        )

        return [self.shell_pipe(
            [genomecov_cmd, awk_cmd],
            output_file=self.output.zero_depth_bed,
            description="Generate zero-depth BED blocks"
        )]

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
            end_pos = AddDeletionstoVCF._extract_info_int(info_field, "END")
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
        reference_sequences = AddDeletionstoVCF._load_reference_sequences(reference_fasta)

        existing_del_keys: set[tuple[str, int, int]] = set()

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
                    AddDeletionstoVCF._existing_deletion_keys(chrom, pos, ref, alt, info_field)
                )

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
            for i, bed_line in enumerate(bed_handle, start=1):
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
                        raise ValueError(
                            f"Cannot encode start-of-contig deletion [{start},{end}) for '{chrom}': "
                            f"requires right anchor at 1-based position {end + 1}, but contig length is {len(ref_seq)}"
                        )

                    pos = 1
                    deleted_seq = ref_seq[:end]
                    anchor_base = ref_seq[end]
                    ref_allele = f"{deleted_seq}{anchor_base}"
                    # TODO: This is a hack to get bcftools consensus to properly apply the deletion
                    alt_allele = "N"
                    end_pos = end
                    svlen = -end
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
                existing_del_keys.add(key)

                info = (
                    "TYPE=INDEL;"
                    f"ZERODEPTH;SVTYPE=DEL;END={end_pos};SVLEN={svlen}"
                )
                base_cols = [chrom, str(pos), f"DEL_{i}", ref_allele, alt_allele, ".", "PASS", info]

                if sample_count > 0:
                    base_cols.extend(["GT", *(["1"] * sample_count)])

                line = "\t".join(base_cols)
                variant_rows.append((chrom, pos, added_index, line))
                added_index += 1

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