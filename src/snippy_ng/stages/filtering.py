from pathlib import Path
from typing import List, Optional
from snippy_ng.stages import BaseStage, BaseOutput, ShellCommand
from snippy_ng.dependencies import samtools, bcftools
from snippy_ng.logging import logger
from pydantic import Field, field_validator


class SamtoolsFilterOutput(BaseOutput):
    cram: Path
    stats: Path = Field(..., description="Flagstat output file for the CRAM file")


class SamtoolsFilter(BaseStage):
    """
    Filter BAM files using Samtools to remove unwanted alignments.
    """
    
    bam: Path = Field(..., description="Input BAM file to filter")
    min_mapq: int = Field(20, description="Minimum mapping quality")
    exclude_flags: int = Field(1796, description="SAM flags to exclude (default: unmapped, secondary, qcfail, duplicate)")
    include_flags: Optional[int] = Field(None, description="SAM flags to include")
    regions: Optional[str] = Field(None, description="Regions to include (BED file or region string)")
    additional_filters: str = Field("", description="Additional samtools view options")
    
    _dependencies = [samtools]
    
    @property
    def output(self) -> SamtoolsFilterOutput:
        filtered_bam = f"{self.prefix}.filtered.cram"
        return SamtoolsFilterOutput(
            cram=filtered_bam,
            stats=f"{filtered_bam}.flagstat.txt"
        )
    
    def build_filter_command(self) -> ShellCommand:
        """Constructs the samtools view command for filtering."""
        cmd_parts = [
            "samtools",
            "view",
            "-O",
            "cram,embed_ref=2",
            "-o",
            str(self.output.cram),
        ]
        
        # Add threading
        if self.cpus > 1:
            cmd_parts.extend(["--threads", str(self.cpus - 1)])
        
        # Add mapping quality filter
        if self.min_mapq > 0:
            cmd_parts.extend(["-q", str(self.min_mapq)])
        
        # Add flag filters
        if self.exclude_flags:
            cmd_parts.extend(["-F", str(self.exclude_flags)])
        
        if self.include_flags is not None:
            cmd_parts.extend(["-f", str(self.include_flags)])
        
        # Add regions if specified as BED file
        if self.regions and Path(self.regions).exists():
            cmd_parts.extend(["-L", str(self.regions)])
        
        # Add additional filters (split if it contains spaces)
        if self.additional_filters:
            import shlex
            cmd_parts.extend(shlex.split(self.additional_filters))
        
        # Add input file
        cmd_parts.append(str(self.bam))
        
        # Add region string if not a file
        if self.regions and not Path(self.regions).exists():
            cmd_parts.append(str(self.regions))
        
        filter_cmd = self.shell_cmd(
            command=cmd_parts,
            description=f"Filter CRAM file with MAPQ>={self.min_mapq}, flags={self.exclude_flags}",
        )
        
        return filter_cmd
    
    def build_index_command(self):
        """Returns the samtools index command."""
        return self.shell_cmd([
            "samtools", "index", str(self.output.cram)
        ], description=f"Index filtered CRAM file: {self.output.cram}")
    
    def build_flagstat_command(self):
        """Returns the samtools flagstat command."""
        return self.shell_cmd(
            ["samtools", "flagstat", str(self.output.cram)],
            description=f"Generate alignment statistics for {self.output.cram}",
            output_file=Path(self.output.stats),
        )

    def test_mapped_reads_not_zero(self):
        """Test that the BAM file contains mapped reads."""
        flagstat_path = self.output.stats
        if not flagstat_path.exists():
            raise FileNotFoundError(
                f"Expected flagstat output file {flagstat_path} was not created"
            )

        with open(flagstat_path) as f:
            for line in f:
                if "mapped (" in line:
                    mapped_reads = int(line.split()[0])
                    if mapped_reads == 0:
                        raise ValueError(
                            "No reads were mapped in BAM file after filtering. Did you use the correct reference?"
                        )
                    break

    def create_commands(self, ctx) -> List:
        """Constructs the filtering commands."""
        filter_cmd = self.build_filter_command()
        index_cmd = self.build_index_command()
        stats_cmd = self.build_flagstat_command()
        return [filter_cmd, index_cmd, stats_cmd]


class SamtoolsFilterByRegion(SamtoolsFilter):
    """
    Filter BAM file to include only alignments in specified regions.
    """
    
    regions: str = Field(..., description="Regions to include (BED file or region string)")
    
    @field_validator("regions")
    @classmethod
    def validate_regions(cls, v):
        if not v:
            raise ValueError("Regions must be specified")
        return v


class SamtoolsFilterByQuality(SamtoolsFilter):
    """
    Filter BAM file based on mapping quality and alignment flags.
    """
    
    min_mapq: int = Field(30, description="Minimum mapping quality (higher than default)")
    exclude_flags: int = Field(3844, description="SAM flags to exclude (default + supplementary)")


class SamtoolsFilterProperPairs(SamtoolsFilter):
    """
    Filter BAM file to include only properly paired reads.
    """
    
    include_flags: int = Field(2, description="Include only properly paired reads")
    exclude_flags: int = Field(1796, description="Exclude unmapped, secondary, qcfail, duplicate")


class VcfFilterOutput(BaseOutput):
    vcf: Path

class VcfFilter(BaseStage):
    vcf: Path = Field(..., description="Input VCF file to filter")
    reference: Path = Field(..., description="Reference FASTA file")
    min_qual: int = Field(100, description="Minimum QUAL score")
    min_depth: int = Field(1, description="Minimum site depth for calling alleles")

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
    min_frac: float = Field(0, description="Minimum proportion for calling alt allele")
    
    # Keep only the tags you want; everything else is dropped.
    _keep_vcf_tags = ",".join(
        [f"INFO/{tag}" for tag in ["TYPE", "DP", "RO", "AO", "AB"]]
        + [f"FORMAT/{tag}" for tag in ["GT", "DP", "RO", "AO", "QR", "QA", "GL"]]
    )

    def create_commands(self, ctx) -> List:
        """Constructs the samtools view command for filtering."""

        # Build the post-norm filter. We filter AFTER splitting/normalizing and after recomputing TYPE.
        base_filter = (
            f'FMT/GT="1/1" && QUAL>={self.min_qual} && FMT/DP>={self.min_depth} '
            f'&& (FMT/AO)/(FMT/DP)>={self.min_frac} && N_ALT=1 && ALT!="*"'
        )
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
                ["bcftools", "view", "--include", f'QUAL>={self.min_qual} && FMT/DP>={self.min_depth}'],
                description="Filter variants based on QUAL, depth, and allele fraction"
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

    
