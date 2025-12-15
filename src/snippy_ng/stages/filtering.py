from pathlib import Path
from typing import List, Optional
from snippy_ng.stages.base import BaseStage
from snippy_ng.dependencies import samtools, bcftools
from pydantic import Field, field_validator, BaseModel


class BamFilterOutput(BaseModel):
    bam: Path
    bam_index: Path


class BamFilter(BaseStage):
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
    def output(self) -> BamFilterOutput:
        filtered_bam = f"{self.prefix}.filtered.bam"
        return BamFilterOutput(
            bam=filtered_bam,
            bam_index=f"{filtered_bam}.bai"
        )
    
    def build_filter_command(self):
        """Constructs the samtools view command for filtering."""
        cmd_parts = ["samtools", "view", "-b"]
        
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
            description=f"Filter BAM file with MAPQ>={self.min_mapq}, flags={self.exclude_flags}"
        )
        
        return self.shell_pipeline(
            commands=[filter_cmd],
            description="Filter BAM alignments",
            output_file=Path(self.output.bam)
        )
    
    def build_index_command(self):
        """Returns the samtools index command."""
        return self.shell_cmd([
            "samtools", "index", str(self.output.bam)
        ], description=f"Index filtered BAM file: {self.output.bam}")
    
    @property
    def commands(self) -> List:
        """Constructs the filtering commands."""
        filter_cmd = self.build_filter_command()
        index_cmd = self.build_index_command()
        return [filter_cmd, index_cmd]


class BamFilterByRegion(BamFilter):
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


class BamFilterByQuality(BamFilter):
    """
    Filter BAM file based on mapping quality and alignment flags.
    """
    
    min_mapq: int = Field(30, description="Minimum mapping quality (higher than default)")
    exclude_flags: int = Field(3844, description="SAM flags to exclude (default + supplementary)")


class BamFilterProperPairs(BamFilter):
    """
    Filter BAM file to include only properly paired reads.
    """
    
    include_flags: int = Field(2, description="Include only properly paired reads")
    exclude_flags: int = Field(1796, description="Exclude unmapped, secondary, qcfail, duplicate")


class VcfFilterOutput(BaseModel):
    vcf: Path


class VcfFilter(BaseStage):
    """
    Filter VCF files using Samtools to remove unwanted variants.
    """

    vcf: Path = Field(..., description="Input VCF file to filter")
    reference: Path = Field(..., description="Reference FASTA file")
    min_qual: int = Field(100, description="Minimum QUAL score")
    min_depth: int = Field(1, description="Minimum site depth for calling alleles")
    min_frac: float = Field(0, description="Minimum proportion for calling alt allele")
    exclude_insertions: bool = Field(
        True,
        description="Exclude insertions from variant calls so the pseudo-alignment remains the same length as the reference",
    )

    _dependencies = [bcftools]

    @property
    def output(self) -> VcfFilterOutput:
        filtered_vcf = f"{self.prefix}.filtered.vcf"
        return VcfFilterOutput(vcf=filtered_vcf)

    @property
    def commands(self) -> List:
        """Constructs the samtools view command for filtering."""

        # Build the post-norm filter. We filter AFTER splitting/normalizing and after recomputing TYPE.
        base_filter = (
            f'FMT/GT="1/1" && QUAL>={self.min_qual} && FMT/DP>={self.min_depth} '
            f'&& (FMT/AO)/(FMT/DP)>={self.min_frac} && N_ALT=1 && ALT!="*"'
        )
        if self.exclude_insertions:
            base_filter += " && strlen(ALT) <= strlen(REF)"

        # Keep only the tags you want; everything else is dropped.
        keep_vcf_tags = ",".join(
            [f"^INFO/{tag}" for tag in ["TYPE", "DP", "RO", "AO", "AB"]]
            + [f"^FORMAT/{tag}" for tag in ["GT", "DP", "RO", "AO", "QR", "QA", "GL"]]
        )

        bcftools_pipeline = self.shell_pipeline(
            commands=[
                self.shell_cmd(
                    [
                        "bcftools",
                        "norm",
                        "-f",
                        str(self.reference),
                        "-m",
                        "-both",
                        str(self.vcf),
                    ],
                    description="Normalize and split multiallelic variants",
                ),
                self.shell_cmd(
                    ["bcftools", "+fill-tags", "-", "--", "-t", "TYPE"],
                    description="Recompute TYPE from REF/ALT",
                ),
                self.shell_cmd(
                    ["bcftools", "view", "--include", base_filter, "-"],
                    description="Filter variants after normalization and TYPE recomputation",
                ),
                self.shell_cmd(
                    ["bcftools", "annotate", "--remove", keep_vcf_tags, "-"],
                    description="Remove unnecessary VCF annotations",
                ),
            ],
            description="Normalize, recompute TYPE, filter, and annotate variants",
            output_file=Path(self.output.vcf),
        )
        return [bcftools_pipeline]


class VcfFilterLong(BaseStage):
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
    
    vcf: Path = Field(..., description="Input VCF file to filter")
    reference: Path = Field(..., description="Reference FASTA file")
    reference_index: Path = Field(..., description="Reference FASTA index file (.fai)")
    
    min_qual: int = Field(100, description="Minimum QUAL score")
    min_depth: int = Field(1, description="Minimum site depth for calling alleles")
    max_indel: int = Field(10000, description="Maximum indel length to keep")
    
    _dependencies = [bcftools]
    
    @property
    def output(self) -> VcfFilterOutput:
        filtered_vcf = f"{self.prefix}.filtered.vcf"
        return VcfFilterOutput(vcf=filtered_vcf)
    
    
    @property
    def commands(self) -> List:
        """Constructs the bcftools pipeline for long-read variant filtering."""
        
        # Create temp files for contigs and header
        contigs_file = f"{self.prefix}.contigs.txt"
        header_file = f"{self.prefix}.header.txt"
        
        # Step 1: Generate contig lines from faidx
        create_contigs = self.shell_cmd(
            ["awk", '{print "##contig=<ID="$1",length="$2">"}', str(self.reference_index)],
            description="Generate contig lines from reference index",
        )
        create_contigs_pipeline = self.shell_pipeline(
            commands=[create_contigs],
            description="Create contig definitions",
            output_file=Path(contigs_file)
        )
        
        # Step 2: Create new header with all contigs
        # This combines: bcftools view -h | grep -v "^##contig=" | sed -e "3r $contigs"
        create_header_cmd = f'bcftools view -h {self.vcf} | grep -v "^##contig=" | sed -e "3r {contigs_file}" > {header_file}'
        create_header = self.shell_cmd(
            ["bash", "-c", create_header_cmd],
            description="Create VCF header with all contigs"
        )
        
        # Step 3: Build the main filtering pipeline
        pipeline_commands = [
            self.shell_cmd(
                ["bcftools", "reheader", "-h", header_file, str(self.vcf)],
                description="Replace VCF header with new header containing all contigs"
            ),
        ]
        
        
        # Continue with the filtering pipeline
        pipeline_commands.extend([
            self.shell_cmd(
                ["bcftools", "view", "-i", 'GT="alt"'],
                description="Remove non-alt alleles"
            ),
            self.shell_cmd(
                ["bcftools", "view", "-e", 'ALT="."'],
                description="Remove sites with no alt allele (NanoCaller bug fix)"
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
                ["bcftools", "sort"],
                description="Sort VCF"
            ),
            self.shell_cmd(
                ["bcftools", "view", "-i", 'GT="A"'],
                description="Remove non-alt alleles and output final VCF"
            ),
        ])
        
        main_pipeline = self.shell_pipeline(
            commands=pipeline_commands,
            description="Long-read variant filtering pipeline",
            output_file=Path(self.output.vcf)
        )
        
        # Step 4: Cleanup temp files
        cleanup_cmd = self.shell_cmd(
            ["rm", "-f", contigs_file, header_file],
            description="Remove temporary files"
        )
        
        return [create_contigs_pipeline, create_header, main_pipeline, cleanup_cmd]

    