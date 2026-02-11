# Concrete Alignment Strategies
from pathlib import Path
from typing import List, Annotated, Optional

from snippy_ng.stages.base import BaseStage, BaseOutput
from snippy_ng.dependencies import freebayes, bcftools, bedtools, paftools, clair3
from snippy_ng.logging import logger

from pydantic import Field, AfterValidator


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
    vcf: Path

class FreebayesCallerOutput(BaseCallerOutput):
    vcf: Path
    regions: str


class FreebayesCaller(Caller):
    """
    Call variants using Freebayes.
    """

    bam: Path = Field(..., description="Input BAM file")
    fbopt: str = Field("", description="Additional Freebayes options")
    mincov: int = Field(10, description="Minimum site depth for calling alleles")
    minfrac: float = Field(
        0.05, description="Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position"
    )
    exclude_insertions: bool = Field(
        True,
        description="Exclude insertions from variant calls so the pseudo-alignment remains the same length as the reference",
    )

    _dependencies = [freebayes, bcftools]

    @property
    def output(self) -> FreebayesCallerOutput:
        return FreebayesCallerOutput(
            vcf=self.prefix + ".raw.vcf",
            regions=str(self.reference) + ".txt",
        )

    @property
    def commands(self) -> List:
        """Constructs the Freebayes variant calling and postprocessing commands."""

        # 1) Regions for parallel FreeBayes
        generate_regions_cmd = self.shell_cmd(
            ["fasta_generate_regions.py", str(self.reference_index), "202106"],
            description="Generate genomic regions for parallel variant calling",
        )
        generate_regions_pipeline = self.shell_pipeline(
            commands=[generate_regions_cmd],
            description="Generate regions file for parallel processing",
            output_file=Path(self.output.regions),
        )

        # 2) FreeBayes parallel call
        freebayes_cmd_parts = [
            "freebayes-parallel",
            str(self.output.regions),
            str(self.cpus),
            "--ploidy", "2",
            "--min-alternate-count", "2",
            "--min-alternate-fraction", str(self.minfrac),
            "--min-coverage", str(self.mincov),
            "--min-repeat-entropy", "1.0",
            "--min-base-quality", "13",
            "--min-mapping-quality", "60",
            "--strict-vcf",
        ]
        if self.fbopt:
            import shlex

            freebayes_cmd_parts.extend(shlex.split(self.fbopt))
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

    @property
    def commands(self) -> List:
        """Constructs the Freebayes variant calling and postprocessing commands."""

        # Regions for parallel FreeBayes
        generate_regions_cmd = self.shell_cmd(
            ["fasta_generate_regions.py", str(self.reference_index), "202106"],
            description="Generate genomic regions for parallel variant calling",
            output_file=Path(self.output.regions),
        )
        # FreeBayes parallel call
        freebayes_cmd_parts = [
            "freebayes-parallel",
            str(self.output.regions),
            str(self.cpus),
            "--haplotype-length", "-1",
            "--min-mapping-quality", "10",
            "--min-base-quality", "10",
            "--ploidy", "2",
            "--min-coverage", str(self.mincov),
            "--min-alternate-fraction", str(self.minfrac),

        ]
        if self.fbopt:
            import shlex

            freebayes_cmd_parts.extend(shlex.split(self.fbopt))
        freebayes_cmd_parts.extend(["-f", str(self.reference), str(self.bam)])

        freebayes_cmd = self.shell_cmd(
            freebayes_cmd_parts,
            description="Call variants with FreeBayes in parallel",
            output_file=Path(self.output.vcf),
        )
        
        return [generate_regions_cmd, freebayes_cmd]

class PAFCallerOutput(BaseCallerOutput):
    vcf: Path
    aln_bed: Path
    missing_bed: Path
    annotations_file: Path
    annotations_file_index: Path


class PAFCaller(Caller):
    """
    Call variants from PAF alignments using paftools.js.
    """

    paf: Path = Field(..., description="Input PAF file")
    ref_dict: Path = Field(..., description="Reference FASTA dictionary file")
    mapq: Optional[int] = Field(
        15, description="Minimum mapping quality for variant calling"
    )
    alen: Optional[int] = Field(
        50, description="Minimum alignment length for variant calling"
    )

    _dependencies = [bedtools, bcftools, paftools]

    @property
    def output(self) -> PAFCallerOutput:
        return PAFCallerOutput(
            vcf=Path(f"{self.prefix}.raw.vcf"),
            aln_bed=Path(f"{self.prefix}.aln.bed"),
            missing_bed=Path(f"{self.prefix}.missing.bed"),
            annotations_file=Path(f"{self.prefix}.annotations.gz"),
            annotations_file_index=Path(f"{self.prefix}.annotations.gz.tbi"),
        )

    @property
    def commands(self) -> List:
        """Constructs the PAF processing and BED generation commands."""

        # 4) Convert PAF to merged aligned reference intervals (BED)
        # Keep primary or pseudo-primary hits: tp:A:P or tp:A:I
        paf_to_pipeline = self.shell_pipeline(
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
        compute_missing_bed_cmd = self.shell_cmd(
            [
                "bedtools",
                "complement",
                "-g",
                str(self.ref_dict),
                "-i",
                str(self.output.aln_bed),
            ],
            description="Compute unaligned (missing) reference regions",
            output_file=self.output.missing_bed,
        )

        # variant calling
        paftools_cmd = self.shell_cmd(
                    [
                        "paftools.js",
                        "call",
                        "-q",
                        str(self.mapq),
                        "-L",
                        str(self.alen),
                        "-l",
                        str(self.alen),
                        "-s",
                        self.prefix,
                        "-f",
                        str(self.reference),
                        str(self.paf),
                    ],
                    description="Call variants from PAF using paftools.js",
                    output_file=self.output.vcf.with_suffix(".tmp"),
                )

        create_annotation_file_pipeline = self.shell_pipeline(
            commands=[
                self.shell_cmd(
                    [
                        "bcftools",
                        "query",
                        str(paftools_cmd.output_file),
                        "-f",
                        "%CHROM\\t%POS\\t1\\t1\\n",
                    ],
                    description="Extract variant positions from VCF",
                ),
                self.shell_cmd(
                    ["bgzip", "-c"],
                    description="Compress annotation file with bgzip",
                ),
            ],
            description="Compress annotation file with bgzip",
            output_file=self.output.annotations_file,
        )
        index_cmd = self.shell_cmd(
            [
                "tabix",
                "-s1",
                "-b2",
                "-e2",
                str(self.output.annotations_file),
            ],
            description="Index annotation file with tabix",
        )
        annotations_cmd = self.shell_cmd(
                    [
                        "bcftools",
                        "annotate",
                        "-H", '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                        "-H", '##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count for each ALT">',
                        "-c", "CHROM,POS,FMT/DP:=FORMAT/DP,FMT/AO:=FORMAT/AO",
                        "-a", str(self.output.annotations_file),
                        str(paftools_cmd.output_file),
                    ],
                    description="Insert FORMAT header lines for DP and AO",
                    output_file=self.output.vcf,
                )

        return [
            paf_to_pipeline,
            compute_missing_bed_cmd,
            paftools_cmd,
            create_annotation_file_pipeline,
            index_cmd,
            annotations_cmd,
        ]

class Clair3CallerOutput(BaseCallerOutput):
    vcf: Path


class Clair3Caller(Caller):
    """
    Call variants using Clair3.
    """
    bam: Path = Field(..., description="Input BAM file")
    clair3_model: Path = Field(..., description="Absolute path to Clair3 model")
    platform: str = Field("ont", description="Sequencing platform (e.g., ont, hifi)")
    fast_mode: bool = Field(True, description="Enable fast mode for Clair3")

    _dependencies = [clair3]

    @property
    def output(self) -> Clair3CallerOutput:
        return Clair3CallerOutput(
            vcf=Path(f"{self.prefix}.raw.vcf"),
        )

    @property
    def commands(self) -> List:
        """Constructs the Clair3 variant calling commands."""

        clair3_cmd = self.shell_cmd(
            [
                "run_clair3.sh",
                f"--model_path={self.clair3_model.resolve()}",
                f"--bam_fn={str(self.bam.resolve())}",
                f"--ref_fn={str(self.reference.resolve())}",
                f"--threads={str(self.cpus)}",
                f"--output={Path(self.prefix + '_clair3_out').resolve()}",
                f"--platform={self.platform}",
                "--include_all_ctgs",
                "--haploid_precise",
                "--no_phasing_for_fa",
                "--enable_long_indel",
            ],
            description="Call variants with Clair3",
        )
        if self.fast_mode:
            clair3_cmd.command.append("--fast_mode")
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
        cleanup_cmd = self.shell_cmd(
            [
                "rm",
                "-rf",
                f"{self.prefix}_clair3_out",
            ],
            description="Clean up Clair3 output directory",
        )

        return [clair3_cmd, unzip_cmd, move_vcf_cmd, cleanup_cmd] 
