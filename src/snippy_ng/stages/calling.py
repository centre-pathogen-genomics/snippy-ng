# Concrete Alignment Strategies
from pathlib import Path
from typing import List, Annotated

from snippy_ng.stages.base import BaseStage, BaseOutput
from snippy_ng.dependencies import freebayes, bcftools

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


class FreebayesCallerOutput(BaseOutput):
    raw_vcf: str
    filter_vcf: str
    regions: str


class FreebayesCaller(Caller):
    """
    Call variants using Freebayes.
    """

    bam: Path = Field(..., description="Input BAM file")
    fbopt: str = Field("", description="Additional Freebayes options")
    mincov: int = Field(10, description="Minimum site depth for calling alleles")
    minfrac: float = Field(
        0.05, description="Minimum proportion for variant evidence (0=AUTO)"
    )
    mindepth: float = Field(0, description="Minimum proportion for calling alt allele")
    minqual: float = Field(100.0, description="Minimum quality in VCF column 6")
    exclude_insertions: bool = Field(
        True,
        description="Exclude insertions from variant calls so the pseudo-alignment remains the same length as the reference",
    )

    _dependencies = [freebayes, bcftools]

    @property
    def output(self) -> FreebayesCallerOutput:
        return FreebayesCallerOutput(
            raw_vcf=self.prefix + ".raw.vcf",
            filter_vcf=self.prefix + ".filt.vcf",
            regions=str(self.reference) + ".txt",
        )

    @property
    def commands(self) -> List:
        """Constructs the Freebayes variant calling commands."""
        bcf_filter = f'FMT/GT="1/1" && QUAL>={self.minqual} && FMT/DP>={self.mincov} && (FMT/AO)/(FMT/DP)>={self.mindepth}'
        if self.exclude_insertions:
            bcf_filter += ' && INFO/TYPE!="ins"'
        keep_vcf_tags = ",".join(
            [f"^INFO/{tag}" for tag in ["TYPE", "DP", "RO", "AO", "AB"]]
            + [f"^FORMAT/{tag}" for tag in ["GT", "DP", "RO", "AO", "QR", "QA", "GL"]]
        )

        # Generate regions command with output redirection
        generate_regions_cmd = self.shell_cmd(
            ["fasta_generate_regions.py", str(self.reference_index), "202106"],
            description="Generate genomic regions for parallel variant calling",
        )

        generate_regions_pipeline = self.shell_pipeline(
            commands=[generate_regions_cmd],
            description="Generate regions file for parallel processing",
            output_file=Path(self.output.regions),
        )

        # Freebayes command with output redirection
        freebayes_cmd_parts = [
            "freebayes-parallel",
            str(self.output.regions),
            str(self.cpus),
            "-p",
            "2",
            "-P",
            "0",
            "-C",
            "2",
            "-F",
            str(self.minfrac),
            "--min-coverage",
            str(self.mincov),
            "--min-repeat-entropy",
            "1.0",
            "-q",
            "13",
            "-m",
            "60",
            "--strict-vcf",
        ]
        if self.fbopt:
            import shlex

            freebayes_cmd_parts.extend(shlex.split(self.fbopt))
        freebayes_cmd_parts.extend(["-f", str(self.reference), str(self.bam)])

        freebayes_cmd = self.shell_cmd(
            freebayes_cmd_parts, description="Call variants with Freebayes in parallel"
        )
        freebayes_pipeline = self.shell_pipeline(
            commands=[freebayes_cmd],
            description="Freebayes variant calling",
            output_file=Path(self.output.raw_vcf),
        )

        # BCFtools pipeline: view | norm | annotate
        bcftools_view_cmd = self.shell_cmd(
            ["bcftools", "view", "--include", bcf_filter, str(self.output.raw_vcf)],
            description="Filter variants by quality and depth",
        )

        bcftools_norm_cmd = self.shell_cmd(
            ["bcftools", "norm", "-a", "-f", str(self.reference), "-"],
            description="Normalize variants",
        )

        bcftools_annotate_cmd = self.shell_cmd(
            ["bcftools", "annotate", "--remove", keep_vcf_tags],
            description="Remove unnecessary VCF annotations",
        )

        bcftools_pipeline = self.shell_pipeline(
            commands=[bcftools_view_cmd, bcftools_norm_cmd, bcftools_annotate_cmd],
            description="Filter, normalize and annotate variants",
            output_file=Path(self.output.filter_vcf),
        )

        return [generate_regions_pipeline, freebayes_pipeline, bcftools_pipeline]
