# Concrete Alignment Strategies
from pathlib import Path
from typing import List, Annotated, Optional

from snippy_ng.stages.base import BaseStage, BaseOutput
from snippy_ng.dependencies import freebayes, bcftools, bedtools, paftools

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
    vcf: str
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
            freebayes_cmd_parts,
            description="Call variants with FreeBayes in parallel",
        )
        freebayes_pipeline = self.shell_pipeline(
            commands=[freebayes_cmd],
            description="FreeBayes variant calling",
            output_file=Path(self.output.vcf),
        )

        return [generate_regions_pipeline, freebayes_pipeline]


class PAFCallerOutput(BaseOutput):
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
        )
        compute_missing_bed_pipeline = self.shell_pipeline(
            commands=[compute_missing_bed_cmd],
            description="Generate missing regions BED file",
            output_file=self.output.missing_bed,
        )

        # variant calling
        paftools_pipeline = self.shell_pipeline(
            commands=[
                self.shell_cmd(
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
                ),
            ],
            description="Variant calling from PAF",
            output_file=self.output.vcf.with_suffix(".tmp"),
        )

        create_annotation_file_pipeline = self.shell_pipeline(
            commands=[
                self.shell_cmd(
                    [
                        "bcftools",
                        "query",
                        str(paftools_pipeline.output_file),
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
        annotations_pipeline = self.shell_pipeline(
            commands=[
                self.shell_cmd(
                    [
                        "bcftools",
                        "annotate",
                        "-H", '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                        "-H", '##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count for each ALT">',
                        "-c", "CHROM,POS,FMT/DP:=FORMAT/DP,FMT/AO:=FORMAT/AO",
                        "-a", str(self.output.annotations_file),
                        str(paftools_pipeline.output_file),
                    ],
                    description="Insert FORMAT header lines for DP and AO",
                ),
            ],
            description="Annotate VCF with DP and AO headers",
            output_file=self.output.vcf,
        )

        return [
            paf_to_pipeline,
            compute_missing_bed_pipeline,
            paftools_pipeline,
            create_annotation_file_pipeline,
            index_cmd,
            annotations_pipeline,
        ]
