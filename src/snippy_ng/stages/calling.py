# Concrete Alignment Strategies
from pathlib import Path
from typing import List

from snippy_ng.stages.base import BaseStage, BaseOutput
from snippy_ng.dependencies import freebayes, bcftools

from pydantic import Field

# Define the base Pydantic model for alignment parameters
class Caller(BaseStage):
    reference: Path = Field(..., description="Reference file",)
    prefix: str = Field(..., description="Output file prefix")

class FreebayesCallerOutput(BaseOutput):
    raw_vcf: str
    filter_vcf: str

class FreebayesCaller(Caller):
    """
    Call variants using Freebayes.
    """
    bam: Path = Field(..., description="Input BAM file")
    fbopt: str = Field("", description="Additional Freebayes options")
    mincov: int = Field(10, description="Minimum site depth for calling alleles")
    minfrac: float = Field(0.0, description="Minimum proportion for variant evidence (0=AUTO)")
    minqual: float = Field(100.0, description="Minimum quality in VCF column 6")

    _dependencies = [
        freebayes,
        bcftools
    ]

    @property
    def output(self) -> FreebayesCallerOutput:
        return FreebayesCallerOutput(
                raw_vcf=self.prefix + ".raw.vcf",
                filter_vcf=self.prefix + ".filt.vcf"
            )

    @property
    def commands(self) -> List[str]:
        """Constructs the Freebayes variant calling commands."""
        bcf_filter = f'FMT/GT="1/1" && QUAL>={self.minqual} && FMT/DP>={self.mincov} && (FMT/AO)/(FMT/DP)>={self.minfrac}'
        keep_vcf_tags = ",".join([
                f"^INFO/{tag}" for tag in ["TYPE", "DP", "RO", "AO", "AB"]
            ] + [
                f"^FORMAT/{tag}" for tag in ["GT", "DP", "RO", "AO", "QR", "QA", "GL"]
            ])
        generate_regions_cmd = f"fasta_generate_regions.py {self.reference}.fai 1000000 > {self.reference}.txt"
        freebayes_cmd = f"freebayes-parallel {self.reference}.txt {self.cpus} {self.fbopt} -f {self.reference} {self.bam} > {self.output.raw_vcf}"
        bcftools_cmd = f"bcftools view --include '{bcf_filter}' {self.output.raw_vcf} | bcftools norm -f {self.reference} - | bcftools annotate --remove '{keep_vcf_tags}' > {self.output.filter_vcf}"
        return [generate_regions_cmd, freebayes_cmd, bcftools_cmd]