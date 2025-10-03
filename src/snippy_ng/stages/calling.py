# Concrete Alignment Strategies
from pathlib import Path
from typing import List
import shlex

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
        reference_fai = shlex.quote(str(self.reference) + ".fai")
        reference_txt = shlex.quote(str(self.reference) + ".txt")
        generate_regions_cmd = f"fasta_generate_regions.py {reference_fai} 1000000 > {reference_txt}"
        freebayes_cmd = f"freebayes-parallel {reference_txt} {self.cpus} {self.fbopt} -f {shlex.quote(str(self.reference))} {shlex.quote(str(self.bam))} > {shlex.quote(self.output.raw_vcf)}"
        bcftools_cmd = f"bcftools view --include '{bcf_filter}' {shlex.quote(self.output.raw_vcf)} | bcftools norm -f {shlex.quote(str(self.reference))} - | bcftools annotate --remove '{keep_vcf_tags}' > {shlex.quote(self.output.filter_vcf)}"
        return [generate_regions_cmd, freebayes_cmd, bcftools_cmd]