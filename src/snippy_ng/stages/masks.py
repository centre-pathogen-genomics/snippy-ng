from pathlib import Path
from typing import List

from snippy_ng.stages.base import BaseStage, BaseOutput
from snippy_ng.dependencies import bedtools, bcftools

from pydantic import Field


class Mask(BaseStage):
    bam: Path = Field(..., description="Input BAM file")
    prefix: str = Field(..., description="Output file prefix")
    filter: str = Field(..., description="Filter for coverage, e.g. '==0' for zero coverage regions")
    reference: Path = Field(..., description="Reference FASTA file")


class DepthMaskOutPut(BaseOutput):
    bed: Path

class DepthMask(Mask):
    """
    Create a BED file of regions from a BAM file based on coverage filter
    """
    _dependencies = [
        bcftools
    ]

    @property
    def output(self) -> DepthMaskOutPut:
        return DepthMaskOutPut(
            bed=Path(f"{self.prefix}.{self.name.lower()}.bed")
        )

    @property
    def commands(self) -> List:
        """Constructs the bedtools genomecov command to create a BED file of regions based on coverage filter."""
        mask_zero_cov_cmd = self.shell_cmd(
            ["bedtools", "genomecov", "-ibam", str(self.bam), "-bga"],
            description="Generate genome coverage in BED format"
        )
        awk_cmd = self.shell_cmd(
            ["awk", f'$4{self.filter} {{print $1"\t"$2"\t"$3}}'],
            description=f"Apply coverage filter ({self.filter}) to create mask file"
        )
        return [self.shell_pipeline(
            [mask_zero_cov_cmd, awk_cmd], 
            output_file=self.output.bed, 
            description=f"Create mask for regions based on coverage filter ({self.filter})")
        ]

class ZeroDepthMask(DepthMask):
    """
    Create a BED file of zero-coverage regions from a BAM file
    """
    filter: str = Field("==0", description="Filter for coverage, e.g. '==0' for zero coverage regions")

class MinDepthMask(DepthMask):
    """
    Create a BED file of low-coverage regions from a BAM file
    """
    filter: str = Field("<10", description="Filter for coverage, e.g. '<10' for regions with coverage less than 10")

class ApplyMaskOutput(BaseOutput):
    fasta: Path

class ApplyMask(BaseStage):
    """
    Apply a BED mask to a reference FASTA file, replacing masked regions with 'N's.
    """
    bed: Path = Field(..., description="Input mask file")
    prefix: str = Field(..., description="Output file prefix")
    reference: Path = Field(..., description="Reference file")
    char: str = Field("N", description="Character to use for masking (default: 'N')")

    _dependencies = [
        bedtools
    ]

    @property
    def output(self) -> ApplyMaskOutput:
        return ApplyMaskOutput(
            fasta=Path(f"{self.prefix}.masked.fasta")
        )

    @property
    def commands(self) -> List:
        """Constructs the bedtools maskfasta command."""
        mask_fasta_cmd = self.shell_cmd([
            "bedtools", "maskfasta",
            "-fi", str(self.reference),
            "-bed", str(self.bed),
            "-fo", str(self.output.fasta),
            "-mc", self.char
        ], description="Apply BED mask to FASTA file")
        return [mask_fasta_cmd] 