# Concrete Alignment Strategies
from pathlib import Path
from typing import List, Optional

from snippy_ng.stages.base import BaseStage, BaseOutput
from snippy_ng.dependencies import bcftools

from pydantic import Field


class PseudoAlignment(BaseStage):
    reference: Path = Field(..., description="Reference file")
    prefix: str = Field(..., description="Output file prefix")

class BcftoolsPseudoAlignmentOutput(BaseOutput):
    fasta: Path

class BcftoolsPseudoAlignment(PseudoAlignment):
    """
    Call pseudo-alignment using Bcftools consensus.
    """
    vcf_gz: Path = Field(..., description="Input VCF.gz file")

    _dependencies = [
        bcftools
    ]

    @property
    def output(self) -> BcftoolsPseudoAlignmentOutput:
        return BcftoolsPseudoAlignmentOutput(
            fasta=Path(f"{self.prefix}.pseudo.afa")
        )

    @property
    def commands(self) -> List:
        """Constructs the bcftools consensus command."""

        bcf_csq_args = ["bcftools", "consensus"]
        bcf_csq_args.extend([
            "-f", str(self.reference),
            "-o", str(self.output.fasta),
            "--mark-del", "-",
            str(self.vcf_gz),
        ]) 
        return [
            self.shell_cmd(["bcftools", "index", str(self.vcf_gz)], description="Indexing VCF file"),
            self.shell_cmd(bcf_csq_args, description="Calling consensus with bcftools"),
            self.shell_cmd(["rm", f"{self.vcf_gz}.csi"], description="Removing VCF index file")
        ]