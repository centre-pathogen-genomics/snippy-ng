# Concrete Alignment Strategies
from pathlib import Path
from typing import List

from snippy_ng.stages.base import BaseStage, BaseOutput
from snippy_ng.dependencies import bcftools

from pydantic import Field


class Caller(BaseStage):
    reference: Path = Field(..., description="Reference file",)
    variants: Path = Field(..., description="Input VCF file",)
    features: Path = Field(..., description="Input features file")
    prefix: str = Field(..., description="Output file prefix")

class BcftoolsConsequencesCallerOutput(BaseOutput):
    annotated_vcf: Path

class BcftoolsConsequencesCaller(Caller):
    """
    Call consequences using Bcftools csq.
    """

    _dependencies = [
        bcftools
    ]

    @property
    def output(self) -> BcftoolsConsequencesCallerOutput:
        return BcftoolsConsequencesCallerOutput(
            annotated_vcf=Path(f"{self.prefix}.vcf")
        )

    @property
    def commands(self) -> List[str]:
        """Constructs the bcftools csq command."""
        # check if features file exists and is not empty
        features_found = True
        with open(self.features, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    break
            else:
                features_found = False
        
        if not features_found:
            cmd = f"cp {self.variants} {self.output.annotated_vcf}"
            return [cmd]
        
        bcf_csq_cmd = (
            f"bcftools csq -f {self.reference} "
            f"-g {self.features} "
            f"-o {self.output.annotated_vcf} "
            f"{self.variants}"
        )
        return [bcf_csq_cmd]