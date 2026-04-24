# Concrete Alignment Strategies
from pathlib import Path
from typing import List

from pydantic import Field

from snippy_ng.stages import BaseStage, BaseOutput
from snippy_ng.dependencies import bcftools
from snippy_ng.envvars import EnvVar


class BcftoolsConsequencesCallerOutput(BaseOutput):
    annotated_vcf: Path = Field(..., description="VCF file with functional consequence annotations from bcftools csq")

class BcftoolsConsequencesCaller(BaseStage):
    """
    Call consequences using Bcftools csq.
    """
    reference: Path = Field(..., description="Reference file",)
    variants: Path = Field(..., description="Input VCF file",)
    features: Path = Field(..., description="Input features file")
    use_local_csq: bool = EnvVar("LOCAL_BCFTOOLS_CSQ", default=True, description="Whether to use bcftools csq's --local-csq mode for consequence annotation. Theres a bug in bcftools csq that means we need to use --local-csq for correct annotation of compound variants https://github.com/samtools/bcftools/issues/2543.")

    _dependencies = [
        bcftools
    ]

    @property
    def output(self) -> BcftoolsConsequencesCallerOutput:
        return BcftoolsConsequencesCallerOutput(
            annotated_vcf=Path(f"{self.prefix}.vcf")
        )

    def create_commands(self, ctx) -> List:
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
            cmd = self.shell_cmd([
                "cp", str(self.variants), str(self.output.annotated_vcf)
            ], description=f"Copy VCF file (no features found): {self.variants} -> {self.output.annotated_vcf}")
            return [cmd]
        
        
        bcf_csq_args = ["bcftools", "csq"]
        if self.use_local_csq:
            bcf_csq_args.append("--local-csq")

        bcf_csq_args.extend([
            "--threads", str(ctx.cpus),
            "--ncsq", "64",
            "-f", str(self.reference),
            "-g", str(self.features),
            "-o", str(self.output.annotated_vcf),
            "--phase", "m",
            str(self.variants),
        ])

        bcf_csq_cmd = self.shell_cmd(
            bcf_csq_args,
            description="Annotate variants with consequences using bcftools csq",
        )
        return [bcf_csq_cmd]