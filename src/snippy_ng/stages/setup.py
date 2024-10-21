from snippy_ng.stages.base import BaseStage
from pydantic import Field, BaseModel
from pathlib import Path


class PrepareReferenceOutput(BaseModel):
    reference: Path
    gff: Path


class PrepareReference(BaseStage):
    input: Path = Field(..., description="Reference file")
    reference_prefix: str = Field("ref", description="Output reference name")
    reference_dir: Path = Field(Path("reference"), description="Reference directory")

    @property
    def output(self) -> PrepareReferenceOutput:
        return PrepareReferenceOutput(
            reference=self.reference_dir / f"{self.reference_prefix}.fa",
            gff=self.reference_dir / f"{self.reference_prefix}.gff",
        )

    @property
    def commands(self):
        
        return [
            f"rm -f {self.output.reference}",
            f"mkdir -p {self.reference_dir}",
            f"ln -s {self.input} {self.output.reference}",
        ]