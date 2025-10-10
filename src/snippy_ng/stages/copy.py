from pathlib import Path
from pydantic import Field

from snippy_ng.stages.base import BaseStage, BaseOutput


class CopyFileOutput(BaseOutput):
    copied_file: Path


class CopyFile(BaseStage):
    """
    Copy a file from input location to output location.
    """
    input: Path = Field(..., description="Input file to copy")
    output_path: Path = Field(..., description="Output file path")

    @property
    def output(self) -> CopyFileOutput:
        return CopyFileOutput(
            copied_file=self.output_path
        )

    @property
    def commands(self):
        return [
            self.shell_cmd([
                "cp", str(self.input), str(self.output_path)
            ], description=f"Copy {self.input} to {self.output_path}")
        ]
