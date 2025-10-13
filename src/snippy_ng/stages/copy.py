from pathlib import Path
from typing import Optional
from pydantic import Field

from snippy_ng.stages.base import BaseStage, BaseOutput
from snippy_ng.dependencies import seqkit


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
        return CopyFileOutput(copied_file=self.output_path)

    @property
    def commands(self):
        return [
            self.shell_cmd(
                ["cp", str(self.input), str(self.output_path)],
                description=f"Copy {self.input} to {self.output_path}",
            )
        ]


class CopyFasta(CopyFile):
    """
    Copy a FASTA file from input location to output location.
    """

    header: Optional[str] = Field(None, description="Header for the FASTA file")

    _dependencies = [seqkit]

    @property
    def commands(self):
        if self.header:
            cmd = self.shell_cmd(
                ["seqkit", "replace", "-p", ".*", "-r", self.header, str(self.input)],
                description=f"Copy and rename FASTA {self.input} to {self.output_path}",
            )
        else:
            cmd = self.shell_cmd(
                ["cp", str(self.input), str(self.output_path)],
                description=f"Copy {self.input} to {self.output_path}",
            )
        return [
            self.shell_pipeline(
                [cmd],
                output_file=self.output.copied_file,
                description="Copy FASTA file to output location",
            )
        ]
