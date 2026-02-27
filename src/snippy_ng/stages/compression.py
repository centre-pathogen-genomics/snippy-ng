# Concrete Alignment Strategies
from pathlib import Path

from snippy_ng.stages import BaseStage, BaseOutput

from pydantic import Field


class Compressor(BaseStage):
    input: Path = Field(..., description="File to compress")
    suffix: str = Field(..., description="Compression suffix")


class BgzipCompressorOutput(BaseOutput):
    compressed: Path

class BgzipCompressor(Compressor):
    """
    Compress a file using bgzip.
    """

    @property
    def output(self) -> BgzipCompressorOutput:
        return BgzipCompressorOutput(
            compressed=self.input.with_name(f"{self.input.name}.{self.suffix}")
        )

    def create_commands(self, ctx):
        """Constructs the gzip compression command."""
        bgzip_cmd = self.shell_cmd([
                "bgzip", "-c", str(self.input)
            ], 
            output_file=self.output.compressed,
            description="Compressing file with bgzip"
        )
        return [bgzip_cmd]

class VcfCompressorOutput(BaseOutput):
    compressed: Path = Field(..., description="Compressed VCF file")

class VcfCompressor(BgzipCompressor):
    """
    Compress a VCF file using bgzip.
    """
    suffix: str = Field("vcf.gz", description="Compression suffix for VCF files")

    @property
    def output(self) -> VcfCompressorOutput:
        return VcfCompressorOutput(
            compressed=self.input.with_name(f"{self.input.stem}.{self.suffix}")
        )