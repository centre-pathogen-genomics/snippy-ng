# Concrete Alignment Strategies
from pathlib import Path

from snippy_ng.stages import BaseStage, BaseOutput
from snippy_ng.dependencies import samtools

from pydantic import Field


class Compressor(BaseStage):
    input: Path = Field(..., description="File to compress")
    suffix: str = Field(..., description="Compression suffix")


class BgzipCompressorOutput(BaseOutput):
    gz: Path = Field(..., description="BGZF-compressed output file")

class BgzipCompressor(Compressor):
    """
    Compress a file using bgzip.
    """

    @property
    def output(self) -> BgzipCompressorOutput:
        return BgzipCompressorOutput(
            gz=self.input.with_name(f"{self.input.name}.{self.suffix}")
        )

    def create_commands(self, ctx):
        """Constructs the gzip compression command."""
        bgzip_cmd = self.shell_cmd([
                "bgzip", "-c", str(self.input)
            ], 
            output_file=self.output.gz,
            description="Compressing file with bgzip"
        )
        return [bgzip_cmd]


class VcfCompressorOutput(BaseOutput):
    gz: Path = Field(..., description="Compressed VCF file")

class VcfCompressor(BgzipCompressor):
    """
    Compress a VCF file using bgzip.
    """
    suffix: str = Field("vcf.gz", description="Compression suffix for VCF files")

    @property
    def output(self) -> VcfCompressorOutput:
        return VcfCompressorOutput(
            gz=self.input.with_name(f"{self.input.stem}.{self.suffix}")
        )


class CramCompressorOutput(BaseOutput):
    cram: Path = Field(..., description="Compressed CRAM file")


class CramCompressor(Compressor):
    """
    Convert a SAM/BAM file to CRAM with embedded reference.
    """

    reference: Path = Field(..., description="Reference FASTA file for CRAM encoding")
    suffix: str = Field("cram", description="Compression suffix for CRAM files")

    _dependencies = [samtools]

    @property
    def output(self) -> CramCompressorOutput:
        return CramCompressorOutput(
            cram=self.input.with_name(f"{self.input.stem}.{self.suffix}")
        )

    def create_commands(self, ctx):
        """Construct the SAM/BAM to CRAM conversion command."""
        cram_cmd = self.shell_cmd(
            [
                "samtools",
                "view",
                "--threads", str(ctx.cpus),
                "-C",
                "-T",
                str(self.reference),
                "--output-fmt-option",
                "embed_ref=1",
                "-o",
                str(self.output.cram),
                str(self.input),
            ],
            description="Convert SAM/BAM to CRAM with embedded reference",
        )
        return [cram_cmd]
