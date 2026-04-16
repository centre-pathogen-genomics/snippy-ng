from pathlib import Path
from typing import Iterable, List, TextIO

from snippy_ng.stages import BaseStage, BaseOutput, TempPath
from snippy_ng.dependencies import bedtools, bcftools

from pydantic import Field


class DepthBedsFromBamOutput(BaseOutput):
    """Output from combined depth BED generation."""
    zero_depth_bed: Path = Field(..., description="BED file with zero-depth regions")
    min_depth_bed: Path = Field(..., description="BED file with regions to mask due to low depth")


class DepthBedsFromBam(BaseStage):
    """Generate zero-depth and minimum-depth BED files from one coverage scan."""

    bam: Path = Field(..., description="Input BAM file")
    min_depth: int = Field(..., description="Minimum depth threshold")

    _dependencies = [
        bedtools
    ]

    @property
    def output(self) -> DepthBedsFromBamOutput:
        return DepthBedsFromBamOutput(
            zero_depth_bed=Path(f"{self.prefix}.zerodepth.bed"),
            min_depth_bed=Path(f"{self.prefix}.mindepth.bed"),
        )

    def create_commands(self, ctx) -> List:
        genomecov_cmd = self.shell_cmd(
            ["bedtools", "genomecov", "-ibam", str(self.bam), "-bga"],
            description="Generate genome coverage in BED format"
        )
        split_cmd = self.shell_cmd(
            [
                "awk",
                "-v", f"zero_depth_bed={self.output.zero_depth_bed}",
                "-v", f"min_depth_bed={self.output.min_depth_bed}",
                "-v", f"min_depth={self.min_depth}",
                'BEGIN { OFS="\\t"; printf "" > zero_depth_bed; printf "" > min_depth_bed } '
                '$4==0 { print $1, $2, $3 > zero_depth_bed; next } '
                '$4>0 && $4<min_depth { print $1, $2, $3 > min_depth_bed } '
                'END { close(zero_depth_bed); close(min_depth_bed) }',
            ],
            description="Split coverage into zero-depth and min-depth BED files"
        )

        return [
            self.shell_pipe(
                [genomecov_cmd, split_cmd],
                description="Generate zero-depth and min-depth BED blocks from one coverage scan",
            )
        ]

    @staticmethod
    def write_depth_beds(
        genomecov_lines: Iterable[str],
        min_depth: int,
        zero_depth_handle: TextIO,
        min_depth_handle: TextIO,
    ) -> None:
        for raw_line in genomecov_lines:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            try:
                depth = float(fields[3])
            except ValueError:
                continue

            bed_line = f"{fields[0]}\t{fields[1]}\t{fields[2]}\n"
            if depth == 0:
                zero_depth_handle.write(bed_line)
            elif 0 < depth < min_depth:
                min_depth_handle.write(bed_line)


class DepthMaskFromBedOutput(BaseOutput):
    """Output from applying a precomputed minimum-depth mask."""
    masked_fasta: Path = Field(..., description="FASTA with low-depth regions (depth < min_depth) masked")


class DepthMaskFromBed(BaseStage):
    """Apply a precomputed minimum-depth BED mask to a FASTA file."""

    fasta: Path = Field(..., description="Input FASTA file to be masked")
    mask_bed: Path = Field(..., description="BED file with regions to mask due to low depth")
    min_depth: int = Field(..., description="Minimum depth threshold")
    mask_char: str = Field("N", description="Character to use for masking")

    _dependencies = [
        bedtools
    ]

    @property
    def output(self) -> DepthMaskFromBedOutput:
        return DepthMaskFromBedOutput(
            masked_fasta=Path(f"{self.prefix}.mindepth_masked.fasta")
        )

    def create_commands(self, ctx) -> List:
        return [
            self.shell_cmd([
                "bedtools", "maskfasta",
                "-fi", str(self.fasta),
                "-bed", str(self.mask_bed),
                "-fo", str(self.output.masked_fasta),
                "-fullHeader",
                "-mc", self.mask_char
            ], description=f"Apply min-depth mask (0 < depth < {self.min_depth})")
        ]


class ApplyMaskOutput(BaseOutput):
    masked_fasta: Path = Field(..., description="FASTA with user-supplied BED mask applied")


class ApplyMask(BaseStage):
    """
    Masking stage that applies a supplied BED mask to a FASTA file.
    """
    fasta: Path = Field(..., description="Input FASTA file to be masked")
    mask_bed: Path = Field(..., description="BED file with regions to mask")
    mask_char: str = Field("X", description="Character to use for masking")

    _dependencies = [
        bedtools
    ]

    @property
    def output(self) -> ApplyMaskOutput:
        return ApplyMaskOutput(
            masked_fasta=Path(f"{self.prefix}.masked.fasta")
        )

    def create_commands(self, ctx) -> List:
        """Apply mask to FASTA file using temporary copy"""
        temp_fasta = self.fasta.with_suffix(".tmp")
        
        return [
            # Copy input FASTA to temporary location
            self.shell_cmd([
                "cp", str(self.fasta), str(temp_fasta)
            ], description=f"Copy input FASTA to temporary location: {temp_fasta}"),
            
            # Apply mask to temporary FASTA
            self.shell_cmd([
                "bedtools", "maskfasta",
                "-fi", str(temp_fasta),
                "-bed", str(self.mask_bed),
                "-fo", str(self.output.masked_fasta),
                "-fullHeader",
                "-mc", self.mask_char
            ], description="Masking FASTA with provided BED file"),
            
            # Clean up temporary file
            self.shell_cmd([
                "rm", str(temp_fasta)
            ], description="Remove temporary FASTA file")
        ]


class HetMaskOutput(BaseOutput):
    """Output from the heterozygous/low quality masking stage"""
    masked_fasta: Path = Field(..., description="FASTA with heterozygous and low-quality sites masked")
    het_sites_bed: TempPath = Field(..., description="BED file of heterozygous and low-quality variant sites")


class QualMask(BaseStage):
    """
    Low quality sites masking stage.
    
    This stage masks low QUAL sites with 'n' characters.
    """
    fasta: Path = Field(..., description="Input FASTA file to be masked")
    vcf: Path = Field(..., description="Input VCF file (raw VCF recommended)")
    min_qual: float = Field(20.0, description="Minimum QUAL threshold for sites (default: 20.0)")
    mask_char: str = Field("n", description="Character to use for masking (default: 'n')")

    _dependencies = [
        bcftools,
        bedtools
    ]

    @property
    def output(self) -> HetMaskOutput:
        return HetMaskOutput(
            masked_fasta=Path(f"{self.prefix}.het_masked.fasta"),
            het_sites_bed=Path(f"{self.prefix}.het_sites.bed")
        )

    def create_commands(self, ctx) -> List:
        """Generate het/low-qual masking commands"""
        commands = []
        
        # Generate BED file of heterozygous and low quality sites
        commands.extend(self._generate_het_sites_bed())
        
        # Apply mask to FASTA
        commands.append(self._apply_het_mask())
        
        return commands
    
    def _generate_het_sites_bed(self) -> List:
        """Generate BED file of heterozygous and low quality sites using bcftools"""
        # Query for het sites and low quality sites
        bcftools_query = self.shell_cmd([
                "bcftools", "query", 
                "-i", f'QUAL<{self.min_qual}',
                "-f", "%CHROM\\t%POS0\\t%POS\\n",
                str(self.vcf)
            ], 
            description=f"Extract heterozygous sites and sites with QUAL < {self.min_qual}",
            output_file=self.output.het_sites_bed,
        )
        
        return [bcftools_query]
    
    def _apply_het_mask(self):
        """Apply het sites mask to FASTA file"""
        return self.shell_cmd([
            "bedtools", "maskfasta",
            "-fi", str(self.fasta),
            "-bed", str(self.output.het_sites_bed),
            "-fo", str(self.output.masked_fasta),
            "-fullHeader",
            "-mc", self.mask_char
        ], description=f"Apply heterozygous sites mask with '{self.mask_char}' character")
