from pathlib import Path
from typing import List

from snippy_ng.stages import BaseStage, BaseOutput, TempPath
from snippy_ng.dependencies import bedtools, bcftools 

from pydantic import Field



class DepthMaskOutput(BaseOutput):
    """Output from the minimum-depth masking stage."""
    masked_fasta: Path = Field(..., description="FASTA with low-depth regions (depth < min_depth) masked")
    min_depth_bed: TempPath = Field(..., description="BED file with regions to mask due to low depth")


class DepthMask(BaseStage):
    """
    Minimum-depth masking stage.

    This stage masks regions with depth strictly below `min_depth` using `N`.
    """
    cram: Path = Field(..., description="Input CRAM file")
    fasta: Path = Field(..., description="Input FASTA file to be masked")
    min_depth: int = Field(..., description="Minimum depth threshold")
    mask_char: str = Field("N", description="Character to use for masking")

    _dependencies = [
        bedtools
    ]

    @property
    def output(self) -> DepthMaskOutput:
        return DepthMaskOutput(
            masked_fasta=Path(f"{self.prefix}.mindepth_masked.fasta"),
            min_depth_bed=Path(f"{self.prefix}.mindepth.bed")
        )

    def create_commands(self, ctx) -> List:
        """Generate commands to create and apply a minimum-depth mask."""
        return [
            *self._generate_depth_mask_commands(
                filter_condition=f"<{self.min_depth}",
                output_bed=self.output.min_depth_bed,
                description=f"Generate min-depth mask (depth < {self.min_depth})",
            ),
            self._apply_mask_command(
                input_fasta=self.fasta,
                mask_bed=self.output.min_depth_bed,
                output_fasta=self.output.masked_fasta,
                mask_char=self.mask_char,
                description=f"Apply min-depth mask (< {self.min_depth})",
            ),
        ]
    
    def _generate_depth_mask_commands(self, filter_condition: str, output_bed: Path, description: str) -> List:
        """Generate commands to create a depth-based mask BED file using bedtools genomecov"""
        view_cram = self.shell_cmd([
            "samtools", "view", "-b", str(self.cram)
        ], description="Convert CRAM to BAM for depth calculation")
        genomecov_cmd = self.shell_cmd(
            ["bedtools", "genomecov", "-ibam", '-', "-bga"],
            description="Generate genome coverage in BED format"
        )
        awk_cmd = self.shell_cmd(
            ["awk", f'$4{filter_condition} {{print $1"\\t"$2"\\t"$3}}'],
            description=f"Filter for regions with depth {filter_condition}"
        )
        
        return [self.shell_pipeline(
            [view_cram, genomecov_cmd, awk_cmd], 
            output_file=output_bed, 
            description=description
        )]
    
    def _apply_mask_command(self, input_fasta: Path, mask_bed: Path, output_fasta: Path, mask_char: str, description: str):
        """Generate command to apply a mask to a FASTA file"""
        return self.shell_cmd([
            "bedtools", "maskfasta",
            "-fi", str(input_fasta),
            "-bed", str(mask_bed),
            "-fo", str(output_fasta),
            "-fullHeader",
            "-mc", mask_char
        ], description=description)


class DelMaskOutput(BaseOutput):
    """Output from the deletion (zero-depth) masking stage."""
    masked_fasta: Path = Field(..., description="FASTA with zero-depth regions masked using the deletion character")
    zero_depth_bed: TempPath = Field(..., description="BED file with zero-depth regions")


class DelMask(BaseStage):
    """
    Zero-depth masking stage.

    This stage masks regions with depth equal to zero using `-`.
    """
    bam: Path = Field(..., description="Input BAM file")
    fasta: Path = Field(..., description="Input FASTA file to be masked")
    mask_char: str = Field("-", description="Character to use for masking")

    _dependencies = [
        bedtools
    ]

    @property
    def output(self) -> DelMaskOutput:
        return DelMaskOutput(
            masked_fasta=Path(f"{self.prefix}.del_masked.fasta"),
            zero_depth_bed=Path(f"{self.prefix}.zerodepth.bed")
        )

    def create_commands(self, ctx) -> List:
        """Generate commands to create and apply a zero-depth mask."""
        return [
            *self._generate_zero_depth_mask_commands(),
            self._apply_zero_depth_mask(),
        ]

    def _generate_zero_depth_mask_commands(self) -> List:
        """Generate commands to create a zero-depth BED mask file."""
        genomecov_cmd = self.shell_cmd(
            ["bedtools", "genomecov", "-ibam", str(self.bam), "-bga"],
            description="Generate genome coverage in BED format"
        )
        awk_cmd = self.shell_cmd(
            ["awk", '$4==0 {print $1"\\t"$2"\\t"$3}'],
            description="Filter for regions with depth == 0"
        )

        return [self.shell_pipe(
            [genomecov_cmd, awk_cmd],
            output_file=self.output.zero_depth_bed,
            description="Generate zero-depth mask"
        )]

    def _apply_zero_depth_mask(self):
        """Apply zero-depth mask to FASTA file."""
        return self.shell_cmd([
            "bedtools", "maskfasta",
            "-fi", str(self.fasta),
            "-bed", str(self.output.zero_depth_bed),
            "-fo", str(self.output.masked_fasta),
            "-fullHeader",
            "-mc", self.mask_char
        ], description="Apply zero-depth mask")


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


class HetMask(BaseStage):
    """
    Heterozygous and low quality sites masking stage.
    
    This stage masks heterozygous sites and low QUAL sites with 'n' characters.
    It identifies sites where:
    - GT="het" (any heterozygous genotype like 0/1, 1/2, etc.)
    - QUAL < min_qual threshold
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
                "-i", f'GT="het" || QUAL<{self.min_qual}',
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

