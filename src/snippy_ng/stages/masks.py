from pathlib import Path
from typing import List, Optional

from snippy_ng.stages.base import BaseStage, BaseOutput
from snippy_ng.dependencies import bedtools, bcftools 

from pydantic import Field



class DepthMaskOutput(BaseOutput):
    """Output from the depth masking stage"""
    masked_fasta: Path
    min_depth_bed: Optional[Path] = None
    zero_depth_bed: Path


class DepthMask(BaseStage):
    """
    Depth masking stage that generates depth-based masks and applies them sequentially.
    
    This stage:
    1. Generates min-depth mask BED file (if min_depth > 0)
    2. Generates zero-depth mask BED file
    3. Applies masks sequentially to the reference FASTA
    """
    bam: Path = Field(..., description="Input BAM file")
    reference: Path = Field(..., description="Reference FASTA file")
    prefix: str = Field(..., description="Output file prefix")
    min_depth: int = Field(0, description="Minimum depth threshold (0 = skip min depth masking)")

    _dependencies = [
        bedtools
    ]

    @property
    def output(self) -> DepthMaskOutput:
        return DepthMaskOutput(
            masked_fasta=Path(f"{self.prefix}.depth_masked.fasta"),
            min_depth_bed=Path(f"{self.prefix}.mindepth.bed") if self.min_depth > 0 else None,
            zero_depth_bed=Path(f"{self.prefix}.zerodepth.bed")
        )

    @property
    def commands(self) -> List:
        """Generate all depth masking commands in sequence"""
        commands = []
        
        # Generate min-depth mask if requested
        min_depth_cmd = self._generate_depth_mask_commands(
            filter_condition=f"<{self.min_depth}",
            output_bed=self.output.min_depth_bed,
            description=f"Generate min-depth mask (depth < {self.min_depth})"
        )
        commands.extend(min_depth_cmd)
        
        # Generate zero-depth mask (always)
        zero_depth_cmd = self._generate_depth_mask_commands(
            filter_condition="==0",
            output_bed=self.output.zero_depth_bed,
            description="Generate zero-depth mask"
        )
        commands.extend(zero_depth_cmd)
        
        # Apply masks sequentially
        current_fasta = self.reference
        
        # Apply min-depth mask first (with 'N')
        min_depth_masked = Path(f"{self.prefix}.mindepth_masked.fasta")
        commands.append(self._apply_mask_command(
            input_fasta=current_fasta,
            mask_bed=self.output.min_depth_bed,
            output_fasta=min_depth_masked,
            mask_char="N",
            description=f"Apply min-depth mask (< {self.min_depth})"
        ))
        current_fasta = min_depth_masked
        
        # Apply zero-depth mask (with '-')
        commands.append(self._apply_mask_command(
            input_fasta=current_fasta,
            mask_bed=self.output.zero_depth_bed,
            output_fasta=self.output.masked_fasta,
            mask_char="-",
            description="Apply zero-depth mask"
        ))
        
        return commands
    
    def _generate_depth_mask_commands(self, filter_condition: str, output_bed: Path, description: str) -> List:
        """Generate commands to create a depth-based mask BED file using bedtools genomecov"""
        genomecov_cmd = self.shell_cmd(
            ["bedtools", "genomecov", "-ibam", str(self.bam), "-bga"],
            description="Generate genome coverage in BED format"
        )
        awk_cmd = self.shell_cmd(
            ["awk", f'$4{filter_condition} {{print $1"\\t"$2"\\t"$3}}'],
            description=f"Filter for regions with depth {filter_condition}"
        )
        
        return [self.shell_pipeline(
            [genomecov_cmd, awk_cmd], 
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


class ApplyMaskOutput(BaseOutput):
    masked_fasta: Path


class ApplyMask(BaseStage):
    """
    Masking stage that applies a supplied BED mask to a FASTA file.
    """
    fasta: Path = Field(..., description="Input FASTA file to be masked")
    mask_bed: Path = Field(..., description="BED file with regions to mask")
    prefix: str = Field(..., description="Output file prefix")
    mask_char: str = Field("X", description="Character to use for masking")

    _dependencies = [
        bedtools
    ]

    @property
    def output(self) -> ApplyMaskOutput:
        return ApplyMaskOutput(
            masked_fasta=Path(f"{self.prefix}.masked.fasta")
        )

    @property
    def commands(self) -> List:
        """Apply mask to FASTA file using temporary copy"""
        temp_fasta = self.tmpdir / f"temp_{self.fasta.name}"
        
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
    masked_fasta: Path
    het_sites_bed: Path


class HetMask(BaseStage):
    """
    Heterozygous and low quality sites masking stage.
    
    This stage masks heterozygous sites and low QUAL sites with 'n' characters.
    It identifies sites where:
    - GT="het" (any heterozygous genotype like 0/1, 1/2, etc.)
    - QUAL < min_qual threshold
    """
    reference: Path = Field(..., description="Reference FASTA file")
    vcf: Path = Field(..., description="Input VCF file (raw VCF recommended)")
    prefix: str = Field(..., description="Output file prefix")
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

    @property
    def commands(self) -> List:
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
        ], description=f"Extract heterozygous sites and sites with QUAL < {self.min_qual}")
        
        return [self.shell_pipeline(
            [bcftools_query],
            output_file=self.output.het_sites_bed,
            description="Create BED file of heterozygous and low quality sites"
        )]
    
    def _apply_het_mask(self):
        """Apply het sites mask to FASTA file"""
        return self.shell_cmd([
            "bedtools", "maskfasta",
            "-fi", str(self.reference),
            "-bed", str(self.output.het_sites_bed),
            "-fo", str(self.output.masked_fasta),
            "-fullHeader",
            "-mc", self.mask_char
        ], description=f"Apply heterozygous sites mask with '{self.mask_char}' character")

