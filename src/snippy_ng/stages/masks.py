from pathlib import Path
from typing import List, Optional

from snippy_ng.stages.base import BaseStage, BaseOutput
from snippy_ng.dependencies import bedtools 

from pydantic import Field



class ComprehensiveMaskOutput(BaseOutput):
    """Output from the comprehensive masking stage"""
    masked_fasta: Path
    min_depth_bed: Optional[Path] = None
    zero_depth_bed: Path
    zero_depth_masked_fasta: Optional[Path] = None
    min_depth_masked_fasta: Optional[Path] = None


class ComprehensiveMask(BaseStage):
    """
    Comprehensive masking stage that generates depth-based masks and applies them sequentially.
    
    This stage:
    1. Generates min-depth mask BED file (if min_depth > 0)
    2. Generates zero-depth mask BED file
    3. Applies masks sequentially to the reference FASTA
    4. Optionally applies user-supplied mask
    """
    bam: Path = Field(..., description="Input BAM file")
    reference: Path = Field(..., description="Reference FASTA file")
    prefix: str = Field(..., description="Output file prefix")
    min_depth: int = Field(0, description="Minimum depth threshold (0 = skip min depth masking)")
    user_mask: Optional[Path] = Field(None, description="Optional user-supplied BED mask file")

    _dependencies = [
        bedtools
    ]

    @property
    def output(self) -> ComprehensiveMaskOutput:
        return ComprehensiveMaskOutput(
            masked_fasta=Path(f"{self.prefix}.masked.afa"),
            min_depth_bed=Path(f"{self.prefix}.mindepth.bed") if self.min_depth > 0 else None,
            zero_depth_bed=Path(f"{self.prefix}.zerodepth.bed"),
            min_depth_masked_fasta=Path(f"{self.prefix}.mindepth_masked.afa") if self.min_depth > 0 else None,
            zero_depth_masked_fasta=Path(f"{self.prefix}.zerodepth_masked.afa")
        )

    @property
    def commands(self) -> List:
        """Generate all masking commands in sequence"""
        commands = []
        
        # Generate min-depth mask if requested
        if self.min_depth > 0:
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
        
        # Apply min-depth mask first (with 'n')
        if self.min_depth > 0:
            commands.append(self._apply_mask_command(
                input_fasta=current_fasta,
                mask_bed=self.output.min_depth_bed,
                output_fasta=self.output.min_depth_masked_fasta,
                mask_char="n",
                description=f"Apply min-depth mask (< {self.min_depth})"
            ))
            current_fasta = self.output.min_depth_masked_fasta
        
        # Apply zero-depth mask (with '-')
        commands.append(self._apply_mask_command(
            input_fasta=current_fasta,
            mask_bed=self.output.zero_depth_bed,
            output_fasta=self.output.zero_depth_masked_fasta,
            mask_char="-",
            description="Apply zero-depth mask"
        ))
        current_fasta = self.output.zero_depth_masked_fasta
        
        # Apply user mask if provided (with 'N')
        if self.user_mask:
            commands.append(self._apply_mask_command(
                input_fasta=current_fasta,
                mask_bed=self.user_mask,
                output_fasta=self.output.masked_fasta,
                mask_char="N",
                description="Apply user-supplied mask"
            ))
        else:
            # Just copy the final result
            commands.append(self.shell_cmd([
                "cp", str(current_fasta), str(self.output.masked_fasta)
            ], description="Copy final masked FASTA"))
        
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