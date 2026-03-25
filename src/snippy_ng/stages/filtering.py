from pathlib import Path
from typing import List, Optional
from snippy_ng.stages import BaseStage, BaseOutput, ShellCommand
from snippy_ng.dependencies import samtools
from pydantic import Field, field_validator


class SamtoolsFilterOutput(BaseOutput):
    bam: Path = Field(..., description="Filtered alignment file in BAM format")
    stats: Path = Field(..., description="Flagstat output file for the BAM file")


class SamtoolsFilter(BaseStage):
    """
    Filter BAM files using Samtools to remove unwanted alignments.
    """
    
    bam: Path = Field(..., description="Input BAM file to filter")
    reference: Path = Field(..., description="Reference FASTA file for BAM output")
    min_mapq: int = Field(20, description="Minimum mapping quality")
    exclude_flags: int = Field(1796, description="SAM flags to exclude (default: unmapped, secondary, qcfail, duplicate)")
    include_flags: Optional[int] = Field(None, description="SAM flags to include")
    regions: Optional[str] = Field(None, description="Regions to include (BED file or region string)")
    additional_filters: str = Field("", description="Additional samtools view options")
    
    _dependencies = [samtools]
    
    @property
    def output(self) -> SamtoolsFilterOutput:
        filtered_bam = f"{self.prefix}.filtered.bam"
        return SamtoolsFilterOutput(
            bam=filtered_bam,
            stats=f"{filtered_bam}.flagstat.txt"
        )
    
    def build_filter_command(self, ctx) -> ShellCommand:
        """Constructs the samtools view command for filtering."""
        cmd_parts = [
            "samtools",
            "view",
            "-O",
            "bam",
            "--reference", str(self.reference),
            "-o",
            str(self.output.bam),
        ]
        
        # Add threading
        if ctx.cpus > 1:
            cmd_parts.extend(["--threads", str(ctx.cpus - 1)])
        
        # Add mapping quality filter
        if self.min_mapq > 0:
            cmd_parts.extend(["-q", str(self.min_mapq)])
        
        # Add flag filters
        if self.exclude_flags:
            cmd_parts.extend(["-F", str(self.exclude_flags)])
        
        if self.include_flags is not None:
            cmd_parts.extend(["-f", str(self.include_flags)])
        
        # Add regions if specified as BED file
        if self.regions and Path(self.regions).exists():
            cmd_parts.extend(["-L", str(self.regions)])
        
        # Add additional filters (split if it contains spaces)
        if self.additional_filters:
            import shlex
            cmd_parts.extend(shlex.split(self.additional_filters))
        
        # Add input file
        cmd_parts.append(str(self.bam))
        
        # Add region string if not a file
        if self.regions and not Path(self.regions).exists():
            cmd_parts.append(str(self.regions))
        
        filter_cmd = self.shell_cmd(
            command=cmd_parts,
            description=f"Filter BAM file with MAPQ>={self.min_mapq}, flags={self.exclude_flags}",
        )
        
        return filter_cmd
    
    def build_index_command(self):
        """Returns the samtools index command."""
        return self.shell_cmd([
            "samtools", "index", str(self.output.bam)
        ], description=f"Index filtered BAM file: {self.output.bam}")
    
    def build_flagstat_command(self):
        """Returns the samtools flagstat command."""
        return self.shell_cmd(
            ["samtools", "flagstat", str(self.output.bam)],
            description=f"Generate alignment statistics for {self.output.bam}",
            output_file=Path(self.output.stats),
        )

    def test_mapped_reads_not_zero(self):
        """Test that the BAM file contains mapped reads."""
        flagstat_path = self.output.stats
        if not flagstat_path.exists():
            raise FileNotFoundError(
                f"Expected flagstat output file {flagstat_path} was not created"
            )

        with open(flagstat_path) as f:
            for line in f:
                if "mapped (" in line:
                    mapped_reads = int(line.split()[0])
                    if mapped_reads == 0:
                        raise ValueError(
                            "No reads were mapped in BAM file after filtering. Did you use the correct reference?"
                        )
                    break

    def create_commands(self, ctx) -> List:
        """Constructs the filtering commands."""
        filter_cmd = self.build_filter_command(ctx)
        index_cmd = self.build_index_command()
        stats_cmd = self.build_flagstat_command()
        return [filter_cmd, index_cmd, stats_cmd]


class SamtoolsFilterByRegion(SamtoolsFilter):
    """
    Filter BAM file to include only alignments in specified regions.
    """
    
    regions: str = Field(..., description="Regions to include (BED file or region string)")
    
    @field_validator("regions")
    @classmethod
    def validate_regions(cls, v):
        if not v:
            raise ValueError("Regions must be specified")
        return v


class SamtoolsFilterByQuality(SamtoolsFilter):
    """
    Filter BAM file based on mapping quality and alignment flags.
    """
    
    min_mapq: int = Field(30, description="Minimum mapping quality (higher than default)")
    exclude_flags: int = Field(3844, description="SAM flags to exclude (default + supplementary)")


class SamtoolsFilterProperPairs(SamtoolsFilter):
    """
    Filter BAM file to include only properly paired reads.
    """
    
    include_flags: int = Field(2, description="Include only properly paired reads")
    exclude_flags: int = Field(1796, description="Exclude unmapped, secondary, qcfail, duplicate")
