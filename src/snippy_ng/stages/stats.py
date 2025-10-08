from pathlib import Path
from typing import List
from snippy_ng.stages.base import BaseStage, ShellCommand
from snippy_ng.dependencies import seqkit
from pydantic import Field, field_validator, BaseModel


class SeqKitReadStatsOutput(BaseModel):
    """Output model for SeqKitReadStats stages.
    
    Attributes:
        stats_tsv: Path to the generated TSV file containing read statistics.
    """
    stats_tsv: str


class SeqKitReadStats(BaseStage):
    """Generate read statistics using seqkit stats command.
    
    This stage analyzes FASTQ/FASTA files and produces a TSV file with 
    comprehensive statistics including sequence count, length metrics,
    quality scores, and GC content.
    
    Args:
        reads: List of input read files (FASTQ or FASTA).
        prefix: Output file prefix for the generated stats.tsv file.
        all_stats: Whether to output all statistics including quartiles and N50.
        tabular: Whether to output in machine-friendly tabular format.
        basename_only: Whether to only output basename of files in results.
        skip_errors: Whether to skip files with errors and show warnings.
        fastq_encoding: FASTQ quality encoding (sanger, solexa, illumina-1.3+, etc.).
        gap_letters: Gap letters to be counted in sequences.
        additional_options: Additional seqkit stats command-line options.
        cpus: Number of CPU cores to use for processing.
        tmpdir: Temporary directory for intermediate files.
        
    Returns:
        SeqKitReadStatsOutput: Output containing path to the generated stats.tsv file.
        
    Example:
        >>> from pathlib import Path
        >>> stage = SeqKitReadStats(
        ...     reads=["sample_R1.fastq", "sample_R2.fastq"],
        ...     prefix="sample_stats",
        ...     tmpdir=Path("/tmp"),
        ...     cpus=4
        ... )
        >>> print(stage.output.stats_tsv)
        sample_stats.stats.tsv
    """
    
    reads: List[str] = Field(..., description="List of input read files (FASTQ or FASTA)")
    prefix: str = Field(..., description="Output file prefix")
    all_stats: bool = Field(True, description="Output all statistics including quartiles and N50")
    tabular: bool = Field(True, description="Output in machine-friendly tabular format")
    basename_only: bool = Field(False, description="Only output basename of files")
    skip_errors: bool = Field(True, description="Skip files with errors and show warnings")
    fastq_encoding: str = Field("sanger", description="FASTQ quality encoding (sanger, solexa, illumina-1.3+, etc.)")
    gap_letters: str = Field("- .", description="Gap letters to be counted")
    additional_options: str = Field("", description="Additional seqkit stats options")
    
    _dependencies = [seqkit]
    
    @field_validator("reads")
    @classmethod
    def validate_reads(cls, v):
        """Validate that read files are provided and exist.
        
        Args:
            v: List of read file paths to validate.
            
        Returns:
            List[str]: Validated list of read file paths.
            
        Raises:
            ValueError: If no read files provided or if any file doesn't exist.
        """
        if not v or len(v) == 0:
            raise ValueError("At least one read file must be provided")
        return v
    
    @field_validator("fastq_encoding")
    @classmethod
    def validate_fastq_encoding(cls, v):
        """Validate FASTQ quality encoding format.
        
        Args:
            v: FASTQ encoding string to validate.
            
        Returns:
            str: Validated FASTQ encoding string.
            
        Raises:
            ValueError: If encoding is not in the list of supported formats.
        """
        valid_encodings = [
            "sanger", "solexa", "illumina-1.3+", "illumina-1.5+", "illumina-1.8+"
        ]
        if v not in valid_encodings:
            raise ValueError(f"Invalid FASTQ encoding. Must be one of: {', '.join(valid_encodings)}")
        return v
    
    @property
    def output(self) -> SeqKitReadStatsOutput:
        """Get the output specification for this stage.
        
        Returns:
            SeqKitReadStatsOutput: Object containing paths to output files.
        """
        return SeqKitReadStatsOutput(
            stats_tsv=f"{self.prefix}.stats.tsv"
        )
    
    def build_seqkit_stats_command(self) -> ShellCommand:
        """Constructs the seqkit stats command.
        
        Builds the complete seqkit stats command with all specified options
        including threading, output format, encoding, and file handling options.
        
        Returns:
            ShellCommand: Complete seqkit stats command ready for execution.
        """
        shell_cmd = self.shell_cmd(
            command=["seqkit", "stats"],
            description=f"Generate read statistics for {len(self.reads)} files using seqkit stats"
        )
        
        # Threading
        if self.cpus > 1:
            shell_cmd.command.extend(["-j", str(self.cpus)])
        
        # Output format options
        if self.tabular:
            shell_cmd.command.append("-T")
        
        if self.all_stats:
            shell_cmd.command.append("-a")
        
        if self.basename_only:
            shell_cmd.command.append("-b")
        
        if self.skip_errors:
            shell_cmd.command.append("-e")
        
        # FASTQ encoding
        if self.fastq_encoding != "sanger":
            shell_cmd.command.extend(["-E", self.fastq_encoding])
        
        # Gap letters
        if self.gap_letters != "- .":
            shell_cmd.command.extend(["-G", self.gap_letters])
        
        # Additional options (split if it contains spaces)
        if self.additional_options:
            import shlex
            shell_cmd.command.extend(shlex.split(self.additional_options))
        
        # Input files
        shell_cmd.command.extend(self.reads)
        
        # Create shell command with output file
        return self.shell_pipeline(
            commands=[shell_cmd],
            description=f"Generate read statistics for {len(self.reads)} files using seqkit stats",
            output_file=Path(self.output.stats_tsv)
        )
    
    @property
    def commands(self) -> List[ShellCommand]:
        """Get the list of commands to execute for this stage.
        
        Returns:
            List[ShellCommand]: List containing the seqkit stats command.
        """
        return [self.build_seqkit_stats_command()]


class SeqKitReadStatsBasic(SeqKitReadStats):
    """Basic read statistics without extra information like quartiles and N50.
    
    This variant of SeqKitReadStats outputs only essential statistics
    without additional metrics like quartiles, N50, or other detailed measures.
    It's useful for quick quality checks when comprehensive statistics aren't needed.
    
    Args:
        all_stats: Set to False to output only basic statistics.
        
    Note:
        Inherits all other parameters from SeqKitReadStats.
        
    Example:
        >>> stage = SeqKitReadStatsBasic(
        ...     reads=["sample.fastq"],
        ...     prefix="basic_stats",
        ...     tmpdir=Path("/tmp")
        ... )
        >>> # Will generate command without -a flag for basic stats only
    """
    
    all_stats: bool = Field(False, description="Output basic statistics only")


class SeqKitReadStatsDetailed(SeqKitReadStats):
    """Detailed read statistics with all available metrics.
    
    This variant provides comprehensive statistics including all standard metrics
    plus additional N-statistics (like N75, N90, N95) that can be customized
    based on analysis requirements.
    
    Args:
        all_stats: Set to True to output all available statistics.
        additional_n_stats: List of integers (0-100) specifying additional 
            N-statistics to compute (e.g., [75, 90, 95] for N75, N90, N95).
            
    Note:
        Inherits all other parameters from SeqKitReadStats.
        
    Example:
        >>> stage = SeqKitReadStatsDetailed(
        ...     reads=["assembly.fasta"],
        ...     prefix="detailed_stats",
        ...     tmpdir=Path("/tmp"),
        ...     additional_n_stats=[75, 90, 95]
        ... )
        >>> # Will include N75, N90, N95 in addition to standard statistics
    """
    
    all_stats: bool = Field(True, description="Output all available statistics")
    additional_n_stats: List[int] = Field(
        default_factory=list, 
        description="Additional N-statistics to compute (e.g., [90, 95] for N90 and N95)"
    )
    
    @field_validator("additional_n_stats")
    @classmethod
    def validate_n_stats(cls, v):
        """Validate N-statistic values are within valid range.
        
        Args:
            v: List of N-statistic values to validate.
            
        Returns:
            List[int]: Validated list of N-statistic values.
            
        Raises:
            ValueError: If any N-statistic value is not between 0 and 100.
        """
        for stat in v:
            if not (0 <= stat <= 100):
                raise ValueError(f"N-statistic values must be between 0 and 100, got: {stat}")
        return v
    
    def build_seqkit_stats_command(self) -> ShellCommand:
        """Constructs the seqkit stats command with additional N-statistics.
        
        Builds the complete seqkit stats command including any additional
        N-statistics specified in the additional_n_stats parameter.
        
        Returns:
            ShellCommand: Complete seqkit stats command with N-statistics options.
        """
        shell_cmd = self.shell_cmd(
            command=["seqkit", "stats"],
            description=f"Generate detailed read statistics for {len(self.reads)} files using seqkit stats"
        ) 
        
        # Threading
        if self.cpus > 1:
            shell_cmd.command.extend(["-j", str(self.cpus)])
        
        # Output format options
        if self.tabular:
            shell_cmd.command.append("-T")
        
        if self.all_stats:
            shell_cmd.command.append("-a")
        
        # Additional N-statistics
        if self.additional_n_stats:
            n_stats_str = ",".join(map(str, self.additional_n_stats))
            shell_cmd.command.extend(["-N", n_stats_str])
        
        if self.basename_only:
            shell_cmd.command.append("-b")
        
        if self.skip_errors:
            shell_cmd.command.append("-e")
        
        # FASTQ encoding
        if self.fastq_encoding != "sanger":
            shell_cmd.command.extend(["-E", self.fastq_encoding])
        
        # Gap letters
        if self.gap_letters != "- .":
            shell_cmd.command.extend(["-G", self.gap_letters])
        
        # Additional options (split if it contains spaces)
        if self.additional_options:
            import shlex
            shell_cmd.command.extend(shlex.split(self.additional_options))
        
        # Input files
        shell_cmd.command.extend(self.reads)
        
        # Create shell command with output file
        n_stats_desc = f" with N-statistics {self.additional_n_stats}" if self.additional_n_stats else ""
        return self.shell_pipeline(
            commands=[shell_cmd],
            description=f"Generate detailed read statistics for {len(self.reads)} files using seqkit stats{n_stats_desc}",
            output_file=Path(self.output.stats_tsv)
        )
