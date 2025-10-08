from pathlib import Path
from typing import List, Optional, Any
from snippy_ng.stages.base import BaseStage
from snippy_ng.dependencies import rasusa
from pydantic import Field, field_validator, model_validator, BaseModel


class RasusaDownsampleReadsOutput(BaseModel):
    """Output model for RasusaDownsampleReads stage.
    
    Attributes:
        downsampled_r1: Path to downsampled R1 read file.
        downsampled_r2: Optional path to downsampled R2 read file (for paired-end reads).
    """
    downsampled_r1: str
    downsampled_r2: Optional[str] = None


class RasusaDownsampleReads(BaseStage):
    """Downsample reads using rasusa based on genome coverage.
    
    This stage uses rasusa to randomly subsample sequencing reads to a specified
    coverage depth based on genome length. It supports both single-end and 
    paired-end reads, and can downsample to either a specific coverage depth 
    or a specific number of reads.
    
    Args:
        reads: List of input read files (FASTQ format, can be gzipped).
        prefix: Output file prefix for downsampled reads.
        genome_length: Genome length in base pairs (required for coverage-based downsampling).
        coverage: Target coverage depth for downsampling (mutually exclusive with num_reads).
        num_reads: Target number of reads for downsampling (mutually exclusive with coverage).
        seed: Random seed for reproducible downsampling.
        output_format: Output format for reads ('fastq' or 'fasta').
        compression_level: Compression level for gzipped output (1-9).
        additional_options: Additional rasusa command-line options.
        tmpdir: Temporary directory for intermediate files.
        cpus: Number of CPU cores to use.
        
    Returns:
        RasusaDownsampleReadsOutput: Object containing paths to downsampled read files.
        
    Example:
        >>> from pathlib import Path
        >>> stage = RasusaDownsampleReads(
        ...     reads=["sample_R1.fastq.gz", "sample_R2.fastq.gz"],
        ...     prefix="downsampled",
        ...     genome_length=197394,
        ...     coverage=50,
        ...     tmpdir=Path("/tmp")
        ... )
        >>> print(stage.output.downsampled_r1)
        'downsampled.downsampled.R1.fastq.gz'
        >>> print(stage.output.downsampled_r2)  
        'downsampled.downsampled.R2.fastq.gz'
    """
    
    reads: List[str] = Field(..., description="List of input read files (FASTQ format)")
    prefix: str = Field(..., description="Output file prefix for downsampled reads")
    genome_length: Optional[Any] = Field(None, description="Genome length in base pairs (required for coverage-based downsampling) - can be int or at_run_time() callable")
    coverage: Optional[float] = Field(None, description="Target coverage depth for downsampling")
    num_reads: Optional[int] = Field(None, description="Target number of reads for downsampling")
    seed: Optional[int] = Field(None, description="Random seed for reproducible downsampling")
    output_format: str = Field("fastq", description="Output format ('fastq' or 'fasta')")
    compression_level: int = Field(6, description="Compression level for gzipped output (1-9)")
    additional_options: str = Field("", description="Additional rasusa command-line options")
    
    _dependencies = [rasusa]
    
    @field_validator("reads")
    @classmethod
    def validate_reads(cls, v):
        """Validate that read files are provided.
        
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
    
    @model_validator(mode='after')
    def validate_coverage_or_reads(self):
        """Validate that either coverage or num_reads is specified, but not both.
        
        Returns:
            Self: The validated model instance.
            
        Raises:
            ValueError: If both or neither coverage and num_reads are specified.
        """
        if self.coverage is not None and self.num_reads is not None:
            raise ValueError("Cannot specify both coverage and num_reads - choose one")
        if self.coverage is None and self.num_reads is None:
            raise ValueError("Must specify either coverage or num_reads")
        
        # Also validate genome_length is provided when coverage is used
        # For AtRunTime objects, we can't check the value without triggering evaluation,
        # so we'll validate this later during command building
        if self.coverage is not None:
            # Check if we have a genome_length value or AtRunTime object
            if self.genome_length is None:
                raise ValueError("genome_length is required when using coverage-based downsampling")
            # If it's an AtRunTime object, skip validation here (will be validated at runtime)
        
        return self
    
    @field_validator("output_format")
    @classmethod
    def validate_output_format(cls, v):
        """Validate output format is supported.
        
        Args:
            v: Output format string.
            
        Returns:
            str: Validated output format.
            
        Raises:
            ValueError: If format is not supported.
        """
        valid_formats = ["fastq", "fasta"]
        if v not in valid_formats:
            raise ValueError(f"Invalid output format. Must be one of: {', '.join(valid_formats)}")
        return v
    
    @field_validator("compression_level")
    @classmethod
    def validate_compression_level(cls, v):
        """Validate compression level is within valid range.
        
        Args:
            v: Compression level value.
            
        Returns:
            int: Validated compression level.
            
        Raises:
            ValueError: If compression level is not between 1 and 9.
        """
        if not (1 <= v <= 9):
            raise ValueError("Compression level must be between 1 and 9")
        return v
    
    @property
    def output(self) -> RasusaDownsampleReadsOutput:
        """Get the output specification for this stage.
        
        Returns:
            RasusaDownsampleReadsOutput: Object containing paths to output files.
        """
        # Use the same style as FastpCleanReadsOutput
        downsampled_r1 = f"{self.prefix}.downsampled.R1.{self.output_format}"
        if self.output_format == "fastq" and any(f.endswith(".gz") for f in self.reads):
            downsampled_r1 += ".gz"
        
        downsampled_r2 = None
        if len(self.reads) == 2:
            downsampled_r2 = f"{self.prefix}.downsampled.R2.{self.output_format}"
            if self.output_format == "fastq" and any(f.endswith(".gz") for f in self.reads):
                downsampled_r2 += ".gz"
        
        return RasusaDownsampleReadsOutput(
            downsampled_r1=downsampled_r1,
            downsampled_r2=downsampled_r2
        )
    
    def build_rasusa_command(self):
        """Build the rasusa command for downsampling reads.
        
        Constructs the complete rasusa command with all specified options
        including coverage/read count, seed, output format, and file handling.
        Uses the correct rasusa command syntax: rasusa reads [OPTIONS] [INPUT_FILES]
        
        Returns:
            ShellCommand: Complete rasusa command ready for execution.
        """
        cmd_parts = ["rasusa", "reads"]
        
        # Coverage or number of reads
        if self.coverage is not None:
            cmd_parts.extend(["--coverage", str(self.coverage)])
            cmd_parts.extend(["--genome-size", str(self.genome_length)])
        elif self.num_reads is not None:
            cmd_parts.extend(["--num", str(self.num_reads)])
        
        # Output files (one -o flag per output file)
        cmd_parts.extend(["-o", str(self.output.downsampled_r1)])
        if self.output.downsampled_r2:
            cmd_parts.extend(["-o", str(self.output.downsampled_r2)])
        
        # Random seed
        if self.seed is not None:
            cmd_parts.extend(["--seed", str(self.seed)])
        
        # Output format
        if self.output_format == "fasta":
            cmd_parts.append("--fasta")
        
        # Compression level
        if (self.output.downsampled_r1.endswith(".gz") or 
            (self.output.downsampled_r2 and self.output.downsampled_r2.endswith(".gz"))):
            cmd_parts.extend(["--compress-level", str(self.compression_level)])
        
        # Additional options (split if it contains spaces)
        if self.additional_options:
            import shlex
            cmd_parts.extend(shlex.split(self.additional_options))
        
        # Input files (at the end)
        cmd_parts.extend([str(read) for read in self.reads])
        
        coverage_desc = f"coverage {self.coverage}x" if self.coverage else f"{self.num_reads} reads"
        return self.shell_cmd(
            command=cmd_parts,
            description=f"Downsample reads to {coverage_desc} using rasusa"
        )
    
    @property
    def commands(self) -> List:
        """Get the list of commands to execute for this stage.
        
        Returns:
            List: List containing the rasusa command.
        """
        return [self.build_rasusa_command()]


class RasusaDownsampleReadsByCoverage(RasusaDownsampleReads):
    """Downsample reads to a specific coverage depth using rasusa.
    
    This variant is specifically configured for coverage-based downsampling
    and requires both coverage and genome size to be specified.
    
    Args:
        coverage: Target coverage depth (required, no default).
        genome_length: Genome length in base pairs (required for coverage calculation).
        
    Note:
        Inherits all other parameters from RasusaDownsampleReads.
        The num_reads parameter is disabled in this variant.
        
    Example:
        >>> stage = RasusaDownsampleReadsByCoverage(
        ...     reads=["sample_R1.fastq.gz", "sample_R2.fastq.gz"],
        ...     prefix="cov50x",
        ...     genome_length=197394,
        ...     coverage=50.0,
        ...     tmpdir=Path("/tmp")
        ... )
    """
    
    coverage: float = Field(..., description="Target coverage depth for downsampling")
    genome_length: Any = Field(..., description="Genome length in base pairs (required for coverage calculation) - can be int or at_run_time() callable")
    num_reads: Optional[int] = Field(None, description="Disabled in coverage-based variant")
    
    @field_validator("num_reads")
    @classmethod
    def disable_num_reads(cls, v):
        """Ensure num_reads is not used in coverage-based variant.
        
        Args:
            v: num_reads value (should be None).
            
        Returns:
            None: Always returns None.
            
        Raises:
            ValueError: If num_reads is specified.
        """
        if v is not None:
            raise ValueError("Cannot specify num_reads in coverage-based downsampling variant")
        return v


class RasusaDownsampleReadsByCount(RasusaDownsampleReads):
    """Downsample reads to a specific number of reads using rasusa.
    
    This variant is specifically configured for count-based downsampling
    and requires the number of reads to be specified.
    
    Args:
        num_reads: Target number of reads (required, no default).
        
    Note:
        Inherits all other parameters from RasusaDownsampleReads.
        The coverage parameter is disabled in this variant.
        
    Example:
        >>> stage = RasusaDownsampleReadsByCount(
        ...     reads=["sample_R1.fastq.gz", "sample_R2.fastq.gz"],
        ...     prefix="1M_reads",
        ...     num_reads=1000000,
        ...     tmpdir=Path("/tmp")
        ... )
    """
    
    num_reads: int = Field(..., description="Target number of reads for downsampling")
    coverage: Optional[float] = Field(None, description="Disabled in count-based variant")
    
    @field_validator("coverage")
    @classmethod
    def disable_coverage(cls, v):
        """Ensure coverage is not used in count-based variant.
        
        Args:
            v: coverage value (should be None).
            
        Returns:
            None: Always returns None.
            
        Raises:
            ValueError: If coverage is specified.
        """
        if v is not None:
            raise ValueError("Cannot specify coverage in count-based downsampling variant")
        return v
