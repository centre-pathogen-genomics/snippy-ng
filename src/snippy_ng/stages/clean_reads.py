from pathlib import Path
from typing import List, Optional
from snippy_ng.stages.base import BaseStage, ShellCommand
from snippy_ng.dependencies import fastp
from pydantic import Field, field_validator, BaseModel


class FastpCleanReadsOutput(BaseModel):
    cleaned_r1: str
    cleaned_r2: Optional[str] = None
    html_report: str
    json_report: str


class FastpCleanReads(BaseStage):
    """
    Clean and filter FASTQ reads using fastp.
    
    This stage removes low-quality reads, trims adapters, and performs
    quality control on paired-end or single-end FASTQ files.
    """
    
    reads: List[str] = Field(..., description="List of input read files (1 or 2 files)")
    prefix: str = Field(..., description="Output file prefix")
    min_length: int = Field(15, description="Minimum read length after trimming")
    quality_cutoff: int = Field(20, description="Quality cutoff for base trimming")
    unqualified_percent_limit: int = Field(20, description="Percentage of unqualified bases allowed")
    n_base_limit: int = Field(5, description="Maximum number of N bases allowed")
    detect_adapter_for_pe: bool = Field(True, description="Auto-detect adapters for paired-end reads")
    correction: bool = Field(True, description="Enable base correction in overlapped regions")
    dedup: bool = Field(False, description="Enable deduplication")
    overrepresentation_analysis: bool = Field(True, description="Enable overrepresentation analysis")
    additional_options: str = Field("", description="Additional fastp options")
    
    _dependencies = [fastp]
    
    @field_validator("reads")
    @classmethod
    def validate_reads(cls, v):
        if not v or len(v) == 0:
            raise ValueError("At least one read file must be provided")
        if len(v) > 2:
            raise ValueError("Maximum of 2 read files (paired-end) supported")
        for read_file in v:
            if not Path(read_file).exists():
                raise ValueError(f"Read file does not exist: {read_file}")
        return v
    
    @property
    def output(self) -> FastpCleanReadsOutput:
        cleaned_r1 = f"{self.prefix}.cleaned.R1.fastq.gz"
        cleaned_r2 = None
        if len(self.reads) == 2:
            cleaned_r2 = f"{self.prefix}.cleaned.R2.fastq.gz"
        
        return FastpCleanReadsOutput(
            cleaned_r1=cleaned_r1,
            cleaned_r2=cleaned_r2,
            html_report=f"{self.prefix}.fastp.html",
            json_report=f"{self.prefix}.fastp.json"
        )
    
    def build_fastp_command(self) -> ShellCommand:
        """Constructs the fastp command for read cleaning."""
        cmd_parts = ["fastp"]
        
        # Input files
        cmd_parts.extend(["-i", str(self.reads[0])])
        if len(self.reads) == 2:
            cmd_parts.extend(["-I", str(self.reads[1])])
        
        # Output files
        cmd_parts.extend(["-o", str(self.output.cleaned_r1)])
        if self.output.cleaned_r2:
            cmd_parts.extend(["-O", str(self.output.cleaned_r2)])
        
        # Reports
        cmd_parts.extend(["-h", str(self.output.html_report)])
        cmd_parts.extend(["-j", str(self.output.json_report)])
        
        # Threading
        if self.cpus > 1:
            cmd_parts.extend(["--thread", str(self.cpus)])
        
        # Quality filtering
        cmd_parts.extend(["--length_required", str(self.min_length)])
        cmd_parts.extend(["--cut_tail_window_size", "4"])
        cmd_parts.extend(["--cut_tail_mean_quality", str(self.quality_cutoff)])
        cmd_parts.extend(["--unqualified_percent_limit", str(self.unqualified_percent_limit)])
        cmd_parts.extend(["--n_base_limit", str(self.n_base_limit)])
        
        # Adapter detection and trimming
        if len(self.reads) == 2 and self.detect_adapter_for_pe:
            cmd_parts.append("--detect_adapter_for_pe")
        
        # Base correction for paired-end overlaps
        if len(self.reads) == 2 and self.correction:
            cmd_parts.append("--correction")
        
        # Deduplication
        if self.dedup:
            cmd_parts.append("--dedup")
        
        # Overrepresentation analysis
        if self.overrepresentation_analysis:
            cmd_parts.append("--overrepresentation_analysis")
        
        # Additional options (split if it contains spaces)
        if self.additional_options:
            import shlex
            cmd_parts.extend(shlex.split(self.additional_options))
        
        read_type = "paired-end" if len(self.reads) == 2 else "single-end"
        return self.shell_cmd(
            command=cmd_parts,
            description=f"Clean and filter {read_type} reads using fastp"
        )
    
    @property
    def commands(self) -> List:
        """Constructs the fastp cleaning command."""
        return [self.build_fastp_command()]


class FastpCleanReadsAggressive(FastpCleanReads):
    """
    Aggressive read cleaning for low-quality samples.
    """
    
    min_length: int = Field(30, description="Longer minimum read length")
    quality_cutoff: int = Field(25, description="Higher quality cutoff")
    unqualified_percent_limit: int = Field(10, description="Lower percentage of unqualified bases")
    n_base_limit: int = Field(2, description="Lower number of N bases allowed")
    dedup: bool = Field(True, description="Enable deduplication by default")


class FastpCleanReadsConservative(FastpCleanReads):
    """
    Conservative read cleaning to retain maximum data.
    """
    
    min_length: int = Field(10, description="Shorter minimum read length")
    quality_cutoff: int = Field(15, description="Lower quality cutoff")
    unqualified_percent_limit: int = Field(30, description="Higher percentage of unqualified bases")
    n_base_limit: int = Field(10, description="Higher number of N bases allowed")
    dedup: bool = Field(False, description="Disable deduplication")
