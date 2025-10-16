from pathlib import Path
from typing import List
from snippy_ng.stages.base import BaseStage, ShellCommandPipeline
from snippy_ng.dependencies import samtools, bwa, samclip, minimap2
from pydantic import Field, field_validator, BaseModel


class AlignerOutput(BaseModel):
    bam: Path

class Aligner(BaseStage):
    reference: Path = Field(..., description="Reference file")
    prefix: str = Field(..., description="Output file prefix")
    maxsoft: int = Field(10, description="Maximum soft clipping to allow")
    aligner_opts: str = Field("", description="Additional options for the aligner")

    @property
    def output(self) -> AlignerOutput:
        return AlignerOutput(bam=self.prefix + ".bam")

    @property
    def common_commands(self) -> List:
        """Common commands for sorting, fixing mates, and marking duplicates."""
        sort_cpus = max(1, int(self.cpus / 2))
        sort_ram = f"{1000 * self.ram // sort_cpus}M"
        sort_threads = str(sort_cpus - 1)
        sort_temp = str(self.tmpdir)
        
        sort_name_cmd = self.shell_cmd([
            "samtools", "sort", "-n", "-l", "0", "-T", sort_temp, 
            "--threads", sort_threads, "-m", sort_ram
        ], description="Sort BAM by read name")
        
        fixmate_cmd = self.shell_cmd([
            "samtools", "fixmate", "-m", "--threads", sort_threads, "-", "-"
        ], description="Fix mate pair information")
        
        sort_coord_cmd = self.shell_cmd([
            "samtools", "sort", "-l", "0", "-T", sort_temp,
            "--threads", sort_threads, "-m", sort_ram
        ], description="Sort BAM by coordinates")
        
        markdup_cmd = self.shell_cmd([
            "samtools", "markdup", "--threads", sort_threads, "-r", "-s", "-", "-"
        ], description="Mark and remove duplicates")
        
        return [sort_name_cmd, fixmate_cmd, sort_coord_cmd, markdup_cmd]

    def build_samclip_command(self):
        """Constructs the samclip command to remove soft-clipped bases."""
        return self.shell_cmd([
            "samclip", "--max", str(self.maxsoft), "--ref", f"{self.reference}.fai"
        ], description="Remove excessive soft-clipped bases")

    def build_alignment_pipeline(self, align_cmd) -> ShellCommandPipeline:
        """Constructs the full alignment pipeline command."""
        samclip_cmd = self.build_samclip_command()
        common_cmds = self.common_commands
        
        # Create pipeline: align_cmd | samclip | sort_name | fixmate | sort_coord | markdup
        pipeline_commands = [align_cmd, samclip_cmd] + common_cmds
        
        return self.shell_pipeline(
            commands=pipeline_commands,
            description="Alignment pipeline: align -> samclip -> sort by name -> fixmate -> sort by coord -> markdup",
            output_file=Path(self.output.bam)
        )

    def build_index_command(self):
        """Returns the samtools index command."""
        return self.shell_cmd([
            "samtools", "index", str(self.output.bam)
        ], description=f"Index BAM file {self.output.bam}")


class BWAMEMReadsAligner(Aligner):
    """
    Align reads to a reference using BWA-MEM.
    """

    reads: List[str] = Field(
        default_factory=list, description="List of input read files"
    )

    @field_validator("reads")
    @classmethod
    def check_reads(cls, v):
        if not v:
            raise ValueError("Reads list must not be empty")
        return v

    _dependencies = [samtools, bwa, samclip]

    @property
    def commands(self) -> List:
        """Constructs the BWA alignment commands."""  
        bwa_index_cmd = self.shell_cmd([
            "bwa", "index", str(self.reference)
        ], description=f"Index reference with BWA: {self.reference}")
        
        # Build BWA mem command
        bwa_cmd_parts = ["bwa", "mem"]
        if self.aligner_opts:
            import shlex
            bwa_cmd_parts.extend(shlex.split(self.aligner_opts))
        bwa_cmd_parts.extend(["-t", str(self.cpus), str(self.reference)])
        bwa_cmd_parts.extend([str(r) for r in self.reads])
        
        bwa_cmd = self.shell_cmd(
            bwa_cmd_parts,
            description=f"Align {len(self.reads)} read files with BWA-MEM"
        )

        alignment_pipeline = self.build_alignment_pipeline(bwa_cmd)
        index_cmd = self.build_index_command()
        
        return [bwa_index_cmd, alignment_pipeline, index_cmd]


class PreAlignedReads(Aligner):
    """
    Use pre-aligned reads in a BAM file.
    """

    bam: Path = Field(..., description="Input BAM file")

    _dependencies = [samtools, samclip]

    @field_validator("bam")
    @classmethod
    def bam_exists(cls, v):
        if not v.exists():
            raise ValueError("BAM file does not exist")
        return v

    @property
    def commands(self) -> List:
        """Constructs the commands to extract reads from a BAM file."""
        view_cmd = self.shell_cmd([
            "samtools", "view", "-h", "-O", "SAM", str(self.bam)
        ], description=f"Extract reads from BAM file: {self.bam}")

        alignment_pipeline = self.build_alignment_pipeline(view_cmd)
        index_cmd = self.build_index_command()
        
        return [alignment_pipeline, index_cmd]
    
class MinimapAligner(Aligner):
    """
    Align reads to a reference using Minimap2.
    """

    reads: List[str] = Field(
        default_factory=list, description="List of input read files"
    )

    @field_validator("reads")
    @classmethod
    def check_reads(cls, v):
        if not v:
            raise ValueError("Reads list must not be empty")
        return v
    
    @property
    def ram_per_thread(self) -> int:
        """Calculate RAM per thread in MB."""
        return max(1, self.ram // self.cpus)

    _dependencies = [minimap2, samtools, samclip]

    @property
    def commands(self) -> List:
        """Constructs the Minimap2 alignment commands."""
        # Build minimap2 command
        minimap_cmd_parts = ["minimap2", "-a"]
        if self.aligner_opts:
            import shlex
            minimap_cmd_parts.extend(shlex.split(self.aligner_opts))
        minimap_cmd_parts.extend(["-t", str(self.cpus), str(self.reference)])
        minimap_cmd_parts.extend([str(r) for r in self.reads])
        
        minimap_cmd = self.shell_cmd(
            minimap_cmd_parts,
            description=f"Align {len(self.reads)} read files with Minimap2"
        )

        alignment_pipeline = self.build_alignment_pipeline(minimap_cmd)
        index_cmd = self.build_index_command()
        
        return [alignment_pipeline, index_cmd]



class AssemblyAlignerOutput(BaseModel):
    paf: Path

class AssemblyAligner(BaseStage):
    """
    Align an assembly to a reference using Minimap2.
    """
    reference: Path = Field(..., description="Reference file")
    assembly: Path = Field(..., description="Input assembly FASTA file")
    prefix: str = Field(..., description="Output file prefix")

    _dependencies = [minimap2]

    @property
    def output(self) -> AssemblyAlignerOutput:
        paf_file = Path(f"{self.prefix}.paf")
        return AssemblyAlignerOutput(
            paf=paf_file
        )

    @property
    def commands(self) -> List:
        """Constructs the Minimap2 alignment commands."""

        minimap_pipeline = self.shell_pipeline(
            commands=[
                self.shell_cmd([
                    "minimap2", "-x", "asm20", "-t", str(self.cpus), "-c", "--cs",
                    str(self.reference), str(self.assembly)
                ], description="Align assembly to reference with Minimap2"),
                self.shell_cmd([
                    "sort", "-k6,6", "-k8,8n"
                ], description="Sort PAF output by reference name and position")
            ],
            description="Align assembly to reference and sort",
            output_file=self.output.paf
        )
        return [minimap_pipeline]