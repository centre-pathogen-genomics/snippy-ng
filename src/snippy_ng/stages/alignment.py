from pathlib import Path
from typing import List
from snippy_ng.stages.base import BaseStage
from snippy_ng.dependencies import samtools, bwa, samclip
from pydantic import Field, field_validator, BaseModel


class AlignerOutput(BaseModel):
    bam: str

class Aligner(BaseStage):
    cpus: int = Field(1, description="Number of CPU cores to use")
    ram: int = Field(4, description="RAM in GB to use")
    reference: Path = Field(..., description="Reference file")
    prefix: str = Field(..., description="Output file prefix")
    tmpdir: Path = Field("/tmp", description="Temporary directory")
    maxsoft: int = Field(10, description="Maximum soft clipping to allow")

    @property
    def output(self) -> AlignerOutput:
        return AlignerOutput(bam=self.prefix + ".bam")

    @property
    def common_commands(self) -> List[str]:
        """Common commands for sorting, fixing mates, and marking duplicates."""
        sort_cpus = max(1, int(self.cpus / 2))
        sort_ram = f"{1000 * self.ram // sort_cpus}M"
        sort_cpus = f"--threads {sort_cpus - 1}"
        sort_temp = f"-T {self.tmpdir}"
        sort_options = f"-l 0 {sort_temp} {sort_cpus} -m {sort_ram}"

        sort_name_cmd = f"samtools sort -n {sort_options}"
        fixmate_cmd = f"samtools fixmate -m {sort_cpus} - -"
        sort_cord_cmd = f"samtools sort {sort_options}"
        markdup_cmd = f"samtools markdup {sort_temp} {sort_cpus} -r -s - -"
        return [sort_name_cmd, fixmate_cmd, sort_cord_cmd, markdup_cmd]

    def build_samclip_command(self) -> str:
        """Constructs the samclip command to remove soft-clipped bases."""
        return f"samclip --max {self.maxsoft} --ref {self.reference}.fai"

    def build_alignment_command(self, align_cmd: str) -> str:
        """Constructs the full alignment pipeline command."""
        samclip_cmd = self.build_samclip_command()
        common_cmds = " | ".join(self.common_commands)
        return f"{align_cmd} | {samclip_cmd} | {common_cmds} > {self.output.bam}"

    def build_index_command(self) -> str:
        """Returns the samtools index command."""
        return f"samtools index {self.output.bam}"


class BWAMEMReadsAligner(Aligner):
    """
    Align reads to a reference using BWA-MEM.
    """

    reads: List[str] = Field(
        default_factory=list, description="List of input read files"
    )
    bwaopt: str = Field("", description="Additional BWA options")

    @field_validator("reads")
    @classmethod
    def check_reads(cls, v):
        if not v:
            raise ValueError("Reads list must not be empty")
        return v

    _dependencies = [samtools, bwa, samclip]

    @property
    def commands(self) -> List[str]:
        """Constructs the BWA alignment commands."""
        fasta_index = f"samtools faidx {self.reference}"
        bwa_index_cmd = f"bwa index {self.reference}"
        bwa_cmd = f"bwa mem {self.bwaopt} -t {self.cpus} {self.reference} {' '.join(self.reads)}"

        full_cmd = self.build_alignment_command(bwa_cmd)
        index_cmd = self.build_index_command()
        return [fasta_index, bwa_index_cmd, full_cmd, index_cmd]


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
    def commands(self) -> List[str]:
        """Constructs the commands to extract reads from a BAM file."""
        fasta_index = f"samtools faidx {self.reference}"
        view_cmd = f"samtools view -h -O SAM {self.bam}"

        full_cmd = self.build_alignment_command(view_cmd)
        index_cmd = self.build_index_command()
        return [fasta_index, full_cmd, index_cmd]
