from pathlib import Path
import sys
from typing import List
from snippy_ng.stages import BaseStage, ShellCommandPipe, BaseOutput
from snippy_ng.dependencies import samtools, bwa, minimap2
from pydantic import Field, field_validator


class AlignerOutput(BaseOutput):
    cram: Path = Field(..., description="Output CRAM file")

class Aligner(BaseStage):
    """
    Base class for read alignment stages. Defines common properties and methods for alignment pipelines.
    """
    reference: Path = Field(..., description="Reference file")
    reads: List[str] = Field(..., default_factory=list, description="List of input read files")
    aligner_opts: str = Field("", description="Additional options for the aligner")

    @property
    def output(self) -> AlignerOutput:
        return AlignerOutput(
            cram=self.prefix + ".cram",
        )

class ShortReadAligner(Aligner):
    """
    Base class for short read alignment pipelines. Implements common steps: sorting, fixing mates, and marking duplicates.
    Optionally filters soft-clipped reads using samclip (recommended only for short reads, as it may discard valid soft clipping in long reads).
    """

    maxsoft: int = Field(10, description="Maximum soft clipping to allow")
    samclip: bool = Field(
        True, description="Whether to run samclip to filter soft-clipped reads"
    )

    _dependencies = [samtools]

    @property
    def common_commands(self) -> List:
        """Common commands for sorting, fixing mates, and marking duplicates."""
        sort_cpus = max(1, int(self.cpus / 2))
        sort_ram = f"{1000 * self.ram // sort_cpus}M"
        sort_threads = str(max(1, sort_cpus - 1))
        sort_temp = str(self.tmpdir)

        commands: List = []

        sort_name_cmd = self.shell_cmd(
            [
                "samtools",
                "sort",
                "-n",
                "-O",
                "sam",
                "-T",
                sort_temp,
                "--threads",
                sort_threads,
                "-m",
                sort_ram,
            ],
            description="Sort BAM by read name",
        )
        commands.append(sort_name_cmd)

        if self.samclip:
            samclip_cmd = self.shell_cmd(
                [
                    sys.executable,
                    "-m",
                    "snippy_ng",
                    "samclip",
                    "--index",
                    f"{self.reference}.fai",
                    "--max",
                    str(self.maxsoft),
                ],
                description="Filter alignments with excessive soft clipping",
            )
            commands.append(samclip_cmd)

        fixmate_cmd = self.shell_cmd(
            ["samtools", "fixmate", "--threads", sort_threads, "-m", "-u", "-", "-"],
            description="Fix mate pair information",
        )
        commands.append(fixmate_cmd)

        sort_coord_cmd = self.shell_cmd(
            [
                "samtools",
                "sort",
                "-u",
                "-T",
                sort_temp,
                "--threads",
                sort_threads,
                "-m",
                sort_ram,
            ],
            description="Sort BAM by coordinates",
        )
        commands.append(sort_coord_cmd)

        markdup_cmd = self.shell_cmd(
            [
                "samtools",
                "markdup",
                "-O",
                "cram,embed_ref=2",
                "--threads",
                sort_threads,
                "-r",
                "-s",
                "-",
                "-",
            ],
            description="Mark and remove duplicates",
        )
        commands.append(markdup_cmd)

        return commands

    def build_alignment_pipeline(self, align_cmd) -> ShellCommandPipe:
        """Constructs the full alignment pipeline command."""
        common_cmds = self.common_commands

        pipeline_commands = [align_cmd] + common_cmds

        return self.shell_pipeline(
            commands=pipeline_commands,
            description="Alignment pipeline: align, filter, fix mates, sort, mark duplicates",
            output_file=Path(self.output.cram),
        )


class BWAMEMShortReadAligner(ShortReadAligner):
    """
    Align reads to a reference using BWA-MEM.
    """

    reference_index: Path = Field(..., description="BWA index file for the reference")

    _dependencies = [samtools, bwa]

    @property
    def commands(self) -> List:
        """Constructs the BWA alignment commands."""
        bwa_index_cmd = self.shell_cmd(
            ["bwa", "index", str(self.reference)],
            description=f"Index reference with BWA: {self.reference}",
        )

        # Build BWA mem command
        bwa_cmd_parts = ["bwa", "mem"]
        if self.aligner_opts:
            import shlex

            bwa_cmd_parts.extend(shlex.split(self.aligner_opts))
        bwa_cmd_parts.extend(["-t", str(self.cpus), str(self.reference)])
        bwa_cmd_parts.extend([str(r) for r in self.reads])

        bwa_cmd = self.shell_cmd(
            bwa_cmd_parts,
            description=f"Align {len(self.reads)} read files with BWA-MEM",
        )

        alignment_pipeline = self.build_alignment_pipeline(bwa_cmd)
        return [bwa_index_cmd, alignment_pipeline]


class Minimap2ShortReadAligner(ShortReadAligner):
    """
    Align reads to a reference using Minimap2.
    """
    _dependencies = [minimap2, samtools]

    @property
    def ram_per_thread(self) -> int:
        """Calculate RAM per thread in MB."""
        return max(1, self.ram // self.cpus)


    @property
    def commands(self) -> List:
        """Constructs the Minimap2 alignment commands."""
        # Build minimap2 command
        minimap_cmd_parts = ["minimap2", "-a", "-x", "sr"]
        if self.aligner_opts:
            import shlex

            minimap_cmd_parts.extend(shlex.split(self.aligner_opts))
        minimap_cmd_parts.extend(["-t", str(self.cpus), str(self.reference)])
        minimap_cmd_parts.extend([str(r) for r in self.reads])

        minimap_cmd = self.shell_cmd(
            minimap_cmd_parts,
            description=f"Align {len(self.reads)} read files with Minimap2",
        )

        alignment_pipeline = self.build_alignment_pipeline(minimap_cmd)
        return [alignment_pipeline]


class Minimap2LongReadAligner(Aligner):
    """
    Align reads to a reference using Minimap2.
    """
    minimap_preset: str = Field(
        "map-ont", description="Minimap2 preset to use for alignment"
    )
    
    _dependencies = [minimap2, samtools]

    @property
    def commands(self) -> List:
        """Constructs the Minimap2 alignment commands."""
        # Build minimap2 command
        minimap_cmd_parts = [
            "minimap2",
            "-a",
            "-L",
            "--cs",
            "--MD",
            "-x",
            self.minimap_preset,
        ]
        if self.aligner_opts:
            import shlex

            minimap_cmd_parts.extend(shlex.split(self.aligner_opts))
        minimap_cmd_parts.extend(["-t", str(self.cpus), str(self.reference)])
        minimap_cmd_parts.extend([str(r) for r in self.reads])

        minimap_pipeline = self.shell_pipeline(
            [
                self.shell_cmd(
                    minimap_cmd_parts,
                    description=f"Align {len(self.reads)} read files with Minimap2",
                ),
                self.shell_cmd(
                    ["samtools", "sort", "--threads", str(self.cpus), "-O",
                "cram,embed_ref=2"],
                    description="Sort and convert to CRAM",
                ),
            ],
            description="Minimap2 alignment pipeline",
            output_file=self.output.cram,
        )

        return [minimap_pipeline]


class AssemblyAlignerOutput(BaseOutput):
    paf: Path


class AssemblyAligner(BaseStage):
    """
    Align an assembly to a reference using Minimap2.
    """

    reference: Path = Field(..., description="Reference file")
    assembly: Path = Field(..., description="Input assembly FASTA file")

    _dependencies = [minimap2]

    @property
    def output(self) -> AssemblyAlignerOutput:
        paf_file = Path(f"{self.prefix}.paf")
        return AssemblyAlignerOutput(paf=paf_file)

    @property
    def commands(self) -> List:
        """Constructs the Minimap2 alignment commands."""

        minimap_pipeline = self.shell_pipeline(
            commands=[
                self.shell_cmd(
                    [
                        "minimap2",
                        "-x",
                        "asm20",
                        "-t",
                        str(self.cpus),
                        "-c",
                        "--cs",
                        str(self.reference),
                        str(self.assembly),
                    ],
                    description="Align assembly to reference with Minimap2",
                ),
                self.shell_cmd(
                    ["sort", "-k6,6", "-k8,8n"],
                    description="Sort PAF output by reference name and position",
                ),
            ],
            description="Align assembly to reference and sort",
            output_file=self.output.paf,
        )
        return [minimap_pipeline]

    def test_paf_output(self):
        """Test that the PAF output file was created and is not empty."""
        paf_path = self.output.paf
        if not paf_path.exists():
            raise FileNotFoundError(
                f"Expected PAF output file {paf_path} was not created"
            )
        if paf_path.stat().st_size == 0:
            raise ValueError(
                f"PAF output file {paf_path} is empty, expected alignment results. Did you use the correct reference?"
            )
