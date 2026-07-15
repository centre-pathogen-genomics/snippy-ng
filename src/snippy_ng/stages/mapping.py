from pathlib import Path
import sys
from typing import List, Optional
from snippy_ng.stages import BaseStage, ShellProcessPipe, BaseOutput
from snippy_ng.dependencies import samtools, bwa, minimap2, nucmer
from snippy_ng.envvars import EnvVarField
from pydantic import Field 


class AlignerOutput(BaseOutput):
    bam: Path = Field(..., description="Coordinate-sorted, deduplicated alignment file in BAM format")

class Aligner(BaseStage):
    """
    Base class for read alignment stages. Defines common properties and methods for alignment pipelines.
    """
    reference: Path = Field(..., description="Reference file")
    reads: List[Path] = Field(..., default_factory=list, description="List of input read files")
    aligner_opts: str = Field("", description="Additional options for the aligner")

    @property
    def output(self) -> AlignerOutput:
        return AlignerOutput(
            bam=self.prefix + ".bam",
        )

class ShortReadAligner(Aligner):
    """
    Base class for short read alignment pipelines. Implements common steps: sorting, fixing mates, and marking duplicates.
    Optionally filters soft-clipped reads using samclip (recommended only for short reads, as it may discard valid soft clipping in long reads).
    """

    samclip: bool = Field(
        True, description="Whether to run samclip to filter soft-clipped reads"
    )
    maxsoft: int = Field(10, description="Maximum soft clipping to allow")

    _dependencies = [samtools]

    def common_commands(self, ctx) -> List:
        """Common commands for sorting, fixing mates, and marking duplicates."""
        sort_cpus = max(1, int(ctx.cpus / 2))
        sort_ram_param = ["-m", f"{1000 * ctx.ram // sort_cpus}M"] if ctx.ram else []
        sort_threads = str(max(1, sort_cpus - 1))
        sort_temp = str(ctx.tmpdir or ".")

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
                *sort_ram_param
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
                    "utils",
                    "aln",
                    "samclip",
                    "--index",
                    f"{self.reference}.fai",
                    "--max",
                    str(self.maxsoft),
                    "--fix-mate",
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
                *sort_ram_param,
            ],
            description="Sort BAM by coordinates",
        )
        commands.append(sort_coord_cmd)

        markdup_cmd = self.shell_cmd(
            [
                "samtools",
                "markdup",
                "-O",
                "bam",
                "--reference",
                str(self.reference),
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

    def build_alignment_pipeline(self, align_cmd, ctx) -> ShellProcessPipe:
        """Constructs the full alignment pipeline command."""
        common_cmds = self.common_commands(ctx)

        pipeline_commands = [align_cmd] + common_cmds

        return self.shell_pipe(
            commands=pipeline_commands,
            description="Alignment pipeline: align, name-sort, filter soft-clipped reads, fix mates, coordinate-sort, mark duplicates",
            output_file=Path(self.output.bam),
        )


class BWAMEMShortReadAligner(ShortReadAligner):
    """
    Align reads to a reference using BWA-MEM.
    """

    reference_index: Path = Field(..., description="BWA index file for the reference")

    _dependencies = [samtools, bwa]

    def create_commands(self, ctx) -> List:
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
        bwa_cmd_parts.extend(["-t", str(ctx.cpus), str(self.reference)])
        bwa_cmd_parts.extend([str(r) for r in self.reads])

        bwa_cmd = self.shell_cmd(
            bwa_cmd_parts,
            description=f"Align {len(self.reads)} read files with BWA-MEM",
        )

        alignment_pipeline = self.build_alignment_pipeline(bwa_cmd, ctx)
        return [bwa_index_cmd, alignment_pipeline]


class Minimap2ShortReadAligner(ShortReadAligner):
    """
    Align reads to a reference using Minimap2.
    """
    _dependencies = [minimap2, samtools]

    def create_commands(self, ctx) -> List:
        """Constructs the Minimap2 alignment commands."""
        # Build minimap2 command
        minimap_cmd_parts = ["minimap2", "-a", "-x", "sr"]
        if self.aligner_opts:
            import shlex

            minimap_cmd_parts.extend(shlex.split(self.aligner_opts))
        minimap_cmd_parts.extend(["-t", str(ctx.cpus), str(self.reference)])
        minimap_cmd_parts.extend([str(r) for r in self.reads])

        minimap_cmd = self.shell_cmd(
            minimap_cmd_parts,
            description=f"Align {len(self.reads)} read files with Minimap2",
        )

        alignment_pipeline = self.build_alignment_pipeline(minimap_cmd, ctx)
        return [alignment_pipeline]


class Minimap2LongReadAligner(Aligner):
    """
    Align reads to a reference using Minimap2.
    """
    minimap_preset: str = Field(
        "map-ont", description="Minimap2 preset to use for alignment"
    )
    reference_index: Path = Field(..., description="Reference FASTA index (.fai)")
    max_clip_fraction: Optional[float] = Field(
        None,
        ge=0,
        le=1,
        description="Optional maximum terminal clipping fraction for samclip",
    )
    
    _dependencies = [minimap2, samtools]

    def create_commands(self, ctx) -> List:
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
        minimap_cmd_parts.extend(["-t", str(ctx.cpus), str(self.reference)])
        minimap_cmd_parts.extend([str(r) for r in self.reads])

        processes = [
            self.shell_cmd(
                minimap_cmd_parts,
                description=f"Align {len(self.reads)} read files with Minimap2",
            ),
        ]
        if self.max_clip_fraction is not None:
            processes.extend(
                [
                    self.shell_cmd(
                        [
                            "samtools",
                            "sort",
                            "-n",
                            "-O",
                            "sam",
                            "--threads",
                            str(ctx.cpus),
                        ],
                        description="Name-sort alignments for samclip",
                    ),
                    self.shell_cmd(
                        [
                            sys.executable,
                            "-m",
                            "snippy_ng",
                            "utils",
                            "aln",
                            "samclip",
                            "--index",
                            str(self.reference_index),
                            "--max-clip-fraction",
                            str(self.max_clip_fraction),
                        ],
                        description="Filter long-read alignments by clipped-read fraction",
                    ),
                ]
            )
        processes.append(
            self.shell_cmd(
                [
                    "samtools",
                    "sort",
                    "--threads",
                    str(ctx.cpus),
                    "-O",
                    "bam",
                    "--reference",
                    str(self.reference),
                ],
                description="Coordinate-sort and convert to BAM",
            )
        )

        minimap_pipeline = self.shell_pipe(
            processes,
            description=(
                "Minimap2 alignment pipeline with fraction-based samclip"
                if self.max_clip_fraction is not None
                else "Minimap2 alignment pipeline"
            ),
            output_file=self.output.bam,
        )

        return [minimap_pipeline]


class AssemblyAlignerOutput(BaseOutput):
    paf: Path = Field(..., description="Sorted pairwise alignment file (PAF) of assembly-to-reference alignments")


class AssemblyAligner(BaseStage):
    """
    Align an assembly to a reference using Minimap2.
    """

    reference: Path = Field(..., description="Reference file")
    assembly: Path = Field(..., description="Input assembly FASTA file")
    minimap_preset: str = EnvVarField(
        "asm20",
        "ASM_MINIMAP_PRESET",
        description="Minimap2 preset to use for assembly-to-reference alignment",
    )

    _dependencies = [minimap2]

    @property
    def output(self) -> AssemblyAlignerOutput:
        paf_file = Path(f"{self.prefix}.paf")
        return AssemblyAlignerOutput(paf=paf_file)

    def create_commands(self, ctx) -> List:
        """Constructs the Minimap2 alignment commands."""

        minimap_pipeline = self.shell_pipe(
            commands=[
                self.shell_cmd(
                    [
                        "minimap2",
                        "-x",
                        self.minimap_preset,
                        "-t",
                        str(ctx.cpus),
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


class AssemblyNucmerAlignerOutput(BaseOutput):
    delta: Path = Field(..., description="Delta file of assembly-to-reference alignments produced by nucmer")
    assembly: Optional[Path] = Field(None, description="Uncompressed assembly FASTA used as the nucmer query")


class AssemblyNucmerAligner(BaseStage):
    """
    Align an assembly to a reference using nucmer from MUMmer.
    """

    reference: Path = Field(..., description="Reference file")
    assembly: Path = Field(..., description="Input assembly FASTA file")
    breaklen: int = EnvVarField(250, "MUMMER_BREAKLEN", description="Maximum poor-scoring extension distance for nucmer")
    mincluster: int = EnvVarField(120, "MUMMER_MINCLUSTER", description="Minimum length of a match cluster for nucmer")
    maxgap: int = EnvVarField(50, "MUMMER_MAXGAP", description="Maximum gap between adjacent matches in a cluster for nucmer")
    minmatch: int = EnvVarField(46, "MUMMER_MINMATCH", description="Minimum exact-match length for nucmer anchors")
    minalign: int = EnvVarField(400, "MUMMER_MINALIGN", description="Minimum alignment length retained by nucmer after extension")

    _dependencies = [nucmer]

    @property
    def output(self) -> AssemblyNucmerAlignerOutput:
        assembly = Path(f"{self.prefix}.assembly.fa") if self.assembly.suffix.lower() == ".gz" else None
        return AssemblyNucmerAlignerOutput(
            delta=Path(f"{self.prefix}.delta"),
            assembly=assembly,
        )

    @property
    def nucmer_assembly(self) -> Path:
        return self.output.assembly or self.assembly

    def create_commands(self, ctx) -> List:
        commands = []
        if self.output.assembly is not None:
            commands.append(
                self.shell_cmd(
                    ["gunzip", "-c", str(self.assembly)],
                    description="Decompress assembly FASTA for nucmer",
                    output_file=self.output.assembly,
                )
            )
        commands.append(
            self.shell_cmd(
                [
                    "nucmer",
                    "--prefix",
                    self.prefix,
                    "--threads",
                    str(ctx.cpus),
                    "--breaklen",
                    str(self.breaklen),
                    "--mincluster",
                    str(self.mincluster),
                    "--maxgap",
                    str(self.maxgap),
                    "--minmatch",
                    str(self.minmatch),
                    "--minalign",
                    str(self.minalign),
                    str(self.reference),
                    str(self.nucmer_assembly),
                ],
                description="Align assembly to reference with nucmer",
            )
        )
        return commands

    def test_delta_output(self):
        delta_path = self.output.delta
        if not delta_path.exists():
            raise FileNotFoundError(
                f"Expected delta output file {delta_path} was not created"
            )
        if delta_path.stat().st_size == 0:
            raise ValueError(
                f"Delta output file {delta_path} is empty, expected alignment results. Did you use the correct reference?"
            )
