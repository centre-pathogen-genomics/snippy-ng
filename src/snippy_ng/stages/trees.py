from snippy_ng.context import Context
from pathlib import Path
from typing import Optional

from Bio import Phylo, SeqIO
from pydantic import Field
from snippy_ng.stages import BaseOutput, BaseStage
from snippy_ng.dependencies import biopython, clonalframeml, iqtree
from snippy_ng.exceptions import StageExecutionError

class TreeBuildingError(StageExecutionError):
    """Custom exception for errors during tree building stage execution."""
    pass


class IQTreeBuildTreeOutput(BaseOutput):
    bionj: Optional[Path] = Field(None, description="BIONJ starting tree generated outside pathogen mode")
    checkpoint: Optional[Path] = Field(None, description="IQ-TREE checkpoint file generated outside pathogen mode")
    contree: Optional[Path] = Field(None, description="Consensus tree from bootstrap replicates (not produced in fast mode)")
    iqtree: Optional[Path] = Field(None, description="IQ-TREE run summary generated outside pathogen mode")
    log: Path = Field(..., description="IQ-TREE execution log file")
    mldist: Optional[Path] = Field(None, description="Pairwise ML distance matrix generated outside pathogen mode")
    splits: Optional[Path] = Field(None, description="Split support file in NEXUS format (not produced in fast mode)")
    tree: Path = Field(..., description="Output tree file in Newick format")


class IQTreeBuildTree(BaseStage):
    """
    Filter a multiple sequence alignment to retain only positions present
    in a specified fraction of samples (soft core).
    """

    aln: Path = Field(..., description="Input multiple sequence alignment file in FASTA format")
    model: str = Field(default="GTR+G", description="Substitution model for IQ-TREE")
    bootstrap: int = Field(default=1000, description="Number of ultrafast bootstrap replicates to perform")
    fconst: Optional[str] = Field(default=None, description="Constant patterns to add into alignment in format a,c,g,t (e.g. 10,5,0,0)")
    fast_mode: bool = Field(default=False, description="Use fast mode for rough tree search (faster but less accurate)")
    cmaple: bool = Field(default=True, description="Use pathogen mode for IQ-TREE with --alrt for SH-aLRT support values")

    _dependencies = [iqtree]

    @property
    def output(self) -> IQTreeBuildTreeOutput:
        prefix = str(self.prefix)
        bionj = None
        checkpoint = None
        contree = None
        iqtree = None
        mldist = None
        splits = None

        if not self.cmaple:
            bionj = Path(f"{prefix}.bionj")
            checkpoint = Path(f"{prefix}.ckp.gz")
            iqtree = Path(f"{prefix}.iqtree")
            mldist = Path(f"{prefix}.mldist")

        if not self.fast_mode and not self.cmaple:
            contree = Path(f"{prefix}.contree")
            splits = Path(f"{prefix}.splits.nex")

        return IQTreeBuildTreeOutput(
            bionj=bionj,
            checkpoint=checkpoint,
            contree=contree,
            iqtree=iqtree,
            log=Path(f"{prefix}.log"),
            mldist=mldist,
            splits=splits,
            tree=Path(f"{prefix}.treefile"),
        )

    def create_commands(self, ctx):
        validate_aln_cmd = self.python_cmd(
            self.check_aln_has_three_samples,
            description="Validating that the input alignment contains at least 3 samples"
        )
        iqtree_cmd = self.shell_cmd(
            [
                "iqtree",
                "-s", str(self.aln),
                "-pre", str(self.prefix),
                "-m", self.model, 
                "-T", str(ctx.cpus),
                "--threads-max", str(ctx.cpus),
                "-st", "DNA",
                "-redo"
            ], 
            description="Building phylogenetic tree with IQ-TREE"
        )
        if ctx.ram:
            iqtree_cmd.command.extend(["-mem", f"{ctx.ram}G"])
        if self.fconst:
            iqtree_cmd.command.extend(["-fconst", self.fconst])
        if self.cmaple:
            iqtree_cmd.command.extend(["--pathogen", "--alrt", str(self.bootstrap)])
        if self.fast_mode:
            iqtree_cmd.command.append("--fast")
        else:
            # Ultrafast bootstrap (-bb) does not work with -fast option
            iqtree_cmd.command.extend(["-bb", str(self.bootstrap)])
        return [validate_aln_cmd, iqtree_cmd]

    def check_aln_has_three_samples(self):
        with open(self.aln) as f:
            num_samples = 0
            for line in f:
                if num_samples >= 3:
                    break
                if line.startswith(">"):
                    num_samples += 1
        if num_samples < 3:
            raise TreeBuildingError(f"Alignment must contain at least 3 samples to build a tree, but found only {num_samples} in {self.aln}") 


class ClonalFrameMLCorrectTreeOutput(BaseOutput):
    recombination_corrected_tree: Path = Field(..., description="ClonalFrameML tree with labelled internal nodes and recombination-corrected branches")
    importation_status: Path = Field(..., description="Recombination events inferred on each tree branch")
    em: Path = Field(..., description="Maximum-likelihood recombination parameter estimates")
    emsim: Optional[Path] = Field(None, description="Simulated uncertainty estimates when EM simulations are requested")
    ml_sequence: Path = Field(..., description="Maximum-likelihood ancestral and imputed sequences")
    position_cross_reference: Path = Field(..., description="Mapping from input alignment sites to reconstructed sequence patterns")

class ClonalFrameMLCorrectTree(BaseStage):
    """Infer recombination and correct the branch lengths of a starting tree."""

    tree: Path = Field(..., description="Starting tree in Newick format")
    aln: Path = Field(..., description="FASTA alignment corresponding to the leaves of the starting tree")
    kappa: float = Field(2.0, gt=0.0, description="Transition/transversion bias")
    emsim: int = Field(0, ge=0, description="Number of simulations used to estimate uncertainty")
    ignore_incomplete_sites: bool = Field(
        default=False,
        description="Ignore alignment sites containing ambiguous or missing bases",
    )

    _dependencies = [clonalframeml]

    @property
    def output(self) -> ClonalFrameMLCorrectTreeOutput:
        prefix = str(self.prefix)
        return ClonalFrameMLCorrectTreeOutput(
            recombination_corrected_tree=Path(f"{prefix}.labelled_tree.newick"),
            importation_status=Path(f"{prefix}.importation_status.txt"),
            em=Path(f"{prefix}.em.txt"),
            emsim=Path(f"{prefix}.emsim.txt") if self.emsim else None,
            ml_sequence=Path(f"{prefix}.ML_sequence.fasta"),
            position_cross_reference=Path(f"{prefix}.position_cross_reference.txt")
        )

    def create_commands(self, ctx: Context):
        command = [
            "ClonalFrameML",
            str(self.tree),
            str(self.aln),
            str(self.prefix),
            "-kappa",
            str(self.kappa),
            "-num_threads",
            str(ctx.cpus),
            "-ignore_incomplete_sites",
            str(self.ignore_incomplete_sites).lower(),
            "-show_progress",
            str(ctx.debug).lower(),
        ]
        if self.emsim:
            command.extend(["-emsim", str(self.emsim)])
        return [
            self.shell_cmd(
                command,
                description="Infer recombination and correct tree branch lengths with ClonalFrameML",
            )
        ]


class ScaleTreeToSNPsOutput(BaseOutput):
    snp_scaled_tree: Path = Field(..., description="Tree with branch lengths scaled to expected SNP counts")


class ScaleTreeToSNPs(BaseStage):
    """Scale Newick branch lengths from substitutions per site to expected SNP counts."""

    tree: Path = Field(..., description="Input tree with branch lengths in substitutions per site")
    aln: Path = Field(..., description="Original alignment used to build the tree")
    fconst: Optional[str] = Field(default=None, description="Constant patterns added to the alignment for IQ-TREE")

    _dependencies = [biopython]

    @property
    def output(self) -> ScaleTreeToSNPsOutput:
        return ScaleTreeToSNPsOutput(snp_scaled_tree=Path(f"{self.prefix}.snps.newick"))

    def create_commands(self, ctx):
        return [
            self.python_cmd(
                func=self.scale_tree_to_snps,
                args=[self.tree, self.aln, self.output.snp_scaled_tree, self.fconst],
                description="Scale tree branch lengths to expected SNP counts",
            )
        ]

    @staticmethod
    def count_constant_sites(fconst: Optional[str]) -> int:
        if not fconst:
            return 0

        fconst_value = str(fconst).strip()

        try:
            return sum(int(count.strip()) for count in fconst_value.split(","))
        except ValueError as exc:
            raise TreeBuildingError(
                f"Invalid constant site counts for fconst: {fconst}"
            ) from exc

    @classmethod
    def count_alignment_sites(cls, aln: Path, fconst: Optional[str] = None) -> int:
        expected_length = None
        sequence_count = 0
        for record in SeqIO.parse(str(aln), "fasta"):
            sequence_count += 1
            sequence_length = len(record.seq)
            if expected_length is None:
                expected_length = sequence_length
            elif sequence_length != expected_length:
                raise TreeBuildingError(
                    f"Alignment sequences must all have the same length in {aln}: "
                    f"expected {expected_length}, found {sequence_length} for {record.id}"
                )

        if sequence_count == 0 or expected_length is None:
            raise TreeBuildingError(f"No sequences found in alignment: {aln}")

        return expected_length + cls.count_constant_sites(fconst)

    @classmethod
    def scale_tree_to_snps(
        cls,
        tree: Path,
        aln: Path,
        output_tree: Path,
        fconst: Optional[str] = None,
    ) -> None:
        alignment_sites = cls.count_alignment_sites(aln, fconst)
        scaled_tree = Phylo.read(str(tree), "newick")

        for clade in scaled_tree.find_clades():
            if clade.branch_length is not None:
                # Scale branch length from substitutions per site to expected SNP count
                # round to nearest integer since we expect SNP counts to be whole numbers
                clade.branch_length = round(clade.branch_length * alignment_sites)

        with output_tree.open("w") as handle:
            Phylo.write(
                scaled_tree,
                handle,
                "newick",
                format_branch_length="%1.0f",
            )


class TreeDistanceMatrixOutput(BaseOutput):
    distance: Path = Field(..., description="Pairwise tree distances in TSV format")


class TreeDistanceMatrix(BaseStage):
    """Write pairwise terminal distances from a Newick tree as long-form TSV."""

    tree: Path = Field(..., description="Input Newick tree")

    _dependencies = [biopython]

    @property
    def output(self) -> TreeDistanceMatrixOutput:
        return TreeDistanceMatrixOutput(distance=Path(f"{self.prefix}.distance.tsv"))

    def create_commands(self, ctx):
        return [
            self.python_cmd(
                func=self.write_distance_matrix,
                args=[self.tree, self.output.distance],
                description="Calculate pairwise distances from the tree",
            )
        ]

    @staticmethod
    def write_distance_matrix(tree: Path, output_distance: Path) -> None:
        phylo_tree = Phylo.read(str(tree), "newick")
        terminals = phylo_tree.get_terminals()
        names = [terminal.name for terminal in terminals]

        if any(name is None for name in names):
            raise TreeBuildingError(f"All tree terminals must have names: {tree}")

        with output_distance.open("w") as handle:
            for index, terminal in enumerate(terminals[1:], start=1):
                for other, other_name in zip(terminals[:index], names[:index]):
                    distance = phylo_tree.distance(terminal, other)
                    distance_text = (
                        str(int(distance)) if distance.is_integer() else str(distance)
                    )
                    handle.write(f"{names[index]}\t{other_name}\t{distance_text}\n")
