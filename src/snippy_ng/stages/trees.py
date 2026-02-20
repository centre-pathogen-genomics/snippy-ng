from pathlib import Path
from typing import Optional

from pydantic import Field
from snippy_ng.stages import BaseOutput, BaseStage
from snippy_ng.dependencies import iqtree


class IQTreeBuildTreeOutput(BaseOutput):
    bionj: Path
    checkpoint: Path
    contree: Optional[Path]
    iqtree: Path
    log: Path
    mldist: Path
    splits: Optional[Path]
    tree: Path


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

    _dependencies = [iqtree]

    @property
    def output(self) -> IQTreeBuildTreeOutput:
        prefix = Path(self.prefix)
        contree = None
        splits = None
        if not self.fast_mode:
            contree=prefix.with_suffix(".contree")
            splits=prefix.with_suffix(".splits.nex")

        return IQTreeBuildTreeOutput(
            bionj=prefix.with_suffix(".bionj"),
            checkpoint=prefix.with_suffix(".ckp.gz"),
            contree=contree,
            iqtree=prefix.with_suffix(".iqtree"),
            log=prefix.with_suffix(".log"),
            mldist=prefix.with_suffix(".mldist"),
            splits=splits,
            tree=prefix.with_suffix(".treefile"),
        )

    def create_commands(self, ctx):
        iqtree_cmd = self.shell_cmd(
            [
                "iqtree",
                "-s", str(self.aln),
                "-pre", str(self.prefix),
                "-m", self.model, 
                "-T", "AUTO",
                "--threads-max", str(self.cpus),
                "-nt", "AUTO",
                "-st", "DNA",
                "-redo"
            ], 
            description="Building phylogenetic tree with IQ-TREE"
        )
        if self.ram:
            iqtree_cmd.command.extend(["-mem", f"{self.ram}G"])
        if self.fconst:
            iqtree_cmd.command.extend(["-fconst", self.fconst])
        if self.fast_mode:
            iqtree_cmd.command.append("--fast")
        else:
            # Ultrafast bootstrap (-bb) does not work with -fast option
            iqtree_cmd.command.extend(["-bb", str(self.bootstrap)])
        return [iqtree_cmd]
    