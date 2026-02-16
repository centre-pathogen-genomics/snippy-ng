from pathlib import Path

from pydantic import Field
from snippy_ng.stages import BaseOutput, BaseStage
from snippy_ng.dependencies import iqtree


class IQTreeBuildTreeOutput(BaseOutput):
    bionj: Path
    checkpoint: Path
    contree: Path
    iqtree: Path
    log: Path
    mldist: Path
    splits: Path
    tree: Path


class IQTreeBuildTree(BaseStage):
    """
    Filter a multiple sequence alignment to retain only positions present
    in a specified fraction of samples (soft core).
    """

    aln: Path = Field(..., description="Input multiple sequence alignment file in FASTA format")
    model: str = "GTR+G4"
    bootstrap: int = 1000
    fconst: str | None = None

    _dependencies = [iqtree]

    @property
    def output(self) -> IQTreeBuildTreeOutput:
        prefix = Path(self.prefix)
        return IQTreeBuildTreeOutput(
            bionj=prefix.with_suffix(".bionj"),
            checkpoint=prefix.with_suffix(".ckp.gz"),
            contree=prefix.with_suffix(".contree"),
            iqtree=prefix.with_suffix(".iqtree"),
            log=prefix.with_suffix(".log"),
            mldist=prefix.with_suffix(".mldist"),
            splits=prefix.with_suffix(".splits.nex"),
            tree=prefix.with_suffix(".treefile"),
        )

    @property
    def commands(self):
        iqtree_cmd = self.shell_cmd(
            [
                "iqtree",
                "-s", str(self.aln),
                "-pre", str(self.prefix),
                "-m", self.model,
                "-bb", str(self.bootstrap),
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
        return [iqtree_cmd]
    