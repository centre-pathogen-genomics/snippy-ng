from pathlib import Path
from snippy_ng.stages.trees import IQTreeBuildTree

def create_tree_pipeline_stages(
    aln: str,
    model: str = "GTR+G",
    bootstrap: int = 1000,
    fconst: str | None = None,
    tmpdir: Path = Path("/tmp"),
    cpus: int = 1,
    ram: int = 8,
) -> list:
    stages = []

    # Stage to build a phylogenetic tree using IQ-TREE
    iqtree_stage = IQTreeBuildTree(
        aln=aln,
        model=model,
        bootstrap=bootstrap,
        fconst=fconst,
        tmpdir=tmpdir,
        cpus=cpus,
        ram=ram,
    )
        
    stages.append(iqtree_stage)

    return stages