from snippy_ng.stages.trees import IQTreeBuildTree
from pathlib import Path

def create_tree_pipeline_stages(
    aln: str,
    model: str,
    bootstrap: int,
    fconst: str | None,
    tmpdir: str | None,
    cpus: int,
    ram: int,
) -> list:
    stages = []

    #if fconst is a path read the content
    if fconst and Path(fconst).is_file():
        with open(fconst, 'r') as f:
            fconst = f.read().strip()

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