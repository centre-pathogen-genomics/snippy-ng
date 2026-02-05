from snippy_ng.stages.trees import IQTreeBuildTree

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