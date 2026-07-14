from pathlib import Path
from typing import Optional
from pydantic import Field
from snippy_ng.pipelines import PipelineBuilder, SnippyPipeline
from snippy_ng.stages.trees import ClonalFrameMLCorrectTree, IQTreeBuildTree, ScaleTreeToSNPs


class TreePipelineBuilder(PipelineBuilder):
    """Builder for phylogenetic tree building pipeline."""
    aln: Path = Field(..., description="Multiple sequence alignment file")
    model: str = Field(default="GTR+G", description="Substitution model for IQ-TREE")
    bootstrap: int = Field(default=1000, description="Number of bootstrap replicates")
    fconst: Optional[str] = Field(default=None, description="Frequency constants for ascertainment bias correction")
    fast_mode: bool = Field(default=False, description="Use fast mode for IQ-TREE (faster but less accurate)")
    cmaple: bool = Field(default=True, description="Use pathogen mode for IQ-TREE with --alrt for SH-aLRT support values")
    clonalframe: bool = Field(default=False, description="Correct the inferred tree for recombination with ClonalFrameML")
    clonalframe_kappa: float = Field(default=2.0, gt=0.0, description="Transition/transversion bias used by ClonalFrameML")
    clonalframe_emsim: int = Field(default=0, ge=0, description="ClonalFrameML simulations used to estimate uncertainty")

    def build(self) -> SnippyPipeline:
        """Build and return the tree building pipeline."""
        stages = []

        # Stage to build a phylogenetic tree using IQ-TREE
        iqtree_stage = IQTreeBuildTree(
            aln=self.aln,
            model=self.model,
            bootstrap=self.bootstrap,
            fconst=self.fconst,
            fast_mode=self.fast_mode,
            prefix=f"{self.prefix}.initial" if self.clonalframe else self.prefix,
            cmaple=self.cmaple,
        )
            
        stages.append(iqtree_stage)

        final_tree = iqtree_stage.output.tree
        clonalframe_stage = None
        if self.clonalframe:
            clonalframe_stage = ClonalFrameMLCorrectTree(
                tree=iqtree_stage.output.tree,
                aln=self.aln,
                kappa=self.clonalframe_kappa,
                emsim=self.clonalframe_emsim,
                prefix=f"{self.prefix}.clonalframe",
            )
            stages.append(clonalframe_stage)
            final_tree = clonalframe_stage.output.labelled_tree

        snp_tree_stage = ScaleTreeToSNPs(
            tree=final_tree,
            aln=self.aln,
            fconst=self.fconst,
            prefix=self.prefix,
        )
        stages.append(snp_tree_stage)

        outputs_to_keep = [snp_tree_stage.output.tree, iqtree_stage.output.tree]
        if clonalframe_stage is not None:
            outputs_to_keep.extend(clonalframe_stage.output.paths)
        return SnippyPipeline(stages=stages, outputs_to_keep=outputs_to_keep)
