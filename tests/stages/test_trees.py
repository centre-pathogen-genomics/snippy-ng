from Bio import Phylo

from snippy_ng.pipelines.tree import TreePipelineBuilder
from snippy_ng.stages.trees import IQTreeBuildTree, ScaleTreeToSNPs


def test_scale_tree_to_snps_multiplies_branch_lengths_by_alignment_sites(tmp_path):
    tree = tmp_path / "tree.treefile"
    tree.write_text("(sample_a:0.1,(sample_b:0.2,sample_c:0.3):0.4);")
    aln = tmp_path / "core.aln"
    aln.write_text(
        ">sample_a\nACGTACGTAA\n"
        ">sample_b\nACGTACGTAT\n"
        ">sample_c\nACGTACGTAG\n"
    )
    output_tree = tmp_path / "tree.snps.treefile"

    ScaleTreeToSNPs.scale_tree_to_snps(tree, aln, output_tree)

    scaled_tree = Phylo.read(str(output_tree), "newick")
    branch_lengths = sorted(
        clade.branch_length
        for clade in scaled_tree.find_clades()
        if clade.branch_length not in (None, 0.0)
    )
    assert branch_lengths == [1.0, 2.0, 3.0, 4.0]


def test_scale_tree_to_snps_includes_fconst_sites(tmp_path):
    tree = tmp_path / "tree.treefile"
    tree.write_text("(sample_a:0.1,sample_b:0.2,sample_c:0.3);")
    aln = tmp_path / "core.aln"
    aln.write_text(
        ">sample_a\nACGT\n"
        ">sample_b\nACGA\n"
        ">sample_c\nACGG\n"
    )
    output_tree = tmp_path / "tree.snps.treefile"

    ScaleTreeToSNPs.scale_tree_to_snps(tree, aln, output_tree, "1,2,3,0")

    scaled_tree = Phylo.read(str(output_tree), "newick")
    branch_lengths = sorted(
        clade.branch_length
        for clade in scaled_tree.find_clades()
        if clade.branch_length not in (None, 0.0)
    )
    assert branch_lengths == [1.0, 2.0, 3.0]


def test_tree_pipeline_adds_snp_scaled_tree_after_iqtree(tmp_path):
    aln = tmp_path / "core.aln"
    aln.write_text(
        ">sample_a\nACGT\n"
        ">sample_b\nACGA\n"
        ">sample_c\nACGG\n"
    )

    pipeline = TreePipelineBuilder(aln=aln, prefix="tree").build()

    assert isinstance(pipeline.stages[0], IQTreeBuildTree)
    assert isinstance(pipeline.stages[1], ScaleTreeToSNPs)
    assert pipeline.stages[1].tree == pipeline.stages[0].output.tree
    assert pipeline.stages[1].aln == aln
    assert pipeline.stages[1].fconst is None
    assert pipeline.stages[1].output.tree.name == "tree.snps.newick"


def test_tree_pipeline_passes_fconst_string(tmp_path):
    aln = tmp_path / "core.aln"
    aln.write_text(
        ">sample_a\nACGT\n"
        ">sample_b\nACGA\n"
        ">sample_c\nACGG\n"
    )

    pipeline = TreePipelineBuilder(aln=aln, fconst="1,2,3,4", prefix="tree").build()

    assert pipeline.stages[0].fconst == "1,2,3,4"
    assert pipeline.stages[1].fconst == "1,2,3,4"
