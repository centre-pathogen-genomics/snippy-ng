from pathlib import Path

from Bio import Phylo

from snippy_ng.context import Context
from snippy_ng.pipelines.tree import TreePipelineBuilder
from snippy_ng.stages.trees import ClonalFrameMLCorrectTree, IQTreeBuildTree, ScaleTreeToSNPs


def test_iqtree_outputs_append_extensions_to_dotted_prefix(tmp_path):
    stage = IQTreeBuildTree(
        aln=tmp_path / "core.aln",
        prefix="snippy.initial",
        cmaple=False,
    )

    assert stage.output.bionj == Path("snippy.initial.bionj")
    assert stage.output.checkpoint == Path("snippy.initial.ckp.gz")
    assert stage.output.contree == Path("snippy.initial.contree")
    assert stage.output.iqtree == Path("snippy.initial.iqtree")
    assert stage.output.log == Path("snippy.initial.log")
    assert stage.output.mldist == Path("snippy.initial.mldist")
    assert stage.output.splits == Path("snippy.initial.splits.nex")
    assert stage.output.tree == Path("snippy.initial.treefile")


def test_iqtree_pathogen_mode_does_not_require_standard_mode_outputs(tmp_path):
    stage = IQTreeBuildTree(
        aln=tmp_path / "core.aln",
        prefix="snippy.initial",
        cmaple=True,
        fast_mode=True,
    )

    assert stage.output.bionj is None
    assert stage.output.checkpoint is None
    assert stage.output.iqtree is None
    assert stage.output.mldist is None
    assert stage.output.paths == [
        Path("snippy.initial.log"),
        Path("snippy.initial.treefile"),
    ]


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


def test_clonalframeml_correct_tree_command(tmp_path):
    stage = ClonalFrameMLCorrectTree(
        tree=tmp_path / "initial.treefile",
        aln=tmp_path / "core.full.aln",
        kappa=3.2,
        emsim=10,
        prefix="tree.clonalframe",
    )

    commands = stage.create_commands(Context(cpus=4))

    assert len(commands) == 1
    assert commands[0].command == [
        "ClonalFrameML",
        str(tmp_path / "initial.treefile"),
        str(tmp_path / "core.full.aln"),
        "tree.clonalframe",
        "-kappa",
        "3.2",
        "-num_threads",
        "4",
        "-ignore_incomplete_sites",
        "false",
        "-show_progress",
        "false",
        "-emsim",
        "10",
    ]
    assert commands[0].output_file is None
    assert stage.output.labelled_tree == Path("tree.clonalframe.labelled_tree.newick")
    assert stage.output.emsim == Path("tree.clonalframe.emsim.txt")


def test_tree_pipeline_uses_clonalframe_corrected_tree(tmp_path):
    aln = tmp_path / "core.full.aln"
    aln.write_text(
        ">sample_a\nACGT\n"
        ">sample_b\nACGA\n"
        ">sample_c\nACGG\n"
    )

    pipeline = TreePipelineBuilder(
        aln=aln,
        fconst="1,2,3,4",
        clonalframe=True,
        clonalframe_kappa=3.0,
        clonalframe_emsim=5,
        prefix="tree",
    ).build()

    assert [type(stage) for stage in pipeline.stages] == [
        IQTreeBuildTree,
        ClonalFrameMLCorrectTree,
        ScaleTreeToSNPs,
    ]
    initial_tree, clonalframe, scaled_tree = pipeline.stages
    assert initial_tree.prefix == "tree.initial"
    assert initial_tree.fconst == "1,2,3,4"
    assert clonalframe.tree == initial_tree.output.tree
    assert clonalframe.aln == aln
    assert clonalframe.kappa == 3.0
    assert clonalframe.emsim == 5
    assert scaled_tree.tree == clonalframe.output.labelled_tree
    assert scaled_tree.aln == aln
    assert scaled_tree.fconst == "1,2,3,4"
    assert clonalframe.output.labelled_tree in pipeline.outputs_to_keep
    assert clonalframe.output.importation_status in pipeline.outputs_to_keep
    assert clonalframe.output.em in pipeline.outputs_to_keep
    assert clonalframe.output.emsim in pipeline.outputs_to_keep
    assert clonalframe.output.ml_sequence in pipeline.outputs_to_keep
    assert clonalframe.output.position_cross_reference in pipeline.outputs_to_keep
