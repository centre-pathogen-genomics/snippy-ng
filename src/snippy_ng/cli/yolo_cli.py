"""
Thanks for looking at the source code! You found a hidden command :D
"""

from pathlib import Path
from typing import Optional
import click

from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options


@click.command(
    cls=CommandWithGlobals, context_settings={"show_default": True}, hidden=True
)
@add_snippy_global_options()
@click.argument(
    "directory",
    required=False,
    type=click.Path(exists=True, readable=True, path_type=AbsolutePath),
)
def yolo(directory: Optional[Path], **config):
    """
    Pipeline that automates everything.

    Not recommended for general use unless you've got no idea what you're doing.

    Example usage:

        snippy-ng yolo
    """
    import os
    from snippy_ng.logging import logger
    from snippy_ng.pipelines.common import load_or_prepare_reference
    from snippy_ng.pipelines.multi import run_multi_pipeline
    from snippy_ng.pipelines import SnippyPipeline
    from snippy_ng.utils.gather import gather_samples_config
    from snippy_ng.exceptions import PipelineExecutionError, InvalidReferenceError

    logger.warning(
        "You are running the YOLO pipeline. This pipeline is not recommended for general use unless you have no idea what you're doing. Please consider using one of the other pipelines with more specific parameters for better results and more control over the analysis."
    )

    # find reference and reads in the directory
    # look for file called reference or ref with fasta, fa, fna, gbk, genbank extension
    directory = directory or Path.cwd()
    reference = None
    for ext in ["fasta", "fa", "fna", "gbk", "genbank"]:
        for name in ["reference", "ref"]:
            candidates = list(directory.rglob(f"{name}.{ext}"))
            if candidates:
                reference = candidates[0].resolve()
                logger.info(f"Found reference file: {reference}")
                break
        if reference:
            break
    if not reference:
        raise InvalidReferenceError(
            "No reference file found! Please ensure you have a file called `reference` with one of the following extensions: fasta, fa, fna, gbk, genbank e.g. reference.fasta or ref.gbk"
        )

    # find all samples in the directory and create config
    samples = gather_samples_config(
        inputs=[directory],
        max_depth=4,
        aggressive_ids=False,
        exclude_name_regex=None,
        exclude_files=[reference],
    )
    logger.info(f"Found {len(samples)} samples: {', '.join(samples.keys())}")
    # use freebayes for long read samples in YOLO mode
    # TODO: need to determine the chemistry of the long reads to choose the best clair3 model
    fb_samples = {}
    for name, sample in samples.items():
        if sample["type"] == "long":
            sample["caller"] = "freebayes"
        fb_samples[name] = sample
    cfg = {
        "reference": str(reference),
        "samples": fb_samples,
    }
    # create reusable reference
    ref_stage = load_or_prepare_reference(
        reference_path=cfg["reference"],
        output_directory=Path(config["outdir"]) / "reference",
    )
    ref_pipeline = SnippyPipeline(stages=[ref_stage])
    ref_pipeline.run(
        skip_check=config["skip_check"],
        check=config["check"],
        cwd=config["outdir"],
        quiet=config["quiet"],
        create_missing=config["create_missing"],
        keep_incomplete=config["keep_incomplete"],
    )

    snippy_reference_dir = ref_stage.output.reference.parent

    config["cpus"] = (
        min(os.cpu_count(), config["cpus"]) if os.cpu_count() else config["cpus"]
    )
    config["cpus_per_sample"] = max(1, config["cpus"] // len(samples))
    try:
        run_multi_pipeline(
            snippy_reference_dir=snippy_reference_dir,
            samples=cfg["samples"],
            config=config,
        )
    except PipelineExecutionError as e:
        logger.horizontal_rule(style="-")
        raise e

    # core alignment
    from snippy_ng.pipelines.core import CorePipelineBuilder

    snippy_dirs = [
        str((Path(config["outdir"]) / "samples" / sample).resolve())
        for sample in cfg["samples"]
    ]
    aln_pipeline = CorePipelineBuilder(
        snippy_dirs=snippy_dirs,
        reference=snippy_reference_dir,
        core=0.95,
        tmpdir=config["tmpdir"],
        cpus=config["cpus"],
        ram=config["ram"],
    ).build()
    outdir = Path(config["outdir"]) / "core"
    outdir.mkdir(parents=True, exist_ok=True)
    aln_pipeline.run(
        skip_check=config["skip_check"],
        check=config["check"],
        cwd=outdir,
        quiet=config["quiet"],
        create_missing=config["create_missing"],
        keep_incomplete=config["keep_incomplete"],
    )

    if len(samples) < 3:
        logger.warning("Less than 3 samples found, skipping tree construction.")
        return 0
    # tree
    from snippy_ng.pipelines.tree import TreePipelineBuilder

    tree_pipeline = TreePipelineBuilder(
        aln=str(outdir / "core.aln"),
        fconst=(outdir / "core.aln.sites").read_text().strip(),
        cpus=config["cpus"],
        ram=config["ram"],
    ).build()
    outdir = Path(config["outdir"]) / "tree"
    outdir.mkdir(parents=True, exist_ok=True)
    tree_pipeline.run(
        skip_check=config["skip_check"],
        check=config["check"],
        cwd=outdir,
        quiet=config["quiet"],
        create_missing=config["create_missing"],
        keep_incomplete=config["keep_incomplete"],
    )
