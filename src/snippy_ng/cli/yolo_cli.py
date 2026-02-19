"""
Thanks for looking at the source code! You found a hidden command :D
"""

from pathlib import Path
from typing import Iterable, List
import click

from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options


@click.command(
    cls=CommandWithGlobals, context_settings={"show_default": True}, hidden=True
)
@add_snippy_global_options()
@click.argument(
    "directory",
    required=False,
    type=click.Path(exists=True, readable=True, path_type=Path),
    nargs=-1
)
def yolo(directory: Iterable[Path], **config):
    """
    Pipeline that automates everything.

    Not recommended for general use unless you've got no idea what you're doing.

    Example usage:

        snippy-ng yolo
    """
    from snippy_ng.logging import logger
    from snippy_ng.pipelines.common import load_or_prepare_reference
    from snippy_ng.pipelines.multi import run_multi_pipeline
    from snippy_ng.pipelines import SnippyPipeline
    from snippy_ng.utils.gather import gather_samples_config
    from snippy_ng.exceptions import PipelineExecutionError, InvalidReferenceError

    logger.warning(
        "You are running the YOLO pipeline. This pipeline is not recommended for general use unless you have no idea what you're doing. Please consider using one of the other pipelines with more specific parameters for better results and more control over the analysis."
    )

    # flatten, normalize, and de-duplicate input directories
    directories: List[Path] = []
    for d in directory or [Path.cwd()]:
        if isinstance(d, (list, tuple, set)):
            directories.extend(Path(item).absolute() for item in d if item)
        elif d:
            directories.append(Path(d).absolute())

    if not directories:
        directories = [Path.cwd()]

    # keep order while removing duplicates
    directories = list(dict.fromkeys(directories))

    # find reference and reads in the directories
    # look for file called reference or ref with fasta, fa, fna, gbk, genbank extension
    reference = None
    for search_dir in directories:
        for ext in ["fasta", "fa", "fna", "gbk", "genbank"]:
            for name in ["reference", "ref"]:
                candidates = list(search_dir.rglob(f"{name}.{ext}"))
                if candidates:
                    reference = candidates[0].absolute()
                    logger.info(f"Found reference file: {reference}")
                    break
            if reference:
                break
        if reference:
            break
    if not reference:
        raise InvalidReferenceError(
            "No reference file found! Please ensure you have a file called `reference` with one of the following extensions: fasta, fa, fna, gbk, genbank e.g. reference.fasta or ref.gbk"
        )

    # find all samples in the directory and create config
    gathered = gather_samples_config(
        inputs=directories,
        max_depth=4,
        aggressive_ids=False,
        exclude_name_regex=None,
        reference=reference,
    )
    samples = gathered["samples"]
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
        outdir=config["outdir"],
        quiet=config["quiet"],
        create_missing=config["create_missing"],
        keep_incomplete=config["keep_incomplete"],
    )

    snippy_reference_dir = ref_stage.output.reference.parent

    # each sample gets 2 CPUs or total_cpus / num_samples, whichever is higher
    config["cpus_per_sample"] = max(4, config["cpus"] // len(samples))
    try:
        run_multi_pipeline(
            snippy_reference_dir=snippy_reference_dir,
            samples=cfg["samples"],
            outdir=config["outdir"],
            prefix=config["prefix"],
            tmpdir=config["tmpdir"],
            cpus=config["cpus"],
            ram=config["ram"],
            skip_check=config["skip_check"],
            check=config["check"],
            quiet=config["quiet"],
            create_missing=config["create_missing"],
            keep_incomplete=config["keep_incomplete"],
            cpus_per_sample=config["cpus_per_sample"],
        )
    except PipelineExecutionError as e:
        logger.horizontal_rule(style="-")
        raise e

    # core alignment
    from snippy_ng.pipelines.core import CorePipelineBuilder

    snippy_dirs = [
        Path(config["outdir"] / "samples" / sample)
        for sample in cfg["samples"]
    ]
    soft_core_threshold = 0.95
    aln_pipeline = CorePipelineBuilder(
        snippy_dirs=snippy_dirs,
        reference=snippy_reference_dir,
        core=soft_core_threshold,
        tmpdir=config["tmpdir"],
        cpus=config["cpus"],
        ram=config["ram"],
    ).build()
    outdir = Path(config["outdir"]) / "core"
    outdir.mkdir(parents=True, exist_ok=True)
    aln_pipeline.run(
        skip_check=config["skip_check"],
        check=config["check"],
        outdir=outdir,
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
        aln=outdir / aln_pipeline.stages[-1].output.aln,
        fconst=(outdir / aln_pipeline.stages[-1].output.constant_sites).read_text().strip(),
        fast_mode=False,
        cpus=config["cpus"],
        ram=config["ram"],
    ).build()
    outdir = Path(config["outdir"]) / "tree"
    outdir.mkdir(parents=True, exist_ok=True)
    tree_pipeline.run(
        skip_check=config["skip_check"],
        check=config["check"],
        outdir=outdir,
        quiet=config["quiet"],
        create_missing=config["create_missing"],
        keep_incomplete=config["keep_incomplete"],
    )
