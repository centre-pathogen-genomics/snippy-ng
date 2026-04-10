"""
Thanks for looking at the source code! You found a hidden command :D
"""

import os

from pathlib import Path
from typing import Any, Iterable, List
import click

from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import CommandWithGlobals, GlobalOption, add_snippy_global_options


@click.command(
    cls=CommandWithGlobals, context_settings={"show_default": True}, hidden=True
)
@click.option("--cpus", "-c", default=None, required=False, type=int, help="Maximum number of CPUs to use. By default, uses all available CPUs.", cls=GlobalOption)
@add_snippy_global_options(['cpus'])
@click.option("--reference", "--ref", required=False, type=AbsolutePath(exists=True, readable=True), help="Reference genome (FASTA or GenBank) or prepared reference directory")
@click.argument(
    "directory",
    required=False,
    type=AbsolutePath(exists=True, readable=True),
    nargs=-1
)
def yolo(directory: Iterable[Path], reference: Path, outdir: Path, prefix: str, **context: Any):
    """
    Pipeline that automates everything.

    Not recommended for general use unless you've got no idea what you're doing.

    Example usage:

        snippy-ng yolo
    """
    import json
    from collections import Counter
    from snippy_ng.context import Context
    from snippy_ng.logging import logger, derive_log_path
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
    if reference is None:
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
            "No reference file found! Please provide `--reference` or ensure you have a file called `reference` with one of the following extensions: fasta, fa, fna, gbk, genbank e.g. reference.fasta or ref.gbk"
        )

    # find all samples in the directory and create config
    logger.info(f"Gathering samples from: {', '.join(str(d) for d in directories)}")
    gathered = gather_samples_config(
        inputs=directories,
        max_depth=4,
        aggressive_ids=False,
        exclude_name_regex=None,
        reference=reference,
    )
    samples = gathered["samples"]
    if not samples:
        logger.error("No samples found. Please ensure you have at least one sample with reads in the input directory.")
        return 1
    type_counts = Counter(sample_data.get("type", "unknown") for sample_data in samples.values())
    type_summary = ", ".join(f"{sample_type}={count}" for sample_type, count in sorted(type_counts.items()))
    logger.info(f"Found {len(samples)} samples ({type_summary})")
    
    # write config to output directory
    outdir.mkdir(parents=True, exist_ok=True)
    with open(Path(outdir) / "samples.json", "w") as f:
        f.write(json.dumps(gathered, indent=2))

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
        output_directory=outdir / "reference",
    )
    ref_pipeline = SnippyPipeline(stages=[ref_stage])
    if context["cpus"] is None:
        context["cpus"] = os.cpu_count() or 1 # YOLO: use all available CPUs by default
    root_log_path = context.get("log_path") or Context.model_fields["log_path"].default
    context["log_path"] = derive_log_path(root_log_path, outdir / "reference")
    context["outdir"] = outdir / 'reference'
    run_ctx = Context(**context)
    ref_pipeline.run(run_ctx)
    run_ctx.outdir = outdir
    run_ctx.log_path = derive_log_path(run_ctx.log_path, outdir)

    snippy_reference_dir = ref_stage.output.reference_directory

    # each sample gets 4 CPUs or total_cpus / num_samples, whichever is higher
    cpus_per_sample = max(4, context["cpus"] // len(samples))
    try:
        successful_samples, failures = run_multi_pipeline(
            snippy_reference_dir=snippy_reference_dir,
            samples=cfg["samples"],
            prefix=prefix,
            run_ctx=run_ctx,
            cpus_per_sample=cpus_per_sample,
            stop_on_failure=False,
        )
    except PipelineExecutionError as e:
        logger.horizontal_rule(style="-")
        raise e

    if len(successful_samples) < 3:
        logger.warning("Less than 3 samples found, skipping tree construction.")
        if failures:
            raise PipelineExecutionError(
                "Some samples failed and were removed:\n"
                + "\n".join(f"Sample '{s}' -> {msg}" for s, msg in failures)
            )
        return 0

    # core alignment
    from snippy_ng.pipelines.core import CorePipelineBuilder

    snippy_dirs = [
        Path(outdir / "samples" / sample)
        for sample in successful_samples
    ]
    soft_core_threshold = 0.95
    inclusion_threshold = 0.3
    aln_pipeline = CorePipelineBuilder(
        snippy_dirs=snippy_dirs,
        reference=snippy_reference_dir,
        core=soft_core_threshold,
        inclusion_threshold=inclusion_threshold,
    ).build()
    core_outdir = Path(outdir) / "core"
    context["log_path"] = derive_log_path(run_ctx.log_path, core_outdir)
    context["outdir"] = core_outdir
    core_run_ctx = Context(**context)
    aln_pipeline.run(core_run_ctx)

    # tree
    from snippy_ng.pipelines.tree import TreePipelineBuilder


    context["ram"] = None # YOLO: disable RAM limiting for this step
    tree_pipeline = TreePipelineBuilder(
        aln=core_outdir / aln_pipeline.stages[-1].output.soft_core,
        fconst=(core_outdir / aln_pipeline.stages[-1].output.constant_sites).read_text().strip(),
        fast_mode=True,
    ).build()
    tree_outdir = Path(outdir) / "tree"
    context["log_path"] = derive_log_path(run_ctx.log_path, tree_outdir)
    context["outdir"] = tree_outdir
    tree_run_ctx = Context(**context)
    tree_pipeline.run(tree_run_ctx)

    # report
    from snippy_ng.pipelines.report import ReportPipelineBuilder

    report_pipeline = ReportPipelineBuilder(
        tree=tree_outdir / tree_pipeline.stages[-1].output.tree,
        title="Snippy-NG Report",
        metadata=Path(outdir) / f"{prefix}.vcf.summary.tsv",
        prefix="report",
    ).build()
    report_outdir = Path(outdir) / "report"
    context["log_path"] = derive_log_path(run_ctx.log_path, report_outdir)
    context["outdir"] = report_outdir
    report_run_ctx = Context(**context)
    result = report_pipeline.run(report_run_ctx)
    if failures:
        raise PipelineExecutionError(
            "Some samples failed and were removed:\n"
            + "\n".join(f"Sample '{s}' -> {msg}" for s, msg in failures)
        )
    return result
