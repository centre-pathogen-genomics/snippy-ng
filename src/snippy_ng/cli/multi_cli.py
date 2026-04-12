from pathlib import Path
from typing import Any
import click

from snippy_ng.cli.utils import AbsolutePath
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options, GlobalOption


@click.command(cls=CommandWithGlobals, context_settings={"show_default": True})
@click.option("--stop-on-failure", is_flag=True, default=False, help="Stop the run when any per-sample analysis fails", cls=GlobalOption)
@click.option("--cpus-per-sample", type=click.INT, default=8, help="Number of CPUs to allocate per sample", cls=GlobalOption)
@add_snippy_global_options()
@click.argument(
    "config",
    required=True,
    type=click.File(mode="r")
)
@click.option(
    "--reference",
    "--ref",
    required=False,
    type=AbsolutePath(exists=True, readable=True),
)
@click.option("--core", type=click.FloatRange(min=0, max=1.0), default=0.95, help="Proportion of samples a site must be present in to be included in the core alignment")
@click.option("--inclusion-threshold", "-i",  type=click.FloatRange(min=0, max=1.0), default=0.1, help="Posterior probability threshold for retaining membership in the main alignment percentage cluster")
def multi(config: click.File, reference: Path | None, cpus_per_sample: int, core: float, inclusion_threshold: float, stop_on_failure: bool, outdir: Path, prefix: str, **context: Any):
    """
    Multi-sample SNP calling pipeline and core alignment construction 

    Example usage:

        $ snippy-ng multi samples.csv --ref reference.fasta

        $ snippy-ng utils gather --ref reference.fasta --json | snippy-ng multi -

    """
    from snippy_ng.pipelines.common import load_or_prepare_reference
    from snippy_ng.pipelines import SnippyPipeline
    from snippy_ng.exceptions import PipelineExecutionError
    from snippy_ng.pipelines.multi import load_multi_config, run_multi_pipeline
    from snippy_ng.logging import derive_log_path
    from snippy_ng.context import Context

    
    try:
        cfg = load_multi_config(config, reference)
    except Exception as e:
        raise click.ClickException(e)
    if reference is not None:
        cfg["reference"] = str(reference)
    if "reference" not in cfg or not cfg["reference"]:
        raise click.ClickException("Reference genome must be specified either in the config file or as a command line option")
    # create reusable reference
    ref_stage = load_or_prepare_reference(
        reference_path=cfg["reference"],
        output_directory=outdir / "reference",
    )
    ref_pipeline = SnippyPipeline(stages=[ref_stage])
    root_log_path = context.get("log_path") or Context.model_fields["log_path"].default
    context["log_path"] = derive_log_path(root_log_path, outdir / "reference")
    context["outdir"] = outdir / 'reference'
    run_ctx = Context(**context)
    ref_pipeline.run(run_ctx)
    run_ctx.outdir = outdir
    run_ctx.log_path = derive_log_path(run_ctx.log_path, outdir)

    snippy_reference_dir = ref_stage.output.reference_directory
    successful_samples, failures = run_multi_pipeline(
        snippy_reference_dir=snippy_reference_dir,
        samples=cfg["samples"],
        prefix=prefix,
        run_ctx=run_ctx,
        cpus_per_sample=cpus_per_sample,
        stop_on_failure=stop_on_failure,
    )
    if failures and stop_on_failure:
        raise PipelineExecutionError(
            "Some samples failed:\n"
            + "\n".join(f"Sample '{s}' -> {msg}" for s, msg in failures)
        )
    
    # core alignment
    from snippy_ng.pipelines.core import CorePipelineBuilder

    aln_pipeline = CorePipelineBuilder(
        snippy_dirs=[str(outdir / 'samples' / sample) for sample in successful_samples],
        reference=snippy_reference_dir,
        core=core,
        inclusion_threshold=inclusion_threshold,
    ).build()
    core_outdir = Path(outdir) / 'core'
    context["log_path"] = derive_log_path(run_ctx.log_path, core_outdir)
    context["outdir"] = core_outdir
    core_run_ctx = Context(**context)
    result = aln_pipeline.run(core_run_ctx)
    if failures:
        raise PipelineExecutionError(
            "Some samples failed and were removed:\n"
            + "\n".join(f"Sample '{s}' -> {msg}" for s, msg in failures)
        )
    return result
