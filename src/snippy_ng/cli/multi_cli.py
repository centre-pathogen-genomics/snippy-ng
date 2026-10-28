from pathlib import Path
from typing import Any, Optional, TextIO
import click

from snippy_ng.cli.utils import AbsolutePath, config_or_url_callback, reference_or_accession_callback
from snippy_ng.cli.utils.globals import CommandWithGlobals, add_snippy_global_options, GlobalOption


def run_multi_config(
    cfg: dict,
    *,
    reference: Optional[Path],
    cpus_per_sample: Optional[int] = None,
    core: float,
    remove_outliers: bool,
    outlier_threshold: float,
    stop_on_failure: bool,
    outdir: Path,
    prefix: str,
    context: dict[str, Any],
    mask: Optional[Path] = None,
):
    from snippy_ng.pipelines.common import (
        download_reference,
        is_assembly_accession,
        load_or_prepare_reference,
    )
    from snippy_ng.pipelines import SnippyPipeline
    from snippy_ng.exceptions import PipelineExecutionError
    from snippy_ng.pipelines.multi import add_core_alignment_qc, run_multi_pipeline
    from snippy_ng.logging import derive_log_path
    from snippy_ng.context import Context

    if reference is not None:
        cfg["reference"] = str(reference)
    if "reference" not in cfg or not cfg["reference"]:
        raise click.ClickException("Reference genome must be specified either in the config file or as a command line option")
    # create reusable reference
    reference_stages = []
    reference_input = cfg["reference"]
    if is_assembly_accession(reference_input):
        reference_input = download_reference(
            reference_input,
            reference_stages,
            output_directory=outdir / "reference",
        )
    ref_stage = load_or_prepare_reference(
        reference_path=reference_input,
        output_directory=outdir / "reference",
    )
    reference_stages.append(ref_stage)
    ref_pipeline = SnippyPipeline(stages=reference_stages)
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
    if not successful_samples:
        raise PipelineExecutionError(
            "All samples failed. No core alignment will be produced.\n"
            + "\n".join(f"Sample '{s}' -> {msg}" for s, msg in failures)
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
        remove_outliers=remove_outliers,
        outlier_threshold=outlier_threshold,
        mask=mask,
    ).build()
    core_outdir = Path(outdir) / 'core'
    context["log_path"] = derive_log_path(run_ctx.log_path, core_outdir)
    context["outdir"] = core_outdir
    core_run_ctx = Context(**context)
    result = aln_pipeline.run(core_run_ctx)
    add_core_alignment_qc(
        qc_tsv=Path(outdir) / f"{prefix}.qc.tsv",
        core_aligned_tsv=core_outdir / "core.aligned.tsv",
    )
    return {
        "result": result,
        "successful_samples": successful_samples,
        "failed_samples": failures,
        "core_outdir": core_outdir,
    }


@click.command(cls=CommandWithGlobals, context_settings={"show_default": True})
@click.option("--stop-on-failure", is_flag=True, default=False, help="Stop the run when any per-sample analysis fails", cls=GlobalOption)
@click.option("--cpus-per-sample", type=click.INT, default=None, help="Number of CPUs to allocate per sample", cls=GlobalOption)
@add_snippy_global_options(exclude=["sample_name"])
@click.argument(
    "config",
    required=True,
    type=click.STRING,
    callback=config_or_url_callback,
    metavar="CONFIG",
)
@click.option(
    "--reference",
    "--ref",
    required=False,
    type=click.STRING,
    callback=reference_or_accession_callback,
)
@click.option("--core", type=click.FloatRange(min=0, max=1.0), default=0.95, help="Proportion of samples a site must be present in to be included in the core alignment")
@click.option("--remove-alignment-outliers", "remove_outliers", is_flag=True, default=False, help="Remove detected low-alignment outlier samples from the alignment before core filtering")
@click.option("--alignment-outlier-threshold", "outlier_threshold", type=click.FloatRange(min=0.0, max=1.0, min_open=True), default=0.20, help="Posterior probability threshold for retaining membership in the retained alignment percentage clusters. Only used with --remove-alignment-outliers and cohorts large enough for mixture modelling")
@click.option(
    "--mask",
    default=None,
    type=AbsolutePath(exists=True, readable=True),
    help="Mask file (BED format, reference coordinates) to mask regions with Ns in every sample of the alignment",
)
def multi(config: TextIO, reference: Optional[Path | str], cpus_per_sample: Optional[int], core: float, remove_outliers: bool, outlier_threshold: float, mask: Optional[Path], stop_on_failure: bool, outdir: Path, prefix: str, **context: Any):
    """
    Multi-sample SNP calling pipeline and core alignment construction

    CONFIG is a CSV, TSV, or JSON file describing the samples. It may be a local
    path, '-' for stdin, or an HTTP(S) URL to download the config from.

    Example usage:

        $ snippy-ng multi samples.csv --ref reference.fasta

        $ snippy-ng multi https://example.com/samples.json --ref reference.fasta

        $ snippy-ng utils gather --ref reference.fasta --json | snippy-ng multi -

    """
    from snippy_ng.pipelines.multi import load_multi_config
    from snippy_ng.exceptions import PipelineExecutionError

    try:
        cfg = load_multi_config(config, reference)
    except Exception as e:
        raise click.ClickException(e)
    result = run_multi_config(
        cfg,
        reference=reference,
        cpus_per_sample=cpus_per_sample,
        core=core,
        remove_outliers=remove_outliers,
        outlier_threshold=outlier_threshold,
        mask=mask,
        stop_on_failure=stop_on_failure,
        outdir=outdir,
        prefix=prefix,
        context=context,
    )
    if result["failed_samples"]:
        raise PipelineExecutionError(
            "Some samples failed and were removed:\n"
            + "\n".join(f"Sample '{s}' -> {msg}" for s, msg in result["failed_samples"])
        )
    return 
