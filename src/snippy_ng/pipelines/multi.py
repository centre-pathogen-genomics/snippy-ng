import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Optional, TextIO, Tuple
from collections import OrderedDict
from snippy_ng.context import Context
from snippy_ng.exceptions import PipelineExecutionError, SnippyError
from snippy_ng.logging import logger, derive_log_path
import csv
import json
import io

MultiPipelineResult = tuple[list[str], list[tuple[str, str]]]


def run_multi_pipeline(
    snippy_reference_dir: Path,
    samples: Dict[str, Any],
    *,
    prefix: str,
    run_ctx: Context,
    cpus_per_sample: int,
    stop_on_failure: bool = False,
) -> MultiPipelineResult:
    """
    Special pipeline runner for multi-sample mode. Runs each sample in parallel using ProcessPoolExecutor.
    """

    total_cpus = int(run_ctx.cpus)
    # cap cpus_per_sample to total_cpus
    cpus_per_sample = min(int(cpus_per_sample), total_cpus)
    # Limit max parallelism to avoid oversubscription
    max_parallel = max(1, total_cpus // cpus_per_sample)


    # Minimal picklable config
    global_config = {
        "reference": str(snippy_reference_dir),
        "outdir": str(run_ctx.outdir),
        "prefix": prefix,
        "run_ctx": run_ctx.model_dump(mode="python"),
        "cpus_per_sample": cpus_per_sample,
    }

    jobs = [
        (sample_name, sample_cfg, global_config)
        for sample_name, sample_cfg in samples.items()
    ]

    ctx = mp.get_context("spawn")  # safest cross-platform
    failures = []
    successful_samples = []

    executor = ProcessPoolExecutor(max_workers=max_parallel, mp_context=ctx)
    fail_fast_shutdown = False
    try:
        futures = {executor.submit(_run_one_sample, job): job[0] for job in jobs}
        pending = set(futures)

        while pending:
            for fut in as_completed(pending):
                pending.remove(fut)
                sample = futures[fut]
                try:
                    fut.result()
                    successful_samples.append(sample)
                except SnippyError as e:
                    failures.append((sample, str(e)))
                    logger.error(f"FAILED: {sample}: {e}")
                except Exception as e:
                    # Catch any unexpected exception types so that multi-sample
                    # runs always complete with a consolidated failure summary.
                    failures.append((sample, str(e)))
                    logger.error(f"FAILED with unexpected error: {sample}: {e}")

                if failures and stop_on_failure:
                    for pending_fut in pending:
                        pending_fut.cancel()
                    executor.shutdown(wait=False, cancel_futures=True)
                    fail_fast_shutdown = True
                    pending.clear()
                    break
    finally:
        if not fail_fast_shutdown:
            executor.shutdown(wait=True)

    if failures:
        failed_samples_tsv = Path(run_ctx.outdir) / f"{prefix}.failed.tsv"
        _write_failed_samples_tsv(failures=failures, output_path=failed_samples_tsv)
        logger.warning(f"Skipping {len(failures)} failed samples; wrote {failed_samples_tsv}")

    if not successful_samples:
        raise PipelineExecutionError("No samples completed successfully.")

    _combine_sample_stats(
        samples_dir=Path(run_ctx.outdir) / "samples",
        sample_names=successful_samples,
        prefix=prefix,
        outdir=Path(run_ctx.outdir),
    )
    return successful_samples, failures
    

def _run_one_sample(job: Tuple[str, Dict[str, Any], Dict[str, Any]]) -> str:
    sample_name, sample_cfg, config = job

    # Import inside worker for clean spawn
    from snippy_ng.pipelines.asm import AsmPipelineBuilder
    from snippy_ng.pipelines.long import LongPipelineBuilder
    from snippy_ng.pipelines.short import ShortPipelineBuilder
    import click
    
    sample_type = sample_cfg.get("type")
    del sample_cfg["type"]

    if sample_type == "short":
        reads = sample_cfg.get("reads")
        if reads and isinstance(reads, str):
            reads = [reads]
        if not reads:
            # if reads not provided, expect left/right
            reads = [str(r) for r in (sample_cfg.get("left"), sample_cfg.get("right")) if r]
        pipeline = ShortPipelineBuilder(
            reference=config["reference"],
            reads=reads,
            bam=str(sample_cfg.get("bam")) if sample_cfg.get("bam") else None,
            prefix=config["prefix"],
            sample_name=sample_name,
            **{k: v for k, v in sample_cfg.items() if k not in ["left", "right", "bam", "reads", "sample_name"]},
        ).build()

    elif sample_type == "long":
        pipeline = LongPipelineBuilder(
            reference=config["reference"],
            reads=str(sample_cfg.get("reads")) if sample_cfg.get("reads") else None,
            bam=str(sample_cfg.get("bam")) if sample_cfg.get("bam") else None,
            prefix=config["prefix"],
            sample_name=sample_name,
            **{k: v for k, v in sample_cfg.items() if k not in ["reads", "bam", "sample_name"]},
        ).build()

    elif sample_type == "asm":
        sample_cfg['assembly'] = str(sample_cfg.get("assembly"))
        pipeline = AsmPipelineBuilder(
            reference=config["reference"],
            prefix=config["prefix"],
            sample_name=sample_name,
            **{k: v for k, v in sample_cfg.items() if k != "sample_name"},
        ).build()

    else:
        raise click.UsageError(
            f"Unknown sample type '{sample_type}' for sample '{sample_name}'"
        )

    # Per-sample output directory to prevent collisions
    outdir = Path(config["outdir"]) / 'samples' / sample_name
    
    # run_snippy_pipeline sets the working dir to outdir
    parent_run_ctx = Context(**config["run_ctx"])
    run_ctx = parent_run_ctx.model_copy(
        update={
            "outdir": outdir,
            "cpus": config["cpus_per_sample"],
            "log_path": derive_log_path(parent_run_ctx.log_path, outdir),
        }
    )
    pipeline.run(run_ctx)

    return sample_name


def _combine_sample_stats(
    samples_dir: Path,
    sample_names: list[str],
    prefix: str,
    outdir: Path,
) -> None:
    tsv_outputs = OrderedDict({
        f"{prefix}.reads.tsv": outdir / f"{prefix}.reads.tsv",
        f"{prefix}.vcf.summary.tsv": outdir / f"{prefix}.vcf.summary.tsv",
        f"{prefix}.vcf.breakdown.tsv": outdir / f"{prefix}.vcf.breakdown.tsv",
    })

    for filename, combined_path in tsv_outputs.items():
        sample_files = [samples_dir / sample_name / filename for sample_name in sample_names]
        _concat_tsv_files(sample_files=sample_files, output_path=combined_path)


def _concat_tsv_files(sample_files: list[Path], output_path: Path) -> None:
    header: list[str] | None = None
    rows_written = 0

    with output_path.open("w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")

        for sample_file in sample_files:
            if not sample_file.exists():
                continue

            with sample_file.open("r", newline="") as in_handle:
                reader = csv.reader(in_handle, delimiter="\t")
                file_header = next(reader, None)
                if file_header is None:
                    continue

                if header is None:
                    header = file_header
                    writer.writerow(header)
                elif file_header != header:
                    raise PipelineExecutionError(
                        f"Cannot combine TSVs with different headers: {sample_file} has {file_header}, expected {header}"
                    )

                for row in reader:
                    writer.writerow(row)
                    rows_written += 1

    if header is None:
        output_path.unlink(missing_ok=True)
    elif rows_written == 0:
        logger.warning(f"Combined TSV has no data rows: {output_path}")


def _write_failed_samples_tsv(
    failures: list[tuple[str, str]],
    output_path: Path,
) -> None:
    with output_path.open("w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["sample", "error"])
        for sample, error in failures:
            writer.writerow([sample, error])


def load_multi_config(config_file: TextIO, reference: Optional[Path]) -> Dict[str, Any]:
    
    # Read all content upfront to avoid seeking issues with stdin
    content = config_file.read()
    
    # Guess format from filename or content
    name = config_file.name
    if name.endswith(".csv"):
        config_type = "csv"
    elif name.endswith(".tsv"):
        config_type = "tsv"
    elif name.endswith(".json"):
        config_type = "json"
    else:
        # Guess from content
        first_line = content.split('\n', 1)[0] if content else ""
        if "," in first_line:
            config_type = "csv"
        elif "\t" in first_line:
            config_type = "tsv"
        elif first_line.strip().startswith("{"):
            config_type = "json"
        else:
            raise ValueError("Could not guess config format from file extension or content. Please use .csv, .tsv, or .json extension.")
    
    # Parse based on type
    if config_type in ("csv", "tsv"):
        if not reference:
            raise ValueError("Reference must be provided when using CSV/TSV config")
        reader = csv.DictReader(io.StringIO(content), delimiter="\t" if config_type == "tsv" else ",")
        cfg = {"samples": {}}
        for row in reader:
            sample = row.get("sample")
            if not sample:
                raise ValueError("Each row in the config file must have a 'sample' column with a unique sample name")
            if sample in cfg["samples"]:
                raise ValueError(f"Duplicate sample name '{sample}' found in config file")
            row.pop("sample", None)
            cfg["samples"][sample] = {k: v for k, v in row.items() if v}
    elif config_type == "json":
        cfg = json.loads(content)
        if not reference and "reference" not in cfg:
            raise ValueError("Reference must be provided in config file or as a command-line option")
    else:
        raise ValueError("Unsupported config file format. Use CSV, TSV, or JSON.")

    for name, row in cfg.get("samples", {}).items():
        if "type" not in row:
            raise ValueError(f"Sample '{name}' is missing required 'type' field")

    return cfg
