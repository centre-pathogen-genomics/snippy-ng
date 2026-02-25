import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Optional, TextIO, Tuple
from snippy_ng.exceptions import PipelineExecutionError, SnippyError
from snippy_ng.logging import logger
import csv
import json
import io


def run_multi_pipeline(
    snippy_reference_dir: Path,
    samples: Dict[str, Any],
    *,
    outdir: Path,
    prefix: str,
    tmpdir: Path,
    cpus: int,
    cpus_per_sample: int,
    ram: int,
    skip_check: bool,
    check: bool,
    quiet: bool,
    create_missing: bool,
    keep_incomplete: bool,
) -> None:
    """
    Special pipeline runner for multi-sample mode. Runs each sample in parallel using ProcessPoolExecutor.
    """

    total_cpus = int(cpus)
    # cap cpus_per_sample to total_cpus
    cpus_per_sample = min(int(cpus_per_sample), total_cpus)
    # Limit max parallelism to avoid oversubscription
    max_parallel = max(1, total_cpus // cpus_per_sample)


    # Minimal picklable config
    global_config = {
        "reference": str(snippy_reference_dir),
        "outdir": str(outdir),
        "prefix": prefix,
        "tmpdir": tmpdir,
        "ram": ram,
        "cpus_per_sample": cpus_per_sample,
        "skip_check": skip_check,
        "check": check,
        "quiet": quiet,
        "create_missing": create_missing,
        "keep_incomplete": keep_incomplete,
    }

    jobs = [
        (sample_name, sample_cfg, global_config)
        for sample_name, sample_cfg in samples.items()
    ]

    ctx = mp.get_context("spawn")  # safest cross-platform
    failures = []

    with ProcessPoolExecutor(max_workers=max_parallel, mp_context=ctx) as executor:
        futures = {executor.submit(_run_one_sample, job): job[0] for job in jobs}

        for fut in as_completed(futures):
            sample = futures[fut]
            try:
                fut.result()
            except SnippyError as e:
                failures.append((sample, str(e)))
                logger.error(f"FAILED: {sample}: {e}")
            except Exception as e:
                # Catch any unexpected exception types so that multi-sample
                # runs always complete with a consolidated failure summary.
                failures.append((sample, str(e)))
                logger.error(f"FAILED with unexpected error: {sample}: {e}")

    if failures:
        raise PipelineExecutionError("Some samples failed (check logs for details):\n" + "\n".join(f"Sample '{s}' -> {msg}" for s, msg in failures))
    

def _run_one_sample(job: Tuple[str, Dict[str, Any], Dict[str, Any]]) -> str:
    sample_name, sample_cfg, config = job

    # Import inside worker for clean spawn
    from snippy_ng.pipelines.asm import AsmPipelineBuilder
    from snippy_ng.pipelines.long import LongPipelineBuilder
    from snippy_ng.pipelines.short import ShortPipelineBuilder
    from snippy_ng.context import Context
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
            **{k: v for k, v in sample_cfg.items() if k not in ["left", "right", "bam", "reads"]},
        ).build()

    elif sample_type == "long":
        pipeline = LongPipelineBuilder(
            reference=config["reference"],
            reads=str(sample_cfg.get("reads")) if sample_cfg.get("reads") else None,
            bam=str(sample_cfg.get("bam")) if sample_cfg.get("bam") else None,
            prefix=config["prefix"],
            **{k: v for k, v in sample_cfg.items() if k not in ["reads", "bam"]},
        ).build()

    elif sample_type == "asm":
        sample_cfg['assembly'] = str(sample_cfg.get("assembly"))
        pipeline = AsmPipelineBuilder(
            reference=config["reference"],
            prefix=config["prefix"],
            **sample_cfg,
        ).build()

    else:
        raise click.UsageError(
            f"Unknown sample type '{sample_type}' for sample '{sample_name}'"
        )

    # Per-sample output directory to prevent collisions
    outdir = Path(config["outdir"]) / 'samples' / sample_name
    outdir.mkdir(parents=True, exist_ok=True)
    
    # run_snippy_pipeline sets the working dir to outdir
    run_ctx = Context(
        outdir=outdir,
        tmpdir=config["tmpdir"],
        cpus=config["cpus_per_sample"],
        ram=config["ram"],
        quiet=config["quiet"],
        create_missing=config["create_missing"],
        keep_incomplete=config["keep_incomplete"],
        skip_check=config["skip_check"],
        check=config["check"],
    )
    pipeline.run(run_ctx)

    return sample_name


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