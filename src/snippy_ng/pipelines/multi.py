import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Optional, TextIO, Tuple
from collections import OrderedDict
from snippy_ng.context import Context
from snippy_ng.exceptions import InvalidReferenceError, PipelineExecutionError, SnippyError
from snippy_ng.logging import logger, derive_log_path
from snippy_ng.pipelines.common import is_ncbi_assembly_accession, is_atb_assembly_accession, is_sra_accession
from snippy_ng.stages.download import DownloadSraReads, DownloadNcbiAssemblyFasta, DownloadAtbAssemblyFasta
from snippy_ng.dependencies import sracha
import csv
import io
import json
import subprocess
from snippy_ng.stages.stats import SAMPLE_QC_COLUMNS

MultiPipelineResult = tuple[list[str], list[tuple[str, str]]]


def _build_jobs(
    samples: dict[str, Any],
    global_config: dict[str, Any],
) -> list[tuple[str, dict[str, Any], dict[str, Any]]]:
    """Build job tuples for *ProcessPoolExecutor*."""
    return [
        (sample_name, sample_cfg, global_config)
        for sample_name, sample_cfg in samples.items()
    ]


def _predownload_multi_inputs(
    samples: dict[str, Any],
    *,
    cache_dir: Path,
    run_ctx: Context,
    failures: list[tuple[str, str]],
) -> dict[str, Any]:
    """Download all unique reads/assembly accessions into *cache_dir*.

    Returns a copy of *samples* with accession keys replaced by absolute local
    file paths.  Failed downloads are appended to *failures* and their samples
    are dropped from the returned dict so they never reach the executor.
    """
    cache_dir = cache_dir.resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)

    dl_ctx = run_ctx.model_copy(
        update={
            "outdir": cache_dir,
            "cpus": 1,
        }
    )

    # ---- 1. collect unique accessions per kind -------------------------------
    reads_acc_reqs: dict[str, list[str]] = {}
    asm_acc_reqs: dict[str, list[str]] = {}

    for sample_name, cfg in samples.items():
        if acc := cfg.get("reads_accession"):
            reads_acc_reqs.setdefault(str(acc), []).append(sample_name)
        if acc := cfg.get("assembly_accession"):
            asm_acc_reqs.setdefault(str(acc), []).append(sample_name)

    downloads: dict[str, list[Path]] = {}
    asm_downloads: dict[str, Path] = {}
    failed_by_sample: dict[str, list[str]] = {}
    checked_deps: set[str] = set()

    def _check_deps(stage) -> None:
        if run_ctx.skip_check:
            return
        for dep in stage._dependencies:
            if dep.citation_only or dep.name in checked_deps:
                continue
            dep.check()
            checked_deps.add(dep.name)

    # ---- 2. download reads (SRA) in a single batch ---------------------------
    if reads_acc_reqs:
        # Validate each accession; reject invalid ones before the batch
        valid_reads: dict[str, list[str]] = {}
        for acc, dependents in reads_acc_reqs.items():
            if not is_sra_accession(acc):
                msg = f"Unsupported SRA accession '{acc}'"
                for sname in dependents:
                    failed_by_sample.setdefault(sname, []).append(msg)
                logger.warning(msg)
            else:
                valid_reads[acc] = dependents

        if valid_reads:
            try:
                # Check sracha dependency once
                if not run_ctx.skip_check:
                    sracha.check()

                acc_list = cache_dir / "sra_accessions.txt"
                acc_list.write_text("\n".join(valid_reads) + "\n")

                subprocess.run(
                    ["sracha", "get", "--accession-list", str(acc_list),
                     "--output-dir", str(cache_dir)],
                    check=True, capture_output=True,
                )

                import shutil
                for acc in valid_reads:
                    single_end = cache_dir / f"{acc}.fastq.gz"
                    r1 = cache_dir / f"{acc}_1.fastq.gz"
                    if single_end.exists() and not r1.exists():
                        shutil.move(str(single_end), str(r1))
                    if not r1.exists():
                        raise FileNotFoundError(f"Expected r1 not found: {r1}")
                    paths = [r1]
                    r2 = cache_dir / f"{acc}_2.fastq.gz"
                    if r2.exists():
                        paths.append(r2)
                    downloads[acc] = paths
                    logger.info(f"Downloaded reads {acc}")
            except Exception as e:
                msg = f"Failed to download reads: {e}"
                for sname in set().union(*valid_reads.values()):
                    failed_by_sample.setdefault(sname, []).append(msg)
                logger.warning(msg)

    # ---- 3. download assembly accessions sequentially ------------------------
    for acc, dependents in asm_acc_reqs.items():
        try:
            if is_ncbi_assembly_accession(acc):
                stage = DownloadNcbiAssemblyFasta(accession=acc, output_directory=cache_dir)
            elif is_atb_assembly_accession(acc):
                stage = DownloadAtbAssemblyFasta(accession=acc, output_directory=cache_dir)
            else:
                raise ValueError(f"Unsupported assembly accession: {acc}")
            _check_deps(stage)
            stage.run(dl_ctx)
            stage.error_if_outputs_missing()

            fasta = stage.output.fasta.resolve()
            asm_downloads[acc] = fasta
            logger.info(f"Downloaded assembly {acc}")
        except Exception as e:
            msg = f"Failed to download assembly {acc}: {e}"
            for sname in dependents:
                failed_by_sample.setdefault(sname, []).append(msg)
            logger.warning(msg)

    # ---- 4. collapse per-sample errors into single entries -------------------
    for sname, msgs in failed_by_sample.items():
        failures.append((sname, "; ".join(msgs)))

    # ---- 5. rebuild sample configs with local paths --------------------------
    resolved: dict[str, Any] = {}
    for sample_name, cfg in samples.items():
        sample_cfg = dict(cfg)
        if acc := sample_cfg.pop("reads_accession", None):
            if str(acc) not in downloads:
                continue
            paths = downloads[str(acc)]
            if sample_cfg.get("type") == "long":
                sample_cfg["reads"] = str(paths[0])
            else:
                sample_cfg["reads"] = [str(p) for p in paths]
        if acc := sample_cfg.pop("assembly_accession", None):
            if str(acc) not in asm_downloads:
                continue
            sample_cfg["assembly"] = str(asm_downloads[str(acc)])
        resolved[sample_name] = sample_cfg

    return resolved


def _allocate_multi_cpus(
    *,
    total_cpus: int,
    sample_count: int,
    cpus_per_sample: Optional[int],
) -> tuple[int, int]:
    total_cpus = max(1, int(total_cpus))
    sample_count = max(1, int(sample_count))

    if cpus_per_sample is not None:
        per_sample = min(max(1, int(cpus_per_sample)), total_cpus)
        return per_sample, max(1, total_cpus // per_sample)

    max_auto_parallel = min(sample_count, max(1, total_cpus // 4))
    while max_auto_parallel > 1 and total_cpus % max_auto_parallel != 0:
        max_auto_parallel -= 1

    per_sample = total_cpus // max_auto_parallel
    return per_sample, max_auto_parallel


def run_multi_pipeline(
    snippy_reference_dir: Path,
    samples: Dict[str, Any],
    *,
    prefix: str,
    run_ctx: Context,
    cpus_per_sample: Optional[int] = None,
    stop_on_failure: bool = False,
) -> MultiPipelineResult:
    """
    Special pipeline runner for multi-sample mode. Runs each sample in parallel using ProcessPoolExecutor.
    """

    cpus_per_sample, max_parallel = _allocate_multi_cpus(
        total_cpus=int(run_ctx.cpus),
        sample_count=len(samples),
        cpus_per_sample=cpus_per_sample,
    )

    # Minimal picklable config
    global_config = {
        "reference": str(snippy_reference_dir),
        "outdir": str(run_ctx.outdir),
        "prefix": prefix,
        "run_ctx": run_ctx.model_dump(mode="python"),
        "cpus_per_sample": cpus_per_sample,
    }

    ctx = mp.get_context("spawn")  # safest cross-platform
    failures = []
    successful_samples = []

    # Pre-download all unique sample accessions centrally so parallel workers
    # never contend on the same network resources.
    cache_dir = Path(run_ctx.outdir) / "downloads"
    samples = _predownload_multi_inputs(
        samples,
        cache_dir=cache_dir,
        run_ctx=run_ctx,
        failures=failures,
    )

    jobs = _build_jobs(samples, global_config) if samples else []

    if failures and stop_on_failure:
        logger.warning(
            f"Skipping pipeline execution due to {len(failures)} download "
            f"failure(s) with stop_on_failure=True."
        )
        jobs = []

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

    failed_samples_tsv = Path(run_ctx.outdir) / f"{prefix}.failed.tsv"
    _write_failed_samples_tsv(failures=failures, output_path=failed_samples_tsv)
    if failures:
        logger.warning(f"Skipping {len(failures)} failed samples. See {failed_samples_tsv} for details.")

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

    def _filter_builder_kwargs(builder_cls, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """Filter out any config keys that are not accepted by the builder's model to avoid validation errors."""
        model_fields = getattr(builder_cls, "model_fields", None)
        if not model_fields:
            return kwargs

        allowed_fields = set(model_fields.keys())
        filtered = {key: value for key, value in kwargs.items() if key in allowed_fields}
        dropped = sorted(set(kwargs) - set(filtered))
        if dropped:
            logger.debug(
                f"Dropping unsupported config keys for sample '{sample_name}' ({builder_cls.__name__}): {', '.join(dropped)}"
            )
        return filtered

    sample_cfg = dict(sample_cfg)
    
    sample_type = sample_cfg.get("type")
    del sample_cfg["type"]

    if sample_type == "short":
        reads = sample_cfg.get("reads")
        if reads and isinstance(reads, str):
            reads = [reads]
        if not reads:
            # if reads not provided, expect left/right
            reads = [str(r) for r in (sample_cfg.get("left"), sample_cfg.get("right")) if r]
        short_builder_kwargs = {
            "reference": config["reference"],
            "reads": reads,
            "bam": str(sample_cfg.get("bam")) if sample_cfg.get("bam") else None,
            "prefix": config["prefix"],
            "sample_name": sample_name,
            **{k: v for k, v in sample_cfg.items() if k not in ["left", "right", "bam", "reads", "sample_name"]},
        }
        pipeline = ShortPipelineBuilder(**_filter_builder_kwargs(ShortPipelineBuilder, short_builder_kwargs)).build()

    elif sample_type == "long":
        long_builder_kwargs = {
            "reference": config["reference"],
            "reads": str(sample_cfg.get("reads")) if sample_cfg.get("reads") else None,
            "bam": str(sample_cfg.get("bam")) if sample_cfg.get("bam") else None,
            "prefix": config["prefix"],
            "sample_name": sample_name,
            **{k: v for k, v in sample_cfg.items() if k not in ["reads", "bam", "sample_name"]},
        }
        pipeline = LongPipelineBuilder(**_filter_builder_kwargs(LongPipelineBuilder, long_builder_kwargs)).build()

    elif sample_type == "asm":
        sample_cfg['assembly'] = str(sample_cfg.get("assembly"))
        asm_builder_kwargs = {
            "reference": config["reference"],
            "prefix": config["prefix"],
            "sample_name": sample_name,
            **{k: v for k, v in sample_cfg.items() if k != "sample_name"},
        }
        pipeline = AsmPipelineBuilder(**_filter_builder_kwargs(AsmPipelineBuilder, asm_builder_kwargs)).build()

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
        f"{prefix}.all.vcf.breakdown.tsv": outdir / f"{prefix}.vcf.breakdown.tsv",
        f"{prefix}.qc.tsv": outdir / f"{prefix}.qc.tsv",
    })

    for filename, combined_path in tsv_outputs.items():
        sample_files = [samples_dir / sample_name / filename for sample_name in sample_names]
        if filename == f"{prefix}.qc.tsv":
            _concat_qc_tsv_files(sample_files=sample_files, output_path=combined_path)
        else:
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


def _concat_qc_tsv_files(sample_files: list[Path], output_path: Path) -> None:
    rows: list[dict[str, str]] = []
    seen_fields: set[str] = set()

    for sample_file in sample_files:
        if not sample_file.exists():
            continue

        with sample_file.open("r", newline="") as in_handle:
            reader = csv.DictReader(in_handle, delimiter="\t")
            seen_fields.update(reader.fieldnames or [])
            rows.extend(reader)

    if not rows:
        output_path.unlink(missing_ok=True)
        return

    fieldnames = list(SAMPLE_QC_COLUMNS)
    fieldnames.extend(sorted(seen_fields - set(fieldnames)))

    with output_path.open("w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def add_core_alignment_qc(
    *,
    qc_tsv: Path,
    core_aligned_tsv: Path,
) -> None:
    if not qc_tsv.exists() or not core_aligned_tsv.exists():
        return

    aligned_by_sample: dict[str, str] = {}
    with core_aligned_tsv.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = row.get("sequence")
            if sample and sample != "reference":
                aligned_by_sample[sample] = row.get("aligned", "")

    with qc_tsv.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    if not rows:
        return
    if "core_aligned_percent" not in fieldnames:
        fieldnames.append("core_aligned_percent")

    for row in rows:
        sample = row.get("sample", "")
        row["core_aligned_percent"] = aligned_by_sample.get(sample, row.get("core_aligned_percent", ""))

    with qc_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


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
