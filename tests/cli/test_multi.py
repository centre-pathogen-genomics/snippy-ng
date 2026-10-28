import json
import csv
from pathlib import Path
import pytest
from click.testing import CliRunner
from unittest.mock import patch
from types import SimpleNamespace

from snippy_ng.cli import snippy_ng
import snippy_ng.pipelines as _pl
from snippy_ng.context import Context
from snippy_ng.pipelines.multi import _run_one_sample
from snippy_ng.stages.setup import LoadReferenceFromMetadataFile
from tests.cli.helpers import get_bad_reference_target, make_prepared_reference, stub_load_or_prepare_reference


@pytest.fixture(autouse=True)
def stub_everything(stub_pipeline, stub_reference_format, stub_common_stages, 
                     stub_short_stages, stub_long_stages, stub_asm_stages):
    """
    Combine all required stubs for multi-sample CLI tests.
    This fixture runs automatically for every test.
    """
    pass


@pytest.fixture
def mock_multi_pipeline():
    """Mock the run_multi_pipeline function to avoid subprocess issues in tests."""
    with patch('snippy_ng.pipelines.multi.run_multi_pipeline') as mock:
        mock.return_value = ([], [])
        yield mock


def _write_multi_sample_inputs(tmp_path, samples_data):
    if not isinstance(samples_data, dict):
        return

    for sample_cfg in samples_data.get("samples", {}).values():
        sample_type = sample_cfg.get("type")
        if sample_type == "short":
            reads = sample_cfg.get("reads", [])
            for read_file in reads:
                (tmp_path / read_file).write_text(">read\nATCG")
            if "left" in sample_cfg:
                (tmp_path / sample_cfg["left"]).write_text(">read\nATCG")
            if "right" in sample_cfg:
                (tmp_path / sample_cfg["right"]).write_text(">read\nATCG")
        elif sample_type == "long":
            reads = sample_cfg.get("reads")
            if reads:
                (tmp_path / reads).write_text(">read\nATCG")
        elif sample_type == "asm":
            assembly = sample_cfg.get("assembly")
            if assembly:
                (tmp_path / assembly).write_text(">contig\nATCG")


def _write_multi_config(tmp_path, config_type, samples_data):
    config_file = tmp_path / f"config.{config_type}"
    if config_type == "json":
        if isinstance(samples_data, dict):
            for sample_cfg in samples_data.get("samples", {}).values():
                if "reads" in sample_cfg and isinstance(sample_cfg["reads"], list):
                    sample_cfg["reads"] = [str(tmp_path / r) for r in sample_cfg["reads"]]
                elif "reads" in sample_cfg and isinstance(sample_cfg["reads"], str):
                    sample_cfg["reads"] = str(tmp_path / sample_cfg["reads"])
                if "left" in sample_cfg:
                    sample_cfg["left"] = str(tmp_path / sample_cfg["left"])
                if "right" in sample_cfg:
                    sample_cfg["right"] = str(tmp_path / sample_cfg["right"])
                if "assembly" in sample_cfg:
                    sample_cfg["assembly"] = str(tmp_path / sample_cfg["assembly"])
            if "reference" in samples_data:
                samples_data["reference"] = str(tmp_path / samples_data["reference"])
        config_file.write_text(json.dumps(samples_data, indent=2))
        return config_file

    content = samples_data
    for filename in ["reads1_R1.fq", "reads1_R2.fq", "reads2_R1.fq", "reads2_R2.fq"]:
        content = content.replace(filename, str(tmp_path / filename))
    config_file.write_text(content)
    return config_file


def _stub_multi_reference_loader(monkeypatch, tmp_path):
    _, prepared_ref = make_prepared_reference(tmp_path, dirname="reference")
    stub_load_or_prepare_reference(monkeypatch, prepared_ref)


def _expected_sample_names(config_type, samples_data):
    if isinstance(samples_data, dict):
        return list(samples_data.get("samples", {}).keys())

    lines = [line for line in samples_data.strip().splitlines() if line]
    if len(lines) < 2:
        return []
    delimiter = "\t" if config_type == "tsv" else ","
    reader = csv.DictReader(lines, delimiter=delimiter)
    return [row["sample"] for row in reader if row.get("sample")]


def _materialize_mock_sample_outputs(outdir, prefix, sample_names):
    variant_bases = ["A", "C", "G", "T"]
    for idx, sample_name in enumerate(sample_names):
        sample_dir = outdir / "samples" / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)
        seq = f"{variant_bases[(idx + 1) % len(variant_bases)]}TCG"
        (sample_dir / f"{prefix}.fna").write_text(f">dummy\n{seq}\n")
        (sample_dir / f"{prefix}.vcf.summary.tsv").write_text(
            "sample\ttotal\tpass\n"
            f"{sample_name}\t1\t1\n"
        )
        (sample_dir / f"{prefix}.all.vcf.breakdown.tsv").write_text(
            "sample\tsection\tlabel\tcount\n"
            f"{sample_name}\tfilter\tPASS\t1\n"
        )
        (sample_dir / f"{prefix}.qc.tsv").write_text(
            "sample\tpipeline_type\tvcf_total\n"
            f"{sample_name}\tshort\t1\n"
        )


def test_allocate_multi_cpus_uses_full_cpu_budget_for_single_auto_worker():
    from snippy_ng.pipelines.multi import _allocate_multi_cpus

    assert _allocate_multi_cpus(
        total_cpus=6,
        sample_count=2,
        cpus_per_sample=None,
    ) == (6, 1)


def test_allocate_multi_cpus_keeps_even_auto_parallelism():
    from snippy_ng.pipelines.multi import _allocate_multi_cpus

    assert _allocate_multi_cpus(
        total_cpus=8,
        sample_count=2,
        cpus_per_sample=None,
    ) == (4, 2)
    assert _allocate_multi_cpus(
        total_cpus=12,
        sample_count=3,
        cpus_per_sample=None,
    ) == (4, 3)


def test_allocate_multi_cpus_reduces_auto_parallelism_to_avoid_idle_cpus():
    from snippy_ng.pipelines.multi import _allocate_multi_cpus

    assert _allocate_multi_cpus(
        total_cpus=10,
        sample_count=3,
        cpus_per_sample=None,
    ) == (5, 2)
    assert _allocate_multi_cpus(
        total_cpus=16,
        sample_count=3,
        cpus_per_sample=None,
    ) == (8, 2)


def test_allocate_multi_cpus_respects_explicit_cpus_per_sample():
    from snippy_ng.pipelines.multi import _allocate_multi_cpus

    assert _allocate_multi_cpus(
        total_cpus=6,
        sample_count=10,
        cpus_per_sample=2,
    ) == (2, 3)
    assert _allocate_multi_cpus(
        total_cpus=6,
        sample_count=10,
        cpus_per_sample=99,
    ) == (6, 1)


def test_add_core_alignment_qc_updates_combined_qc(tmp_path):
    from snippy_ng.pipelines.multi import add_core_alignment_qc

    qc_tsv = tmp_path / "snippy.qc.tsv"
    core_aligned_tsv = tmp_path / "core.aligned.tsv"
    qc_tsv.write_text(
        "sample\tpipeline_type\tvcf_total\tcore_aligned_percent\n"
        "sample1\tshort\t1\t\n"
        "sample2\tlong\t2\t\n"
    )
    core_aligned_tsv.write_text(
        "sequence\taligned\n"
        "reference\t100.00\n"
        "sample1\t99.50\n"
        "sample2\t97.25\n"
    )

    add_core_alignment_qc(qc_tsv=qc_tsv, core_aligned_tsv=core_aligned_tsv)

    rows = list(csv.DictReader(qc_tsv.open("r", newline=""), delimiter="\t"))
    assert rows[0]["core_aligned_percent"] == "99.50"
    assert rows[1]["core_aligned_percent"] == "97.25"


def test_combine_sample_stats_only_writes_multi_breakdown_and_qc(tmp_path):
    from snippy_ng.pipelines.multi import _combine_sample_stats

    sample_names = ["sample1", "sample2"]
    _materialize_mock_sample_outputs(tmp_path, "snippy", sample_names)

    _combine_sample_stats(
        samples_dir=tmp_path / "samples",
        sample_names=sample_names,
        prefix="snippy",
        outdir=tmp_path,
    )

    assert (tmp_path / "snippy.vcf.breakdown.tsv").exists()
    assert (tmp_path / "snippy.qc.tsv").exists()


def test_combine_sample_stats_unions_variable_qc_headers(tmp_path):
    from snippy_ng.pipelines.multi import _combine_sample_stats

    samples_dir = tmp_path / "samples"
    (samples_dir / "sample1").mkdir(parents=True)
    (samples_dir / "sample2").mkdir(parents=True)
    for sample in ("sample1", "sample2"):
        (samples_dir / sample / "snippy.vcf.breakdown.tsv").write_text(
            "sample\tsection\tlabel\tcount\n"
            f"{sample}\tfilter\tPASS\t1\n"
        )
    (samples_dir / "sample1" / "snippy.qc.tsv").write_text(
        "sample\tpipeline_type\tvcf_total\talignment_mapped_percent\n"
        "sample1\tshort\t1\t99.0\n"
    )
    (samples_dir / "sample2" / "snippy.qc.tsv").write_text(
        "sample\tpipeline_type\tvcf_total\tfinal_gap_fraction\n"
        "sample2\tasm\t2\t0.01\n"
    )

    _combine_sample_stats(
        samples_dir=samples_dir,
        sample_names=["sample1", "sample2"],
        prefix="snippy",
        outdir=tmp_path,
    )

    rows = list(csv.DictReader((tmp_path / "snippy.qc.tsv").open("r", newline=""), delimiter="\t"))
    assert "alignment_mean_depth" in rows[0]
    assert rows[0]["alignment_mapped_percent"] == "99.0"
    assert rows[0]["final_gap_fraction"] == ""
    assert rows[1]["alignment_mapped_percent"] == ""
    assert rows[1]["alignment_mean_depth"] == ""
    assert rows[1]["final_gap_fraction"] == "0.01"


##############################################################################
#                               TEST DATA                                     #
##############################################################################

@pytest.mark.parametrize(
    "case_name, config_type, samples_data, expect_exit, expect_run",
    [
        (
            "json_short_reads",
            "json",
            {
                "samples": {
                    "sample1": {
                        "type": "short",
                        "reads": ["reads1_R1.fq", "reads1_R2.fq"]
                    },
                    "sample2": {
                        "type": "short",
                        "reads": ["reads2_R1.fq", "reads2_R2.fq"]
                    }
                }
            },
            0,
            True,
        ),
        (
            "json_mixed_types",
            "json",
            {
                "samples": {
                    "short_sample": {
                        "type": "short",
                        "reads": ["short_R1.fq", "short_R2.fq"]
                    },
                    "long_sample": {
                        "type": "long",
                        "reads": "long_reads.fq"
                    },
                    "asm_sample": {
                        "type": "asm",
                        "assembly": "assembly.fasta"
                    }
                }
            },
            0,
            True,
        ),
        (
            "json_with_reference",
            "json",
            {
                "reference": "ref.fa",
                "samples": {
                    "sample1": {
                        "type": "short",
                        "reads": ["reads_R1.fq", "reads_R2.fq"]
                    }
                }
            },
            0,
            True,
        ),
        (
            "csv_format",
            "csv",
            "sample,type,left,right\nsample1,short,reads1_R1.fq,reads1_R2.fq\nsample2,short,reads2_R1.fq,reads2_R2.fq\n",
            0,
            True,
        ),
        (
            "tsv_format",
            "tsv",
            "sample\ttype\tleft\tright\nsample1\tshort\treads1_R1.fq\treads1_R2.fq\nsample2\tshort\treads2_R1.fq\treads2_R2.fq\n",
            0,
            True,
        ),
        (
            "outdir_exists",
            "json",
            {
                "samples": {
                    "sample1": {
                        "type": "short",
                        "reads": ["reads_R1.fq", "reads_R2.fq"]
                    }
                }
            },
            2,
            False,
        ),
        (
            "missing_type_field",
            "json",
            {
                "samples": {
                    "sample1": {
                        "reads": ["reads_R1.fq", "reads_R2.fq"]
                    }
                }
            },
            1,
            False,
        ),
        (
            "bad_reference",
            "json",
            {
                "samples": {
                    "sample1": {
                        "type": "short",
                        "reads": ["reads_R1.fq", "reads_R2.fq"]
                    }
                }
            },
            1,
            False,
        ),
    ],
)
def test_multi_cli(monkeypatch, tmp_path, mock_multi_pipeline, case_name, config_type, samples_data, expect_exit, expect_run):
    """
    Parameterised test for the `multi` command.
    """

    # --------------- Arrange --------------------------------------------------
    ref_file = tmp_path / "ref.fa"
    ref_file.write_text(">dummy\nATCG")
    
    _write_multi_sample_inputs(tmp_path, samples_data)
    config_file = _write_multi_config(tmp_path, config_type, samples_data)
    
    outdir = tmp_path / "output"
    
    # Setup test conditions
    if case_name == "outdir_exists":
        outdir.mkdir(parents=True, exist_ok=True)
    
    if case_name == "bad_reference":
        monkeypatch.setattr("snippy_ng.pipelines.common.guess_reference_format", lambda _: None)
        monkeypatch.setattr(
            get_bad_reference_target("multi"),
            lambda *args, **kwargs: (_ for _ in ()).throw(ValueError("Could not determine reference format")),
        )
    
    # Build command args
    args = ["multi", str(config_file)]
    
    # Add reference if not in JSON
    if config_type != "json" or (isinstance(samples_data, dict) and "reference" not in samples_data):
        args.extend(["--reference", str(ref_file)])

    args.extend(["--outdir", str(outdir)])
    args.append("--skip-check")

    runner = CliRunner()

    class DummyCorePipeline:
        def run(self, ctx):
            if not ctx.check:
                Path(ctx.outdir).mkdir(parents=True, exist_ok=True)
            return 0

    class DummyCorePipelineBuilder:
        def __init__(self, **_kwargs):
            pass

        def build(self):
            return DummyCorePipeline()

    monkeypatch.setattr("snippy_ng.pipelines.core.CorePipelineBuilder", DummyCorePipelineBuilder)

    if expect_run:
        sample_names = _expected_sample_names(config_type, samples_data)

        def _fake_run_multi_pipeline(**kwargs):
            _materialize_mock_sample_outputs(outdir, kwargs["prefix"], sample_names)
            return sample_names, []

        mock_multi_pipeline.side_effect = _fake_run_multi_pipeline

    # --------------- Act ------------------------------------------------------
    result = runner.invoke(snippy_ng, args)

    # --------------- Assert ---------------------------------------------------
    assert result.exit_code == expect_exit, result.output
    
    if expect_run:
        # Check that run_multi_pipeline was called
        assert mock_multi_pipeline.called, "run_multi_pipeline should have been called"
        
        # Check that pipeline was actually run (DummyPipeline.last.ran)
        # For multi, multiple pipelines are created but we can check the alignment pipeline
        assert _pl.SnippyPipeline.last is not None
        assert _pl.SnippyPipeline.last.ran
        
        # Check that core alignment directory was created
        core_dir = outdir / "core"
        assert core_dir.exists(), "Core alignment directory was not created"
    else:
        if expect_exit == 2:
            # outdir_exists case - fails before pipeline creation
            pass
        else:
            # Error cases - may or may not have created pipeline
            pass


def test_multi_cli_duplicate_sample_names(tmp_path):
    """Test that duplicate sample names in CSV/TSV are rejected."""
    ref_file = tmp_path / "ref.fa"
    ref_file.write_text(">dummy\nATCG")
    
    # Create CSV with duplicate sample names
    config_file = tmp_path / "config.csv"
    config_file.write_text("sample,type,left,right\n" +
                          "sample1,short,r1.fq,r2.fq\n" +
                          "sample1,short,r3.fq,r4.fq\n")
    
    # Create dummy files
    for f in ["r1.fq", "r2.fq", "r3.fq", "r4.fq"]:
        (tmp_path / f).write_text(">read\nATCG")
    
    outdir = tmp_path / "output"
    
    args = ["multi", str(config_file), "--reference", str(ref_file), 
            "--outdir", str(outdir), "--skip-check"]
    
    runner = CliRunner()
    result = runner.invoke(snippy_ng, args)
    
    # Should fail with error about duplicate sample
    assert result.exit_code != 0
    assert "Duplicate sample name" in result.output


def test_multi_cli_csv_without_reference(tmp_path):
    """Test that CSV/TSV format requires --reference option."""
    config_file = tmp_path / "config.csv"
    config_file.write_text("sample,type,left,right\n" +
                          "sample1,short,r1.fq,r2.fq\n")
    
    outdir = tmp_path / "output"
    
    args = ["multi", str(config_file), "--outdir", str(outdir), "--skip-check"]
    
    runner = CliRunner()
    result = runner.invoke(snippy_ng, args)
    
    # Should fail with error about missing reference
    assert result.exit_code != 0
    assert "Reference must be provided" in result.output


def test_multi_cli_downloads_config_from_url(tmp_path):
    """The CONFIG argument accepts an HTTP(S) URL and downloads it with urllib."""
    import http.server
    import socketserver
    import threading

    (tmp_path / "config.csv").write_text(
        "sample,type,left,right\n" + "sample1,short,r1.fq,r2.fq\n"
    )

    class Handler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=str(tmp_path), **kwargs)

        def log_message(self, *args):
            pass

    runner = CliRunner()
    with socketserver.TCPServer(("127.0.0.1", 0), Handler) as httpd:
        port = httpd.server_address[1]
        threading.Thread(target=httpd.serve_forever, daemon=True).start()
        try:
            result = runner.invoke(
                snippy_ng,
                [
                    "multi",
                    f"http://127.0.0.1:{port}/config.csv",
                    "--outdir",
                    str(tmp_path / "output"),
                    "--skip-check",
                ],
            )
            missing = runner.invoke(
                snippy_ng,
                [
                    "multi",
                    f"http://127.0.0.1:{port}/missing.csv",
                    "--outdir",
                    str(tmp_path / "output_missing"),
                    "--skip-check",
                ],
            )
        finally:
            httpd.shutdown()
            httpd.server_close()

    # The remote config was downloaded and parsed, so the CSV reference check fires
    assert result.exit_code != 0
    assert "Reference must be provided" in result.output
    # Download failures are reported against the CONFIG argument
    assert missing.exit_code == 2
    assert "Could not download config" in missing.output


def test_multi_cli_reports_missing_local_config(tmp_path):
    """A missing local config path is reported against the CONFIG argument."""
    result = CliRunner().invoke(
        snippy_ng,
        [
            "multi",
            str(tmp_path / "nope.csv"),
            "--outdir",
            str(tmp_path / "output"),
            "--skip-check",
        ],
    )

    assert result.exit_code == 2
    assert "does not exist" in result.output


def test_run_one_sample_passes_disambiguated_sample_name_to_builder(monkeypatch, tmp_path):
    captured = {}

    class DummyPipeline:
        def run(self, _ctx):
            return 0

    class DummyShortPipelineBuilder:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def build(self):
            return DummyPipeline()

    monkeypatch.setattr("snippy_ng.pipelines.short.ShortPipelineBuilder", DummyShortPipelineBuilder)

    job = (
        "outbreak_b-JKD6159",
        {
            "type": "short",
            "left": str(tmp_path / "JKD6159_R1.fastq.gz"),
            "right": str(tmp_path / "JKD6159_R2.fastq.gz"),
        },
        {
            "reference": str(tmp_path / "reference.fa"),
            "outdir": str(tmp_path / "out"),
            "prefix": "snippy",
            "run_ctx": Context(outdir=tmp_path / "out", cpus=4).model_dump(mode="python"),
            "cpus_per_sample": 2,
        },
    )

    result = _run_one_sample(job)

    assert result == "outbreak_b-JKD6159"
    assert captured["sample_name"] == "outbreak_b-JKD6159"
    assert captured["reads"] == [
        str(tmp_path / "JKD6159_R1.fastq.gz"),
        str(tmp_path / "JKD6159_R2.fastq.gz"),
    ]


def test_run_one_sample_ignores_unsupported_keys_for_asm_builder(monkeypatch, tmp_path):
    captured = {}

    class DummyPipeline:
        def run(self, _ctx):
            return 0

    class DummyAsmPipelineBuilder:
        model_fields = {"reference": None, "prefix": None, "sample_name": None, "assembly": None}

        def __init__(self, **kwargs):
            captured.update(kwargs)

        def build(self):
            return DummyPipeline()

    monkeypatch.setattr("snippy_ng.pipelines.asm.AsmPipelineBuilder", DummyAsmPipelineBuilder)

    job = (
        "example",
        {
            "type": "asm",
            "assembly": str(tmp_path / "assembly.fa"),
            "clair3_model": "clair3_models/r1041_e82_400bps_sup_v520",
        },
        {
            "reference": str(tmp_path / "reference.fa"),
            "outdir": str(tmp_path / "out"),
            "prefix": "snippy",
            "run_ctx": Context(outdir=tmp_path / "out", cpus=4).model_dump(mode="python"),
            "cpus_per_sample": 2,
        },
    )

    result = _run_one_sample(job)

    assert result == "example"
    assert captured["sample_name"] == "example"
    assert captured["assembly"] == str(tmp_path / "assembly.fa")
    assert "clair3_model" not in captured


def test_multi_cli_default_keeps_going_and_uses_only_successful_samples_for_core(monkeypatch, tmp_path):
    ref_file = tmp_path / "ref.fa"
    ref_file.write_text(">dummy\nATCG")

    r1 = tmp_path / "sample1_R1.fq"
    r2 = tmp_path / "sample1_R2.fq"
    r1.write_text("@read\nACGT\n+\nIIII\n")
    r2.write_text("@read\nTGCA\n+\nIIII\n")

    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps({
        "samples": {
            "sample1": {"type": "short", "reads": [str(r1), str(r2)]},
            "sample2": {"type": "short", "reads": [str(r1), str(r2)]},
        }
    }))

    captured_multi = {}
    captured_core = {}

    def fake_run_multi_pipeline(**kwargs):
        captured_multi.update(kwargs)
        return ["sample1"], []

    monkeypatch.setattr("snippy_ng.pipelines.multi.run_multi_pipeline", fake_run_multi_pipeline)

    _stub_multi_reference_loader(monkeypatch, tmp_path)

    class DummyCorePipeline:
        def build(self):
            return self

        def run(self, _ctx):
            return 0

    class DummyCorePipelineBuilder:
        def __init__(self, **kwargs):
            captured_core.update(kwargs)

        def build(self):
            return DummyCorePipeline()

    monkeypatch.setattr("snippy_ng.pipelines.core.CorePipelineBuilder", DummyCorePipelineBuilder)

    outdir = tmp_path / "output"
    runner = CliRunner()
    result = runner.invoke(
        snippy_ng,
        [
            "multi",
            str(config_file),
            "--reference",
            str(ref_file),
            "--outdir",
            str(outdir),
            "--skip-check",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured_core["snippy_dirs"] == [str(outdir / "samples" / "sample1")]
    assert captured_core["mask"] is None
    assert captured_multi["run_ctx"].log_path == (outdir / "LOG.txt").absolute()


def test_multi_cli_passes_mask_to_core(monkeypatch, tmp_path):
    ref_file = tmp_path / "ref.fa"
    ref_file.write_text(">dummy\nATCG")

    mask_file = tmp_path / "mask.bed"
    mask_file.write_text("dummy\t0\t2\n")

    r1 = tmp_path / "sample1_R1.fq"
    r2 = tmp_path / "sample1_R2.fq"
    r1.write_text("@read\nACGT\n+\nIIII\n")
    r2.write_text("@read\nTGCA\n+\nIIII\n")

    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps({
        "samples": {
            "sample1": {"type": "short", "reads": [str(r1), str(r2)]},
        }
    }))

    captured_core = {}

    monkeypatch.setattr(
        "snippy_ng.pipelines.multi.run_multi_pipeline",
        lambda **kwargs: (["sample1"], []),
    )
    _stub_multi_reference_loader(monkeypatch, tmp_path)

    class DummyCorePipeline:
        def run(self, _ctx):
            return 0

    class DummyCorePipelineBuilder:
        def __init__(self, **kwargs):
            captured_core.update(kwargs)

        def build(self):
            return DummyCorePipeline()

    monkeypatch.setattr("snippy_ng.pipelines.core.CorePipelineBuilder", DummyCorePipelineBuilder)

    outdir = tmp_path / "output"
    result = CliRunner().invoke(
        snippy_ng,
        [
            "multi",
            str(config_file),
            "--reference",
            str(ref_file),
            "--mask",
            str(mask_file),
            "--outdir",
            str(outdir),
            "--skip-check",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured_core["mask"] == mask_file


def test_multi_cli_partial_success_exits_nonzero_after_core(monkeypatch, tmp_path):
    ref_file = tmp_path / "ref.fa"
    ref_file.write_text(">dummy\nATCG")

    r1 = tmp_path / "sample1_R1.fq"
    r2 = tmp_path / "sample1_R2.fq"
    r1.write_text("@read\nACGT\n+\nIIII\n")
    r2.write_text("@read\nTGCA\n+\nIIII\n")

    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps({
        "samples": {
            "sample1": {"type": "short", "reads": [str(r1), str(r2)]},
            "sample2": {"type": "short", "reads": [str(r1), str(r2)]},
        }
    }))

    def fake_run_multi_pipeline(**kwargs):
        return ["sample1"], [("sample2", "boom")]

    monkeypatch.setattr("snippy_ng.pipelines.multi.run_multi_pipeline", fake_run_multi_pipeline)

    _stub_multi_reference_loader(monkeypatch, tmp_path)

    class DummyCorePipeline:
        def run(self, _ctx):
            return 0

    class DummyCorePipelineBuilder:
        def __init__(self, **kwargs):
            pass

        def build(self):
            return DummyCorePipeline()

    monkeypatch.setattr("snippy_ng.pipelines.core.CorePipelineBuilder", DummyCorePipelineBuilder)

    outdir = tmp_path / "output"
    runner = CliRunner()
    result = runner.invoke(
        snippy_ng,
        [
            "multi",
            str(config_file),
            "--reference",
            str(ref_file),
            "--outdir",
            str(outdir),
            "--skip-check",
        ],
    )

    assert result.exit_code != 0
    assert "Some samples failed" in result.output


def test_multi_cli_stop_on_failure_passes_fail_fast_mode_to_runner(monkeypatch, tmp_path):
    ref_file = tmp_path / "ref.fa"
    ref_file.write_text(">dummy\nATCG")

    r1 = tmp_path / "sample1_R1.fq"
    r2 = tmp_path / "sample1_R2.fq"
    r1.write_text("@read\nACGT\n+\nIIII\n")
    r2.write_text("@read\nTGCA\n+\nIIII\n")

    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps({
        "samples": {
            "sample1": {"type": "short", "reads": [str(r1), str(r2)]},
        }
    }))

    captured_multi = {}

    def fake_run_multi_pipeline(**kwargs):
        captured_multi.update(kwargs)
        return ["sample1"], []

    monkeypatch.setattr("snippy_ng.pipelines.multi.run_multi_pipeline", fake_run_multi_pipeline)

    _stub_multi_reference_loader(monkeypatch, tmp_path)

    class DummyCorePipeline:
        def run(self, _ctx):
            return 0

    class DummyCorePipelineBuilder:
        def __init__(self, **kwargs):
            pass

        def build(self):
            return DummyCorePipeline()

    monkeypatch.setattr("snippy_ng.pipelines.core.CorePipelineBuilder", DummyCorePipelineBuilder)

    outdir = tmp_path / "output"
    runner = CliRunner()
    result = runner.invoke(
        snippy_ng,
        [
            "multi",
            str(config_file),
            "--reference",
            str(ref_file),
            "--outdir",
            str(outdir),
            "--stop-on-failure",
            "--skip-check",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured_multi["stop_on_failure"] is True
    assert captured_multi["run_ctx"].log_path == (outdir / "LOG.txt").absolute()


def test_multi_cli_uses_existing_prepared_reference_directory(monkeypatch, tmp_path):
    prepared_ref_dir, ref_file = make_prepared_reference(tmp_path, dirname="prepared_reference")
    monkeypatch.setattr(
        "snippy_ng.pipelines.common.load_or_prepare_reference",
        lambda *args, **kwargs: LoadReferenceFromMetadataFile(metadata=prepared_ref_dir / "metadata.json"),
    )

    r1 = tmp_path / "sample1_R1.fq"
    r2 = tmp_path / "sample1_R2.fq"
    r1.write_text("@read\nACGT\n+\nIIII\n")
    r2.write_text("@read\nTGCA\n+\nIIII\n")

    config_file = tmp_path / "config.csv"
    config_file.write_text(
        "sample,type,left,right\n"
        f"sample1,short,{r1},{r2}\n"
    )

    captured_multi = {}

    def fake_run_multi_pipeline(**kwargs):
        captured_multi.update(kwargs)
        return ["sample1"], []

    monkeypatch.setattr("snippy_ng.pipelines.multi.run_multi_pipeline", fake_run_multi_pipeline)

    class DummyCorePipeline:
        def run(self, _ctx):
            return 0

    class DummyCorePipelineBuilder:
        def __init__(self, **kwargs):
            pass

        def build(self):
            return DummyCorePipeline()

    monkeypatch.setattr("snippy_ng.pipelines.core.CorePipelineBuilder", DummyCorePipelineBuilder)

    outdir = tmp_path / "output"
    runner = CliRunner()
    result = runner.invoke(
        snippy_ng,
        [
            "multi",
            str(config_file),
            "--reference",
            str(prepared_ref_dir),
            "--outdir",
            str(outdir),
            "--skip-check",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured_multi["snippy_reference_dir"] == prepared_ref_dir


def test_multi_cli_accepts_reference_accession(monkeypatch, tmp_path):
    r1 = tmp_path / "sample1_R1.fq"
    r2 = tmp_path / "sample1_R2.fq"
    r1.write_text("@read\nACGT\n+\nIIII\n")
    r2.write_text("@read\nTGCA\n+\nIIII\n")

    config_file = tmp_path / "config.csv"
    config_file.write_text(
        "sample,type,left,right\n"
        f"sample1,short,{r1},{r2}\n"
    )

    _stub_multi_reference_loader(monkeypatch, tmp_path)

    captured_download = {}

    def fake_download_assembly(reference_accession, stages, output_directory):
        captured_download["reference_accession"] = reference_accession
        captured_download["output_directory"] = output_directory
        downloaded_reference = Path(output_directory) / f"{reference_accession}.fa"
        stages.append(SimpleNamespace(output=SimpleNamespace(fasta=downloaded_reference)))
        return downloaded_reference

    monkeypatch.setattr("snippy_ng.pipelines.common.download_reference", fake_download_assembly)

    captured_multi = {}

    def fake_run_multi_pipeline(**kwargs):
        captured_multi.update(kwargs)
        return ["sample1"], []

    monkeypatch.setattr("snippy_ng.pipelines.multi.run_multi_pipeline", fake_run_multi_pipeline)

    class DummyCorePipeline:
        def run(self, _ctx):
            return 0

    class DummyCorePipelineBuilder:
        def __init__(self, **kwargs):
            pass

        def build(self):
            return DummyCorePipeline()

    monkeypatch.setattr("snippy_ng.pipelines.core.CorePipelineBuilder", DummyCorePipelineBuilder)

    outdir = tmp_path / "output"
    result = CliRunner().invoke(
        snippy_ng,
        [
            "multi",
            str(config_file),
            "--reference",
            "SAMN123456",
            "--outdir",
            str(outdir),
            "--skip-check",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured_download == {
        "reference_accession": "SAMN123456",
        "output_directory": outdir / "reference",
    }
    assert captured_multi["snippy_reference_dir"] == tmp_path / "reference"


@pytest.fixture
def monkeypatch_download_stages(monkeypatch):
    """Replace sracha subprocess with a fake that creates output files."""
    from pathlib import Path
    import subprocess

    call_counter = {"reads": 0, "asm": 0}

    # Make sracha dependency check pass silently
    monkeypatch.setattr("snippy_ng.pipelines.multi.sracha.check", lambda: ("/usr/bin/sracha", "1.0"))

    orig_subprocess_run = subprocess.run

    def fake_subprocess_run(cmd, *args, **kwargs):
        asm_cmd = cmd if isinstance(cmd, (list, tuple)) else None
        if asm_cmd and asm_cmd[:2] == ["sracha", "get"] and "--accession-list" in asm_cmd:
            call_counter["reads"] += 1
            ai = asm_cmd.index("--accession-list")
            acc_list_path = Path(asm_cmd[ai + 1])
            out_dir = Path(asm_cmd[ai + 3]) if "--output-dir" in asm_cmd else Path(".")
            for line in acc_list_path.read_text().splitlines():
                acc = line.strip()
                if not acc:
                    continue
                (out_dir / f"{acc}_1.fastq.gz").write_text("fq1\n")
            return subprocess.CompletedProcess(cmd, 0)
        return orig_subprocess_run(cmd, *args, **kwargs)

    monkeypatch.setattr("subprocess.run", fake_subprocess_run)

    # Also replace assembly download stages with fakes (they still use stage.run)
    from types import SimpleNamespace

    class FakeDownloadNcbi:
        _dependencies = []
        def __init__(self, **kwargs):
            self.accession = kwargs["accession"]
            outdir = kwargs.get("output_directory", Path("."))
            self.output = SimpleNamespace(fasta=outdir / f"{self.accession}.fa")
        def run(self, _ctx):
            call_counter["asm"] += 1
            self.output.fasta.parent.mkdir(parents=True, exist_ok=True)
            self.output.fasta.write_text(">seq\nATCG\n")
        def error_if_outputs_missing(self):
            pass

    monkeypatch.setattr("snippy_ng.pipelines.multi.DownloadNcbiAssemblyFasta", FakeDownloadNcbi)
    monkeypatch.setattr("snippy_ng.pipelines.multi.DownloadAtbAssemblyFasta", FakeDownloadNcbi)
    return call_counter


def test_predownload_multi_inputs_replaces_accessions_with_abs_paths(
    tmp_path, monkeypatch_download_stages
):
    from snippy_ng.pipelines.multi import _predownload_multi_inputs
    from pathlib import Path

    call_counter = monkeypatch_download_stages
    original = {
        "S1": {"type": "short", "reads_accession": "SRR0000001"},
        "S2": {"type": "long", "reads_accession": "SRR0000002"},
        "S3": {"type": "short", "reads_accession": "SRR0000003"},
        "S4": {"type": "asm", "assembly_accession": "GCF_000000001.1"},
    }
    cache = tmp_path / "cache"
    failures = []
    result = _predownload_multi_inputs(
        original, cache_dir=cache, run_ctx=Context(cpus=1, quiet=True), failures=failures,
    )

    # original samples dict is not mutated
    assert "reads_accession" in original["S1"]

    # keys removed, local paths in place
    assert "reads_accession" not in result["S1"]
    assert "reads_accession" not in result["S2"]
    assert "assembly_accession" not in result["S4"]

    # short single-end
    s1_reads = result["S1"]["reads"]
    assert isinstance(s1_reads, list) and len(s1_reads) == 1
    assert Path(s1_reads[0]).is_absolute() and Path(s1_reads[0]).exists()

    # long
    s2_reads = result["S2"]["reads"]
    assert isinstance(s2_reads, str)
    assert Path(s2_reads).is_absolute() and Path(s2_reads).exists()

    # asm
    s4_asm = result["S4"]["assembly"]
    assert isinstance(s4_asm, str)
    assert Path(s4_asm).is_absolute() and Path(s4_asm).exists()
    assert failures == []
    # 3 unique SRA calls, 1 assembly call
    assert call_counter["reads"] == 1  # single batch call for all SRA
    assert call_counter["asm"] == 1


def test_predownload_multi_inputs_paired_end(tmp_path, monkeypatch):
    """Short paired-end: create both r1 and r2; reads list has both paths."""
    import subprocess
    from pathlib import Path
    from snippy_ng.pipelines.multi import _predownload_multi_inputs

    monkeypatch.setattr("snippy_ng.pipelines.multi.sracha.check", lambda: ("/usr/bin/sracha", "1.0"))

    def paired_subprocess(cmd, *args, **kwargs):
        if isinstance(cmd, (list, tuple)) and cmd[:2] == ["sracha", "get"] and "--accession-list" in cmd:
            ai = cmd.index("--accession-list")
            acc_list_path = Path(cmd[ai + 1])
            out_dir = Path(cmd[ai + 3]) if "--output-dir" in cmd else Path(".")
            for line in acc_list_path.read_text().splitlines():
                acc = line.strip()
                if not acc:
                    continue
                (out_dir / f"{acc}_1.fastq.gz").write_text("r1\n")
                (out_dir / f"{acc}_2.fastq.gz").write_text("r2\n")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr("subprocess.run", paired_subprocess)

    samples = {"PE": {"type": "short", "reads_accession": "SRR0000004"}}
    result = _predownload_multi_inputs(
        samples, cache_dir=tmp_path / "c", run_ctx=Context(cpus=1, quiet=True), failures=[],
    )
    reads = result["PE"]["reads"]
    assert len(reads) == 2
    assert all(Path(p).exists() for p in reads)


def test_predownload_multi_inputs_deduplicates_accessions(
    tmp_path, monkeypatch_download_stages
):
    from snippy_ng.pipelines.multi import _predownload_multi_inputs

    call_counter = monkeypatch_download_stages
    samples = {
        "A": {"type": "short", "reads_accession": "SRR0000001"},
        "B": {"type": "short", "reads_accession": "SRR0000001"},  # duplicate
    }
    failures = []
    result = _predownload_multi_inputs(
        samples, cache_dir=tmp_path / "cache", run_ctx=Context(cpus=1, quiet=True), failures=failures,
    )

    assert result["A"]["reads"] == result["B"]["reads"]
    assert call_counter["reads"] == 1  # one download for two samples
    assert failures == []


def test_predownload_multi_inputs_conditional_failure(tmp_path, monkeypatch):
    """A failing download drops all batch accessions; none succeed."""
    import subprocess
    from pathlib import Path
    from snippy_ng.pipelines.multi import _predownload_multi_inputs

    monkeypatch.setattr("snippy_ng.pipelines.multi.sracha.check", lambda: ("/usr/bin/sracha", "1.0"))

    def raising_subprocess(cmd, *args, **kwargs):
        if isinstance(cmd, (list, tuple)) and cmd[:2] == ["sracha", "get"] and "--accession-list" in cmd:
            raise RuntimeError("batch failed")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr("subprocess.run", raising_subprocess)

    samples = {
        "Good": {"type": "short", "reads_accession": "SRR0000010"},
        "Bad":  {"type": "short", "reads_accession": "SRR9999990"},
    }
    failures = []
    result = _predownload_multi_inputs(
        samples, cache_dir=tmp_path / "c", run_ctx=Context(cpus=1, quiet=True), failures=failures,
    )
    # Batch: one call; failure drops all
    assert result == {}
    assert len(failures) == 2


def test_predownload_multi_inputs_leaves_original_unchanged(
    tmp_path, monkeypatch_download_stages
):
    from snippy_ng.pipelines.multi import _predownload_multi_inputs

    original = {"X": {"type": "short", "reads_accession": "SRR0000001"}}
    failures = []
    _ = _predownload_multi_inputs(
        original, cache_dir=tmp_path / "c", run_ctx=Context(cpus=1, quiet=True), failures=failures,
    )
    assert "reads_accession" in original["X"]


def test_run_multi_pipeline_jobs_have_no_accessions(tmp_path, monkeypatch_download_stages):
    """_build_jobs receives resolved configs with no accession keys."""
    from pathlib import Path
    from copy import deepcopy
    from snippy_ng.pipelines.multi import _predownload_multi_inputs, _build_jobs

    original = {"X": {"type": "short", "reads_accession": "SRR0000001"}}
    failures = []
    resolved = _predownload_multi_inputs(
        deepcopy(original),
        cache_dir=tmp_path / "c",
        run_ctx=Context(cpus=1, quiet=True),
        failures=failures,
    )
    jobs = _build_jobs(resolved, {"reference": str(tmp_path)})
    assert len(jobs) == 1
    _, job_cfg, *_ = jobs[0]
    assert "reads_accession" not in job_cfg
    assert "reads" in job_cfg
    assert isinstance(job_cfg["reads"], list)


def test_predownload_stop_on_failure_writes_failed_tsv(tmp_path, monkeypatch):
    """stop_on_failure with predownload failures sets jobs=[] and writes failed.tsv."""
    from pathlib import Path
    from copy import deepcopy

    # Inject a failure and a resolved sample into _predownload_multi_inputs
    import snippy_ng.pipelines.multi as m

    def fake_predownload(samples, **kw):
        kw["failures"].append(("Bad", "download failed"))
        return {}

    monkeypatch.setattr(m, "_predownload_multi_inputs", fake_predownload)

    ref_dir = tmp_path / "ref"
    ref_dir.mkdir(); (ref_dir / "metadata.json").write_text("{}")

    samples = {"Bad": {"type": "short", "reads_accession": "SRR0000001"}}
    successes, failures = m.run_multi_pipeline(
        snippy_reference_dir=ref_dir,
        samples=deepcopy(samples),
        prefix="test",
        run_ctx=Context(outdir=tmp_path, cpus=1, quiet=True, skip_check=True),
        stop_on_failure=True,
    )
    assert successes == []
    assert len(failures) == 1
    assert "download" in failures[0][1].lower()
    failed_tsv = tmp_path / "test.failed.tsv"
    assert failed_tsv.exists()
    assert failed_tsv.read_text().strip()


def test_predownload_dependency_check(monkeypatch, tmp_path):
    """Dependencies are validated once per unique dependency name."""
    from pathlib import Path
    from snippy_ng.pipelines.multi import _predownload_multi_inputs

    calls = {"sracha": 0}

    def counted_check():
        calls["sracha"] += 1
        return ("/usr/bin/sracha", "1.0")

    # Patch the module-level sracha imported in multi.py
    import snippy_ng.pipelines.multi as m

    monkeypatch.setattr(m.sracha, "check", counted_check)

    # Two samples sharing the same accession → one batch call
    samples = {
        "A": {"type": "short", "reads_accession": "SRR0000001"},
        "B": {"type": "short", "reads_accession": "SRR0000001"},
    }
    # subprocess.run must succeed so batch proceeds
    import subprocess

    def ok_subprocess(cmd, *a, **kw):
        if isinstance(cmd, (list, tuple)) and cmd[:2] == ["sracha", "get"]:
            ai = cmd.index("--accession-list")
            acc_list_path = Path(cmd[ai + 1])
            out_dir = Path(cmd[ai + 3]) if "--output-dir" in cmd else Path(".")
            for line in acc_list_path.read_text().splitlines():
                acc = line.strip()
                if acc:
                    (out_dir / f"{acc}_1.fastq.gz").write_text("ok\n")
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr("subprocess.run", ok_subprocess)

    _ = _predownload_multi_inputs(
        samples, cache_dir=tmp_path / "c",
        run_ctx=Context(cpus=1, quiet=True),
        failures=[],
    )
    # checked once (deduplicated), not twice
    assert calls["sracha"] == 1

    # Skip-check mode
    calls["sracha"] = 0
    _ = _predownload_multi_inputs(
        samples, cache_dir=tmp_path / "c2",
        run_ctx=Context(cpus=1, quiet=True, skip_check=True),
        failures=[],
    )
    assert calls["sracha"] == 0
