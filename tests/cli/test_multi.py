import json
import csv
from pathlib import Path
import pytest
from click.testing import CliRunner
from unittest.mock import patch

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
        (sample_dir / f"{prefix}.pseudo.fna").write_text(f">dummy\n{seq}\n")
        (sample_dir / f"{prefix}.vcf.summary.tsv").write_text(
            "sample\ttotal\tpass\n"
            f"{sample_name}\t1\t1\n"
        )
        (sample_dir / f"{prefix}.vcf.breakdown.tsv").write_text(
            "sample\tsection\tlabel\tcount\n"
            f"{sample_name}\tfilter\tPASS\t1\n"
        )


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
    assert captured_multi["run_ctx"].log_path == (outdir / "LOG.txt").absolute()


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
