import json
import pytest
from click.testing import CliRunner
from pathlib import Path
from unittest.mock import Mock, patch

from snippy_ng.cli import snippy_ng
import snippy_ng.pipelines as _pl


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
        yield mock


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
            "check_only",
            "json",
            {
                "samples": {
                    "sample1": {
                        "type": "short",
                        "reads": ["reads_R1.fq", "reads_R2.fq"]
                    }
                }
            },
            0,
            False,
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
    
    # Create sample files referenced in config
    if isinstance(samples_data, dict):
        for sample_name, sample_cfg in samples_data.get("samples", {}).items():
            sample_type = sample_cfg.get("type")
            
            if sample_type == "short":
                reads = sample_cfg.get("reads", [])
                for read_file in reads:
                    (tmp_path / read_file).write_text(">read\nATCG")
                # Also create left/right if they exist
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
    
    # Create config file
    config_file = tmp_path / f"config.{config_type}"
    if config_type == "json":
        # Update paths in JSON to be absolute
        if isinstance(samples_data, dict):
            for sample_name, sample_cfg in samples_data.get("samples", {}).items():
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
    else:  # CSV or TSV
        # Update paths in CSV/TSV content
        content = samples_data
        for filename in ["reads1_R1.fq", "reads1_R2.fq", "reads2_R1.fq", "reads2_R2.fq"]:
            content = content.replace(filename, str(tmp_path / filename))
        config_file.write_text(content)
    
    outdir = tmp_path / "output"
    
    # Setup test conditions
    if case_name == "outdir_exists":
        outdir.mkdir(parents=True, exist_ok=True)
    
    if case_name == "bad_reference":
        monkeypatch.setattr("snippy_ng.pipelines.common.guess_reference_format", lambda _: None)
    
    # Build command args
    args = ["multi", str(config_file)]
    
    # Add reference if not in JSON
    if config_type != "json" or (isinstance(samples_data, dict) and "reference" not in samples_data):
        args.extend(["--reference", str(ref_file)])
    
    args.extend(["--outdir", str(outdir)])

    # For check-only mode, use --check WITHOUT --skip-check so validation happens
    if case_name == "check_only":
        args.append("--check")
    else:
        # For other tests, skip dependency checks to avoid test dependencies
        args.append("--skip-check")

    runner = CliRunner()

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
        if case_name == "check_only":
            # In check mode, validation happens but pipeline doesn't run
            assert _pl.SnippyPipeline.last is not None
            assert _pl.SnippyPipeline.last.validated
            assert not _pl.SnippyPipeline.last.ran
        elif expect_exit == 2:
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
