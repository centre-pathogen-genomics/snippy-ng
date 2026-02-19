"""
Shared fixtures and stubs for CLI tests.
"""
import types

import pytest

import snippy_ng.pipelines as _pl


class DummyPipeline:
    """Light-weight stand-in that records what happened."""
    last = None                        # will hold most‑recent instance

    def __init__(self, *_, **__):
        self.validated = False
        self.ran       = False
        DummyPipeline.last = self      # remember myself

    def run(self, quiet=False, create_missing=False, keep_incomplete=False, skip_check=False, check=False, outdir=None):
        """Match the new pipeline structure with run method."""
        self.welcome()
        
        if not skip_check:
            self.validate_dependencies()
        
        if check:
            return None
        
        self.set_working_directory(outdir)
        self.ran = True
        self.cleanup(None)
        self.goodbye()

    def welcome(self):                 pass
    def validate_dependencies(self):   self.validated = True
    def set_working_directory(self, *_): pass
    def cleanup(self, dir):                 pass
    def goodbye(self):                 pass
    def error(self, *_):               pass


def stage_factory(output):
    """Factory function to create dummy stage classes with specified output."""
    class _Stage:
        def __init__(self, *_, **__):
            self.output = types.SimpleNamespace(
                **{out_key: out_val for out_key, out_val in output.items()}
            )
    return _Stage


@pytest.fixture
def stub_pipeline(monkeypatch):
    """
    Replace SnippyPipeline constructor to return DummyPipeline.
    This allows builder.build() to execute (triggering Pydantic validation)
    while still mocking the actual pipeline execution.
    """
    def mock_snippy_pipeline(stages=None):
        # Builder instantiates stages (Pydantic validation happens here)
        # But we return a DummyPipeline for execution
        return DummyPipeline(stages=stages)
    
    # Preserve the 'last' attribute for tests to check
    mock_snippy_pipeline.last = None
    
    # After DummyPipeline is created, update the function's last reference
    original_init = DummyPipeline.__init__
    def patched_init(self, *args, **kwargs):
        original_init(self, *args, **kwargs)
        mock_snippy_pipeline.last = self
    
    monkeypatch.setattr(DummyPipeline, "__init__", patched_init)
    monkeypatch.setattr(_pl, "SnippyPipeline", mock_snippy_pipeline)
    
    return DummyPipeline


@pytest.fixture
def stub_reference_format(monkeypatch):
    """Always recognize the reference format as FASTA."""
    monkeypatch.setattr("snippy_ng.pipelines.common.guess_reference_format", lambda _: "fasta")


@pytest.fixture
def stub_common_stages(monkeypatch, tmp_path):
    """Stub out the common stages used across all pipelines."""
    monkeypatch.setattr(
        "snippy_ng.stages.setup.PrepareReference",
        stage_factory({
            "reference": tmp_path / "ref.fa",
            "gff": tmp_path / "ref.gff",
            "reference_index": tmp_path / "ref.fa.fai",
            "reference_dict": tmp_path / "ref.dict",
            "metadata": tmp_path / "metadata.json"
        }),
    )


@pytest.fixture
def stub_asm_stages(monkeypatch, tmp_path):
    """Stub out assembly-specific stages."""
    monkeypatch.setattr(
        "snippy_ng.stages.alignment.AssemblyAligner",
        stage_factory({"paf": tmp_path / "align.paf"}),
    )
    monkeypatch.setattr(
        "snippy_ng.stages.calling.PAFCaller",
        stage_factory({
            "vcf": tmp_path / "calls.vcf",
            "missing_bed": tmp_path / "missing.bed"
        }),
    )


@pytest.fixture
def stub_long_stages(monkeypatch, tmp_path):
    """Stub out long-read-specific stages."""
    monkeypatch.setattr(
        "snippy_ng.stages.alignment.Minimap2LongReadAligner",
        stage_factory({"cram": tmp_path / "align.cram"}),
    )
    monkeypatch.setattr(
        "snippy_ng.stages.calling.FreebayesCallerLong",
        stage_factory({"vcf": tmp_path / "calls.vcf"}),
    )
    monkeypatch.setattr(
        "snippy_ng.stages.clean_reads.SeqkitCleanLongReads",
        stage_factory({"cleaned_reads": tmp_path / "cleaned.fq"}),
    )
    monkeypatch.setattr(
        "snippy_ng.stages.stats.SeqKitReadStatsBasic",
        stage_factory({"stats": tmp_path / "stats.txt"}),
    )


@pytest.fixture
def stub_short_stages(monkeypatch, tmp_path):
    """Stub out short-read-specific stages."""
    monkeypatch.setattr(
        "snippy_ng.stages.alignment.BWAMEMShortReadAligner",
        stage_factory({"cram": tmp_path / "align.cram"}),
    )
    monkeypatch.setattr(
        "snippy_ng.stages.calling.FreebayesCaller",
        stage_factory({"vcf": tmp_path / "calls.vcf"}),
    )
