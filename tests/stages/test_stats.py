"""
Test module for SeqKit read statistics stages
"""

import csv

import pytest
from pydantic import ValidationError

from snippy_ng.stages.stats import (
    FastaCompositionStats,
    SampleQcSummary,
    SamtoolsAlignmentQcStats,
    SeqKitReadStats,
    SeqKitReadStatsBasic,
    SeqKitReadStatsDetailed,
)
from snippy_ng.context import Context


def test_samtools_alignment_qc_parses_outputs(tmp_path):
    flagstat = tmp_path / "sample.flagstat.txt"
    stats = tmp_path / "sample.stats.txt"
    coverage = tmp_path / "sample.coverage.tsv"
    output = tmp_path / "sample.alignment.tsv"
    flagstat.write_text(
        "100 + 0 in total (QC-passed reads + QC-failed reads)\n"
        "90 + 0 primary\n"
        "5 + 0 secondary\n"
        "5 + 0 supplementary\n"
        "2 + 0 duplicates\n"
        "95 + 0 mapped (95.00% : N/A)\n"
        "80 + 0 properly paired (80.00% : N/A)\n"
    )
    stats.write_text(
        "SN\traw total sequences:\t100\n"
        "SN\treads mapped:\t95\n"
        "SN\treads unmapped:\t5\n"
        "SN\tbases mapped:\t4000\n"
        "SN\terror rate:\t0.001\n"
        "SN\taverage length:\t100\n"
        "SN\taverage quality:\t35\n"
        "SN\tinsert size average:\t450\n"
    )
    coverage.write_text(
        "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
        "chr1\t1\t100\t10\t90\t90\t20\t30\t60\n"
        "chr2\t1\t300\t20\t150\t50\t10\t30\t60\n"
    )

    SamtoolsAlignmentQcStats.write_alignment_qc(flagstat, stats, coverage, output, "sample")

    with output.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["sample"] == "sample"
    assert row["alignment_total"] == "100"
    assert row["alignment_mapped_percent"] == "95.0"
    assert row["alignment_reads_mapped"] == "95"
    assert row["alignment_mean_depth"] == "12.5"
    assert row["alignment_mean_breadth"] == "60.0"


def test_fasta_composition_counts_missing_and_gap_characters(tmp_path):
    fasta = tmp_path / "sample.fna"
    fasta.write_text(">sample\nACGTNNnn--X\n")

    metrics = FastaCompositionStats.count_fasta_composition(fasta)

    assert metrics["final_length"] == 11
    assert metrics["final_acgt"] == 4
    assert metrics["final_N"] == 2
    assert metrics["final_n"] == 2
    assert metrics["final_gap"] == 2
    assert metrics["final_other"] == 1
    assert metrics["final_gap_fraction"] == round(2 / 11, 6)


def test_sample_qc_summary_merges_component_tsvs(tmp_path):
    reads_tsv = tmp_path / "sample.reads.tsv"
    alignment_tsv = tmp_path / "sample.alignment.tsv"
    vcf_tsv = tmp_path / "sample.vcf.summary.tsv"
    fasta_tsv = tmp_path / "sample.fasta.tsv"
    output = tmp_path / "sample.qc.tsv"
    reads_tsv.write_text(
        "sample\tfile\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\n"
        "sample\tR1.fq\tFASTQ\tDNA\t10\t1000\t90\t100\t110\n"
        "sample\tR2.fq\tFASTQ\tDNA\t10\t1000\t90\t100\t110\n"
    )
    alignment_tsv.write_text(
        "sample\talignment_mapped_percent\talignment_mean_depth\n"
        "sample\t95.5\t30.25\n"
    )
    vcf_tsv.write_text(
        "sample\ttotal\tpass\tsnp\n"
        "sample\t7\t6\t5\n"
    )
    fasta_tsv.write_text(
        "sample\tfinal_length\tfinal_gap_fraction\tfinal_N_fraction\n"
        "sample\t100\t0.01\t0.02\n"
    )

    SampleQcSummary.write_sample_qc_summary(
        output,
        "sample",
        "short",
        reads_tsv,
        alignment_tsv,
        vcf_tsv,
        fasta_tsv,
    )

    with output.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["sample"] == "sample"
    assert row["pipeline_type"] == "short"
    assert row["read_num_seqs"] == "20"
    assert row["vcf_total"] == "7"
    assert row["alignment_mapped_percent"] == "95.5"
    assert row["final_gap_fraction"] == "0.01"
    assert "core_aligned_percent" not in row


def test_sample_qc_summary_omits_blank_component_columns(tmp_path):
    vcf_tsv = tmp_path / "sample.vcf.summary.tsv"
    output = tmp_path / "sample.qc.tsv"
    vcf_tsv.write_text(
        "sample\ttotal\tpass\n"
        "sample\t0\t0\n"
    )

    SampleQcSummary.write_sample_qc_summary(
        output,
        "sample",
        "asm",
        None,
        None,
        vcf_tsv,
        None,
    )

    with output.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row == {
        "sample": "sample",
        "pipeline_type": "asm",
        "vcf_total": "0",
        "vcf_pass": "0",
    }


class TestSeqKitReadStats:
    """Test SeqKitReadStats stage"""
    
    def test_init_valid_inputs(self, tmp_path):
        """Test initialization with valid inputs"""
        # Create test files
        read1 = tmp_path / "reads1.fastq"
        read2 = tmp_path / "reads2.fastq" 
        read1.touch()
        read2.touch()
        
        stage = SeqKitReadStats(
            reads=[str(read1), str(read2)],
            prefix="test_stats",
        )
        
        assert stage.reads == [read1, read2]
        assert stage.prefix == "test_stats"
        assert stage.all_stats is True
        assert stage.tabular is True
        assert stage.basename_only is False
        assert stage.skip_errors is True
        assert stage.fastq_encoding == "sanger"
        assert stage.gap_letters == "- ."
        
    def test_init_empty_reads_list(self, tmp_path):
        """Test initialization with empty reads list should fail"""
        with pytest.raises(ValidationError) as excinfo:
            SeqKitReadStats(
                reads=[],
                prefix="test_stats",
                tmpdir=tmp_path
            )
        assert "At least one read file must be provided" in str(excinfo.value)
        
        
    def test_invalid_fastq_encoding(self, tmp_path):
        """Test initialization with invalid FASTQ encoding should fail"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        with pytest.raises(ValidationError) as excinfo:
            SeqKitReadStats(
                reads=[str(read_file)],
                prefix="test_stats",
                tmpdir=tmp_path,
                fastq_encoding="invalid_encoding"
            )
        assert "Invalid FASTQ encoding" in str(excinfo.value)
        
    def test_output_property(self, tmp_path):
        """Test output property returns correct paths"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        stage = SeqKitReadStats(
            reads=[str(read_file)],
            prefix="test_stats",
        )
        
        output = stage.output
        assert str(output.stats_tsv) == "test_stats.reads.tsv"
        assert str(output.raw_tsv) == "test_stats.reads.raw.tsv"
        
    def test_basic_command(self, tmp_path):
        """Test basic seqkit stats command generation"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        stage = SeqKitReadStats(
            reads=[str(read_file)],
            prefix="test_stats",
        )
        ctx = Context(cpus=2)
        commands = stage.create_commands(ctx)
        assert len(commands) == 2
        
        cmd = str(commands[0])
        assert "seqkit stats" in cmd
        assert "-j 2" in cmd
        assert "-T" in cmd  # tabular
        assert "-a" in cmd  # all stats
        assert "-e" in cmd  # skip errors
        assert str(read_file) in cmd
        assert "> test_stats.reads.raw.tsv" in cmd

        postprocess_cmd = str(commands[1])
        assert "add_sample_column" in postprocess_cmd
        assert "test_stats.reads.raw.tsv" in postprocess_cmd
        assert "test_stats.reads.tsv" in postprocess_cmd
        
    def test_command_with_custom_options(self, tmp_path):
        """Test command generation with custom options"""
        read1 = tmp_path / "reads1.fastq"
        read2 = tmp_path / "reads2.fastq"
        read1.touch()
        read2.touch()
        
        stage = SeqKitReadStats(
            reads=[str(read1), str(read2)],
            prefix="custom_stats",
            all_stats=False,
            tabular=False,
            basename_only=True,
            skip_errors=False,
            fastq_encoding="illumina-1.3+",
            gap_letters="N -",
            additional_options="--some-option"
        )
        
        ctx = Context(cpus=4)
        commands = stage.create_commands(ctx)
        cmd = str(commands[0])
        
        assert "seqkit stats" in cmd
        assert "-j 4" in cmd
        assert "-T" not in cmd  # tabular disabled
        assert "-a" not in cmd  # all_stats disabled
        assert "-b" in cmd  # basename_only
        assert "-e" not in cmd  # skip_errors disabled
        assert "-E illumina-1.3+" in cmd
        assert "'N -'" in cmd
        assert "--some-option" in cmd
        assert str(read1) in cmd
        assert str(read2) in cmd


class TestSeqKitReadStatsBasic:
    """Test SeqKitReadStatsBasic stage"""
    
    def test_basic_defaults(self, tmp_path):
        """Test that basic variant has all_stats set to False by default"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        stage = SeqKitReadStatsBasic(
            reads=[str(read_file)],
            prefix="basic_stats"
        )
        
        assert stage.all_stats is False
        ctx = Context()
        commands = stage.create_commands(ctx)
        cmd = str(commands[0])
        assert "seqkit stats" in cmd
        assert "-a" not in cmd  # all_stats disabled


class TestSeqKitReadStatsDetailed:
    """Test SeqKitReadStatsDetailed stage"""
    
    def test_detailed_defaults(self, tmp_path):
        """Test that detailed variant has all_stats enabled and additional N-stats"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        stage = SeqKitReadStatsDetailed(
            reads=[str(read_file)],
            prefix="detailed_stats",
            additional_n_stats=[90, 95]
        )
        
        assert stage.all_stats is True
        assert stage.additional_n_stats == [90, 95]
        
        ctx = Context()
        commands = stage.create_commands(ctx)
        cmd = str(commands[0])
        assert "seqkit stats" in cmd
        assert "-a" in cmd  # all_stats enabled
        assert "-N 90,95" in cmd  # additional N-stats
        
    def test_invalid_n_stats(self, tmp_path):
        """Test that invalid N-statistics raise validation error"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        with pytest.raises(ValidationError) as excinfo:
            SeqKitReadStatsDetailed(
                reads=[str(read_file)],
                prefix="detailed_stats",
                tmpdir=tmp_path,
                additional_n_stats=[150]  # Invalid: > 100
            )
        assert "N-statistic values must be between 0 and 100" in str(excinfo.value)
        
        with pytest.raises(ValidationError) as excinfo:
            SeqKitReadStatsDetailed(
                reads=[str(read_file)],
                prefix="detailed_stats",
                tmpdir=tmp_path,
                additional_n_stats=[-10]  # Invalid: < 0
            )
        assert "N-statistic values must be between 0 and 100" in str(excinfo.value)
        
    def test_empty_n_stats(self, tmp_path):
        """Test that empty N-statistics list is valid"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        stage = SeqKitReadStatsDetailed(
            reads=[str(read_file)],
            prefix="detailed_stats",
            additional_n_stats=[]
        )
        ctx = Context()
        
        commands = stage.create_commands(ctx)
        cmd = commands[0]
        assert "-N" not in cmd  # No N-stats option when empty
