"""
Test module for SeqKit read statistics stages
"""

import pytest
from pydantic import ValidationError

from snippy_ng.stages.stats import (
    SeqKitReadStats,
    SeqKitReadStatsBasic,
    SeqKitReadStatsDetailed
)


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
            tmpdir=tmp_path,
            cpus=4
        )
        
        assert stage.reads == [str(read1), str(read2)]
        assert stage.prefix == "test_stats"
        assert stage.cpus == 4
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
            tmpdir=tmp_path
        )
        
        output = stage.output
        assert output.stats_tsv == "test_stats.stats.tsv"
        
    def test_basic_command(self, tmp_path):
        """Test basic seqkit stats command generation"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        stage = SeqKitReadStats(
            reads=[str(read_file)],
            prefix="test_stats",
            tmpdir=tmp_path,
            cpus=2
        )
        
        commands = stage.commands
        assert len(commands) == 1
        
        cmd = str(commands[0])
        assert "seqkit stats" in cmd
        assert "-j 2" in cmd
        assert "-T" in cmd  # tabular
        assert "-a" in cmd  # all stats
        assert "-e" in cmd  # skip errors
        assert str(read_file) in cmd
        assert "> test_stats.stats.tsv" in cmd
        
    def test_command_with_custom_options(self, tmp_path):
        """Test command generation with custom options"""
        read1 = tmp_path / "reads1.fastq"
        read2 = tmp_path / "reads2.fastq"
        read1.touch()
        read2.touch()
        
        stage = SeqKitReadStats(
            reads=[str(read1), str(read2)],
            prefix="custom_stats",
            tmpdir=tmp_path,
            cpus=4,
            all_stats=False,
            tabular=False,
            basename_only=True,
            skip_errors=False,
            fastq_encoding="illumina-1.3+",
            gap_letters="'N -'",
            additional_options="--some-option"
        )
        
        commands = stage.commands
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
            prefix="basic_stats",
            tmpdir=tmp_path
        )
        
        assert stage.all_stats is False
        
        commands = stage.commands
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
            tmpdir=tmp_path,
            additional_n_stats=[90, 95]
        )
        
        assert stage.all_stats is True
        assert stage.additional_n_stats == [90, 95]
        
        commands = stage.commands
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
            tmpdir=tmp_path,
            additional_n_stats=[]
        )
        
        commands = stage.commands
        cmd = commands[0]
        assert "-N" not in cmd  # No N-stats option when empty
