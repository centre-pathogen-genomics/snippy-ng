"""
Test module for Rasusa read downsampling stages
"""

import pytest
from pydantic import ValidationError

from snippy_ng.stages.downsample_reads import (
    RasusaDownsampleReads,
    RasusaDownsampleReadsByCoverage,
    RasusaDownsampleReadsByCount
)


class TestRasusaDownsampleReads:
    """Test RasusaDownsampleReads stage"""
    
    def test_init_valid_inputs_coverage(self, tmp_path):
        """Test initialization with valid inputs using coverage"""
        # Create test files
        read1 = tmp_path / "reads1.fastq.gz"
        read2 = tmp_path / "reads2.fastq.gz"
        read1.touch()
        read2.touch()
        
        stage = RasusaDownsampleReads(
            reads=[str(read1), str(read2)],
            prefix="downsampled",
            genome_length=197394,
            coverage=50.0,
            tmpdir=tmp_path,
            cpus=4
        )
        
        assert stage.reads == [str(read1), str(read2)]
        assert stage.prefix == "downsampled"
        assert stage.coverage == 50.0
        assert stage.num_reads is None
        assert stage.genome_length == 197394
        assert stage.output_format == "fastq"
        assert stage.compression_level == 6
    
    def test_init_valid_inputs_num_reads(self, tmp_path):
        """Test initialization with valid inputs using num_reads"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        stage = RasusaDownsampleReads(
            reads=[str(read_file)],
            prefix="downsampled",
            num_reads=1000000,
            tmpdir=tmp_path
        )
        
        assert stage.num_reads == 1000000
        assert stage.coverage is None
        assert stage.genome_length is None
    
    def test_init_empty_reads_list(self, tmp_path):
        """Test initialization with empty reads list should fail"""
        with pytest.raises(ValidationError) as excinfo:
            RasusaDownsampleReads(
                reads=[],
                prefix="downsampled",
                genome_length=197394,
                coverage=50.0,
                tmpdir=tmp_path
            )
        assert "At least one read file must be provided" in str(excinfo.value)
    
    def test_coverage_without_genome_length_fails(self, tmp_path):
        """Test initialization with coverage but no genome length should fail"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        with pytest.raises(ValidationError) as excinfo:
            RasusaDownsampleReads(
                reads=[str(read_file)],
                prefix="downsampled",
                coverage=50.0,
                tmpdir=tmp_path
            )
        assert "genome_length is required when using coverage-based downsampling" in str(excinfo.value)
    
    def test_both_coverage_and_num_reads_fails(self, tmp_path):
        """Test initialization with both coverage and num_reads should fail"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        with pytest.raises(ValidationError) as excinfo:
            RasusaDownsampleReads(
                reads=[str(read_file)],
                prefix="downsampled",
                genome_length=197394,
                coverage=50.0,
                num_reads=1000000,
                tmpdir=tmp_path
            )
        assert "Cannot specify both coverage and num_reads" in str(excinfo.value)
    
    def test_neither_coverage_nor_num_reads_fails(self, tmp_path):
        """Test initialization with neither coverage nor num_reads should fail"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        with pytest.raises(ValidationError) as excinfo:
            RasusaDownsampleReads(
                reads=[str(read_file)],
                prefix="downsampled",
                tmpdir=tmp_path
            )
        assert "Must specify either coverage or num_reads" in str(excinfo.value)
    
    def test_invalid_output_format(self, tmp_path):
        """Test initialization with invalid output format should fail"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        with pytest.raises(ValidationError) as excinfo:
            RasusaDownsampleReads(
                reads=[str(read_file)],
                prefix="downsampled",
                genome_length=197394,
                coverage=50.0,
                output_format="invalid",
                tmpdir=tmp_path
            )
        assert "Invalid output format" in str(excinfo.value)
    
    def test_invalid_compression_level(self, tmp_path):
        """Test initialization with invalid compression level should fail"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        with pytest.raises(ValidationError) as excinfo:
            RasusaDownsampleReads(
                reads=[str(read_file)],
                prefix="downsampled",
                genome_length=197394,
                coverage=50.0,
                compression_level=10,  # Invalid: > 9
                tmpdir=tmp_path
            )
        assert "Compression level must be between 1 and 9" in str(excinfo.value)
    
    def test_output_property_paired_reads(self, tmp_path):
        """Test output property with paired reads"""
        read1 = tmp_path / "sample_R1.fastq.gz"
        read2 = tmp_path / "sample_R2.fastq.gz"
        read1.touch()
        read2.touch()
        
        stage = RasusaDownsampleReads(
            reads=[str(read1), str(read2)],
            prefix="downsampled",
            genome_length=197394,
            coverage=50.0,
            tmpdir=tmp_path
        )
        
        output = stage.output
        assert output.downsampled_r1 == "downsampled.downsampled.R1.fastq.gz"
        assert output.downsampled_r2 == "downsampled.downsampled.R2.fastq.gz"
    
    def test_output_property_single_read_uncompressed(self, tmp_path):
        """Test output property with single uncompressed read"""
        read_file = tmp_path / "sample.fastq"
        read_file.touch()
        
        stage = RasusaDownsampleReads(
            reads=[str(read_file)],
            prefix="ds",
            num_reads=1000000,
            output_format="fasta",
            tmpdir=tmp_path
        )
        
        output = stage.output
        assert output.downsampled_r1 == "ds.downsampled.R1.fasta"
        assert output.downsampled_r2 is None
    
    def test_basic_command_coverage(self, tmp_path):
        """Test basic rasusa command generation with coverage"""
        read1 = tmp_path / "sample_R1.fastq.gz"
        read2 = tmp_path / "sample_R2.fastq.gz"
        read1.touch()
        read2.touch()
        
        stage = RasusaDownsampleReads(
            reads=[str(read1), str(read2)],
            prefix="ds",
            genome_length=197394,
            coverage=50.0,
            seed=42,
            tmpdir=tmp_path
        )
        
        commands = stage.commands
        assert len(commands) == 1
        
        cmd = str(commands[0])
        assert cmd.startswith("rasusa reads")
        assert "--coverage 50.0" in cmd
        assert "--genome-size 197394" in cmd
        assert "-o ds.downsampled.R1.fastq.gz" in cmd
        assert "-o ds.downsampled.R2.fastq.gz" in cmd
        assert "--seed 42" in cmd
        assert "--compress-level 6" in cmd
        # Input files should be at the end
        assert cmd.endswith(f"{read1} {read2}")
    
    def test_basic_command_num_reads(self, tmp_path):
        """Test basic rasusa command generation with num_reads"""
        read_file = tmp_path / "sample.fastq"
        read_file.touch()
        
        stage = RasusaDownsampleReads(
            reads=[str(read_file)],
            prefix="ds",
            genome_length=197394,  # Still need genome_length for the class
            num_reads=1000000,
            output_format="fasta",
            tmpdir=tmp_path
        )
        
        commands = stage.commands
        cmd = str(commands[0])
        
        assert cmd.startswith("rasusa reads")
        assert "--num 1000000" in cmd
        assert "--fasta" in cmd
        assert "-o ds.downsampled.R1.fasta" in cmd
        assert "--genome-size" not in cmd  # Not used with num_reads
        assert cmd.endswith(str(read_file))
    
    def test_command_with_custom_options(self, tmp_path):
        """Test command generation with custom options"""
        read_file = tmp_path / "sample.fastq.gz"
        read_file.touch()
        
        stage = RasusaDownsampleReads(
            reads=[str(read_file)],
            prefix="custom",
            genome_length=150000,
            coverage=75.0,
            seed=123,
            compression_level=9,
            additional_options="--verbose",
            tmpdir=tmp_path
        )
        
        commands = stage.commands
        cmd = str(commands[0])
        
        assert "--coverage 75.0" in cmd
        assert "--genome-size 150000" in cmd
        assert "--seed 123" in cmd
        assert "--compress-level 9" in cmd
        assert "--verbose" in cmd
    
    def test_at_run_time_genome_length(self, tmp_path):
        """Test that at_run_time works correctly with genome_length"""
        from snippy_ng.stages import at_run_time
        
        # Create test files
        read1 = tmp_path / "reads1.fastq.gz"
        read2 = tmp_path / "reads2.fastq.gz"
        read1.touch()
        read2.touch()
        
        # Track function calls
        call_count = 0
        def get_genome_length():
            nonlocal call_count
            call_count += 1
            return 197394
        
        # Create stage with at_run_time genome_length
        stage = RasusaDownsampleReads(
            reads=[str(read1), str(read2)],
            prefix="test",
            genome_length=at_run_time(get_genome_length),
            coverage=50.0,
            tmpdir=tmp_path
        )
        
        # Function should not be called during initialization
        assert call_count == 0
        
        # Function should be called when building command
        command = str(stage.build_rasusa_command())
        assert call_count == 1
        assert "--genome-size 197394" in command
        
        # Function should not be called again (cached)
        command2 = str(stage.build_rasusa_command())
        assert call_count == 1
        assert command == command2

class TestRasusaDownsampleReadsByCoverage:
    """Test RasusaDownsampleReadsByCoverage stage"""
    
    def test_coverage_required(self, tmp_path):
        """Test that coverage is required in coverage variant"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        stage = RasusaDownsampleReadsByCoverage(
            reads=[str(read_file)],
            prefix="cov_ds",
            genome_length=150000,
            coverage=30.0,
            tmpdir=tmp_path
        )
        
        assert stage.coverage == 30.0
        assert stage.genome_length == 150000
        assert stage.num_reads is None
    
    def test_num_reads_disabled(self, tmp_path):
        """Test that num_reads is disabled in coverage variant"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        with pytest.raises(ValidationError) as excinfo:
            RasusaDownsampleReadsByCoverage(
                reads=[str(read_file)],
                prefix="cov_ds",
                genome_length=150000,
                coverage=30.0,
                num_reads=1000000,  # Should fail
                tmpdir=tmp_path
            )
        assert "Cannot specify num_reads in coverage-based downsampling" in str(excinfo.value)
    
    def test_at_run_time_coverage_variant(self, tmp_path):
        """Test that at_run_time works with RasusaDownsampleReadsByCoverage"""
        from snippy_ng.stages import at_run_time
        
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        # Track function calls
        call_count = 0
        def get_genome_length():
            nonlocal call_count
            call_count += 1
            return 150000
        
        stage = RasusaDownsampleReadsByCoverage(
            reads=[str(read_file)],
            prefix="cov_ds",
            genome_length=at_run_time(get_genome_length),
            coverage=30.0,
            tmpdir=tmp_path
        )
        
        # Function should not be called during initialization
        assert call_count == 0
        
        # Function should be called when building command
        command = str(stage.build_rasusa_command())
        assert call_count == 1
        assert "--genome-size 150000" in command
        assert "--coverage 30.0" in command

class TestRasusaDownsampleReadsByCount:
    """Test RasusaDownsampleReadsByCount stage"""
    
    def test_num_reads_required(self, tmp_path):
        """Test that num_reads is required in count variant"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        stage = RasusaDownsampleReadsByCount(
            reads=[str(read_file)],
            prefix="count_ds",
            num_reads=500000,
            tmpdir=tmp_path
        )
        
        assert stage.num_reads == 500000
        assert stage.coverage is None
        assert stage.genome_length is None
    
    def test_coverage_disabled(self, tmp_path):
        """Test that coverage is disabled in count variant"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        with pytest.raises(ValidationError) as excinfo:
            RasusaDownsampleReadsByCount(
                reads=[str(read_file)],
                prefix="count_ds",
                num_reads=500000,
                coverage=25.0,  # Should fail
                tmpdir=tmp_path
            )
        assert "Cannot specify coverage in count-based downsampling" in str(excinfo.value)
