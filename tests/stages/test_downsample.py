"""
Test module for Rasusa read downsampling stages
"""

import pytest
from pydantic import ValidationError

from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.stages import Context
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
            ref_metadata=ReferenceMetadata(total_length=197394),
            reads=[str(read1), str(read2)],
            prefix="downsampled",
            coverage=50.0,
        )
        
        assert stage.reads == [read1, read2]
        assert stage.prefix == "downsampled"
        assert stage.coverage == 50.0
        assert stage.num_reads is None
        assert stage.ref_metadata.total_length == 197394
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
        )
        
        assert stage.num_reads == 1000000
        assert stage.coverage is None
    
    def test_init_empty_reads_list(self, tmp_path):
        """Test initialization with empty reads list should fail"""
        with pytest.raises(ValidationError) as excinfo:
            RasusaDownsampleReads(
                metadata=ReferenceMetadata(total_length=197394),
                reads=[],
                prefix="downsampled",
                coverage=50.0,
                tmpdir=tmp_path
            )
        assert "At least one read file must be provided" in str(excinfo.value)
    
    def test_coverage_without_metadata_fails(self, tmp_path):
        """Test initialization with coverage but no metadata should fail"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        with pytest.raises(ValidationError) as excinfo:
            RasusaDownsampleReads(
                reads=[str(read_file)],
                prefix="downsampled",
                coverage=50.0,
            )
        assert "is required when using coverage-based downsampling" in str(excinfo.value)
    
    def test_both_coverage_and_num_reads_fails(self, tmp_path):
        """Test initialization with both coverage and num_reads should fail"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        with pytest.raises(ValidationError) as excinfo:
            RasusaDownsampleReads(
                reads=[str(read_file)],
                prefix="downsampled",
                ref_metadata=ReferenceMetadata(total_length=197394),
                coverage=50.0,
                num_reads=1000000,
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
                metadata=ReferenceMetadata(total_length=197394),
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
                metadata=ReferenceMetadata(total_length=197394),
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
            ref_metadata=ReferenceMetadata(total_length=197394),
            coverage=50.0,
        )
        
        output = stage.output
        assert str(output.downsampled_r1) == "downsampled.downsampled.R1.fastq.gz"
        assert str(output.downsampled_r2) == "downsampled.downsampled.R2.fastq.gz"
    
    def test_output_property_single_read_uncompressed(self, tmp_path):
        """Test output property with single uncompressed read"""
        read_file = tmp_path / "sample.fastq"
        read_file.touch()
        
        stage = RasusaDownsampleReads(
            reads=[str(read_file)],
            prefix="ds",
            num_reads=1000000,
            output_format="fasta",
        )
        
        output = stage.output
        assert str(output.downsampled_r1) == "ds.downsampled.R1.fasta"
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
            ref_metadata=ReferenceMetadata(total_length=197394),
            coverage=50.0,
            seed=42,
        )
        ctx = Context()
        commands = stage.create_commands(ctx)
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
            ref_metadata=ReferenceMetadata(total_length=197394),
            num_reads=1000000,
            output_format="fasta",
        )
        ctx = Context()
        commands = stage.create_commands(ctx)
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
            ref_metadata=ReferenceMetadata(total_length=150000),
            coverage=75.0,
            seed=123,
            compression_level=9,
            additional_options="--verbose",
        )
        ctx = Context()
        commands = stage.create_commands(ctx)
        cmd = str(commands[0])
        
        assert "--coverage 75.0" in cmd
        assert "--genome-size 150000" in cmd
        assert "--seed 123" in cmd
        assert "--compress-level 9" in cmd
        assert "--verbose" in cmd
    
class TestRasusaDownsampleReadsByCoverage:
    """Test RasusaDownsampleReadsByCoverage stage"""
    
    def test_coverage_required(self, tmp_path):
        """Test that coverage is required in coverage variant"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        stage = RasusaDownsampleReadsByCoverage(
            reads=[str(read_file)],
            prefix="cov_ds",
            ref_metadata=ReferenceMetadata(total_length=150000),
            coverage=30.0,
        )
        
        assert stage.coverage == 30.0
        assert stage.ref_metadata.total_length == 150000
        assert stage.num_reads is None
    
    def test_num_reads_disabled(self, tmp_path):
        """Test that num_reads is disabled in coverage variant"""
        read_file = tmp_path / "reads.fastq"
        read_file.touch()
        
        with pytest.raises(ValidationError) as excinfo:
            RasusaDownsampleReadsByCoverage(
                reads=[str(read_file)],
                prefix="cov_ds",
                metadata=ReferenceMetadata(total_length=150000),
                coverage=30.0,
                num_reads=1000000,  # Should fail
                tmpdir=tmp_path
            )
        assert "Cannot specify num_reads in coverage-based downsampling" in str(excinfo.value)
    
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
        )
        
        assert stage.num_reads == 500000
        assert stage.coverage is None
        assert stage.ref_metadata is None
    
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
