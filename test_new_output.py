#!/usr/bin/env python3

"""
Test script to verify the new FastpCleanReadsOutput-style structure works correctly
"""

import tempfile
from pathlib import Path
from snippy_ng.stages.downsample_reads import RasusaDownsampleReadsByCoverage
from snippy_ng.stages import at_run_time

def test_new_output_structure():
    """Test the new output structure matches FastpCleanReadsOutput style"""
    
    def get_genome_length():
        return 197394
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        
        # Test paired-end reads
        read1 = tmp_path / "sample_R1.fastq.gz"
        read2 = tmp_path / "sample_R2.fastq.gz"
        read1.touch()
        read2.touch()
        
        print("Testing paired-end reads...")
        stage = RasusaDownsampleReadsByCoverage(
            reads=[str(read1), str(read2)],
            prefix="test",
            genome_length=at_run_time(get_genome_length),
            coverage=30.0,
            tmpdir=tmp_path
        )
        
        output = stage.output
        print(f"R1 output: {output.downsampled_r1}")
        print(f"R2 output: {output.downsampled_r2}")
        
        assert output.downsampled_r1 == "test.downsampled.R1.fastq.gz"
        assert output.downsampled_r2 == "test.downsampled.R2.fastq.gz"
        
        # Test single-end reads
        single_read = tmp_path / "sample.fastq"
        single_read.touch()
        
        print("\nTesting single-end reads...")
        stage_single = RasusaDownsampleReadsByCoverage(
            reads=[str(single_read)],
            prefix="test_single",
            genome_length=197394,
            coverage=50.0,
            output_format="fasta",
            tmpdir=tmp_path
        )
        
        output_single = stage_single.output
        print(f"R1 output: {output_single.downsampled_r1}")
        print(f"R2 output: {output_single.downsampled_r2}")
        
        assert output_single.downsampled_r1 == "test_single.downsampled.R1.fasta"
        assert output_single.downsampled_r2 is None
        
        # Test command generation
        print("\nTesting command generation...")
        command = stage.build_rasusa_command()
        print(f"Command: {command}")
        
        assert "-o test.downsampled.R1.fastq.gz" in command
        assert "-o test.downsampled.R2.fastq.gz" in command
        
        print("\n✅ Success! The new FastpCleanReadsOutput-style structure works correctly.")

if __name__ == "__main__":
    test_new_output_structure()
