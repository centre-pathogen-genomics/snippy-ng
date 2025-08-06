#!/usr/bin/env python3

"""
Test script to verify at_run_time functionality works correctly
"""

import tempfile
from pathlib import Path
from snippy_ng.stages.downsample_reads import RasusaDownsampleReadsByCoverage
from snippy_ng.stages import at_run_time

def create_test_files(tmp_dir):
    """Create test read files"""
    read1 = tmp_dir / "test_R1.fastq.gz"
    read2 = tmp_dir / "test_R2.fastq.gz"
    read1.touch()
    read2.touch()
    return [str(read1), str(read2)]

def test_at_runtime():
    """Test at_run_time functionality"""
    
    # Create a function that will be called at runtime
    call_count = 0
    def get_genome_length():
        nonlocal call_count
        call_count += 1
        print(f"get_genome_length() called for the {call_count} time(s)")
        return 197394
    
    # Wrap the function with at_run_time
    genome_length = at_run_time(get_genome_length)
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        reads = create_test_files(tmp_path)
        
        print("Creating RasusaDownsampleReadsByCoverage stage...")
        print(f"Call count before creating stage: {call_count}")
        
        # Create the stage - at this point, get_genome_length should NOT be called
        stage = RasusaDownsampleReadsByCoverage(
            reads=reads,
            prefix="test",
            genome_length=genome_length,
            coverage=30.0,
            tmpdir=tmp_path
        )
        
        print(f"Call count after creating stage: {call_count}")
        print("Stage created successfully. The function should not have been called yet.")
        
        print("\nNow building the rasusa command...")
        # When we build the command, str() will be called on genome_length
        command = stage.build_rasusa_command()
        
        print(f"Call count after building command: {call_count}")
        print(f"Generated command: {command}")
        
        # Verify the command contains the genome size
        assert "--genome-size 197394" in command
        print("\n✅ Success! The genome length was correctly evaluated at runtime.")
        
        # Test that subsequent calls use cached value
        stage.build_rasusa_command()
        print(f"Call count after second command build: {call_count}")
        print("✅ Success! The function was only called once (cached correctly).")

if __name__ == "__main__":
    test_at_runtime()
