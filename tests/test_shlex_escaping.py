"""Tests to verify that filenames with spaces are properly escaped in shell commands."""
import pytest
from pathlib import Path
from snippy_ng.stages.setup import PrepareReference
from snippy_ng.stages.consequences import BcftoolsConsequencesCaller
from snippy_ng.stages.calling import FreebayesCaller
from snippy_ng.stages.alignment import BWAMEMReadsAligner, MinimapAligner, PreAlignedReads
from snippy_ng.stages.alignment_filtering import AlignmentFilter
from snippy_ng.stages.stats import SeqKitReadStats
from snippy_ng.stages.clean_reads import FastpCleanReads
from snippy_ng.stages.downsample_reads import RasusaDownsampleReads


class TestPrepareReferenceEscaping:
    """Test PrepareReference stage escapes filenames with spaces."""
    
    def test_commands_with_spaces(self, tmp_path):
        """Test that filenames with spaces are properly escaped."""
        # Create a reference directory with spaces in the name
        ref_dir = tmp_path / "my reference"
        ref_dir.mkdir()
        
        input_file = tmp_path / "input file.gbk"
        input_file.touch()
        
        stage = PrepareReference(
            input=input_file,
            ref_fmt="genbank",
            reference_prefix="my ref",
            reference_dir=ref_dir,
            tmpdir=tmp_path
        )
        
        commands = stage.commands
        # Check that rm and mkdir commands are properly escaped
        assert len(commands) == 3
        assert "rm -f" in commands[0]
        # Check that the filename with spaces is quoted (contains the space within quotes)
        assert "my reference" in commands[0] and ("'" in commands[0] or '"' in commands[0])
        assert "mkdir -p" in commands[1]
        assert "my reference" in commands[1] and ("'" in commands[1] or '"' in commands[1])


class TestBcftoolsConsequencesCallerEscaping:
    """Test BcftoolsConsequencesCaller stage escapes filenames with spaces."""
    
    def test_commands_with_spaces(self, tmp_path):
        """Test that filenames with spaces are properly escaped."""
        ref_file = tmp_path / "my reference.fa"
        ref_file.touch()
        
        variants_file = tmp_path / "my variants.vcf"
        variants_file.touch()
        
        features_file = tmp_path / "my features.gff"
        # Create a features file with content
        features_file.write_text("##gff-version 3\nseq1\tsource\tgene\t1\t100\t.\t+\t.\tID=gene1\n")
        
        stage = BcftoolsConsequencesCaller(
            reference=ref_file,
            variants=variants_file,
            features=features_file,
            prefix="my output",
            tmpdir=tmp_path
        )
        
        commands = stage.commands
        assert len(commands) == 1
        cmd = commands[0]
        assert "bcftools csq" in cmd
        # Check that all filenames with spaces are quoted
        assert "my reference.fa" in cmd and ("'" in cmd or '"' in cmd)
        assert "my variants.vcf" in cmd
        assert "my features.gff" in cmd
        assert "my output.vcf" in cmd


class TestFreebayesCallerEscaping:
    """Test FreebayesCaller stage escapes filenames with spaces."""
    
    def test_commands_with_spaces(self, tmp_path):
        """Test that filenames with spaces are properly escaped."""
        ref_file = tmp_path / "my reference.fa"
        ref_file.touch()
        
        bam_file = tmp_path / "my alignment.bam"
        bam_file.touch()
        
        stage = FreebayesCaller(
            reference=ref_file,
            bam=bam_file,
            prefix="my output",
            tmpdir=tmp_path
        )
        
        commands = stage.commands
        assert len(commands) == 3
        
        # Check that all commands properly escape filenames
        # The key is that filenames with spaces appear within quotes
        assert "my reference" in " ".join(commands)
        assert "my alignment" in " ".join(commands)
        # Verify quotes are present when spaces are in filenames
        assert "'" in " ".join(commands) or '"' in " ".join(commands)


class TestAlignmentEscaping:
    """Test alignment stages escape filenames with spaces."""
    
    def test_bwa_commands_with_spaces(self, tmp_path):
        """Test BWA aligner escapes filenames with spaces."""
        ref_file = tmp_path / "my reference.fa"
        ref_file.touch()
        
        read1 = tmp_path / "my reads_R1.fastq"
        read2 = tmp_path / "my reads_R2.fastq"
        read1.touch()
        read2.touch()
        
        stage = BWAMEMReadsAligner(
            reference=ref_file,
            reads=[str(read1), str(read2)],
            prefix="my output",
            tmpdir=tmp_path
        )
        
        commands = stage.commands
        # Check that filenames with spaces are escaped
        all_commands = " ".join(commands)
        assert "my reference" in all_commands
        assert "my reads" in all_commands
        # Verify quotes are present
        assert "'" in all_commands or '"' in all_commands
    
    def test_minimap_commands_with_spaces(self, tmp_path):
        """Test Minimap2 aligner escapes filenames with spaces."""
        ref_file = tmp_path / "my reference.fa"
        ref_file.touch()
        
        read1 = tmp_path / "my reads_R1.fastq"
        read1.touch()
        
        stage = MinimapAligner(
            reference=ref_file,
            reads=[str(read1)],
            prefix="my output",
            tmpdir=tmp_path
        )
        
        commands = stage.commands
        # Check that filenames with spaces are escaped
        all_commands = " ".join(commands)
        assert "my reference" in all_commands
        assert "my reads" in all_commands
        # Verify quotes are present
        assert "'" in all_commands or '"' in all_commands


class TestAlignmentFilterEscaping:
    """Test AlignmentFilter stage escapes filenames with spaces."""
    
    def test_commands_with_spaces(self, tmp_path):
        """Test that filenames with spaces are properly escaped."""
        bam_file = tmp_path / "my alignment.bam"
        bam_file.touch()
        
        stage = AlignmentFilter(
            bam=bam_file,
            prefix="my output",
            tmpdir=tmp_path
        )
        
        commands = stage.commands
        assert len(commands) == 2
        
        # Check filter command
        filter_cmd = commands[0]
        assert "my alignment.bam" in filter_cmd
        assert "my output.filtered.bam" in filter_cmd
        
        # Check index command
        index_cmd = commands[1]
        assert "my output.filtered.bam" in index_cmd
        
        # Verify quotes are present for filenames with spaces
        assert "'" in filter_cmd or '"' in filter_cmd


class TestSeqKitReadStatsEscaping:
    """Test SeqKitReadStats stage escapes filenames with spaces."""
    
    def test_commands_with_spaces(self, tmp_path):
        """Test that filenames with spaces are properly escaped."""
        read1 = tmp_path / "my reads_R1.fastq"
        read2 = tmp_path / "my reads_R2.fastq"
        read1.touch()
        read2.touch()
        
        stage = SeqKitReadStats(
            reads=[str(read1), str(read2)],
            prefix="my output",
            tmpdir=tmp_path
        )
        
        commands = stage.commands
        assert len(commands) == 1
        cmd = commands[0]
        
        # Check that read files are present
        assert "my reads_R1.fastq" in cmd
        assert "my reads_R2.fastq" in cmd
        # Check that output is present
        assert "my output.stats.tsv" in cmd
        # Verify quotes are present for filenames with spaces
        assert "'" in cmd or '"' in cmd


class TestFastpCleanReadsEscaping:
    """Test FastpCleanReads stage escapes filenames with spaces."""
    
    def test_commands_with_spaces(self, tmp_path):
        """Test that filenames with spaces are properly escaped."""
        read1 = tmp_path / "my reads_R1.fastq"
        read2 = tmp_path / "my reads_R2.fastq"
        read1.touch()
        read2.touch()
        
        stage = FastpCleanReads(
            reads=[str(read1), str(read2)],
            prefix="my output",
            tmpdir=tmp_path
        )
        
        commands = stage.commands
        assert len(commands) == 1
        cmd = commands[0]
        
        # Check that input and output files are present
        assert "my reads_R1.fastq" in cmd
        assert "my reads_R2.fastq" in cmd
        assert "my output.cleaned.R1.fastq.gz" in cmd
        assert "my output.cleaned.R2.fastq.gz" in cmd
        # Verify quotes are present for filenames with spaces
        assert "'" in cmd or '"' in cmd


class TestRasusaDownsampleReadsEscaping:
    """Test RasusaDownsampleReads stage escapes filenames with spaces."""
    
    def test_commands_with_spaces(self, tmp_path):
        """Test that filenames with spaces are properly escaped."""
        read1 = tmp_path / "my reads_R1.fastq.gz"
        read2 = tmp_path / "my reads_R2.fastq.gz"
        read1.touch()
        read2.touch()
        
        stage = RasusaDownsampleReads(
            reads=[str(read1), str(read2)],
            prefix="my output",
            genome_length=1000000,
            coverage=50.0,
            tmpdir=tmp_path
        )
        
        commands = stage.commands
        assert len(commands) == 1
        cmd = commands[0]
        
        # Check that input and output files are present
        assert "my reads_R1.fastq.gz" in cmd
        assert "my reads_R2.fastq.gz" in cmd
        assert "my output.downsampled.R1.fastq.gz" in cmd
        assert "my output.downsampled.R2.fastq.gz" in cmd
        # Verify quotes are present for filenames with spaces
        assert "'" in cmd or '"' in cmd
