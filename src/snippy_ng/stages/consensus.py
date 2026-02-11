from pathlib import Path
from typing import List

from snippy_ng.exceptions import StageExecutionError
from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.stages.base import BaseStage, BaseOutput
from snippy_ng.dependencies import bcftools

from pydantic import Field

class ConsensusValidationError(StageExecutionError):
    pass

class PseudoAlignment(BaseStage):
    reference: Path = Field(..., description="Reference file")

class BcftoolsPseudoAlignmentOutput(BaseOutput):
    fasta: Path

class BcftoolsPseudoAlignment(PseudoAlignment):
    """
    Call pseudo-alignment using Bcftools consensus.
    """
    vcf_gz: Path = Field(..., description="Input VCF.gz file")
    no_insertions: bool = Field(True, description="Do not apply insertions to the consensus sequence")
    ref_metadata: ReferenceMetadata = Field(..., description="Metadata for the run")

    _dependencies = [
        bcftools
    ]
    
    @property
    def output(self) -> BcftoolsPseudoAlignmentOutput:
        return BcftoolsPseudoAlignmentOutput(
            fasta=Path(f"{self.prefix}.pseudo.raw.fna")
        )
    
    def test_output_matches_reference(self):
        """Asserts that the FASTA file has the expected length and contig count as per the reference metadata."""
        fasta_file = self.output.fasta
        ref_length = self.ref_metadata.total_length
        ref_contigs = self.ref_metadata.num_sequences
        assert ref_length is not None, "Reference length metadata is missing."
        assert ref_contigs is not None, "Reference contig count metadata is missing."
        actual_length = 0
        num_contigs = 0 
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                   num_contigs += 1
                   continue
                actual_length += len(line.strip())
        if actual_length != ref_length:
            raise ConsensusValidationError(f"FASTA file {fasta_file} length mismatch: expected {ref_length}, got {actual_length}")
        if num_contigs != ref_contigs:
            raise ConsensusValidationError(f"FASTA file {fasta_file} contig count mismatch: expected {ref_contigs}, got {num_contigs}")

    @property
    def commands(self) -> List:
        """Constructs the bcftools consensus command."""

        bcf_csq_args = ["bcftools", "consensus"]

        if self.no_insertions:
            # Ensure that only indels where ALT is shorter than REF are applied i.e. no insertions
            bcf_csq_args.extend(["-i", "strlen(ALT)<=strlen(REF)"])
        bcf_csq_args.extend([
            "-f", str(self.reference),
            "-o", str(self.output.fasta),
            "--mark-del", "-",
            str(self.vcf_gz),
        ])
        return [
            self.shell_cmd(["bcftools", "index", str(self.vcf_gz)], description="Indexing VCF file"),
            self.shell_cmd(bcf_csq_args, description="Calling consensus with bcftools"),
            self.shell_cmd(["rm", f"{self.vcf_gz}.csi"], description="Removing VCF index file") 
        ]