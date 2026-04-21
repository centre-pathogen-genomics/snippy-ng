from pathlib import Path
from typing import List

from snippy_ng.exceptions import StageExecutionError
from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.stages import BaseStage, BaseOutput, TempPath
from snippy_ng.dependencies import bcftools

from pydantic import Field

class ConsensusValidationError(StageExecutionError):
    pass

class PseudoAlignment(BaseStage):
    reference: Path = Field(..., description="Reference file")

class BcftoolsPseudoAlignmentOutput(BaseOutput):
    fasta: Path = Field(..., description="Pseudo-alignment consensus FASTA generated from reference + variants")
    vcf_gz: TempPath = Field(..., description="BGZF-compressed VCF file used for consensus calling")
    vcf_index: TempPath = Field(..., description="Index file for the input VCF.gz")

class BcftoolsPseudoAlignment(PseudoAlignment):
    """
    Call pseudo-alignment using Bcftools consensus.
    """
    vcf: Path = Field(..., description="Input VCF file")
    iupac_ambiguity_codes: bool = Field(True, description="Use IUPAC ambiguity codes for heterozygous sites")
    ref_metadata: ReferenceMetadata = Field(..., description="Metadata for the run")

    _dependencies = [
        bcftools
    ]
    
    @property
    def output(self) -> BcftoolsPseudoAlignmentOutput:
        return BcftoolsPseudoAlignmentOutput(
            fasta=Path(f"{self.prefix}.pseudo.raw.fna"),
            vcf_gz=self.vcf.with_suffix(".vcf.gz"),
            vcf_index=self.vcf.with_suffix(".vcf.gz.csi"),
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

    def create_commands(self, ctx) -> List:
        """Constructs the bcftools consensus command."""

        bcf_csq_args = ["bcftools", "consensus"]

        if self.iupac_ambiguity_codes:
            bcf_csq_args.append("--iupac-codes")
        bcf_csq_args.extend([
            "-f", str(self.reference),
            "-o", str(self.output.fasta),
            "--mark-del", "-",
            "--mark-snv", "lc",
            str(self.output.vcf_gz),
        ])
        return [
            self.shell_cmd(["bgzip", "-o", str(self.output.vcf_gz), str(self.vcf)], description="Compressing file with bgzip"),
            self.shell_cmd(["bcftools", "index", "-f", str(self.output.vcf_gz)], description="Indexing VCF file"),
            self.shell_cmd(bcf_csq_args, description="Calling consensus with bcftools"),
        ]
