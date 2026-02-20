from snippy_ng.metadata import ReferenceMetadata
from snippy_ng.stages import BaseStage, BaseOutput
from snippy_ng.__about__ import __version__
from pydantic import Field, computed_field, field_validator
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

from snippy_ng.dependencies import biopython
from snippy_ng.exceptions import StageExecutionError


METADATA_FILE_NAME = Path("metadata.json")


class ReferenceValidationError(StageExecutionError):
    pass


class ReferenceOutput(BaseOutput):
    reference: Path
    reference_index: Path
    reference_dict: Path
    gff: Path
    metadata: Path


class PrepareReference(BaseStage):
    input: Path = Field(..., description="Reference file")
    ref_fmt: str = Field("genbank", description="Reference format")
    directory: Path = Field(Path("reference"), description="Output directory for reference files")

    _dependencies = [
        biopython,
    ]

    @computed_field
    @property
    def reference_prefix(self) -> str:
        """Derives the reference prefix from the input filename stem."""
        return self.input.stem

    @property
    def output(self) -> ReferenceOutput:
        return ReferenceOutput(
            reference=self.directory / f"{self.reference_prefix}.fa",
            reference_index=self.directory / f"{self.reference_prefix}.fa.fai",
            reference_dict=self.directory / f"{self.reference_prefix}.dict",
            gff=self.directory / f"{self.reference_prefix}.gff",
            metadata=self.directory / METADATA_FILE_NAME,
        )

    @property
    def commands(self):
        process_reference_cmd = self.python_cmd(
            func=self.process_reference,
            args=(self.input, self.ref_fmt, self.output.reference, self.output.gff),
            description=f"Extract FASTA and GFF from reference ({self.ref_fmt})",
        )
        return [
            self.shell_cmd(
                ["rm", "-f", str(self.output.reference)],
                description=f"Remove existing reference FASTA: {self.output.reference}",
            ),
            self.shell_cmd(
                ["mkdir", "-p", str(self.directory)],
                description=f"Create reference directory: {self.directory}",
            ),
            process_reference_cmd,
            self.shell_cmd(
                ["samtools", "faidx", str(self.output.reference)],
                description=f"Index reference FASTA with samtools faidx: {self.output.reference}",
            ),
            self.shell_pipe(
                commands=[
                    self.shell_cmd(
                        ["cut", "-f1,2", str(self.output.reference_index)],
                        description="Extract sequence names and lengths from FASTA index",
                    ),
                    self.shell_cmd(
                        ["sort"], description="Sort sequence names and lengths"
                    ),
                ],
                output_file=self.output.reference_dict,
                description=f"Create reference dictionary: {self.output.reference_dict}",
            ),
        ]

    def process_reference(
        self,
        reference_path: Path,
        ref_fmt: str,
        output_fasta_path: Path,
        output_gff_path: Path,
    ):
        """
        Extracts FASTA and GFF3 from a reference file.
        Determines input format and writes Ensembl-style GFF3 only if features exist.

        Args:
            reference_path (Path): Path to the reference file.
            ref_fmt (str): Input format (e.g., 'genbank', 'embl').
            output_fasta_path (Path): Path to save the extracted FASTA file.
            output_gff_path (Path): Path to save the extracted GFF3 file.
        """
        import gzip

        try:
            # Open gzipped or plain text reference
            open_func = open
            try:
                with open(reference_path, "rt") as test_fh:
                    test_fh.read(1)
            except UnicodeDecodeError:
                open_func = gzip.open
            with open_func(reference_path, "rt") as ref_fh:
                seq_records = list(SeqIO.parse(ref_fh, ref_fmt))
        except Exception as e:
            raise ValueError(
                f"Failed to parse {reference_path} with format {ref_fmt}: {e}"
            )

        # Prepare outputs
        ref_seq_dict = {}
        gene_counter = 0
        nseq = 0
        nfeat = 0
        total_length = 0

        with open(output_fasta_path, "w") as fasta_out, open(
            output_gff_path, "w"
        ) as gff_out:
            # Write GFF3 header
            gff_out.write("##gff-version 3\n")

            for seq_record in seq_records:
                # Check for duplicate sequences
                if seq_record.id in ref_seq_dict:
                    raise ValueError(
                        f"Duplicate sequence {seq_record.id} in {reference_path}"
                    )

                # Clean sequence: uppercase and replace non-standard bases with 'N'
                dna = Seq(str(seq_record.seq).upper().replace("U", "T"))
                dna = Seq("".join([base if base in "AGTCN" else "N" for base in dna]))
                seq_record.seq = dna
                ref_seq_dict[seq_record.id] = dna

                # Write to FASTA
                SeqIO.write(seq_record, fasta_out, "fasta")
                nseq += 1
                total_length += len(dna)

                # Write sequence region directive for GFF3
                gff_out.write(f"##sequence-region {seq_record.id} 1 {len(dna)}\n")

                # Group features by gene/transcript hierarchy
                genes = {}
                standalone_features = []

                for feature in seq_record.features:
                    ftype = feature.type
                    if ftype in ("source", "misc_feature"):
                        continue  # Skip unwanted features

                    # Determine gene and transcript IDs
                    gene_id = None
                    transcript_id = None

                    if "locus_tag" in feature.qualifiers:
                        gene_id = feature.qualifiers["locus_tag"][0]
                    elif "gene" in feature.qualifiers:
                        gene_id = feature.qualifiers["gene"][0]
                    else:
                        gene_counter += 1
                        gene_id = f"gene_{gene_counter}"

                    # For features that need transcripts (CDS, exon, UTRs)
                    if ftype in ("CDS", "exon", "five_prime_UTR", "three_prime_UTR"):
                        transcript_id = f"{gene_id}_transcript_1"

                        if gene_id not in genes:
                            genes[gene_id] = {
                                "feature": feature,
                                "transcripts": {},
                                "start": feature.location.start
                                + 1,  # Convert to 1-based
                                "end": feature.location.end,
                                "strand": "+"
                                if feature.location.strand == 1
                                else "-"
                                if feature.location.strand == -1
                                else ".",
                            }
                        else:
                            # Extend gene boundaries
                            genes[gene_id]["start"] = min(
                                genes[gene_id]["start"], feature.location.start + 1
                            )
                            genes[gene_id]["end"] = max(
                                genes[gene_id]["end"], feature.location.end
                            )

                        if transcript_id not in genes[gene_id]["transcripts"]:
                            genes[gene_id]["transcripts"][transcript_id] = {
                                "features": [],
                                "start": feature.location.start + 1,
                                "end": feature.location.end,
                                "strand": "+"
                                if feature.location.strand == 1
                                else "-"
                                if feature.location.strand == -1
                                else ".",
                            }
                        else:
                            # Extend transcript boundaries
                            transcript = genes[gene_id]["transcripts"][transcript_id]
                            transcript["start"] = min(
                                transcript["start"], feature.location.start + 1
                            )
                            transcript["end"] = max(
                                transcript["end"], feature.location.end
                            )

                        genes[gene_id]["transcripts"][transcript_id]["features"].append(
                            feature
                        )

                    elif ftype == "gene":
                        if gene_id not in genes:
                            genes[gene_id] = {
                                "feature": feature,
                                "transcripts": {},
                                "start": feature.location.start + 1,
                                "end": feature.location.end,
                                "strand": "+"
                                if feature.location.strand == 1
                                else "-"
                                if feature.location.strand == -1
                                else ".",
                            }
                    else:
                        # Standalone features (tRNA, rRNA, etc.)
                        standalone_features.append(feature)

                # Write genes and their transcripts
                for gene_id, gene_data in genes.items():
                    gene_feature = gene_data["feature"]

                    # Determine biotype
                    biotype = "protein_coding"  # Default
                    if gene_feature.type == "tRNA":
                        biotype = "tRNA"
                    elif gene_feature.type == "rRNA":
                        biotype = "rRNA"
                    elif gene_feature.type in ("ncRNA", "misc_RNA"):
                        biotype = "misc_RNA"

                    # Get gene name
                    gene_name = ""
                    if "gene" in gene_feature.qualifiers:
                        gene_name = f";Name={gene_feature.qualifiers['gene'][0]}"
                    elif "product" in gene_feature.qualifiers:
                        gene_name = f";Name={gene_feature.qualifiers['product'][0]}"

                    # Write gene line
                    gff_out.write(
                        f"{seq_record.id}\tsnipy-ng\tgene\t{gene_data['start']}\t{gene_data['end']}\t.\t{gene_data['strand']}\t.\tID=gene:{gene_id};biotype={biotype}{gene_name}\n"
                    )
                    nfeat += 1

                    # Write transcripts and their features
                    for transcript_id, transcript_data in gene_data[
                        "transcripts"
                    ].items():
                        # Write transcript line
                        gff_out.write(
                            f"{seq_record.id}\tsnipy-ng\ttranscript\t{transcript_data['start']}\t{transcript_data['end']}\t.\t{transcript_data['strand']}\t.\tID=transcript:{transcript_id};Parent=gene:{gene_id};biotype={biotype}\n"
                        )
                        nfeat += 1

                        # Write transcript features (CDS, exon, UTRs)
                        for feature in transcript_data["features"]:
                            start = feature.location.start + 1  # Convert to 1-based
                            end = feature.location.end
                            strand = (
                                "+"
                                if feature.location.strand == 1
                                else "-"
                                if feature.location.strand == -1
                                else "."
                            )
                            phase = "0" if feature.type == "CDS" else "."

                            gff_out.write(
                                f"{seq_record.id}\tsnipy-ng\t{feature.type}\t{start}\t{end}\t.\t{strand}\t{phase}\tParent=transcript:{transcript_id}\n"
                            )
                            nfeat += 1

                # Write standalone features as genes without transcripts
                for feature in standalone_features:
                    gene_counter += 1
                    gene_id = f"gene_{gene_counter}"

                    if "locus_tag" in feature.qualifiers:
                        gene_id = feature.qualifiers["locus_tag"][0]
                    elif "gene" in feature.qualifiers:
                        gene_id = feature.qualifiers["gene"][0]

                    start = feature.location.start + 1  # Convert to 1-based
                    end = feature.location.end
                    strand = (
                        "+"
                        if feature.location.strand == 1
                        else "-"
                        if feature.location.strand == -1
                        else "."
                    )

                    # Determine biotype
                    biotype = "protein_coding"  # Default
                    if feature.type == "tRNA":
                        biotype = "tRNA"
                    elif feature.type == "rRNA":
                        biotype = "rRNA"
                    elif feature.type in ("ncRNA", "misc_RNA"):
                        biotype = "misc_RNA"

                    # Get gene name
                    gene_name = ""
                    if "gene" in feature.qualifiers:
                        gene_name = f";Name={feature.qualifiers['gene'][0]}"
                    elif "product" in feature.qualifiers:
                        gene_name = f";Name={feature.qualifiers['product'][0]}"

                    gff_out.write(
                        f"{seq_record.id}\tsnipy-ng\tgene\t{start}\t{end}\t.\t{strand}\t.\tID=gene:{gene_id};biotype={biotype}{gene_name}\n"
                    )
                    nfeat += 1
        

        # Remove existing metadata
        if self.output.metadata.exists():
            self.output.metadata.unlink()
        # Write JSON metadata
        metadata = ReferenceMetadata(
            self.output.metadata,
            reference=reference_path.name,
            format=ref_fmt,
            num_sequences=nseq,
            total_length=total_length,
            num_features=nfeat,
            prefix=str(self.reference_prefix),
            version=__version__,
            datetime=__import__("datetime").datetime.now().isoformat(),
        )
        metadata.validate()
        metadata.write()

        print(f"Wrote {nseq} sequences to {output_fasta_path}")
        print(
            f"Wrote {nfeat} features to {output_gff_path}"
            if nfeat > 0
            else f"No features found in {reference_path}"
        )


class ReferenceOutputImmutable(ReferenceOutput):
    _immutable = True


class LoadReferenceFromMetadataFile(BaseStage):
    metadata: Path = Field(METADATA_FILE_NAME, description="Path to reference metadata file")

    _dependencies = [
        biopython,
    ]

    @field_validator("metadata")
    @classmethod
    def must_not_be_none(cls, v):
        if v is None:
            raise ValueError("Metadata path must be provided")
        return v

    @property
    def output(self) -> ReferenceOutputImmutable:
        ref_metadata = ReferenceMetadata(self.metadata)
        # assume reference files are in the same directory as metadata
        reference_prefix = self.metadata.parent / ref_metadata.prefix
        return ReferenceOutputImmutable(
            reference=f"{reference_prefix}.fa",
            reference_index=f"{reference_prefix}.fa.fai",
            reference_dict=f"{reference_prefix}.dict",
            gff=f"{reference_prefix}.gff",
            metadata=self.metadata,
        )

    @property
    def commands(self):
        validate_reference_cmd = self.python_cmd(
            func=self.validate_reference,
            description="Validate reference files and metadata",
        )
        return [validate_reference_cmd]

    def validate_reference(self):
        """
        Validates that reference files exist and metadata is consistent.
        """

        # Check if all required files exist
        required_files = {
            "metadata": self.output.metadata,
            "reference FASTA": self.output.reference,
            "reference index": self.output.reference_index,
            "reference dict": self.output.reference_dict,
            "GFF": self.output.gff,
        }

        missing_files = []
        for file_type, file_path in required_files.items():
            if not file_path.exists():
                missing_files.append(f"{file_type}: {file_path}")

        if missing_files:
            raise ReferenceValidationError(
                "Missing required reference files:\n" + "\n".join(missing_files)
            )

        # Load and validate metadata
        metadata = ReferenceMetadata(self.metadata)
        metadata.validate()
        print(f"Reference metadata validated: {self.metadata}")
