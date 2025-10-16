from snippy_ng.stages.base import BaseStage, BaseOutput
from pydantic import Field, field_validator
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

from snippy_ng.dependencies import biopython

class ReferenceOutput(BaseOutput):
    reference: Path
    reference_index: Path
    reference_dict: Path
    gff: Path
    meta: Path


class PrepareReference(BaseStage):
    input: Path = Field(..., description="Reference file")
    ref_fmt: str = Field("genbank", description="Reference format")
    reference_prefix: str = Field("ref", description="Output reference name")
    reference_dir: Path = Field(Path("reference"), description="Reference directory")

    _dependencies = [
        biopython,
    ]

    @property
    def output(self) -> ReferenceOutput:
        return ReferenceOutput(
            reference=self.reference_dir / f"{self.reference_prefix}.fa",
            reference_index=self.reference_dir / f"{self.reference_prefix}.fa.fai",
            reference_dict=self.reference_dir / f"{self.reference_prefix}.dict",
            gff=self.reference_dir / f"{self.reference_prefix}.gff",
            meta=self.reference_dir / "metadata.json",
            
        )

    @property
    def commands(self):
            process_reference_cmd = self.python_cmd(
                func=self.process_reference,
                args=(self.input, self.ref_fmt, self.output.reference, self.output.gff),
                description=f"Extract FASTA and GFF from reference ({self.ref_fmt})"
            )
            return [
                self.shell_cmd([
                    "rm", "-f", str(self.output.reference)
                ], description=f"Remove existing reference FASTA: {self.output.reference}"),
                self.shell_cmd([
                    "mkdir", "-p", str(self.reference_dir)
                ], description=f"Create reference directory: {self.reference_dir}"),
                self.shell_cmd([
                     "cp", str(self.input), str(self.reference_dir / self.input.name)
                ], description=f"Copy reference file to reference directory: {self.input} -> {self.reference_dir / self.input.name}"),
                process_reference_cmd,
                self.shell_cmd([
                    "samtools", "faidx", str(self.output.reference)
                ], description=f"Index reference FASTA with samtools faidx: {self.output.reference}"),
                self.shell_pipeline(
                    commands=[
                        self.shell_cmd([
                            "cut", "-f1,2", str(self.output.reference_index)
                        ], description="Extract sequence names and lengths from FASTA index"),
                        self.shell_cmd([
                            "sort"
                        ], description="Sort sequence names and lengths"),
                    ],
                    output_file=self.output.reference_dict,
                    description=f"Create reference dictionary: {self.output.reference_dict}"
                )
            ]
    
    def process_reference(self, reference_path: Path, ref_fmt: str, output_fasta_path: Path, output_gff_path: Path):
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
                with open(reference_path, 'rt') as test_fh:
                    test_fh.read(1)
            except UnicodeDecodeError:
                open_func = gzip.open
            with open_func(reference_path, 'rt') as ref_fh:
                seq_records = list(SeqIO.parse(ref_fh, ref_fmt))
        except Exception as e:
            raise ValueError(f"Failed to parse {reference_path} with format {ref_fmt}: {e}")

        # Prepare outputs
        ref_seq_dict = {}
        gene_counter = 0
        nseq = 0
        nfeat = 0
        total_length = 0
        
        with open(output_fasta_path, "w") as fasta_out, open(output_gff_path, "w") as gff_out:
            # Write GFF3 header
            gff_out.write("##gff-version 3\n")
            
            for seq_record in seq_records:
                # Check for duplicate sequences
                if seq_record.id in ref_seq_dict:
                    raise ValueError(f"Duplicate sequence {seq_record.id} in {reference_path}")

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
                                'feature': feature,
                                'transcripts': {},
                                'start': feature.location.start + 1,  # Convert to 1-based
                                'end': feature.location.end,
                                'strand': '+' if feature.location.strand == 1 else '-' if feature.location.strand == -1 else '.'
                            }
                        else:
                            # Extend gene boundaries
                            genes[gene_id]['start'] = min(genes[gene_id]['start'], feature.location.start + 1)
                            genes[gene_id]['end'] = max(genes[gene_id]['end'], feature.location.end)
                        
                        if transcript_id not in genes[gene_id]['transcripts']:
                            genes[gene_id]['transcripts'][transcript_id] = {
                                'features': [],
                                'start': feature.location.start + 1,
                                'end': feature.location.end,
                                'strand': '+' if feature.location.strand == 1 else '-' if feature.location.strand == -1 else '.'
                            }
                        else:
                            # Extend transcript boundaries
                            transcript = genes[gene_id]['transcripts'][transcript_id]
                            transcript['start'] = min(transcript['start'], feature.location.start + 1)
                            transcript['end'] = max(transcript['end'], feature.location.end)
                        
                        genes[gene_id]['transcripts'][transcript_id]['features'].append(feature)
                    
                    elif ftype == "gene":
                        if gene_id not in genes:
                            genes[gene_id] = {
                                'feature': feature,
                                'transcripts': {},
                                'start': feature.location.start + 1,
                                'end': feature.location.end,
                                'strand': '+' if feature.location.strand == 1 else '-' if feature.location.strand == -1 else '.'
                            }
                    else:
                        # Standalone features (tRNA, rRNA, etc.)
                        standalone_features.append(feature)

                # Write genes and their transcripts
                for gene_id, gene_data in genes.items():
                    gene_feature = gene_data['feature']
                    
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
                    gff_out.write(f"{seq_record.id}\tsnipy-ng\tgene\t{gene_data['start']}\t{gene_data['end']}\t.\t{gene_data['strand']}\t.\tID=gene:{gene_id};biotype={biotype}{gene_name}\n")
                    nfeat += 1
                    
                    # Write transcripts and their features
                    for transcript_id, transcript_data in gene_data['transcripts'].items():
                        # Write transcript line
                        gff_out.write(f"{seq_record.id}\tsnipy-ng\ttranscript\t{transcript_data['start']}\t{transcript_data['end']}\t.\t{transcript_data['strand']}\t.\tID=transcript:{transcript_id};Parent=gene:{gene_id};biotype={biotype}\n")
                        nfeat += 1
                        
                        # Write transcript features (CDS, exon, UTRs)
                        for feature in transcript_data['features']:
                            start = feature.location.start + 1  # Convert to 1-based
                            end = feature.location.end
                            strand = '+' if feature.location.strand == 1 else '-' if feature.location.strand == -1 else '.'
                            phase = "0" if feature.type == "CDS" else "."
                            
                            gff_out.write(f"{seq_record.id}\tsnipy-ng\t{feature.type}\t{start}\t{end}\t.\t{strand}\t{phase}\tParent=transcript:{transcript_id}\n")
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
                    strand = '+' if feature.location.strand == 1 else '-' if feature.location.strand == -1 else '.'
                    
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
                    
                    gff_out.write(f"{seq_record.id}\tsnipy-ng\tgene\t{start}\t{end}\t.\t{strand}\t.\tID=gene:{gene_id};biotype={biotype}{gene_name}\n")
                    nfeat += 1
                    
        # Write JSON metadata
        metadata = {
            "reference": reference_path.name,
            "format": ref_fmt,
            "num_sequences": nseq,
            "total_length": total_length,
            "num_features": nfeat,
            "prefix": self.reference_prefix,
            "datetime": __import__("datetime").datetime.now().isoformat(),
        }
        with open(self.output.meta, "w") as json_out:
            import json
            json.dump(metadata, json_out, indent=4)

        print(f"Wrote {nseq} sequences to {output_fasta_path}")
        print(f"Wrote {nfeat} features to {output_gff_path}" if nfeat > 0 else f"No features found in {reference_path}")


class LoadReference(BaseStage):
    reference_dir: Path = Field(..., description="Path to reference directory")
    reference_prefix: str = Field("ref", description="Reference prefix")

    _dependencies = [
        biopython,
    ]

    @field_validator("reference_dir")
    @classmethod
    def must_contain_metadata(cls, v: Path):
        meta_path = v / "metadata.json"
        if not meta_path.exists():
            raise ValueError(f"Reference was given as a directory ({v}) but does not contain the required metadata.json file!")
        return v

    @property
    def output(self) -> ReferenceOutput:
        return ReferenceOutput(
            reference=self.reference_dir / f"{self.reference_prefix}.fa",
            reference_index=self.reference_dir / f"{self.reference_prefix}.fa.fai",
            reference_dict=self.reference_dir / f"{self.reference_prefix}.dict",
            gff=self.reference_dir / f"{self.reference_prefix}.gff",
            meta=self.reference_dir / "metadata.json",
        )

    @property
    def commands(self):
        validate_reference_cmd = self.python_cmd(
            func=self.validate_reference,
            description="Validate reference files and metadata"
        )
        return [validate_reference_cmd]
    
    def validate_reference(self):
        """
        Validates that reference files exist and metadata is consistent.
        Optionally compares against an expected reference file.

        Args:
            meta_path (Path): Path to the metadata.json file.
            reference_path (Path): Path to the reference FASTA file.
            reference_index_path (Path): Path to the reference index file.
            gff_path (Path): Path to the GFF file.
        """
        import json

        # Check if all required files exist
        required_files = {
            "metadata": self.output.meta,
            "reference FASTA": self.output.reference,
            "reference index": self.output.reference_index,
            "reference dict": self.output.reference_dict,
            "GFF": self.output.gff
        }
        
        missing_files = []
        for file_type, file_path in required_files.items():
            if not file_path.exists():
                missing_files.append(f"{file_type}: {file_path}")
        
        if missing_files:
            raise FileNotFoundError("Missing required reference files:\n" + "\n".join(missing_files))

        # Load and validate metadata
        try:
            with open(self.output.meta, 'r') as f:
                metadata = json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            raise ValueError(f"Failed to read metadata file {self.output.meta}: {e}")

        required_metadata_fields = ["reference", "format", "num_sequences", "total_length", "num_features", "prefix", "datetime"]
        missing_fields = [field for field in required_metadata_fields if field not in metadata]
        if missing_fields:
            raise ValueError(f"Missing required metadata fields: {missing_fields}")