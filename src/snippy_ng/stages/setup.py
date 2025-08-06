from snippy_ng.stages.base import BaseStage, BaseOutput, Command
from pydantic import Field
from pathlib import Path
from Bio import SeqIO
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation

from snippy_ng.dependencies import biopython

class PrepareReferenceOutput(BaseOutput):
    reference: Path
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
    def output(self) -> PrepareReferenceOutput:
        return PrepareReferenceOutput(
            reference=self.reference_dir / f"{self.reference_prefix}.fa",
            gff=self.reference_dir / f"{self.reference_prefix}.gff",
            meta=self.reference_dir / f"{self.reference_prefix}.json",
        )

    @property
    def commands(self):
        process_reference_cmd = Command(
            func=self.process_reference, 
            args=(self.input, self.ref_fmt, self.output.reference, self.output.gff), 
            description=f"Extract FASTA and GFF from reference ({self.ref_fmt})"
        ) 
        return [
            f"rm -f {self.output.reference}",
            f"mkdir -p {self.reference_dir}",
            process_reference_cmd,
        ]
    
    def process_reference(self, reference_path: Path, ref_fmt: str, output_fasta_path: Path, output_gff_path: Path):
        """
        Extracts FASTA and GFF3 from a reference file.
        Determines input format and writes GFF only if features exist.

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
        feature_id_counter = {}
        nseq = 0
        nfeat = 0
        total_length = 0
        with open(output_fasta_path, "w") as fasta_out, open(output_gff_path, "w") as gff_out:
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

                # Process features for GFF
                new_features = []
                for feature in seq_record.features:
                    ftype = feature.type
                    if ftype in ("source", "gene", "misc_feature"):
                        continue  # Skip unwanted features

                    # Count features by type
                    if ftype not in feature_id_counter:
                        feature_id_counter[ftype] = 0
                    feature_id_counter[ftype] += 1

                    # Add ID to qualifiers
                    if "locus_tag" in feature.qualifiers:
                        feature_id = feature.qualifiers["locus_tag"][0]
                    else:
                        feature_id = f"{ftype}_{feature_id_counter[ftype]}"
                    feature.qualifiers["ID"] = feature_id

                    # Add Name to qualifiers if gene tag is present
                    if "gene" in feature.qualifiers:
                        feature.qualifiers["Name"] = feature.qualifiers["gene"][0]

                    # Assign source
                    feature.qualifiers["source"] = "snippy-ng"

                    # Set phase for CDS features; '.' for others
                    phase = "0" if ftype == "CDS" else "."
                    feature.qualifiers["phase"] = phase

                    new_feature = SeqFeature(
                        location=SimpleLocation(feature.location.start, feature.location.end, strand=feature.location.strand),
                        type=ftype,
                        qualifiers=feature.qualifiers,
                    )
                    new_features.append(new_feature)
                    nfeat += 1

                # Update record features
                seq_record.features = new_features

                # Write GFF features if any
                if new_features:
                    GFF.write([seq_record], gff_out)
                
        # Write JSON metadata
        metadata = {
            "reference": str(reference_path),
            "format": ref_fmt,
            "num_sequences": nseq,
            "total_length": total_length,
            "num_features": nfeat,
            "feature_counts": {ftype: count for ftype, count in feature_id_counter.items()},
        }
        with open(self.output.meta, "w") as json_out:
            import json
            json.dump(metadata, json_out, indent=4)


        print(f"Wrote {nseq} sequences to {output_fasta_path}")
        print(f"Wrote {nfeat} features to {output_gff_path}" if nfeat > 0 else f"No features found in {reference_path}")
