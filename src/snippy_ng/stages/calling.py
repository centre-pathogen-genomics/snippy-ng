import json
import os

# Concrete Alignment Strategies
from pathlib import Path
from typing import List, Annotated, Optional

from snippy_ng.stages import BaseStage, BaseOutput, TempPath
from snippy_ng.dependencies import freebayes, bcftools, bedtools, paftools, clair3, longbow
from snippy_ng.exceptions import StageExecutionError
from snippy_ng.logging import logger

from pydantic import Field, AfterValidator


MIN_FREEBAYES_CHUNK_SIZE = 1000
MIN_CLAIR3_CHUNK_SIZE = 10000


def estimate_reference_bases(reference: Path, reference_index: Path) -> int:
    """Estimate reference size in bases, preferring the FASTA index if available."""
    if reference_index.exists():
        total = 0
        with open(reference_index, "r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    continue
                try:
                    total += int(fields[1])
                except ValueError:
                    continue
        if total > 0:
            return total
    return max(1, reference.stat().st_size)


def get_calling_chunk_size(reference: Path, reference_index: Path, cpus: int, min_chunk_size: int) -> tuple[int, int]:
    """Choose a chunk size that oversamples slightly relative to the available CPUs."""
    refsize = estimate_reference_bases(reference, reference_index)
    num_chunks = 1 + 2 * (max(1, cpus) - 1)
    chunk_size = max(min_chunk_size, int(refsize / num_chunks))
    return num_chunks, chunk_size


def get_short_chunk_size(reference: Path, reference_index: Path, cpus: int) -> tuple[int, int]:
    """Determine Freebayes chunk size based on reference size and available CPUs, with a minimum threshold."""
    return get_calling_chunk_size(reference, reference_index, cpus, MIN_FREEBAYES_CHUNK_SIZE)

def get_long_chunk_size(reference: Path, reference_index: Path, cpus: int) -> tuple[int, int]:
    """Determine Clair3 chunk size based on reference size and available CPUs, with a minimum threshold."""
    return get_calling_chunk_size(reference, reference_index, cpus, MIN_CLAIR3_CHUNK_SIZE)

def no_spaces(v: str) -> str:
    """Ensure that a string contains no spaces."""
    if " " in v:
        raise ValueError(
            "Prefix must not contain spaces, please use underscores or hyphens instead."
        )
    return v


# Define the base Pydantic model for alignment parameters
class Caller(BaseStage):
    reference: Path = Field(
        ...,
        description="Reference file",
    )
    reference_index: Path = Field(..., description="Reference index file")
    prefix: Annotated[str, AfterValidator(no_spaces)] = Field(
        ..., description="Output file prefix"
    )

    def test_check_if_vcf_has_variants(self):
        """Test that the output VCF file is not empty."""
        vcf_path = self.output.vcf
        with open(vcf_path, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    return  # Found a non-header line, VCF is not empty
        # give warning instead of error as some callers may produce empty VCFs if no variants are found
        logger.warning(f"Output VCF file {vcf_path} has no variants (only header lines). Please check if this is expected based on your data and parameters.")

class BaseCallerOutput(BaseOutput):
    vcf: Path = Field(..., description="VCF file containing raw variant calls")

class FreebayesCallerOutput(BaseCallerOutput):
    vcf: Path = Field(..., description="VCF file containing raw variant calls")
    regions: TempPath = Field(..., description="BED file with regions used for parallel calling")


class FreebayesCaller(Caller):
    """
    Call variants using Freebayes.
    """

    bam: Path = Field(..., description="Input BAM file")
    bam_index: Path = Field(..., description="Index file for the input BAM")
    fbopt: str = Field("", description="Additional Freebayes options")
    ploidy: int = Field(2, description="Ploidy for variant calling")
    exclude_insertions: bool = Field(
        True,
        description="Exclude insertions from variant calls so the pseudo-alignment remains the same length as the reference",
    )

    _dependencies = [freebayes, bcftools]

    @property
    def output(self) -> FreebayesCallerOutput:
        return FreebayesCallerOutput(
            vcf=self.prefix + ".raw.vcf",
            regions=self.prefix + ".regions.txt",
        )
    
    def create_commands(self, ctx) -> List:
        """Constructs the Freebayes variant calling and postprocessing commands."""
        num_chunks, chunk_size = get_short_chunk_size(self.reference, self.reference_index, ctx.cpus)
        logger.info(
            f"Freebayes will process {num_chunks} chunks of {chunk_size} bp, {ctx.cpus} chunks at a time."
        )

        # 1) Regions for parallel FreeBayes
        generate_regions_cmd = self.shell_cmd(
            ["fasta_generate_regions.py", str(self.reference_index), str(chunk_size)],
            description="Generate genomic regions for parallel variant calling",
        )
        generate_regions_pipeline = self.shell_pipe(
            commands=[generate_regions_cmd],
            description="Generate regions file for parallel processing",
            output_file=Path(self.output.regions),
        )

        # 2) FreeBayes parallel call
        freebayes_cmd_parts = [
            "freebayes-parallel",
            str(self.output.regions),
            str(ctx.cpus),
            "--ploidy", str(self.ploidy),
            "--genotype-qualities",
            "--min-alternate-count", "2",
            "--min-repeat-entropy", "1.0",
            "--min-base-quality", "13",
            "--min-mapping-quality", "60",
            "--strict-vcf",
        ]
        if self.fbopt:
            import shlex

            freebayes_cmd_parts.extend(shlex.split(self.fbopt))
        freebayes_cmd_parts.extend(["-f", str(self.reference), str(self.bam)])

        freebayes_cmd = self.shell_cmd(
            freebayes_cmd_parts,
            description="Call variants with FreeBayes in parallel",
            output_file=Path(self.output.vcf),
        )
        return [generate_regions_pipeline, freebayes_cmd]


class FreebayesCallerLong(FreebayesCaller):
    """
    Call variants using Freebayes for long-read data.
    """
    # def test_freebayes_ran_successfully(self):
        

    def create_commands(self, ctx) -> List:
        """Constructs the Freebayes variant calling and postprocessing commands."""
        num_chunks, chunk_size = get_long_chunk_size(self.reference, self.reference_index, ctx.cpus)
        logger.info(
            f"Freebayes will process {num_chunks} chunks of {chunk_size} bp, {ctx.cpus} chunks at a time."
        )

        # Regions for parallel FreeBayes
        generate_regions_cmd = self.shell_cmd(
            ["fasta_generate_regions.py", str(self.reference_index), str(chunk_size)],
            description="Generate genomic regions for parallel variant calling",
            output_file=Path(self.output.regions),
        )
        # FreeBayes parallel call
        freebayes_cmd_parts = [
            "freebayes-parallel",
            str(self.output.regions),
            str(ctx.cpus),
            "--ploidy", str(self.ploidy),
            "--genotype-qualities",
            "--haplotype-length", "-1",
            "--min-mapping-quality", "10",
            "--min-base-quality", "10",
        ]
        if self.fbopt:
            import shlex

            freebayes_cmd_parts.extend(shlex.split(self.fbopt))
        freebayes_cmd_parts.extend(["-f", str(self.reference), str(self.bam)])

        freebayes_cmd = self.shell_cmd(
            freebayes_cmd_parts,
            description="Call variants with FreeBayes in parallel",
            output_file=Path(self.output.vcf),
        )
        
        return [generate_regions_cmd, freebayes_cmd]

class PAFCallerOutput(BaseCallerOutput):
    vcf: Path = Field(..., description="VCF file with raw PAF-derived variant calls and annotations")
    tmp_vcf: TempPath = Field(..., description="Temporary intermediate VCF generated by paftools.js")
    missing_bed: Path = Field(..., description="BED file of unaligned (missing) reference regions")
    aln_bed: TempPath = Field(..., description="Temporary BED file of merged aligned reference intervals")
    annotations_file: TempPath = Field(..., description="Temporary bgzip-compressed annotation table used to populate VCF FORMAT fields")
    annotations_file_index: TempPath = Field(..., description="Temporary tabix index for the compressed annotation table")


class PAFCaller(Caller):
    """
    Call variants from PAF alignments using paftools.js.
    """

    paf: Path = Field(..., description="Input PAF file")
    ref_dict: Path = Field(..., description="Reference FASTA dictionary file")
    mapq: Optional[int] = Field(
        15, description="Minimum mapping quality for variant calling"
    )
    alen: Optional[int] = Field(
        50, description="Minimum alignment length for variant calling"
    )

    _dependencies = [bedtools, bcftools, paftools]

    @property
    def output(self) -> PAFCallerOutput:
        return PAFCallerOutput(
            vcf=Path(f"{self.prefix}.raw.vcf"),
            aln_bed=Path(f"{self.prefix}.aln.bed"),
            missing_bed=Path(f"{self.prefix}.missing.bed"),
            annotations_file=Path(f"{self.prefix}.annotations.gz"),
            annotations_file_index=Path(f"{self.prefix}.annotations.gz.tbi"),
            tmp_vcf=Path(f"{self.prefix}.tmp.vcf"),
        )

    def create_commands(self, ctx) -> List:
        """Constructs the PAF processing and BED generation commands."""

        # 4) Convert PAF to merged aligned reference intervals (BED)
        # Keep primary or pseudo-primary hits: tp:A:P or tp:A:I
        paf_to_pipeline = self.shell_pipe(
            commands=[
                self.shell_cmd(
                    ["grep", "-E", "tp:A:[PI]", str(self.paf)],
                    description="Filter PAF for primary or pseudo-primary alignments",
                ),
                self.shell_cmd(
                    ["cut", "-f6,8,9"],
                    description="Extract relevant PAF fields (reference name, start, end)",
                ),
                self.shell_cmd(
                    ["bedtools", "sort", "-i", "-"],
                    description="Sort BED entries",
                ),
                self.shell_cmd(
                    ["bedtools", "merge", "-i", "-"],
                    description="Merge overlapping BED intervals",
                ),
            ],
            description="Convert PAF to merged aligned reference intervals (BED)",
            output_file=self.output.aln_bed,
        )

        # 5) Compute unaligned (missing) reference regions
        compute_missing_bed_cmd = self.shell_cmd(
            [
                "bedtools",
                "complement",
                "-g",
                str(self.ref_dict),
                "-i",
                str(self.output.aln_bed),
            ],
            description="Compute unaligned (missing) reference regions",
            output_file=self.output.missing_bed,
        )

        # variant calling
        paftools_cmd = self.shell_cmd(
                    [
                        "paftools.js",
                        "call",
                        "-q",
                        str(self.mapq),
                        "-L",
                        str(self.alen),
                        "-l",
                        str(self.alen),
                        "-s",
                        self.prefix,
                        "-f",
                        str(self.reference),
                        str(self.paf),
                    ],
                    description="Call variants from PAF using paftools.js",
                    output_file=self.output.tmp_vcf,
                )

        create_annotation_file_pipeline = self.shell_pipe(
            commands=[
                self.shell_cmd(
                    [
                        "bcftools",
                        "query",
                        str(paftools_cmd.output_file),
                        "-f",
                        "%CHROM\\t%POS\\t1\\t1\\n",
                    ],
                    description="Extract variant positions from VCF",
                ),
                self.shell_cmd(
                    ["bgzip", "-c"],
                    description="Compress annotation file with bgzip",
                ),
            ],
            description="Compress annotation file with bgzip",
            output_file=self.output.annotations_file,
        )
        index_cmd = self.shell_cmd(
            [
                "tabix",
                "-f", # Force overwrite if index already exists
                "-s1",
                "-b2",
                "-e2",
                str(self.output.annotations_file),
            ],
            description="Index annotation file with tabix",
        )
        annotations_cmd = self.shell_cmd(
                    [
                        "bcftools",
                        "annotate",
                        "-H", '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                        "-H", '##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count for each ALT">',
                        "-c", "CHROM,POS,FMT/DP:=FORMAT/DP,FMT/AO:=FORMAT/AO",
                        "-a", str(self.output.annotations_file),
                        str(paftools_cmd.output_file),
                    ],
                    description="Insert FORMAT header lines for DP and AO",
                    output_file=self.output.vcf,
                )
        # TODO: Combine vcf with missing bed using <DEL> blocks so we have a single 
        # VCF output with both variant calls and missing regions annotated as deletions.
        return [
            paf_to_pipeline,
            compute_missing_bed_cmd,
            paftools_cmd,
            create_annotation_file_pipeline,
            index_cmd,
            annotations_cmd,
        ]

class Clair3CallerOutput(BaseCallerOutput):
    vcf: Path = Field(..., description="VCF file containing raw variant calls produced by Clair3")


class Clair3ModelSelectorError(StageExecutionError):
    """Raised when the Clair3 model cannot be resolved based on Longbow predictions."""
    pass

class LongbowClair3ModelOutput(BaseOutput):
    prediction_json: TempPath = Field(..., description="Temporary Longbow prediction report")
    clair3_model: Path = Field(..., description="Text file containing the resolved Clair3 model path for the current run")


class LongbowClair3ModelSelector(BaseStage):
    reads: Path = Field(..., description="Input FASTQ file used for Longbow prediction")

    _dependencies = [longbow]

    @property
    def output(self) -> LongbowClair3ModelOutput:
        return LongbowClair3ModelOutput(
            prediction_json=Path(f"{self.prefix}.longbow.json"),
            clair3_model=Path(f"{self.prefix}.clair3_model.txt"),
        )

    def create_commands(self, ctx) -> List:
        logger.warning("Running Longbow to predict ONT basecalling configuration for Clair3 model selection. This may take a few minutes... Specify --clair3-model to skip this step and use a specific model.")
        return [
            self.shell_cmd(
                [
                    "longbow",
                    "-i", str(self.reads),
                    "-o", str(self.output.prediction_json),
                    "-t", str(ctx.cpus),
                ],
                description="Predict ONT basecalling configuration with Longbow",
            ),
            self.python_cmd(
                self.resolve_clair3_model,
                args=[self.output.prediction_json, self.output.clair3_model],
                description="Resolve the Clair3 model from Longbow output",
            ),
        ]

    @staticmethod
    def _normalize_string(value) -> Optional[str]:
        if value is None:
            return None
        if isinstance(value, (int, float)):
            value = str(value)
        if not isinstance(value, str):
            return None
        normalized = value.strip().lower()
        return normalized or None

    @classmethod
    def _get_nested_value(cls, data, *paths) -> Optional[str]:
        for path in paths:
            current = data
            for key in path:
                if not isinstance(current, dict):
                    current = None
                    break
                current = current.get(key)
            normalized = cls._normalize_string(current)
            if normalized is not None:
                return normalized
        return None

    @classmethod
    def _parse_longbow_prediction(cls, prediction_json: Path) -> dict[str, Optional[str]]:
        with open(prediction_json, "r", encoding="utf-8") as handle:
            data = json.load(handle)

        if not isinstance(data, dict):
            raise ValueError(f"Unexpected Longbow JSON format in '{prediction_json}'.")

        return {
            "basecaller": cls._get_nested_value(
                data,
                ("basecaller",),
                ("Software",),
                ("basecalling_software",),
                ("prediction", "basecaller"),
                ("result", "basecaller"),
            ),
            "flowcell": cls._get_nested_value(
                data,
                ("flowcell",),
                ("Flowcell",),
                ("flowcell_version",),
                ("prediction", "flowcell"),
                ("prediction", "flowcell_version"),
                ("result", "flowcell"),
            ),
            "major_version": cls._get_nested_value(
                data,
                ("major_version",),
                ("Version",),
                ("basecaller_version",),
                ("basecaller_major_version",),
                ("prediction", "major_version"),
                ("prediction", "basecaller_version"),
                ("result", "major_version"),
            ),
            "mode": cls._get_nested_value(
                data,
                ("mode",),
                ("Mode",),
                ("basecalling_mode",),
                ("prediction", "mode"),
                ("prediction", "basecalling_mode"),
                ("result", "mode"),
            ),
            "dorado_model_version": cls._get_nested_value(
                data,
                ("dorado_model_version",),
                ("prediction", "dorado_model_version"),
                ("result", "dorado_model_version"),
            ),
        }

    @staticmethod
    def _longbow_summary(prediction: dict[str, Optional[str]]) -> str:
        return ", ".join(
            [
                f"flowcell={prediction.get('flowcell') or 'unknown'}",
                f"basecaller={prediction.get('basecaller') or 'unknown'}",
                f"major_version={prediction.get('major_version') or 'unknown'}",
                f"mode={prediction.get('mode') or 'unknown'}",
            ]
        )

    @staticmethod
    def _candidate_roots() -> list[Path]:
        launch_dir = Path(os.environ.get("PWD", str(Path.cwd()))).expanduser()
        roots: list[Path] = []
        for env_var in ("CLAIR3_MODELS", "CLAIR3_MODEL_ROOT"):
            value = os.environ.get(env_var)
            if value:
                root = Path(value).expanduser()
                if not root.is_absolute():
                    root = launch_dir / root
                roots.append(root)

        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix:
            roots.extend(
                [
                    Path(conda_prefix) / "bin" / "models",
                    Path(conda_prefix) / "models",
                ]
            )

        roots.extend(
            [
                Path("/opt/models"),
                Path("/opt/models/clair3_models"),
            ]
        )

        unique_roots: list[Path] = []
        seen: set[Path] = set()
        for root in roots:
            root = root.expanduser()
            if root in seen:
                continue
            seen.add(root)
            unique_roots.append(root)
        return unique_roots

    @staticmethod
    def _container_roots() -> list[Path]:
        roots: list[Path] = []
        for env_var in ("CLAIR3_MODEL_CONTAINER_ROOT", "CLAIR3_MODELS_CONTAINER_ROOT"):
            value = os.environ.get(env_var)
            if value:
                roots.append(Path(value))
        return roots

    @staticmethod
    def _unique_preserving_order(items: list[str]) -> list[str]:
        seen: set[str] = set()
        unique: list[str] = []
        for item in items:
            if item in seen:
                continue
            seen.add(item)
            unique.append(item)
        return unique

    @staticmethod
    def _r10_model_series(mode: str) -> list[str]:
        return [
            f"r1041_e82_400bps_{mode}_v520",
            f"r1041_e82_400bps_{mode}_v500",
            f"r1041_e82_400bps_{mode}_v430",
            f"r1041_e82_400bps_{mode}_v410",
        ]

    @staticmethod
    def _r10_guppy56_model_series(mode: str) -> list[str]:
        if mode == "sup":
            return [
                "r1041_e82_400bps_sup_g615",
                "r1041_e82_260bps_sup_g632",
            ]
        if mode == "hac":
            return [
                "r1041_e82_400bps_hac_g632",
                "r1041_e82_400bps_hac_g615",
            ]
        if mode == "fast":
            return [
                "r1041_e82_400bps_fast_g632",
                "r1041_e82_400bps_fast_g615",
                "r1041_e82_260bps_fast_g632",
            ]
        return []

    @classmethod
    def _candidate_model_names(cls, prediction: dict[str, Optional[str]]) -> list[str]:
        basecaller = prediction.get("basecaller") or ""
        flowcell = prediction.get("flowcell") or ""
        major_version = prediction.get("major_version") or ""
        mode = prediction.get("mode") or ""
        is_guppy56 = "guppy" in basecaller and (
            "5or6" in major_version or "guppy5" in major_version or "guppy6" in major_version
        )

        if "r10" in flowcell:
            candidates: list[str] = []
            if "sup" in mode:
                candidates.extend(cls._r10_model_series("sup")[:2])
                if is_guppy56:
                    candidates.extend(cls._r10_guppy56_model_series("sup"))
                candidates.extend(cls._r10_model_series("sup")[2:])
            elif "hac" in mode:
                candidates.extend(cls._r10_model_series("hac")[:2])
                if is_guppy56:
                    candidates.extend(cls._r10_guppy56_model_series("hac"))
                candidates.extend(cls._r10_model_series("hac")[2:])
            elif "fast" in mode:
                candidates.extend(cls._r10_guppy56_model_series("fast"))
            else:
                candidates.extend(
                    [
                        "r1041_e82_400bps_sup_v520",
                        "r1041_e82_400bps_sup_v500",
                        "r1041_e82_400bps_hac_v520",
                        "r1041_e82_400bps_hac_v500",
                        "r1041_e82_400bps_sup_v410",
                        "r1041_e82_400bps_hac_v410",
                    ]
                )
            return cls._unique_preserving_order(candidates)

        if "r9" in flowcell and "guppy" in basecaller:
            if "guppy2" in major_version:
                return ["ont_guppy2", "r941_prom_hac_g238"]
            if "3or4" in major_version or "guppy3" in major_version or "guppy4" in major_version:
                return ["ont", "r941_prom_hac_g360+g422", "r941_prom_hac_g360+g422_1235"]
            if "5or6" in major_version or "guppy5" in major_version or "guppy6" in major_version:
                if "sup" in mode:
                    return ["ont_guppy5", "r941_prom_sup_g5014", "ont"]
                return ["ont_guppy5", "ont", "r941_prom_sup_g5014"]

        if "r9" in flowcell and ("dorado" in basecaller or "dorado0" in major_version or major_version == "0"):
            return ["ont_guppy5", "ont", "r941_prom_sup_g5014"]

        if "guppy2" in major_version:
            return ["ont_guppy2", "r941_prom_hac_g238"]
        if "guppy5" in major_version or "guppy6" in major_version or "dorado0" in major_version:
            return ["ont_guppy5", "r941_prom_sup_g5014", "ont"]
        if "guppy3" in major_version or "guppy4" in major_version:
            return ["ont", "r941_prom_hac_g360+g422", "r941_prom_hac_g360+g422_1235"]

        return ["ont_guppy5", "ont", "ont_guppy2", "r941_prom_sup_g5014", "r941_prom_hac_g360+g422", "r941_prom_hac_g238"]

    @staticmethod
    def _find_model_directory(root: Path, model_name: str) -> Optional[Path]:
        direct = root / model_name
        if direct.is_dir():
            return direct
        nested = root / "clair3_models" / model_name
        if nested.is_dir():
            return nested
        return None

    @classmethod
    def _prefer_v500_for_container(cls, candidates: list[str]) -> list[str]:
        def sort_key(model_name: str) -> tuple[int, int]:
            if "_v500" in model_name:
                return (0, 0)
            if "_v520" in model_name:
                return (0, 1)
            return (1, 0)

        head = sorted(
            [name for name in candidates if "_v500" in name or "_v520" in name],
            key=sort_key,
        )
        tail = [name for name in candidates if "_v500" not in name and "_v520" not in name]
        return cls._unique_preserving_order(head + tail)

    @classmethod
    def _resolve_local_model(cls, prediction: dict[str, Optional[str]]) -> Path:
        roots = [root for root in cls._candidate_roots() if root.exists()]
        candidates = cls._candidate_model_names(prediction)
        for model_name in candidates:
            for root in roots:
                resolved = cls._find_model_directory(root, model_name)
                if resolved is not None:
                    return resolved.resolve()

        container_roots = cls._container_roots()
        if container_roots:
            # TODO: the containers currently ship with only the v500 models
            # so we reorder the candidates to prefer the v500 models if 
            # container roots are being used. Remove once v520 is standard in containers.
            container_candidates = cls._prefer_v500_for_container(candidates)
            return container_roots[0] / container_candidates[0]

        if not roots:
            raise Clair3ModelSelectorError(
                "Could not find any Clair3 model roots. Set CLAIR3_MODELS for local models or CLAIR3_MODEL_CONTAINER_ROOT for container-only model paths."
            )

        raise Clair3ModelSelectorError(
            "Could not resolve a Clair3 model from Longbow prediction "
            f"({cls._longbow_summary(prediction)}). Searched roots: {', '.join(str(root) for root in roots)}. "
            f"Tried models: {', '.join(candidates)}. Provide --clair3-model to override."
        )

    @classmethod
    def resolve_clair3_model(cls, prediction_json: Path, output_path: Path) -> None:
        prediction = cls._parse_longbow_prediction(prediction_json)
        model_dir = cls._resolve_local_model(prediction)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(f"{model_dir}\n", encoding="utf-8")

        logger.info(
            f"Resolved Clair3 model '{model_dir}' from Longbow prediction ({cls._longbow_summary(prediction)})."
        )


class Clair3Caller(Caller):
    """
    Call variants using Clair3.
    """
    bam: Path = Field(..., description="Input BAM file")
    clair3_model: Path = Field(..., description="Path to Clair3 model")
    platform: str = Field("ont", description="Sequencing platform (e.g., ont, hifi)")
    fast_mode: bool = Field(True, description="Enable fast mode for Clair3")

    _dependencies = [clair3]

    @property
    def output(self) -> Clair3CallerOutput:
        return Clair3CallerOutput(
            vcf=Path(f"{self.prefix}.raw.vcf"),
        )

    def _resolve_model_path_arg(self) -> Path:
        model_path = Path(self.clair3_model)
        if model_path.exists() and model_path.is_file():
            resolved = model_path.read_text(encoding="utf-8").strip()
            if not resolved:
                raise Clair3ModelSelectorError(f"Resolved Clair3 model file '{model_path}' is empty.")
            return Path(resolved)
        return model_path

    def create_commands(self, ctx) -> List:
        """Constructs the Clair3 variant calling commands."""
        model_path = self._resolve_model_path_arg()
        _, chunk_size = get_long_chunk_size(self.reference, self.reference_index, ctx.cpus)

        clair3_cmd = self.shell_cmd(
            [
                "run_clair3.sh",
                f"--model_path={model_path.absolute()}",
                f"--bam_fn={str(self.bam.absolute())}",
                f"--ref_fn={str(self.reference.absolute())}",
                f"--threads={str(ctx.cpus)}",
                f"--output={Path(self.prefix + '_clair3_out').absolute()}",
                f"--platform={self.platform}",
                f"--chunk_size={chunk_size}",
                "--include_all_ctgs",
                "--no_phasing_for_fa",
                "--enable_long_indel",
            ],
            description="Call variants with Clair3",
        )
        if self.fast_mode:
            clair3_cmd.command.append("--fast_mode")
        unzip_cmd = self.shell_cmd(
            [
                "gunzip",
                f"{self.prefix}_clair3_out/merge_output.vcf.gz",
            ],
            description="Unzip Clair3 VCF output",
        )
        move_vcf_cmd = self.shell_cmd(
            [
                "mv",
                f"{self.prefix}_clair3_out/merge_output.vcf",
                str(self.output.vcf),
            ],
            description="Move Clair3 VCF to final output location",
        )
        cleanup_cmd = self.shell_cmd(
            [
                "rm",
                "-rf",
                f"{self.prefix}_clair3_out",
            ],
            description="Clean up Clair3 output directory",
        )

        return [clair3_cmd, unzip_cmd, move_vcf_cmd, cleanup_cmd] 
