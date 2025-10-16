from typing import Optional
from dataclasses import dataclass
from shutil import which
import subprocess
import re

from snippy_ng.exceptions import (
    InvalidDependencyError,
    MissingDependencyError,
    InvalidDependencyVersionError,
)

from packaging.version import parse, InvalidVersion, Version


@dataclass
class Dependency:
    name: str
    citation: Optional[str] = None
    version_pattern: str = r"(\d+\.\d+\.\d+)"  # Regex pattern to extract version
    version_arg: Optional[str] = "--version"
    version: Optional[str] = None
    min_version: Optional[str] = None
    max_version: Optional[str] = None
    less_then: Optional[str] = None

    def check(self):
        if not which(self.name):
            raise MissingDependencyError(
                f"Could not find dependency {self.name}! Please install it."
            )
        version = self.get_version_from_cli()
        return self._base_validator(version)

    def format_version_requirements(self):
        requirements = []
        if self.version:
            requirements.append(f"={self.version}")
        if self.min_version:
            requirements.append(f">={self.min_version}")
        if self.max_version:
            requirements.append(f"<={self.max_version}")
        if self.less_then:
            requirements.append(f"<{self.less_then}")
        if not requirements:
            return self.name
        return f'{self.name} {",".join(requirements)}'
    
    def get_version_from_cli(self):
        cmd = [self.name]
        if self.version_arg:
            cmd.append(self.version_arg)

        result = subprocess.run(
            cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, text=True
        ).stdout.strip()

        # Apply regex pattern from version_extractor to extract version
        match = re.search(self.version_pattern, result)
        if not match:
            raise InvalidDependencyVersionError(
                f"Could not extract version from '{result}' for {self.name}."
            )

        return match.group(0)
    
    def _base_validator(self, version_str) -> Version:
        try:
            parsed_version = parse(version_str)
        except InvalidVersion:
            raise InvalidDependencyVersionError(
                f"Could not parse version '{version_str}' for {self.name}."
            )

        # Validate the extracted version against the given constraints
        if self.version and parsed_version != parse(self.version):
            raise InvalidDependencyError(
                f"{self.name} version must be {self.version} (found {parsed_version})"
            )
        if self.min_version and parsed_version < parse(self.min_version):
            raise InvalidDependencyError(
                f"{self.name} minimum version allowed is {self.min_version} (found {parsed_version})"
            )
        if self.max_version and parsed_version > parse(self.max_version):
            raise InvalidDependencyError(
                f"{self.name} maximum version allowed is {self.max_version} (found {parsed_version})"
            )
        if self.less_then and parsed_version >= parse(self.less_then):
            raise InvalidDependencyError(
                f"{self.name} version must be less than {self.less_then} (found {parsed_version})"
            )

        return parsed_version
 
class PythonDependency(Dependency):
    def check(self):
        from importlib.metadata import version
        from importlib.metadata import PackageNotFoundError

        try:
            version = version(self.name)
        except PackageNotFoundError:
            raise MissingDependencyError(
                f"Could not find dependency {self.name}! Please install it."
            )
        return self._base_validator(version)

# Python Dependencies
biopython = PythonDependency(
    "biopython",
    citation="Peter J. A. Cock, Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg, Thomas Hamelryck, Frank Kauff, Bartek Wilczynski, Michiel J. L. de Hoon, Biopython: freely available Python tools for computational molecular biology and bioinformatics, Bioinformatics, Volume 25, Issue 11, June 2009, Pages 1422-1423, https://doi.org/10.1093/bioinformatics/btp163"
)

# Alignment
samtools = Dependency(
    "samtools",
    citation="Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li, Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008",
    less_then="1.21",
    version_pattern=r"(\d+\.\d+)",
)
samclip = Dependency("samclip", citation="", min_version="0.4.0")
bwa = Dependency(
    "bwa",
    citation="Heng Li, Richard Durbin, Fast and accurate short read alignment with Burrows-Wheeler transform, Bioinformatics, Volume 25, Issue 14, July 2009, Pages 1754-1760, https://doi.org/10.1093/bioinformatics/btp324",
    version_arg=None,
    version_pattern=r"(\d+\.\d+\.\d+)",
)
minimap2 = Dependency(
    "minimap2",
    citation="Heng Li, Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34, Issue 18, 15 September 2018, Pages 3099–3103, https://doi.org/10.1093/bioinformatics/bty191",
    version_pattern=r"(\d+\.\d+)(?:-r\d+)?",
    min_version="2.17",
)
paftools = Dependency(
    "paftools.js",
    min_version="2.30",
    version_arg="version",
    version_pattern=r"(\d+\.\d+)(?:-r\d+)?",
)

# Calling
freebayes = Dependency(
    "freebayes",
    citation="Garrison, E. and Marth, G., 2012. Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907.",
    min_version="1.3.2",
)
bcftools = Dependency(
    "bcftools",
    citation="Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li, Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008",
    version_pattern=r"(\d+\.\d+)",
)

# Read cleaning
fastp = Dependency(
    "fastp",
    citation="Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu, fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560",
    min_version="0.20.0",
)

# Read statistics
seqkit = Dependency(
    "seqkit",
    citation="Wei Shen, Shuai Le, Yan Li, Fuquan Hu, SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation, PLoS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962",
    min_version="2.0.0",
    version_arg="version"
)

# Read downsampling
rasusa = Dependency(
    "rasusa",
    citation="Hall, M. B., (2022). Rasusa: Randomly subsample sequencing reads to a specified coverage. Journal of Open Source Software, 7(69), 3941, https://doi.org/10.21105/joss.03941",
    min_version="0.7.0",
    version_arg="--version"
)

# Masking and BED operations
bedtools = Dependency(
    "bedtools",
    citation="Aaron R. Quinlan, Ira M. Hall, BEDTools: a flexible suite of utilities for comparing genomic features, Bioinformatics, Volume 26, Issue 6, March 2010, Pages 841–842, https://doi.org/10.1093/bioinformatics/btq033",
    min_version="2.29.0",
    version_pattern=r"(\d+\.\d+\.\d+)",
)