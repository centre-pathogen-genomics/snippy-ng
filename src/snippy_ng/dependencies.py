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
        ).stdout

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
vt = Dependency(
    "vt",
    citation="Adrian Tan, Gonçalo R. Abecasis, Hyun Min Kang, Unified representation of genetic variants, Bioinformatics, Volume 31, Issue 13, July 2015, Pages 2202–2204, https://doi.org/10.1093/bioinformatics/btv112",
    version_pattern=r"v(\d+\.\d+)",
)
