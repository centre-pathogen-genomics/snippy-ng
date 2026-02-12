# Snippy-NG


[![CZI's Essential Open Source Software for Science](https://img.shields.io/badge/funded%20by-EOSS-FF414B)](https://czi.co/EOSS)
[![PyPI - Version](https://img.shields.io/pypi/v/snippy-nextgen.svg)](https://pypi.org/project/snippy-nextgen)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/snippy-nextgen.svg)](https://pypi.org/project/snippy-nextgen)
[![Test Coverage](https://raw.githubusercontent.com/centre-pathogen-genomics/snippy-ng/refs/heads/_xml_coverage_reports/data/tests/badge.svg)](https://app.codecov.io/github/centre-pathogen-genomics/snippy-ng)
![Downloads](https://img.shields.io/github/downloads/centre-pathogen-genomics/snippy-ng/total)
[![Benchmark](https://byob.yarr.is/centre-pathogen-genomics/snippy-ng/benchmark)](https://cpg.org.au/snippy-blog/posts/2024-10-10-snappy-snippy)
-----

🚨 Snippy-NG is under construction and should not replace Snippy 🚨
----

Check out our progress in the [Snippy-NG Development Blog](https://snippy.cpg.org.au/)!

**Table of Contents**

- [Installation](#installation)
- [License](#license)

## Installation

`snippy` is available on [PyPI](https://pypi.org/project/snippy-nextgen/) and can be installed using `pip` (without external dependencies):
```console
pip install snippy-nextgen
```

The complete snippy-ng environment ([including most dependencies](https://github.com/centre-pathogen-genomics/snippy-ng/issues/110#issuecomment-3857834560)) can be installed using the `snippy-ng` installer script. This script will install the latest version of `snippy-ng` and all its dependencies.

```console
curl -sSL https://cpg.org.au/snippy-ng/install.sh | bash -s -- --force
```

## Development

To set up a development environment, clone the repository and install `pixi` and `hatch`. Pixi is used to manage external dependencies, and Hatch is used to manage the Python package development.

```console
git clone git@github.com:centre-pathogen-genomics/snippy-ng.git && cd snippy-ng
# install pixi if not already installed
curl -fsSL https://pixi.sh/install.sh | bash
# install hatch if not already installed
pixi global install hatch
```

Activate the pixi environment and launch a hatch shell. This will install all dependencies and set up the development environment.

```console
pixi shell
hatch shell
```

```console
snippy-ng --help
```

## Examples 

```bash
snippy-ng asm --ref tests/data/wildtype.gbk --asm tests/data/wildtype.contigs.fa
```

```bash
snippy-ng short --ref tests/data/wildtype.gbk --R1 tests/data/mutant_R1.fastq --R2 tests/data/mutant_R2.fastq
```

```bash
snippy-ng long --clair3-model "/opt/models/r1041_e82_400bps_sup_v430_bacteria_finetuned" --clair3-fast-mode --ref tests/data/JKD6159.fasta --reads tests/data/JKD6159.fastq.gz
```

```bash
snippy-ng gather --ref tests/data/wildtype.gbk tests/data/{wildtype,mutant}* > tests/data/samples.csv 
snippy-ng multi tests/data/samples.csv \
  --ref tests/data/wildtype.gbk \
  -o multi \
  --cpus 6
```

## License

`snippy-ng` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.
