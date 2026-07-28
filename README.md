# Snippy-NG

[![CZI's Essential Open Source Software for Science](https://img.shields.io/badge/funded%20by-EOSS-FF414B)](https://czi.co/EOSS)
[![PyPI - Version](https://img.shields.io/pypi/v/snippy-nextgen.svg)](https://pypi.org/project/snippy-nextgen)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/snippy-nextgen.svg)](https://pypi.org/project/snippy-nextgen)
[![Test Coverage](https://raw.githubusercontent.com/centre-pathogen-genomics/snippy-ng/refs/heads/_xml_coverage_reports/data/tests/badge.svg)](https://app.codecov.io/github/centre-pathogen-genomics/snippy-ng)
![Downloads](https://img.shields.io/github/downloads/centre-pathogen-genomics/snippy-ng/total)
[![Benchmark](https://byob.yarr.is/centre-pathogen-genomics/snippy-ng/benchmark)](https://cpg.org.au/snippy-blog/posts/2024-10-10-snappy-snippy)

## Synopsis

🚨 This software is still under development 🚨

`snippy-ng` is a modern version of the classic `snippy` tool.
It finds variants in haploid microbial
genomes relative to a reference genome.
It can accept short reads, long reads, or
assembled genomes.
It can also combine the variants from multiple
samples to make a multipe sequence alignment
and a phylogenomic tree.
Final results can be viewed from an interactive
HTML report.

### Documentation

* [Snippy-NG Manual](https://cpg.org.au/snippy-ng/) 
* [Snippy-NG Development Blog](https://snippy.cpg.org.au/)!

## Installation

### Pixi
```console
curl -sSL https://cpg.org.au/snippy-ng/install.sh | bash -s -- --force
```

### Conda
Will be available on first official release.

## Quick Start

```bash
snippy-ng short --ref tests/data/reference.gbk --R1 tests/data/mutant_R1.fastq.gz --R2 tests/data/mutant_R2.fastq.gz
```

```bash
snippy-ng asm --ref tests/data/reference.gbk tests/data/wildtype.contigs.fa
```

```bash
export CLAIR3_MODELS=./clair3_models # try to find appropriate models in this directory
snippy-ng long --ref tests/data/JKD6159.fasta tests/data/JKD6159.fastq.gz
```

```bash
snippy-ng utils gather --json --ref tests/data/reference.gbk tests/data/{wildtype,mutant}* > samples.json 
snippy-ng multi samples.json --cpus 6 -o multi
snippy-ng tree --fast multi/core/core.095.aln -o multi/tree
snippy-ng utils report tree multi/tree/tree.snps.newick --metadata multi/snippy.vcf.summary.tsv -o multi/report
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

## Testing

Run the default fast test suite with:

```console
pixi run test
```

Run the simulation-backed integration tests with:

```console
pixi run -e integration test-integration
```

Simulated inputs are generated at test time and cached under `.cache/integration-sim/`.
See [docs/integration-tests.md](/Users/wwirth/programming/snippy-ng/docs/integration-tests.md) for the scenario format and harness contract.

## License

`snippy-ng` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.

## Funding

This work was funded by the
CZI EOSS scheme from 2024-2026.

## Authors

* Wytammma Wirth
* Torsten Seemann


