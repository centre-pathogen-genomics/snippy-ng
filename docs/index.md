---
title: Home
---

# Snippy-NG

[![CZI's Essential Open Source Software for Science](https://img.shields.io/badge/funded%20by-EOSS-FF414B)](https://czi.co/EOSS)
[![PyPI - Version](https://img.shields.io/pypi/v/snippy-nextgen.svg)](https://pypi.org/project/snippy-nextgen)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/snippy-nextgen.svg)](https://pypi.org/project/snippy-nextgen)
[![Test Coverage](https://raw.githubusercontent.com/centre-pathogen-genomics/snippy-ng/refs/heads/_xml_coverage_reports/data/tests/badge.svg)](https://app.codecov.io/github/centre-pathogen-genomics/snippy-ng)
[![Benchmark](https://byob.yarr.is/centre-pathogen-genomics/snippy-ng/benchmark)](https://cpg.org.au/snippy-blog/posts/2024-10-10-snappy-snippy)

Snippy-NG is the next generation of microbial variant calling with Snippy.
It provides command-line pipelines for assembly, short-read, long-read, and
multi-sample SNP analysis, plus tools for core alignment construction,
phylogenetic tree building, and interactive report generation.

Snippy-NG is under active development and should not replace Snippy for
production analyses without validation.

Follow progress on the [Snippy-NG Development Blog](https://snippy.cpg.org.au/).

## Installation

For most users, install the complete Snippy-NG environment with:

```console
curl -sSL https://cpg.org.au/snippy-ng/install.sh | bash -s -- --force
```

The Python package is also available from PyPI:

```console
pip install snippy-nextgen
```

The installed command is:

```console
snippy-ng --help
```

See the [installation guide](installation/index.md) for dependency details,
[Clair3 setup](installation/clair3.md), and development installs.

## Commands

| Command | Purpose |
| --- | --- |
| `snippy-ng asm` | Assembly-based SNP calling from a reference and assembly FASTA. |
| `snippy-ng short` | Short-read SNP calling from paired-end reads or an existing BAM. |
| `snippy-ng long` | Long-read SNP calling from FASTQ reads or an existing BAM. |
| `snippy-ng multi` | Run multi-sample analysis and construct a core alignment. |
| `snippy-ng core` | Create a core alignment from multiple Snippy-NG run directories. |
| `snippy-ng tree` | Build a phylogenetic tree from an alignment. |
| `snippy-ng gui` | Launch the optional Gradio graphical interface. |
| `snippy-ng utils ref` | Prepare a reference genome for Snippy-NG. |
| `snippy-ng utils gather` | Discover sample files and write TSV or JSON input for `multi`. |
| `snippy-ng utils samclip` | Filter clipped reads from SAM alignments. |
| `snippy-ng utils report-tree` | Create an interactive HTML report from a Newick tree. |

## Examples

Run assembly-based SNP calling:

```console
snippy-ng asm --ref tests/data/reference.gbk --asm tests/data/wildtype.contigs.fa
```

Run short-read SNP calling:

```console
snippy-ng short \
  --ref tests/data/reference.gbk \
  --R1 tests/data/mutant_R1.fastq.gz \
  --R2 tests/data/mutant_R2.fastq.gz
```

Run long-read SNP calling:

```console
export CLAIR3_MODELS=./clair3_models
snippy-ng long --ref tests/data/JKD6159.fasta --reads tests/data/JKD6159.fastq.gz
```

Run a multi-sample workflow and build a report:

```console
snippy-ng utils gather --json --ref tests/data/reference.gbk tests/data/{wildtype,mutant}* > samples.json
snippy-ng multi samples.json --cpus 6 -o multi
snippy-ng tree --fast multi/core/core.095.aln -o multi/tree
snippy-ng utils report-tree multi/tree/tree.treefile --metadata multi/snippy.vcf.summary.tsv -o multi/report
```

## Documentation

- [Installation](installation/index.md)
- [Clair3 Setup](installation/clair3.md)
- [Interactive Report](report/index.md)
- [Development Guide](development/index.md)
- [CLI Development](development/cli.md)
- [Pipeline Development](development/pipelines.md)
- [Stage Development](development/stages.md)
- [Integration Tests](development/integration-tests.md)

## Development

To set up a development environment, clone the repository and install `pixi` and
`hatch`.

```console
git clone git@github.com:centre-pathogen-genomics/snippy-ng.git
cd snippy-ng
curl -fsSL https://pixi.sh/install.sh | bash
pixi global install hatch
```

Activate the pixi environment and launch a hatch shell:

```console
pixi shell
hatch shell
```

Run the default fast test suite:

```console
pixi run test
```

Run simulation-backed integration tests:

```console
pixi run -e integration test-integration
```

## License

Snippy-NG is distributed under the terms of the
[MIT license](https://spdx.org/licenses/MIT.html).
