---
title: Installation
---

# Installation

Snippy-NG is published on PyPI as `snippy-nextgen`. The command installed by
the package is `snippy-ng`.

For most users, the installer script is the best starting point because the SNP
calling pipelines depend on external bioinformatics tools as well as the Python
package.

## Recommended Install

Install the complete Snippy-NG environment with:

```console
curl -sSL https://cpg.org.au/snippy-ng/install.sh | bash -s -- --force
```

The published `docs/install.sh` shim downloads and runs the latest release
installer from GitHub:

```console
curl -sSL https://github.com/centre-pathogen-genomics/snippy-ng/releases/latest/download/install.sh | bash -s -- --force
```

After installation, check that the command is available:

```console
snippy-ng --help
```

## Python-Only Install

Install only the Python package with:

```console
pip install snippy-nextgen
```

This is enough to install the `snippy-ng` command and Python dependencies such
as Click, Pydantic, Biopython, and NumPy. It does not install the external
pipeline tools used for alignment, variant calling, masking, tree building, or
read processing.

Use this option when you already manage the external tools yourself, for example
with conda, micromamba, pixi, containers, or a shared HPC module system.

## GUI Install

The `snippy-ng gui` command uses Gradio, which is an optional pip extra. Install Snippy-NG with the GUI extra:

```console
pip install 'snippy-nextgen[gui]'
```

Or add Gradio to an existing Snippy-NG environment:

```console
pip install gradio
```

If Gradio is not installed, `snippy-ng gui` will print an install message and exit without affecting the rest of the command-line interface.

## External Tools

Snippy-NG checks stage dependencies before running a pipeline. The main external
tools used by the current pipelines are:

| Tool | Used for |
| --- | --- |
| `samtools` | BAM/SAM processing, sorting, indexing, reference indexing. |
| `bcftools` | VCF filtering, consensus generation, consequence annotation. |
| `freebayes` | Short-read and optional long-read variant calling. |
| `minimap2` | Short-read, long-read, and assembly alignment. |
| `bwa` | Optional short-read alignment. |
| `bedtools` | BED operations and FASTA masking. |
| `fastp` | Optional short-read cleaning. |
| `seqkit` | FASTA/FASTQ statistics and sequence utilities. |
| `rasusa` | Optional read downsampling. |
| `core-snp-filter` | Core SNP alignment filtering. |
| `iqtree` | Phylogenetic tree construction. |
| `longbow` | Clair3 model prediction for long-read runs. |
| `run_clair3.sh` | Clair3 long-read variant calling. |

The versions used by the development environment are recorded in
`pyproject.toml`.

## Checking An Install

Use `--check` on a pipeline command to validate dependencies without running the
analysis:

```console
snippy-ng short \
  --ref reference.gbk \
  --R1 reads_1.fastq.gz \
  --R2 reads_2.fastq.gz \
  --check
```

For long-read runs with the default Clair3 caller, also make sure a Clair3 model
can be found or pass one explicitly:

```console
snippy-ng long \
  --ref reference.fasta \
  --reads reads.fastq.gz \
  --clair3-model /path/to/clair3/model \
  --check
```

## Clair3

Long-read SNP calling uses Clair3 by default, but the Snippy-NG installer does
not include Clair3 because it is large. Use `--caller freebayes` for an
out-of-the-box long-read run, or install Clair3 separately.

See the [Clair3 installation guide](clair3.md) for conda, Docker, and manual
environment options.

## Development Install

For repository development, use `pixi` and `hatch`.

```console
git clone git@github.com:centre-pathogen-genomics/snippy-ng.git
cd snippy-ng
curl -fsSL https://pixi.sh/install.sh | bash
pixi global install hatch
pixi shell
hatch shell
```

Run the fast test suite with:

```console
pixi run test
```

Run simulation-backed integration tests with:

```console
pixi run -e integration test-integration
```
