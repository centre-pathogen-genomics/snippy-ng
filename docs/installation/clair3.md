---
title: Clair3
---

# Clair3

The Snippy-NG installer does not include Clair3 by default because it is large.
The Clair3 environment and models are roughly 800 MB, so long-read users need to
choose how they want to provide it.

Snippy-NG uses Clair3 by calling `run_clair3.sh`. Any installation method is
fine as long as a working `run_clair3.sh` command is available on `PATH`.

## Option 1: Use FreeBayes

This is the simplest option and should run out of the box with the Snippy-NG
installer:

```console
snippy-ng long --ref reference.fasta --reads reads.fastq.gz --caller freebayes
```

This is not the preferred long-read caller, but it avoids the Clair3 dependency.

## Option 2: Install Clair3 Separately

Create a conda environment with Clair3 installed:

```console
conda create -n clair3 -c bioconda clair3
```

Add a wrapper called `run_clair3.sh` to a directory on your `PATH`. This example
assumes `~/.local/bin` is already on your `PATH`:

```console
cat << 'EOF' > ~/.local/bin/run_clair3.sh
#!/usr/bin/env bash

conda run -n clair3 run_clair3.sh "$@"
EOF

chmod u+x ~/.local/bin/run_clair3.sh
```

Check that Snippy-NG can find it:

```console
which run_clair3.sh
run_clair3.sh --help
```

Then run Snippy-NG with the default long-read caller:

```console
snippy-ng long --ref reference.fasta --reads reads.fastq.gz
```

## Option 3: Use Docker

The repository includes a [Docker wrapper](run_clair3_docker.sh) for Clair3:

```console
install -m 755 docs/installation/run_clair3_docker.sh ~/.local/bin/run_clair3.sh
```

The wrapper runs the Clair3 container and mounts directories for common Clair3
input and output arguments. By default it uses:

```text
hkubal/clair3:latest
```

Override the image if needed:

```console
export CLAIR3_DOCKER_IMAGE=hkubal/clair3:latest
```

On some platforms, you may also need to set a Docker platform:

```console
export CLAIR3_DOCKER_PLATFORM=linux/amd64
```

## Option 4: Install Snippy-NG With Clair3

You can manually create an environment that contains both Clair3 and the
Snippy-NG Python package:

```console
conda create -n clair3 -c bioconda clair3 pip
conda activate clair3
pip install snippy-nextgen
```

Then install any other external tools you need into the same environment. For
example:

```console
conda install -c bioconda minimap2 samtools bcftools bedtools freebayes
```

This works, but it is more manual than the Snippy-NG installer because you are
responsible for filling in missing pipeline dependencies.

## Clair3 Models

Snippy-NG can use Longbow to predict a Clair3 model for read input, but the
model files still need to exist somewhere accessible.

Set a model root with either:

```console
export CLAIR3_MODELS=/path/to/clair3_models
export CLAIR3_MODEL_ROOT=/path/to/clair3_models
```

Or pass a model directly:

```console
snippy-ng long \
  --ref reference.fasta \
  --reads reads.fastq.gz \
  --clair3-model /path/to/model
```

When using BAM or CRAM input with Clair3, pass `--clair3-model` explicitly:

```console
snippy-ng long \
  --ref reference.fasta \
  --bam reads.bam \
  --clair3-model /path/to/model
```

## Future Conda Package

Once Snippy-NG is available as a conda package, the intended install path will be
closer to:

```console
conda create -n snippy-ng -c bioconda snippy-ng clair3
```

That package is not available yet.
