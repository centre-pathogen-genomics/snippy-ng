---
title: Clair3 models
---

# Clair3 model utilities

Download individual Clair3 models without installing the full model bundle.

```console
snippy-ng utils model download r1041_e82_400bps_hac_v500
```

The command accepts either a Clair3 model directory name or a basecall model id:

```console
snippy-ng utils model download dna_r10.4.1_e8.2_400bps_hac@v5.0.0
```

Both download `r1041_e82_400bps_hac_v500`.

## Check which model would be used

```console
snippy-ng utils model check reads.fastq.gz
```

The command reads the first read header and resolves the matching model from the
configured model roots. If the header has no `basecall_model_version_id` or
Dorado `RG:Z` tag, it runs Longbow and resolves the predicted model instead.
URL reads are downloaded to a temporary file first.

```console
snippy-ng utils model check https://example.org/reads.fastq.gz --outdir check
```

`--outdir` (default: current directory) stores the Longbow prediction
(`PREFIX.longbow.json`) and selected model (`PREFIX.clair3_model.txt`).
`--cpus` sets the Longbow thread count. SRA accessions are not supported by
this command because the read file must be available to inspect its header.

## Destination

Models are written to:

1. `--outdir` / `-o` when provided.
2. `$CLAIR3_MODELS` when set.
3. `$CLAIR3_MODEL_ROOT` when set.
4. `./clair3_models` otherwise.

```console
snippy-ng utils model download r1041_e82_400bps_sup_v520 -o /data/clair3_models
```

Use `--force` to replace an existing model directory.

## Model source

Models are downloaded from the Clair3 PyTorch model index at
`https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/`. Each model directory
contains `pileup.pt` and `full_alignment.pt`.

If a model name is not present in that index, the command exits with an error.
The index is browsable at
<https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/>.
