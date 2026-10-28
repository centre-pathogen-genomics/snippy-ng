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
