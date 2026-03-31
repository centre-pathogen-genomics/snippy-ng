---
title: Integration Tests
---
# Simulation-based integration tests

Snippy-NG has a simulation-backed integration test layer under [`tests/integration`](/Users/wwirth/programming/snippy-ng/tests/integration). These tests build synthetic datasets at test time, run the real pipeline, and check that the final VCF recovers the variants that were injected into the reference.

## Test contract

Each scenario defines a small set of truth variants as Python `VariantRecord`s. The harness then:

1. writes a truth VCF from those records
2. applies the variants to the reference FASTA
3. generates the requested input type
4. runs `snippy-ng`
5. compares the final VCF against the ground-truth variants

For example, this scenario:

```python
{
    "name": "short_indel",
    "input_type": "short",
    "variants": (
        VariantRecord("Wildtype", 70, "CA", "C"),
    ),
    "untouched_regions": (("Wildtype", 200, 260),),
    "strict_region": ("Wildtype", 50, 180),
}
```

means:

- apply `Wildtype 70 CA C` to the reference
- simulate short reads from the mutated reference
- run `snippy-ng short` on those reads
- assert that the final VCF contains the expected deletion

This keeps each test definition small and makes the recovery requirement explicit.

## Scenario shape

The integration tests in [`test_simulated_pipelines.py`](/Users/wwirth/programming/snippy-ng/tests/integration/test_simulated_pipelines.py) are written as a simple parametrized config:

```python
SCENARIOS = [
    {
        "name": "short_snp",
        "input_type": "short",
        "variants": (
            VariantRecord("Wildtype", 20, "A", "C"),
        ),
        "untouched_regions": (("Wildtype", 200, 260),),
        "strict_region": None,
    },
]
```

Fields:

- `name`: scenario name used in test ids and cache paths
- `input_type`: one of `short`, `long`, or `asm`
- `variants`: truth variants to inject into the reference
- `untouched_regions`: regions that should stay free of called variants
- `strict_region`: optional region where only the expected truth variants are allowed

## Harness API

The main fixture is [`simulated_dataset`](/Users/wwirth/programming/snippy-ng/tests/integration/conftest.py). Each test passes:

- an `Iterable[VariantRecord]`
- an input type
- a scenario name
- optional untouched regions

Example:

```python
dataset = simulated_dataset(
    (
        VariantRecord("Wildtype", 70, "CA", "C"),
    ),
    "short",
    name="short_indel",
    untouched_regions=(("Wildtype", 200, 260),),
)
```

The returned `IntegrationDataset` exposes:

- the generated truth VCF
- the mutated reference
- the simulated reads or assembly
- the output directory
- the final called VCF
- assertion helpers for expected and unexpected variants

## Input generation

The simulation code lives in [`simulation.py`](/Users/wwirth/programming/snippy-ng/tests/integration/simulation.py).

- `short`: applies truth variants, then generates Illumina-style reads with `art_illumina`
- `long`: applies truth variants, then generates long reads with `badread`
- `asm`: applies truth variants and uses the mutated FASTA directly as the assembly input

The current scope is deliberately small:

- SNPs
- small insertions
- small deletions

The helper rejects unsupported symbolic alleles and overlapping truth edits.

## Caching

Generated artifacts are not checked into git. They are created on demand and cached under [`.cache/integration-sim`](/Users/wwirth/programming/snippy-ng/.cache/integration-sim).

Cache keys include:

- scenario name
- input type
- reference checksum
- truth variant content
- simulation parameters
- simulator versions for `art_illumina` and `badread`

If a matching cache entry already exists, the harness reuses the generated inputs instead of regenerating them.

## Running the tests

Run the integration suite with:

```console
pixi run -e integration hatch run pytest tests/integration
```

Run only the simulation-marked tests with:

```console
pixi run -e integration hatch run pytest -m integration_sim tests/integration
```
