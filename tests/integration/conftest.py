from __future__ import annotations

from pathlib import Path
import subprocess
import sys
from collections.abc import Iterable

import pytest

from tests.integration.simulation import (
    DEFAULT_CACHE_ROOT,
    DEFAULT_REFERENCE,
    PROJECT_ROOT,
    IntegrationDataset,
    SimulationRequest,
    VariantRecord,
    materialize_scenario,
)


@pytest.fixture(scope="session")
def integration_cache_root() -> Path:
    DEFAULT_CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    return DEFAULT_CACHE_ROOT


@pytest.fixture
def simulated_dataset(tmp_path: Path, integration_cache_root: Path):
    def factory(
        variants: Iterable[VariantRecord],
        input_type: str,
        *,
        name: str,
        reference: Path = DEFAULT_REFERENCE,
        untouched_regions: tuple[tuple[str, int, int], ...] = (),
    ) -> IntegrationDataset:
        request = SimulationRequest(
            name=name,
            reference=reference,
            truth_variants=tuple(variants),
            untouched_regions=untouched_regions,
        )
        materialized = materialize_scenario(
            request=request,
            input_type=input_type,
            cache_root=integration_cache_root,
        )
        outdir = tmp_path / f"{request.name}-{input_type}"
        prefix = "snippy"

        if input_type == "short":
            assert materialized.reads_r1 is not None
            assert materialized.reads_r2 is not None
            args = [
                "short",
                "--reference", str(request.reference),
                "--R1", str(materialized.reads_r1),
                "--R2", str(materialized.reads_r2),
                "--outdir", str(outdir),
                "--prefix", prefix,
                "--skip-check",
                "--cpus", "1",
                "--ram", "4",
                "--min-qual", "0",
                "--depth-mask", "1",
            ]
        elif input_type == "long":
            assert materialized.long_reads is not None
            args = [
                "long",
                "--reference", str(request.reference),
                "--reads", str(materialized.long_reads),
                "--outdir", str(outdir),
                "--prefix", prefix,
                "--skip-check",
                "--cpus", "1",
                "--ram", "4",
                "--caller", "freebayes",
                "--min-qual", "0",
                "--depth-mask", "1",
                "--min-read-len", "100",
                "--min-read-qual", "1",
            ]
        else:
            assert materialized.assembly is not None
            args = [
                "asm",
                "--reference", str(request.reference),
                "--assembly", str(materialized.assembly),
                "--outdir", str(outdir),
                "--prefix", prefix,
                "--skip-check",
                "--cpus", "1",
                "--ram", "4",
            ]

        result = subprocess.run(
            [sys.executable, "-m", "snippy_ng", *args],
            check=False,
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0, result.stdout + "\n" + result.stderr

        return IntegrationDataset(
            request=request,
            input_type=input_type,
            reference=request.reference,
            truth_vcf=materialized.truth_vcf,
            mutated_reference=materialized.mutated_reference,
            outdir=outdir,
            called_vcf=outdir / f"{prefix}.vcf",
            assembly=materialized.assembly,
            reads_r1=materialized.reads_r1,
            reads_r2=materialized.reads_r2,
            long_reads=materialized.long_reads,
        )

    return factory
