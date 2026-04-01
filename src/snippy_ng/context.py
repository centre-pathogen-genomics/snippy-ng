from __future__ import annotations

from typing import Optional
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field

class Context(BaseModel):
    outdir: Optional[Path] = Field(None, description="Output directory for the pipeline run")
    tmpdir: Optional[Path] = Field(default=None, description="Temporary directory")
    cpus: int = Field(default=1, description="Max number of CPUs to use")
    ram: Optional[int] = Field(default=None, description="RAM in GB")
    quiet: bool = Field(False, description="Suppress command stdout/stderr output during stage execution")
    create_missing: bool = Field(False, description="Skip stages whose expected outputs already exist")
    keep_incomplete: bool = Field(False, description="Keep partially generated outputs if a stage fails")
    skip_check: bool = Field(False, description="Skip dependency validation before running the pipeline")
    check: bool = Field(False, description="Only validate dependencies and exit without executing stages")
    force: bool = Field(False, description="Allow overwriting existing output directories or files where applicable")
    debug: bool = Field(False, description="Enable debug-oriented runtime behavior and verbose diagnostics")
    break_points: bool = Field(False, description="Break at the start of each stage for debugging")
    no_cleanup: bool = Field(False, description="Disable cleanup of temporary/intermediate files after pipeline execution")

    model_config = ConfigDict(extra='forbid', arbitrary_types_allowed=True)