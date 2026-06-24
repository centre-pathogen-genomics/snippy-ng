from __future__ import annotations

from typing import Optional
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field

class Context(BaseModel):
    outdir: Optional[Path] = Field(default=None, description="Output directory for the pipeline run")
    log_path: Optional[Path] = Field(default=Path("LOG.txt"), description="Optional log file path for this pipeline run")
    tmpdir: Optional[Path] = Field(default=None, description="Temporary directory")
    cpus: int = Field(default=1, description="Max number of CPUs to use")
    ram: Optional[int] = Field(default=None, description="RAM in GB")
    quiet: bool = Field(default=False, description="Suppress command stdout/stderr output during stage execution")
    create_missing: bool = Field(default=False, description="Skip stages whose expected outputs already exist")
    keep_incomplete: bool = Field(default=False, description="Keep partially generated outputs if a stage fails")
    skip_check: bool = Field(default=False, description="Skip dependency validation before running the pipeline")
    check: bool = Field(default=False, description="Only validate dependencies and exit without executing stages")
    force: bool = Field(default=False, description="Allow overwriting existing output directories or files where applicable")
    debug: bool = Field(default=False, description="Enable debug-oriented runtime behavior and verbose diagnostics")
    break_points: bool = Field(default=False, description="Break at the start of each stage for debugging")
    no_cleanup: bool = Field(default=False, description="Disable cleanup of temporary/intermediate files after pipeline execution")

    model_config = ConfigDict(extra='forbid', arbitrary_types_allowed=True)
