"""
Shared pipeline execution logic for Snippy CLI commands.
"""
from pathlib import Path
from snippy_ng.snippy import Snippy
from snippy_ng.exceptions import DependencyError, MissingOutputError, StageExecutionError


def run_snippy_pipeline(
    stages: list,
    skip_check: bool = False,
    check: bool = False,
    outdir: Path = Path("."),
    quiet: bool = False,
    continue_last_run: bool = False,
    keep_incomplete: bool = False
) -> int:
    """
    Common pipeline execution logic for all Snippy CLI commands.
    
    Args:
        stages: List of pipeline stages to execute
        skip_check: If True, skip dependency validation
        check: If True, only check dependencies and exit
        outdir: Output directory for pipeline results
        quiet: If True, suppress output
        continue_last_run: If True, continue from last run
        keep_incomplete: If True, keep incomplete results
    Returns:
        Exit code (0 for success, 1 for failure)
    """
    snippy = Snippy(stages=stages)
    snippy.welcome()

    if not skip_check:
        try:
            snippy.validate_dependencies()
        except DependencyError as e:
            snippy.error(f"Invalid dependencies! Please install '{e}' or use --skip-check to ignore.")
            return 1
    
    if check:
        return 0

    # Set working directory to output folder
    snippy.set_working_directory(outdir)
    try:
        snippy.run(
            quiet=quiet,
            continue_last_run=continue_last_run,
            keep_incomplete=keep_incomplete
        )
    except MissingOutputError as e:
        snippy.error(e)
        return 1
    except StageExecutionError as e:
        snippy.error(e)
        return 1
    
    snippy.cleanup()
    snippy.goodbye()
    
    return 0
