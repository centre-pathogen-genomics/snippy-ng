"""
Shared pipeline execution logic for Snippy CLI commands.
"""
from pathlib import Path
from snippy_ng.snippy import Snippy
from snippy_ng.exceptions import DependencyError, SnippyError


def run_snippy_pipeline(
    stages: list,
    skip_check: bool = False,
    check: bool = False,
    outdir: Path = Path("."),
    quiet: bool = False,
    create_missing: bool = False,
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
        create_missing: If True, continue from last run by creating missing outputs
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
    current_dir = Path.cwd()
    snippy.set_working_directory(outdir)
    try:
        snippy.run(
            quiet=quiet,
            create_missing=create_missing,
            keep_incomplete=keep_incomplete
        )
    except SnippyError as e:
        snippy.error(e)
        return 1
    
    snippy.cleanup(current_dir)
    snippy.goodbye()
    
    return 0
