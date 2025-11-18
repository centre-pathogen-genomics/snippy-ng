"""
Shared pipeline execution logic for Snippy CLI commands.
"""
from snippy_ng.snippy import Snippy
from snippy_ng.exceptions import DependencyError, MissingOutputError


def run_snippy_pipeline(config: dict, stages: list) -> int:
    """
    Common pipeline execution logic for all Snippy CLI commands.
    
    Args:
        config: Configuration dictionary from CLI options
        stages: List of pipeline stages to execute
        
    Returns:
        Exit code (0 for success, 1 for failure)
    """
    snippy = Snippy(stages=stages)
    snippy.welcome()

    if not config.get("skip_check", False):
        try:
            snippy.validate_dependencies()
        except DependencyError as e:
            snippy.error(f"Invalid dependencies! Please install '{e}' or use --skip-check to ignore.")
            return 1
    
    if config["check"]:
        return 0

    # Set working directory to output folder
    snippy.set_working_directory(config["outdir"])
    try:
        snippy.run(
            quiet=config["quiet"],
            continue_last_run=config["continue_last_run"],
            keep_incomplete=config["keep_incomplete"]
        )
    except MissingOutputError as e:
        snippy.error(e)
        return 1
    except RuntimeError as e:
        snippy.error(e)
        return 1
    
    snippy.cleanup()
    snippy.goodbye()
    
    return 0
