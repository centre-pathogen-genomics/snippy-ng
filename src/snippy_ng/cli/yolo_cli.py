"""
Thanks for looking at the source code! You found a hidden command :D
"""

from pathlib import Path
from typing import Optional
import click


@click.command(context_settings={'show_default': True}, hidden=True)
@click.argument("directory", required=False, type=click.Path(exists=True, resolve_path=True, readable=True))
def yolo(directory: Optional[Path]):
    """
    Pipeline that automates everything.
    
    Not recommended for general use unless you've got no idea what you're doing.
    """
    from snippy_ng.logging import logger
    logger.warning("You are running the YOLO pipeline. This pipeline is not recommended for general use unless you have no idea what you're doing. Please consider using one of the other pipelines with more specific parameters for better results and more control over the analysis.")
    raise NotImplementedError("YOLO pipeline is not implemented yet. This is a placeholder for a future pipeline that will attempt to automate the entire process from raw data to tree. Stay tuned!")
    
