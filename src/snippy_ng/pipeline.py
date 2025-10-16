import random
from typing import List
import os
from pathlib import Path
import time

from snippy_ng.exceptions import DependencyError, SkipStageError, MissingOutputError
from snippy_ng.logging import logger
from snippy_ng.__about__ import __version__, DOCS_URL, GITHUB_URL
from snippy_ng.stages.base import BaseStage

class Pipeline:
    """
    Main class for creating Snippy-NG Pipelines.
    """
    stages: List[BaseStage]

    def __init__(self, stages: List[BaseStage] = None):
        if stages is None:
            stages = []
        self.stages = stages
    
    def add_stage(self, stage):
        self.stages.append(stage)

    @property
    def dependencies(self):
        return [dep for stage in self.stages for dep in stage._dependencies]

    def hr(self, msg="", style='-', color='blue'):
        """Print a horizontal rule."""
        logger.horizontal_rule(msg, style=style, color=color)

    def echo(self, message):
        logger.echo(message, err=True)

    def log(self, msg):
        logger.info(msg)
    
    def warning(self, msg):
        logger.warning(msg)
    
    def debug(self, msg):
        logger.debug(msg)

    def error(self, msg):
        logger.error(msg)

    def validate_dependencies(self):
        invalid = []
        checked = set()
        self.hr("CheckDependencies")
        for stage in self.stages:
            self.log(f"Checking dependencies for {stage.name}...")
            for dependency in stage._dependencies:
                if dependency.name in checked:
                    self.log(f"Skipping {dependency.name}, already checked.")
                    continue
                checked.add(dependency.name)
                try:
                    version = dependency.check()
                    self.log(f"Found {dependency.name} v{version}")
                except DependencyError as e:
                    # Capture general dependency error
                    self.error(f"{e}")
                    invalid.append(dependency)
        if invalid:
            raise DependencyError(f"{', '.join([d.format_version_requirements() for d in invalid])}")
        self.log("Dependencies look good!")
        
    def set_working_directory(self, directory):
        # Set the working directory
        self.hr("SetWorkingDirectory")
        self.log(f"Setting working directory to '{directory}'")
        os.chdir(directory)

    def run(self, quiet=False, continue_last_run=False, keep_incomplete=False):
        # Run pipeline sequentially
        self.start_time = time.perf_counter()
        for stage in self.stages:
            self.hr(f"{stage.name}")
            self.debug(stage)
            try:
                try:
                    if continue_last_run and stage.check_outputs():
                        self.log(f"{stage.name} already completed, skipping...")
                        continue
                    stage.run(quiet)
                except (RuntimeError, KeyboardInterrupt) as e:
                    # remove outputs if stage fails
                    if not keep_incomplete:
                        for name, path in stage.output:
                            if path and Path(path).exists():
                                self.warning(f"Removing incomplete output '{name}' ({path}) due to error.")
                                Path(path).unlink()
                    raise e
                for name, path in stage.output:
                    if not path:  # Skip empty outputs
                        continue
                    if Path(path).exists():
                        continue
                    raise MissingOutputError(f"Expected output '{name}' ({path}) not found after running '{stage.name}'")
            except SkipStageError:
                self.stages.remove(stage)
                self.warning(f"STAGE {stage.name} SKIPPED!")
        self.end_time = time.perf_counter()

    def cleanup(self):
        # Clean up unnecessary files
        pass
    
    @property
    def citations(self):
        citations = []
        for stage in self.stages:
            for dependency in stage._dependencies:
                if dependency.citation:
                    citations.append(dependency.citation)
        return sorted(set(citations))

    def welcome(self):
        self.hr()
        self.hr("Running Snippy-NG", style=" ", color="green")
        self.hr()
        self.log(f"Version: {__version__}")

        self.log("Stages:")
        for i, stage in enumerate(self.stages, 1):
            self.log(f"  {i}. {stage.name}")


    def goodbye(self):
        messages = [
            "May the SNPs be with you.",
            "Wishing you a life free of homopolymer errors.",
            f"Found a bug? Post it at {GITHUB_URL}/issues",
            f"Have a suggestion? Tell me at {GITHUB_URL}/issues",
            f"The Snippy manual is at {GITHUB_URL}/blob/master/README.md",
            "Questionable SNP? Try the --report option to see the alignments.",
            "Did you know? Snippy is a combination of SNP, Skippy, and snappy.",
            "Set phasers to align… your SNPs.",
            "To boldly SNP where no one has SNPped before.",
            "Resistance to accurate SNP calling is futile.",
            "Wishing you a genome that's logically consistent…",
            "The final frontier of variant calling.",
            f"Make it so: Report your issues at {GITHUB_URL}/issues.",
            "Highly logical and warp-speed fast.",
            "Live long and SNP accurately.",
            f"Looking for guidance? The Snippy manual is at {DOCS_URL}.",
            "Beam me up, Snippy! Your SNPs are ready.",
            "SNP analysis at warp factor 9!",
            "Keep calm and trust Snippy to get your variants right.",
            "By Grabthar's Hammer… oh wait, wrong reference. Check your SNPs!",
            "Assimilate accurate SNP data with Snippy.",
            "There's no such thing as a no-win SNP scenario with Snippy.",
            "Do your SNPs feel out of phase? Realign them with Snippy.",
            "Snippy: The only logical choice for variant detection.",
            "SNPs detected, Captain! Ready for the next mission.",
        ]
        self.hr()
        self.echo('')
        self.hr("Snippy-NG completed!", style=" ", color="green")
        self.echo('')
        self.hr()
        self.echo(f"Total runtime: {self.end_time - self.start_time:.2f} seconds")
        self.echo(f"Documentation: {DOCS_URL}")
        self.echo(f"GitHub: {GITHUB_URL}")
        self.hr("Citations")
        self.echo('\n'.join(f"{i}. {cite}" for i, cite in enumerate(self.citations, 1)))
        self.hr()
        # Print a random goodbye message
        self.hr(f"{random.choice(messages)}", style=" ", color=None)
