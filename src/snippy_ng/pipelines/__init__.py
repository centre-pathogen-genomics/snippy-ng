import random
from typing import List
import os
from pathlib import Path
import time

from snippy_ng.exceptions import DependencyError
from snippy_ng.logging import logger
from snippy_ng.__about__ import __version__, DOCS_URL, GITHUB_URL
from snippy_ng.stages import BaseStage, Context
from pydantic import BaseModel, ConfigDict, Field


class PipelineBuilder(BaseModel):
    """Base class for building Snippy pipelines."""
    model_config = ConfigDict(extra='forbid')
    cpus: int = Field(default=1, description="Number of CPUs to use")
    ram: int = Field(default=8, description="RAM in GB")
    prefix: str = Field(default="tree", description="Output file prefix")

    def build(self) -> 'SnippyPipeline':
        """Build and return the SnippyPipeline."""
        raise NotImplementedError("Subclasses must implement the build method.")


class SnippyPipeline:
    """
    Main class for creating Snippy-NG Pipelines.
    """
    stages: List[BaseStage]

    def __init__(self, stages: List[BaseStage] = None):
        if stages is None:
            stages = []
        self.stages = stages

    def run(self, outdir, tmpdir, cpus, ram, quiet=False, create_missing=False, keep_incomplete=False, skip_check=False, check=False, force=False, debug=False, no_cleanup=False):
        self.welcome()

        if not skip_check:
            try:
                self.validate_dependencies()
            except DependencyError as e:
                raise DependencyError(f"Invalid dependencies! Please install '{e}' or use --skip-check to ignore.") from e
        
        if check:
            return None

        # Set working directory to output folder
        current_dir = Path.cwd()
        run_ctx = Context(
            outdir=outdir,
            tmpdir=tmpdir,
            cpus=cpus,
            ram=ram,
            quiet_mode=quiet,
            create_missing=create_missing,
            keep_incomplete=keep_incomplete
        )
        try:
            self.set_working_directory(outdir)
            self._execute_pipeline_stages_in_order(
                run_ctx,
            )
        finally:
            self.cleanup(current_dir)
        self.goodbye()
    
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
            self.debug(f"Checking dependencies for {stage.name}...")
            for dependency in stage._dependencies:
                if dependency.name in checked:
                    self.debug(f"Skipping {dependency.name}, already checked.")
                    continue
                checked.add(dependency.name)
                try:
                    version = dependency.check()
                    self.log(f"{dependency.name} v{version}")
                except DependencyError as e:
                    # Capture general dependency error
                    self.error(f"{e}")
                    invalid.append(dependency)
        if invalid:
            raise DependencyError(f"{', '.join([d.format_version_requirements() for d in invalid])}")
        
    def set_working_directory(self, directory):
        # Set the working directory
        self.hr("SetWorkingDirectory")
        self.log(f"Setting working directory to '{directory}'")
        os.chdir(directory)

    def _execute_pipeline_stages_in_order(self, run_ctx: Context):
        if not self.stages:
            raise ValueError("No stages to run in the pipeline!")
        # Run pipeline sequentially
        self.start_time = time.perf_counter()
        for stage in self.stages:
            self.hr(f"{stage.name}")
            self.debug(stage)
            start = time.perf_counter()
            stage.run(run_ctx)
            end = time.perf_counter()
            self.debug(f"Runtime: {(end - start):.2f} seconds")
            # After running each stage,
            # check all expected outputs were produced
            stage.error_if_outputs_missing()
            # run any tests defined for the stage
            stage.run_tests()
        self.end_time = time.perf_counter()

    
    def cleanup(self, directory: Path = None):
        # reset working directory
        if directory:
            os.chdir(directory)
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
        self.hr(f"Running Snippy-NG v{__version__}", style=" ", color="green")
        self.hr()

        self.log("Stages:")
        for i, stage in enumerate(self.stages, 1):
            self.log(f"{i:3}. {stage.name}")


    def goodbye(self):
        messages = [
            "May the SNPs be with you.",
            "Wishing you a life free of homopolymer errors.",
            f"Found a bug? Post it at {GITHUB_URL}/issues",
            f"Have a suggestion? Tell me at {GITHUB_URL}/issues",
            f"The Snippy manual is at {GITHUB_URL}/blob/master/README.md",
            "Questionable SNP? Try the --report option to see the alignments.",
            "Did you know? Snippy is a combination of SNP, Skippy, and snappy.",
            "Set phasers to align… your reads.",
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
        total_run_time = self.end_time - self.start_time
        self.echo(f"Total runtime: {total_run_time:.2f} seconds")
        self.echo(f"Documentation: {DOCS_URL}")
        self.echo(f"GitHub: {GITHUB_URL}")
        self.hr("Citations")
        self.echo('\n'.join(f"{i}. {cite}" for i, cite in enumerate(self.citations, 1)))
        self.hr()
        # Print a random goodbye message
        self.hr(f"{random.choice(messages)}", style=" ", color='green')
