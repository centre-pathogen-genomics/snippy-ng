import random
from collections import OrderedDict
from typing import Dict, List, Optional
import os
from pathlib import Path
import time

from snippy_ng.exceptions import DependencyError, PipelineExecutionError
from snippy_ng.logging import logger
from snippy_ng.__about__ import __version__, DOCS_URL, GITHUB_URL
from snippy_ng.stages import BaseStage, Context
from snippy_ng.utils.files import human_readable_size
from pydantic import BaseModel, ConfigDict, Field


class PipelineBuilder(BaseModel):
    """Base class for building Snippy pipelines."""
    model_config = ConfigDict(extra='forbid')
    prefix: str = Field(default="snps", description="Output file prefix")

    def build(self) -> 'SnippyPipeline':
        """Build and return the SnippyPipeline."""
        raise NotImplementedError("Subclasses must implement the build method.")


class SnippyPipeline:
    """
    Main class for creating Snippy-NG Pipelines.

    If outputs_to_keep is not specified, all stage outputs will be kept by default. If outputs_to_keep is provided, only the specified outputs will be kept and all other stage outputs will be deleted during cleanup (unless --no-cleanup is used).
    """
    stages: List[BaseStage]
    outputs_to_keep: List[Path]
    outdir: Optional[Path] = None

    def __init__(self, stages: List[BaseStage] = None, outputs_to_keep: List[Path] = None):
        self.stages = stages or []
        if outputs_to_keep is None:
            # keep all outputs by default if not specified
            outputs_to_keep = []
            for stage in self.stages:
                for _, output in stage.output.non_temporary_outputs():
                    if isinstance(output, Path):
                        outputs_to_keep.append(output)
        self.outputs_to_keep = outputs_to_keep
        # check that outputs_to_keep are Paths else error
        for output in self.outputs_to_keep:
            if not isinstance(output, Path):
                raise ValueError(f"All outputs to keep must be of type Path. Invalid output: {output} ({type(output)})")


    def run(self, context: Context):
        self.welcome()

        if not context.skip_check:
            try:
                self.validate_dependencies()
            except DependencyError as e:
                raise DependencyError(f"Invalid dependencies! Please install '{e}' or use --skip-check to ignore.") from e
        
        if context.check:
            return None

        # Set working directory to output folder
        current_dir = Path.cwd()
        self.outdir = context.outdir
        try:
            self.set_working_directory(context.outdir)
            self._execute_pipeline_stages_in_order(
                context,
            )
            self.cleanup(context, outputs_to_keep=self.outputs_to_keep)
        finally:
            # Ensure we always return to the original working directory
            os.chdir(current_dir)
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
    
    def is_debug(self):
        return logger.is_debug()

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
                if dependency.citation_only:
                    self.debug(f"Skipping {dependency.name}, citation only.")
                    continue
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
        try:
            Path(directory).mkdir(parents=True, exist_ok=True)
            os.chdir(directory)
        except PermissionError as e:
            raise PipelineExecutionError(
                f"Permission denied while creating or entering output directory '{directory}'"
            ) from e

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
            # Clean up stage tmp outputs
            if not run_ctx.no_cleanup:
                stage.cleanup_tmp_outputs()
        self.end_time = time.perf_counter()

    def cleanup(self, context: Context = None, outputs_to_keep: List[Path] | None = None):
        """Delete all declared stage outputs except those explicitly kept."""
        if context and context.no_cleanup:
            self.debug("Skipping output cleanup (--no-cleanup).")
            return

        keep_keys = {os.path.abspath(str(p)) for p in (outputs_to_keep or [])}
        if not keep_keys:
            self.debug("No outputs specified to keep, skipping cleanup.")
            return
        self.debug(f"Keeping {len(keep_keys)} output file(s) from cleanup.")
        self.debug(f"Keep list: {outputs_to_keep}")
        removed = 0
        seen = set()
        dirs_to_remove = set()

        for stage in self.stages:
            if stage.output._immutable:
                self.debug(f"Skipping cleanup for stage '{stage.name}' with immutable output.")
                continue
            # Assume temporary outputs are already cleaned up by the stage
            for _, output in stage.output.non_temporary_outputs():
                key = os.path.abspath(str(output))
                if key in seen:
                    continue
                seen.add(key)

                if key in keep_keys:
                    continue

                if output.exists():
                    if output.is_dir():
                        dirs_to_remove.add(output)
                    else:
                        self.debug(f"Removing output file: {output}")
                        output.unlink()
                        removed += 1

        # Remove directories only if they are empty, after files are deleted.
        for output_dir in sorted(dirs_to_remove, reverse=True):
            key = os.path.abspath(str(output_dir))
            if key in keep_keys:
                continue
            if output_dir.exists() and output_dir.is_dir():
                try:
                    output_dir.rmdir()
                    self.debug(f"Removed empty output directory: {output_dir}")
                    removed += 1
                except OSError:
                    self.warning(f"Could not remove output directory {output_dir} (not empty).")
                    continue

        if removed:
            self.debug(f"Removed {removed} output file(s) not in keep list.")
    
    @property
    def citations(self):
        citations = []
        for stage in self.stages:
            for dependency in stage._dependencies:
                if dependency.citation:
                    citations.append(dependency.citation)
        return sorted(set(citations))
    
    def output_descriptions(self) -> OrderedDict[str, Dict[str, str]]:
        """Return kept outputs as a dict keyed by path with output metadata."""
        keep_keys = {os.path.abspath(str(p)) for p in self.outputs_to_keep}
        descriptions: OrderedDict[str, Dict[str, str]] = OrderedDict()
        seen = set()

        for stage in self.stages:
            output_model = stage.output
            if output_model._immutable:
                continue
            for field_name in output_model.__class__.model_fields:
                output_value = getattr(output_model, field_name, None)
                if not isinstance(output_value, Path):
                    continue

                key = os.path.abspath(str(output_value))
                if key not in keep_keys or key in seen:
                    continue
                seen.add(key)

                description = stage.get_output_description(field_name)
                description_formatted = f"{description}" if description else ""
                size_formatted = human_readable_size(Path(self.outdir) / output_value if self.outdir else output_value)
                descriptions[str(output_value)] = {
                    "output": field_name,
                    "path": str(output_value),
                    "size": size_formatted,
                    "description": description_formatted,
                }

        return descriptions

    def write_output_descriptions(self, output_file: Path):
        descriptions = self.output_descriptions()
        if descriptions:
            header = "output\tpath\tsize\tdescription"
            rows = [
                f"{d['output']}\t{d['path']}\t{d['size']}\t{d['description']}"
                for d in descriptions.values()
            ]
            output_file.write_text("\n".join([header, *rows]))
            self.debug(f"Wrote output descriptions to {output_file}")
        else:
            self.debug("No output descriptions to write.")
     
    def write_output_citations(self, output_file: Path):
        citations = self.citations
        if citations:
            output_file.write_text("\n".join(f"{i}. {cite}" for i, cite in enumerate(self.citations, 1)))
            self.debug(f"Wrote {len(citations)} citations to {output_file}")
        else:
            self.debug("No citations to write.")

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
        outdir = Path.cwd() if self.outdir is None else Path(self.outdir)
        self.write_output_descriptions(outdir / "FILES.tsv")
        self.write_output_citations(outdir / "CITATIONS.txt")
        self.hr()
        # Print a random goodbye message
        self.hr(f"{random.choice(messages)}", style=" ", color='green')
