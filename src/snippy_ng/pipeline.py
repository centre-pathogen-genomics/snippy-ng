import random
from typing import List
import os
from pathlib import Path

from snippy_ng.exceptions import DependencyError, SkipStageError, MissingOutputError
from snippy_ng.logging import logger
from snippy_ng.__about__ import __version__, URL
from snippy_ng.stages.base import BaseStage

class Pipeline:
    """
    Main class for creating Snippy-NG Pipelines.
    """
    stages: List[BaseStage]
    line: str = "-" * 39

    def __init__(self, stages: List[BaseStage] = None):
        if stages is None:
            stages = []
        self.stages = stages
    
    def add_stage(self, stage):
        self.stages.append(stage)

    @property
    def dependencies(self):
        return [dep for stage in self.stages for dep in stage._dependencies]

    def log(self, msg):
        logger.info(msg)
    
    def warning(self, msg):
        logger.warning(msg)
    
    def debug(self, msg):
        logger.debug(msg)

    def error(self, msg):
        logger.error(msg)

    def welcome(self):
        self.log(self.line)
        self.log("  🦘 ⚡✂️   Running Snippy-NG  ✂️  ⚡🦘")
        self.log(self.line)
        self.log(f"Version: {__version__}")
        self.log(self.line)
        self.log(f"Stages: {' -> '.join([stage.name for stage in self.stages])}")
        self.log(self.line)

    def validate_dependencies(self):
        invalid = []
        self.log("CHECKING DEPENDENCIES...")
        for stage in self.stages:
            self.log(f"Checking dependencies for {stage.name}...")
            for dependency in stage._dependencies:
                # TODO: skip if already checked
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
        self.log(f"Setting working directory to '{directory}'")
        os.chdir(directory)

    def run(self, quiet=False):
        # Run pipeline sequentially
        for stage in self.stages:
            self.log(f"RUNNING {stage.name} STAGE...")
            self.log(stage)
            try:
                stage.run(quiet)
                for name, output in stage.output:
                    if not Path(output).exists():
                        self.error(f"Output file {output} not found!")
                        raise MissingOutputError("Output file not found!")
                self.log(f"STAGE {stage.name} COMPLETE!")
            except SkipStageError:
                self.stages.remove(stage)
                self.warning(f"STAGE {stage.name} SKIPPED!")

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

    def goodbye(self):
        messages = [
            "May the SNPs be with you.",
            "Wishing you a life free of homopolymer errors.",
            f"Found a bug? Post it at {URL}/issues",
            f"Have a suggestion? Tell me at {URL}/issues",
            f"The Snippy manual is at {URL}/blob/master/README.md",
            "Questionable SNP? Try the --report option to see the alignments.",
            "Did you know? Snippy is a combination of SNP, Skippy, and snappy.",
            "Set phasers to align… your SNPs.",
            "To boldly SNP where no one has SNPped before.",
            "Resistance to accurate SNP calling is futile.",
            "Wishing you a genome that's logically consistent…",
            "The final frontier of variant calling.",
            f"Make it so: Report your issues at {URL}/issues.",
            "Highly logical and warp-speed fast.",
            "Live long and SNP accurately.",
            f"Looking for guidance? The Snippy manual is at {URL}.",
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
        # Print a random goodbye message
        self.log(self.line)
        self.log("  🦘 ⚡ ✂️ Snippy-NG complete! ✂️ ⚡ 🦘")
        self.log(self.line)
        self.log(f"Please cite the following:\n{'- ' + '\n- '.join(self.citations)}")
        self.log(self.line)
        self.log(f"{random.choice(messages)}")
