from abc import abstractmethod
from typing import List
import subprocess

from snippy_ng.logging import logger
from snippy_ng.dependencies import Dependency

from pydantic import BaseModel

class BaseOutput(BaseModel):
    pass

class BaseStage(BaseModel):
    _dependencies: List[Dependency] = []
    
    @property
    def name(self) -> str:
        """Returns the name of the stage."""
        return self.__class__.__name__

    @property
    @abstractmethod
    def output(self) -> BaseOutput:
        """Defines the output of the stage."""
        pass

    @property
    @abstractmethod
    def commands(self) -> List[str]:
        """Constructs the shell commands."""
        pass

    def run(self, quiet=False):
        """Runs the commands in the shell."""
        for cmd in self.commands:
            logger.info(f"Running: {cmd}")
            stdout = None
            stderr = subprocess.STDOUT
            if quiet:
                stdout = subprocess.DEVNULL
                stderr = subprocess.DEVNULL
            try:
                subprocess.run(cmd, shell=True, check=True, stdout=stdout, stderr=stderr, text=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"Command failed with exit code {e.returncode}")
                raise RuntimeError(f"Failed to run command: {cmd}")