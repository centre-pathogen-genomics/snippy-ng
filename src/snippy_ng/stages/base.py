from contextlib import contextmanager
from abc import abstractmethod
from typing import List, Callable, Optional
import subprocess
import sys
from io import StringIO
from pathlib import Path

from snippy_ng.logging import logger
from snippy_ng.dependencies import Dependency
from snippy_ng.exceptions import SkipStageError

from pydantic import BaseModel, Field

class BaseOutput(BaseModel):
    pass

class Command(BaseModel):
    func: Callable
    args: List = []
    description: str
    
    def __str__(self):
        return self.description

class BaseStage(BaseModel):
    cpus: int = Field(1, description="Number of CPU cores to use")
    ram: Optional[int] = Field(4, description="RAM in GB to use")
    tmpdir: Path = Field(description="Temporary directory")

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
    def commands(self) -> List[str | Command]:
        """Constructs the commands."""
        pass

    def run(self, quiet=False):
        """Runs the commands in the shell or calls the function."""
        for cmd in self.commands:
            logger.info(f"Running: {cmd}") 
            stdout = sys.stderr
            stderr = sys.stderr
            if quiet:
                stdout = subprocess.DEVNULL
                stderr = subprocess.DEVNULL
            try:
                if isinstance(cmd, Command):
                    if quiet:
                        # Capture output if quiet mode is enabled
                        with StringIO() as out, StringIO() as err, \
                                self.redirect_output(out), self.redirect_error(err):
                            cmd.func(*cmd.args)
                    else:
                        with self.redirect_stdout_to_err():
                            cmd.func(*cmd.args)
                elif isinstance(cmd, str):
                    subprocess.run(cmd, shell=True, check=True, stdout=stdout, stderr=stderr, text=True)
                else:
                    raise ValueError(f"Invalid command type: {type(cmd)}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Command failed with exit code {e.returncode}")
                raise RuntimeError(f"Failed to run command: {cmd}")
            except SkipStageError as e:
                logger.warning(f"Skipping stage: {e}")
                raise e
            except Exception as e:
                logger.error(f"Function call failed: {e}")
                raise RuntimeError(f"Failed to run function: {cmd}")
    
    @contextmanager
    def redirect_stdout_to_err(self):
        original_stdout = sys.stdout
        sys.stdout = sys.stderr
        try:
            yield
        finally:
            sys.stdout = original_stdout

    @contextmanager
    def redirect_output(self, output_stream):
        original_stdout = sys.stdout
        sys.stdout = output_stream
        try:
            yield
        finally:
            sys.stdout = original_stdout

    @contextmanager
    def redirect_error(self, error_stream):
        original_stderr = sys.stderr
        sys.stderr = error_stream
        try:
            yield
        finally:
            sys.stderr = original_stderr
        