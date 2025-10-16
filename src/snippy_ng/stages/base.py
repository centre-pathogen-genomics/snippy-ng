from contextlib import contextmanager, nullcontext
from abc import abstractmethod
from typing import List, Callable, Optional
import subprocess
import sys
import tempfile
from io import StringIO
from pathlib import Path

from snippy_ng.logging import logger
from snippy_ng.dependencies import Dependency
from snippy_ng.exceptions import InvalidCommandTypeError, SkipStageError

from pydantic import BaseModel, Field
from shlex import quote


class BaseOutput(BaseModel):
    pass

class PythonCommand(BaseModel):
    func: Callable
    args: List = []
    description: str

    def __str__(self):
        return f"{self.func.__name__}({', '.join(map(str, self.args))})"

class ShellCommand(BaseModel):
    command: List[str]
    description: str

    def __str__(self):
        return " ".join(quote(arg) for arg in self.command)

class ShellCommandPipeline(BaseModel):
    commands: List[ShellCommand]
    description: str
    output_file: Optional[Path] = None

    def __str__(self):
        output_file = ''
        if self.output_file:
            output_file = f" > {self.output_file}"
        return f"{' | '.join(str(cmd) for cmd in self.commands)}{output_file}"
    
    def __iter__(self):
        return iter(self.commands)

class BaseStage(BaseModel):
    cpus: int = Field(1, description="Number of CPU cores to use")
    ram: Optional[int] = Field(4, description="RAM in GB to use")
    tmpdir: Optional[Path] = Field(default_factory=lambda: Path(tempfile.gettempdir()), description="Temporary directory")

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
    def commands(self) -> List[ShellCommandPipeline | ShellCommand | PythonCommand]:
        """Constructs the commands."""
        pass

    def escape(self, user_value) -> str:
        """Returns an escaped string for shell commands."""
        return quote(str(user_value))
    
    def python_cmd(self, func: Callable, args: List = [], description: Optional[str] = None) -> PythonCommand:
        """Creates a Python command."""
        if description is None:
            description = f"{func.__name__} with arguments {', '.join(args)}"
        return PythonCommand(func=func, args=args, description=description)
    
    def shell_cmd(self, command: List[str], description: str) -> ShellCommand:
        """Creates a shell command."""
        return ShellCommand(command=command, description=description)
    
    def shell_pipeline(self, commands: List[ShellCommand], description: str, output_file: Optional[Path] = None) -> ShellCommandPipeline:
        """Creates a shell pipeline."""
        # Validate that all commands are ShellCommand objects
        for i, cmd in enumerate(commands):
            if not isinstance(cmd, ShellCommand):
                raise InvalidCommandTypeError(
                    f"Pipeline command at index {i} must be a ShellCommand, got {type(cmd).__name__}. "
                    f"Use self.shell_cmd() to create ShellCommand objects."
                )
        return ShellCommandPipeline(commands=commands, description=description, output_file=output_file)


    def run(self, quiet=False):
        """Runs the commands in the shell or calls the function."""
        for cmd in self.commands:
            logger.info(cmd.description)
            logger.info(f"❯ {cmd}") 
            stdout = sys.stderr
            stderr = sys.stderr
            if quiet:
                stdout = subprocess.DEVNULL
                stderr = subprocess.DEVNULL
            try:
                if isinstance(cmd, PythonCommand):
                    if quiet:
                        # Capture output if quiet mode is enabled
                        with StringIO() as out, StringIO() as err, \
                                self.redirect_output(out), self.redirect_error(err):
                            cmd.func(*cmd.args)
                    else:
                        with self.redirect_stdout_to_err():
                            cmd.func(*cmd.args)
                elif isinstance(cmd, ShellCommand):
                    subprocess.run(cmd.command, check=True, stdout=stdout, stderr=stderr, text=True)
                elif isinstance(cmd, ShellCommandPipeline):
                    processes = []
                    prev_stdout = None
                    output_file = cmd.output_file
                    
                    if not cmd.commands:
                        raise ValueError("No commands to run in the pipeline")
                    
                    with (open(output_file, 'w') if output_file else nullcontext(stdout)) as final_stdout:
                        for i, pipeline_part in enumerate(cmd.commands):
                            # Determine stdout for this process
                            if i == len(cmd.commands) - 1:
                                # Last process - output to final destination
                                process_stdout = final_stdout
                            else:
                                # Intermediate process - output to pipe
                                process_stdout = subprocess.PIPE
                            
                            p = subprocess.Popen(
                                pipeline_part.command, 
                                stdin=prev_stdout, 
                                stdout=process_stdout, 
                                stderr=stderr,
                                text=True
                            )
                            
                            # Close the previous stdout pipe (we're done with it)
                            if prev_stdout is not None:
                                prev_stdout.close()
                            
                            # Set up for next iteration
                            if i < len(cmd.commands) - 1:
                                prev_stdout = p.stdout
                            
                            processes.append(p)
                        
                        # Wait for all processes to complete
                        for p in processes:
                            p.wait()
                        
                        # Check for failures
                        failed_processes = [p for p in processes if p.returncode != 0]
                        if failed_processes:
                            failed_process = failed_processes[0]  # Report first failure
                            raise subprocess.CalledProcessError(
                                returncode=failed_process.returncode, 
                                cmd=failed_process.args
                            )
                else:
                    raise InvalidCommandTypeError(f"Command must be of type List or PythonCommand, got {type(cmd)}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Command failed with exit code {e.returncode}")
                cmd = " ".join(quote(arg) for arg in e.cmd) 
                raise RuntimeError(f"Failed to run command: {cmd}")
            except SkipStageError as e:
                logger.warning(f"Skipping stage: {e}")
                raise e
            except InvalidCommandTypeError as e:
                raise e
            except Exception as e:
                logger.error(f"Function call failed: {e}")
                raise RuntimeError(f"Failed to run function: {cmd}")

    def check_outputs(self) -> bool:
        """Check if all expected output files exist."""
        for _, path in self.output:
            if not path:
                continue
            if not Path(path).exists():
                return False
        return True

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
