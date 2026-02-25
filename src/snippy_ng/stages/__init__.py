from __future__ import annotations

from contextlib import contextmanager, nullcontext
from abc import abstractmethod
from typing import ClassVar, Iterable, List, Callable, Optional
import subprocess
import sys
from io import StringIO
from pathlib import Path
import signal

from snippy_ng.logging import logger
from snippy_ng.dependencies import Dependency
from snippy_ng.exceptions import InvalidCommandTypeError, MissingOutputError, StageExecutionError, StageTestFailure
from snippy_ng.context import Context

from pydantic import BaseModel, ConfigDict, Field
from shlex import quote


class BaseOutput(BaseModel):
    model_config = ConfigDict(extra='forbid')
    _immutable: bool = False

class PythonCommand(BaseModel):
    func: Callable
    args: List = []
    description: str

    def __str__(self):
        return f"{self.func.__name__}({', '.join(map(str, self.args))})"

class ShellCommand(BaseModel):
    command: List[str]
    description: str
    output_file: Optional[Path] = None

    def __str__(self):
        output_file = ''
        if self.output_file:
            output_file = f" > {self.output_file}"
        return f"{' '.join(quote(str(arg)) for arg in self.command)}{output_file}"

class ShellProcessPipe(BaseModel):
    processes: List[ShellCommand] # TODO rename to parts
    description: str
    output_file: Optional[Path] = None

    def __str__(self):
        output_file = ''
        if self.output_file:
            output_file = f" > {self.output_file}"
        return f"{' | '.join(str(cmd) for cmd in self.processes)}{output_file}"
    
    def __iter__(self):
        return iter(self.processes)

TestFn = Callable[["BaseStage"], None]
class BaseStage(BaseModel):
    model_config = ConfigDict(extra='forbid', arbitrary_types_allowed=True)
    prefix: str = Field("snps", description="Prefix for output files") # TODO maybe this should be in the context instead?

    _dependencies: List[Dependency] = []
    _tests: ClassVar[List[TestFn]] = []

    @property
    def name(self) -> str:
        """Returns the name of the stage."""
        return self.__class__.__name__

    @property
    @abstractmethod
    def output(self) -> BaseOutput:
        """Defines the output of the stage."""
        pass

    @abstractmethod
    def create_commands(self, ctx: Context) -> List[ShellProcessPipe | ShellCommand | PythonCommand]:
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
    
    def shell_cmd(self, command: List[str], description: str, output_file: Optional[Path] = None) -> ShellCommand:
        """Creates a shell command."""
        assert isinstance(command, list), f"Command must be a list of strings, got {command}"
        assert all(isinstance(arg, str) for arg in command), f"All command arguments must be strings, got {command}"
        assert isinstance(description, str), f"Description must be a string, got {description}"
        return ShellCommand(command=command, description=description, output_file=output_file)
    
    def shell_pipe(self, commands: List[ShellCommand], description: str, output_file: Optional[Path] = None) -> ShellProcessPipe:
        """Creates a shell pipeline."""
        # Validate that all commands are ShellCommand objects
        for i, cmd in enumerate(commands):
            if not isinstance(cmd, ShellCommand):
                raise InvalidCommandTypeError(
                    f"Pipeline command at index {i} must be a ShellCommand, got {type(cmd).__name__}. "
                    f"Use self.shell_cmd() to create ShellCommand objects."
                )
        if not commands:
            raise InvalidCommandTypeError("Pipeline must contain at least one ShellCommand.")
        for i, cmd in enumerate(commands[:-1]):
            if cmd.output_file:
                raise InvalidCommandTypeError(
                    f"Pipeline command at index {i} cannot set output_file. "
                    "Only the final command may set output_file."
                )
        if output_file and commands[-1].output_file:
            raise InvalidCommandTypeError(
                "Pipeline output_file conflicts with final command output_file. "
                "Set output_file on the pipeline or the final command, not both."
            )
        return ShellProcessPipe(processes=commands, description=description, output_file=output_file)


    def run(self, ctx: Context):
        """Runs the commands in the shell or calls the function."""
        if ctx.create_missing:
            try:
                self.error_if_outputs_missing()
                logger.info(f"{self.name} already completed, skipping...")
                return
            except MissingOutputError:
                pass
        try:
            for cmd in self.create_commands(ctx):
                assert isinstance(cmd, (ShellCommand, PythonCommand, ShellProcessPipe)), f"Invalid command type: {type(cmd)} in stage {self.name}"
                logger.info(cmd.description)
                logger.info(f" ❯ {cmd}") 
                stdout = sys.stderr
                stderr = sys.stderr
                if ctx.quiet:
                    stdout = subprocess.DEVNULL
                    stderr = subprocess.DEVNULL
                try:
                    if isinstance(cmd, PythonCommand):
                        if ctx.quiet:
                            # Capture output if quiet mode is enabled
                            with StringIO() as out, StringIO() as err, \
                                    self.redirect_output(out), self.redirect_error(err):
                                cmd.func(*cmd.args)
                        else:
                            with self.redirect_stdout_to_err():
                                cmd.func(*cmd.args)
                    elif isinstance(cmd, ShellCommand):
                        if cmd.output_file:
                            with open(cmd.output_file, 'w') as out:
                                subprocess.run(cmd.command, check=True, stdout=out, stderr=stderr, text=True)
                        else:
                            subprocess.run(cmd.command, check=True, stdout=stdout, stderr=stderr, text=True)
                    elif isinstance(cmd, ShellProcessPipe):
                        if not cmd.processes:
                            raise ValueError("No commands to run in the pipeline")

                        processes: List[subprocess.Popen] = []
                        last_output_file = cmd.output_file or (cmd.processes[-1].output_file if cmd.processes else None)

                        try:
                            with (open(last_output_file, 'wb') if last_output_file else nullcontext(stdout)) as final_stdout:
                                prev_proc: Optional[subprocess.Popen] = None
                                for i, pipeline_part in enumerate(cmd.processes):
                                    is_last = i == len(cmd.processes) - 1
                                    process_stdout = final_stdout if is_last else subprocess.PIPE
                                    process_stdin = None if prev_proc is None else prev_proc.stdout

                                    p = subprocess.Popen(
                                        pipeline_part.command,
                                        stdin=process_stdin,
                                        stdout=process_stdout,
                                        stderr=stderr,
                                        text=False, # pipes should be in binary mode
                                    )
                                    processes.append(p)

                                    # Allow upstream processes to receive SIGPIPE if downstream exits.
                                    if prev_proc is not None and prev_proc.stdout is not None:
                                        prev_proc.stdout.close()

                                    prev_proc = p

                                # Wait for all processes to complete
                                for p in processes:
                                    p.wait()

                                failures = [p for p in processes if p.returncode not in (0, None)]
                                if failures:
                                    sigpipe_rc = -int(getattr(signal, "SIGPIPE", 13))
                                    non_sigpipe = [p for p in failures if p.returncode != sigpipe_rc]
                                    primary = (non_sigpipe[-1] if non_sigpipe else failures[0])
                                    raise subprocess.CalledProcessError(
                                        returncode=primary.returncode,
                                        cmd=primary.args,
                                    )
                        except Exception:
                            # Ensure we don't leave a half-running pipeline behind.
                            for p in processes:
                                try:
                                    if p.poll() is None:
                                        p.terminate()
                                except Exception:
                                    # Ignore errors during best-effort process termination in cleanup.
                                    pass
                            for p in processes:
                                try:
                                    p.wait(timeout=2)
                                except Exception:
                                    # Ignore errors while waiting during best-effort cleanup.
                                    # Ignore errors while waiting for processes to exit during cleanup.
                                    pass
                            raise
                    else:
                        raise InvalidCommandTypeError(f"Command must be of type List or PythonCommand, got {type(cmd)}")
                except subprocess.CalledProcessError as e:
                    logger.error(f"Command failed with exit code {e.returncode}")
                    cmd = " ".join(quote(arg) for arg in e.cmd) 
                    raise StageExecutionError(f"Failed to run command: {cmd}")
                except InvalidCommandTypeError as e:
                    raise e
        except (Exception, KeyboardInterrupt) as e:
                # remove outputs if stage fails
                if ctx.keep_incomplete:
                    raise e 
                if self.output._immutable:
                    raise e
                output_removed = False
                for name, path in self.output:
                    if path and Path(path).exists():
                        output_removed = True
                        logger.warning(f"Removing incomplete output '{name}' ({path}).")
                        Path(path).unlink()
                if output_removed:
                    logger.warning("Set `keep_incomplete=True` to retain incomplete outputs on error.")
                raise e

    def _discover_test_methods(self) -> Iterable[Callable[[], None]]:
        """
        Finds bound instance methods named test_*.
        """
        discovered = []
        for name in sorted(dir(self)):
            if not name.startswith("test_"):
                continue
            attr = getattr(self, name, None)
            if not callable(attr):
                continue

            discovered.append(attr) 
        return discovered 

    def run_tests(self) -> None:
        """
        Runs all test methods defined on the stage. Test methods should be instance methods that take no arguments and are named with a "test_" prefix.
        """
        discovered_tests = self._discover_test_methods()
        if not discovered_tests:
            logger.debug(f"No tests found for {self.name}")
            return
        logger.debug("Running tests...")
        for t in discovered_tests:
            logger.debug(f"{t.__name__}...")
            try:
                t()
            except Exception as e:
                name = getattr(t, "__name__", repr(t))
                raise StageTestFailure(f"{name.upper()}: {e}") from e

    def error_if_outputs_missing(self):
        """Raises an error if any expected output files are missing."""
        missing = []
        for name, path in self.output:
            if not path:
                continue
            if not Path(path).exists():
                missing.append((name, path))
        if missing:
            missing_str = ", ".join(f"{name} ({path})" for name, path in missing)
            raise MissingOutputError(f"Expected output files are missing: {missing_str}")

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
