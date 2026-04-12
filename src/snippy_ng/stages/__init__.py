from __future__ import annotations

from contextlib import ExitStack, contextmanager, nullcontext
from abc import abstractmethod
from types import UnionType
from typing import Annotated, Any, ClassVar, Iterable, List, Callable, Optional, Union, get_args, get_origin
import subprocess
import sys
from io import StringIO
from pathlib import Path
import signal
import shutil
import threading

from snippy_ng.logging import logger
from snippy_ng.dependencies import Dependency
from snippy_ng.exceptions import InvalidCommandTypeError, MissingOutputError, StageExecutionError, StageTestFailure
from snippy_ng.context import Context

from pydantic import BaseModel, ConfigDict, Field, PrivateAttr
from shlex import quote


class TemporaryOutput:
    """Marker metadata for temporary output fields."""


TempPath = Annotated[Path, TemporaryOutput()]


class BaseOutput(BaseModel):
    model_config = ConfigDict(extra='forbid')
    _immutable: bool = False
    _description_overrides: dict[str, str] = PrivateAttr(default_factory=dict)

    def get_description(self, field_name: str) -> str:
        """Get an output field description, honoring instance overrides."""
        if field_name in self._description_overrides:
            return self._description_overrides[field_name]
        field_info = self.__class__.model_fields[field_name]
        return field_info.description or ""

    @classmethod
    def _annotation_has_temporary_marker(cls, annotation: Any) -> bool:
        """Recursively detect TemporaryOutput marker in nested type annotations."""
        if annotation is None:
            return False

        origin = get_origin(annotation)

        if origin is Annotated:
            args = get_args(annotation)
            if not args:
                return False
            base_type, *metadata = args
            if any(isinstance(meta, TemporaryOutput) for meta in metadata):
                return True
            return cls._annotation_has_temporary_marker(base_type)

        if origin in (Union, UnionType):
            return any(cls._annotation_has_temporary_marker(arg) for arg in get_args(annotation))

        return False

    @classmethod
    def _is_temporary_field(cls, field_name: str) -> bool:
        field_info = cls.model_fields[field_name]
        if any(isinstance(meta, TemporaryOutput) for meta in field_info.metadata):
            return True
        return cls._annotation_has_temporary_marker(field_info.annotation)

    def temporary_outputs(self) -> List[tuple[str, Path]]:
        """Return temporary output fields and their resolved path values."""
        tmp_outputs = []
        for name in self.__class__.model_fields:
            if not self._is_temporary_field(name):
                continue
            value = getattr(self, name, None)
            if isinstance(value, Path):
                tmp_outputs.append((name, value))
        return tmp_outputs

    def non_temporary_outputs(self) -> List[tuple[str, Path]]:
        """Return non-temporary output fields and their resolved path values."""
        outputs = []
        for name in self.__class__.model_fields:
            if self._is_temporary_field(name):
                continue
            value = getattr(self, name, None)
            if isinstance(value, Path):
                outputs.append((name, value))
        return outputs

    def all_outputs(self) -> List[tuple[str, Path]]:
        """Return all output fields and their resolved path values."""
        outputs = []
        for name in self.__class__.model_fields:
            value = getattr(self, name, None)
            if isinstance(value, Path):
                outputs.append((name, value))
        return outputs

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
    prefix: str = Field("snippy", description="Prefix for output files") # TODO maybe this should be in the context instead?

    _dependencies: List[Dependency] = []
    _tests: ClassVar[List[TestFn]] = []
    _output_description_overrides: dict[str, str] = PrivateAttr(default_factory=dict)

    @property
    def name(self) -> str:
        """Returns the name of the stage."""
        return self.__class__.__name__
    
    @property
    @abstractmethod
    def output(self) -> BaseOutput:
        """Defines the output of the stage."""
        pass

    def update_output_description(self, field_name: str, description: str) -> None:
        """Override an output field description for this stage instance."""
        if field_name not in self.output.__class__.model_fields:
            raise ValueError(f"Unknown output field '{field_name}' for {self.name}")
        self._output_description_overrides[field_name] = description

    def get_output_description(self, field_name: str) -> str:
        """Get the effective output description for this stage instance."""
        if field_name in self._output_description_overrides:
            return self._output_description_overrides[field_name]
        return self.output.get_description(field_name)

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


    def _log_command_output(self, output: str | bytes | None, *, quiet: bool) -> None:
        if output in (None, b"", ""):
            return
        if isinstance(output, bytes):
            output = output.decode("utf-8", errors="replace")
        logger.echo(output, err=True, nl=False, console=not quiet)

    def _stream_stderr(self, stream, *, quiet: bool, stderr_chunks: list[str]) -> None:
        try:
            for chunk in iter(stream.readline, b""):
                if not chunk:
                    break
                decoded = chunk.decode("utf-8", errors="replace")
                stderr_chunks.append(decoded)
                logger.echo(decoded, err=True, nl=False, console=not quiet)
        finally:
            stream.close()

    def _run_streaming_shell_command(self, command: list[str], *, quiet: bool, output_file: Optional[Path] = None) -> None:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=False,
            bufsize=-1,
        )
        stdout_chunks: list[bytes] = []
        stderr_chunks: list[str] = []

        def consume_stdout() -> None:
            if process.stdout is None:
                return
            try:
                with (open(output_file, "wb") if output_file is not None else nullcontext()) as output_handle:
                    for chunk in iter(process.stdout.readline, b""):
                        if not chunk:
                            break
                        stdout_chunks.append(chunk)
                        if output_handle is not None:
                            output_handle.write(chunk)
                            output_handle.flush()
                        else:
                            self._log_command_output(chunk, quiet=quiet)
            finally:
                process.stdout.close()

        stdout_thread = threading.Thread(target=consume_stdout, daemon=True)
        stderr_thread = threading.Thread(
            target=self._stream_stderr,
            args=(process.stderr,),
            kwargs={"quiet": quiet, "stderr_chunks": stderr_chunks},
            daemon=True,
        )
        stdout_thread.start()
        stderr_thread.start()
        returncode = process.wait()
        stdout_thread.join()
        stderr_thread.join()
        if returncode != 0:
            raise subprocess.CalledProcessError(
                returncode=returncode,
                cmd=command,
                output=b"".join(stdout_chunks),
                stderr="".join(stderr_chunks),
            )

    def _run_python_command(self, cmd: PythonCommand, ctx: Context) -> None:
        with StringIO() as out, StringIO() as err, self.redirect_output(out), self.redirect_error(err):
            cmd.func(*cmd.args)
            stdout_text = out.getvalue()
            stderr_text = err.getvalue()
        self._log_command_output(stdout_text, quiet=ctx.quiet)
        self._log_command_output(stderr_text, quiet=ctx.quiet)

    def _run_shell_command(self, cmd: ShellCommand, ctx: Context) -> None:
        self._run_streaming_shell_command(cmd.command, quiet=ctx.quiet, output_file=cmd.output_file)

    def _run_shell_pipeline(self, cmd: ShellProcessPipe, ctx: Context) -> None:
        if not cmd.processes:
            raise ValueError("No commands to run in the pipeline")
        processes: List[subprocess.Popen] = []
        last_output_file = cmd.output_file or (cmd.processes[-1].output_file if cmd.processes else None)
        stdout_chunks: list[bytes] = []
        stderr_chunks: dict[int, list[str]] = {}
        stderr_threads: list[threading.Thread] = []
        try:
            with ExitStack() as stack:
                final_stdout_handle = stack.enter_context(open(last_output_file, "wb")) if last_output_file else None
                prev_proc: Optional[subprocess.Popen] = None
                for pipeline_part in cmd.processes:
                    process_stdin = None if prev_proc is None else prev_proc.stdout
                    p = subprocess.Popen(
                        pipeline_part.command,
                        stdin=process_stdin,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=False,
                        bufsize=-1,
                    )
                    processes.append(p)
                    stderr_chunks[id(p)] = []
                    if p.stderr is not None:
                        thread = threading.Thread(
                            target=self._stream_stderr,
                            args=(p.stderr,),
                            kwargs={"quiet": ctx.quiet, "stderr_chunks": stderr_chunks[id(p)]},
                            daemon=True,
                        )
                        stderr_threads.append(thread)
                        thread.start()
                    if prev_proc is not None and prev_proc.stdout is not None:
                        prev_proc.stdout.close()
                    prev_proc = p

                if processes and processes[-1].stdout is not None:
                    for chunk in iter(processes[-1].stdout.readline, b""):
                        if not chunk:
                            break
                        stdout_chunks.append(chunk)
                        if final_stdout_handle is not None:
                            final_stdout_handle.write(chunk)
                            final_stdout_handle.flush()
                        else:
                            self._log_command_output(chunk, quiet=ctx.quiet)
                    processes[-1].stdout.close()

                for p in processes:
                    p.wait()
                for thread in stderr_threads:
                    thread.join()

                failures = [p for p in processes if p.returncode not in (0, None)]
                if failures:
                    sigpipe_rc = -int(getattr(signal, "SIGPIPE", 13))
                    non_sigpipe = [p for p in failures if p.returncode != sigpipe_rc]
                    primary = (non_sigpipe[-1] if non_sigpipe else failures[0])
                    raise subprocess.CalledProcessError(
                        returncode=primary.returncode,
                        cmd=primary.args,
                        output=b"".join(stdout_chunks),
                        stderr="".join(
                            chunk
                            for proc in processes
                            for chunk in stderr_chunks.get(id(proc), [])
                        ),
                    )
        except Exception:
            for p in processes:
                try:
                    if p.poll() is None:
                        p.terminate()
                except Exception:
                    pass
            for p in processes:
                try:
                    p.wait(timeout=2)
                except Exception:
                    pass
            raise

    def _run_command(self, cmd: ShellProcessPipe | ShellCommand | PythonCommand, ctx: Context) -> None:
        if isinstance(cmd, PythonCommand):
            self._run_python_command(cmd, ctx)
        elif isinstance(cmd, ShellCommand):
            self._run_shell_command(cmd, ctx)
        elif isinstance(cmd, ShellProcessPipe):
            self._run_shell_pipeline(cmd, ctx)
        else:
            raise InvalidCommandTypeError(f"Command must be of type List or PythonCommand, got {type(cmd)}")

    def _should_skip(self, ctx: Context) -> bool:
        if not ctx.create_missing:
            return False
        persistent_outputs = self.output.non_temporary_outputs()
        if persistent_outputs:
            try:
                self.error_if_outputs_missing()
                logger.info(f"{self.name} already completed, skipping...")
                return True
            except MissingOutputError:
                return False
        logger.debug(f"{self.name} has no persistent outputs (all temporary or none), rerunning with --create-missing.")
        return False

    def _cleanup_incomplete_outputs(self, ctx: Context) -> None:
        if ctx.keep_incomplete or self.output._immutable:
            return
        output_removed = False
        for name, path in self.output.all_outputs():
            output_path = Path(path)
            if path and output_path.exists():
                resolved_output = output_path.resolve()
                if resolved_output == Path.cwd():
                    logger.warning(f"Skipping cleanup of current working directory output '{name}' ({path}).")
                    continue
                output_removed = True
                logger.warning(f"Removing incomplete output '{name}' ({path}).")
                if output_path.is_dir() and not output_path.is_symlink():
                    shutil.rmtree(output_path)
                else:
                    output_path.unlink()
        if output_removed:
            logger.warning("Set `keep_incomplete=True` to retain incomplete outputs on error.")

    def run(self, ctx: Context):
        """Runs the commands in the shell or calls the function."""
        if self._should_skip(ctx):
            return
        if ctx.break_points:
            self.breakpoint()
        try:
            for cmd in self.create_commands(ctx):
                assert isinstance(cmd, (ShellCommand, PythonCommand, ShellProcessPipe)), f"Invalid command type: {type(cmd)} in stage {self.name}"
                logger.info(cmd.description)
                logger.info(f" ❯ {cmd}")
                try:
                    self._run_command(cmd, ctx)
                except subprocess.CalledProcessError as e:
                    logger.error(f"Command failed with exit code {e.returncode}")
                    cmd = " ".join(quote(arg) for arg in e.cmd)
                    raise StageExecutionError(f"Failed to run command: {cmd}")
                except InvalidCommandTypeError as e:
                    raise e
        except (Exception, KeyboardInterrupt) as e:
            self._cleanup_incomplete_outputs(ctx)
            raise e

    def cleanup_tmp_outputs(self):
        """Delete temporary outputs."""
        dirs_to_remove = set()
        for name, path in self.output.temporary_outputs():
            if path and Path(path).exists():
                if path.is_dir():
                    dirs_to_remove.add(path)
                else:   
                    logger.debug(f"Removing temporary output '{name}' ({path}).")
                    Path(path).unlink()
        for d in sorted(dirs_to_remove, reverse=True):
            if d.exists():
                logger.debug(f"Removing temporary output directory ({d}).")
                try:
                    d.rmdir()
                except OSError:
                    logger.warning(f"Could not remove temporary directory {d} (not empty).")
                    continue

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

    def error_if_outputs_missing(self, include_temporary: bool = False):
        """Raises an error if any expected output files are missing."""
        missing = []
        outputs = self.output if include_temporary else self.output.non_temporary_outputs()
        for name, path in outputs:
            if not path:
                continue
            if not Path(path).exists():
                missing.append((name, path))
            else:
                logger.debug(f"Output '{name}' exists at {path}.")
        if missing:
            missing_str = ", ".join(f"{name} ({path})" for name, path in missing)
            raise MissingOutputError(f"Expected output files are missing: {missing_str}")

    def breakpoint(self):
        input("Press Enter to continue...")

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
