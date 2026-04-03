import sys
from pathlib import Path

from snippy_ng.context import Context
from snippy_ng.logging import logger
from snippy_ng.stages import BaseOutput, BaseStage


class NoOutput(BaseOutput):
    pass


class EchoStage(BaseStage):
    @property
    def output(self) -> BaseOutput:
        return NoOutput()

    def create_commands(self, ctx):
        return [
            self.shell_cmd(
                [
                    sys.executable,
                    "-c",
                    "import sys; print('hello from stdout'); print('hello from stderr', file=sys.stderr)",
                ],
                description="Emit test output",
            )
        ]


class PythonEchoStage(BaseStage):
    @property
    def output(self) -> BaseOutput:
        return NoOutput()

    def create_commands(self, ctx):
        return [
            self.python_cmd(
                func=self._emit_output,
                description="Emit Python command output",
            )
        ]

    @staticmethod
    def _emit_output():
        print("python stdout")
        print("python stderr", file=sys.stderr)


class BinaryPipelineStage(BaseStage):
    @property
    def output(self) -> BaseOutput:
        return NoOutput()

    def create_commands(self, ctx):
        return [
            self.shell_pipe(
                commands=[
                    self.shell_cmd(
                        [
                            sys.executable,
                            "-c",
                            "import gzip, sys; sys.stdout.buffer.write(gzip.compress(b'hello'))",
                        ],
                        description="Emit gzipped bytes",
                    )
                ],
                description="Write binary pipeline output",
                output_file=Path("binary.out.gz"),
            )
        ]


def test_stage_output_is_written_to_log_file(tmp_path):
    log_path = tmp_path / "LOG.txt"
    previous_log_path = logger.get_log_path()
    try:
        logger.set_log_path(log_path)
        logger.reset_log_file()

        stage = EchoStage(prefix="snippy")
        stage.run(Context(outdir=tmp_path, log_path=log_path))

        contents = log_path.read_text()
        assert "hello from stdout" in contents
        assert "hello from stderr" in contents
    finally:
        logger.set_log_path(previous_log_path)


def test_python_command_output_is_written_to_log_file(tmp_path):
    log_path = tmp_path / "LOG.txt"
    previous_log_path = logger.get_log_path()
    try:
        logger.set_log_path(log_path)
        logger.reset_log_file()

        stage = PythonEchoStage(prefix="snippy")
        stage.run(Context(outdir=tmp_path, log_path=log_path))

        contents = log_path.read_text()
        assert "python stdout" in contents
        assert "python stderr" in contents
    finally:
        logger.set_log_path(previous_log_path)


def test_binary_pipeline_output_file_does_not_decode_stdout_as_utf8(tmp_path):
    log_path = tmp_path / "LOG.txt"
    previous_log_path = logger.get_log_path()
    cwd = Path.cwd()
    try:
        logger.set_log_path(log_path)
        logger.reset_log_file()
        import os
        os.chdir(tmp_path)

        stage = BinaryPipelineStage(prefix="snippy")
        stage.run(Context(outdir=tmp_path, log_path=log_path))

        output_path = tmp_path / "binary.out.gz"
        assert output_path.exists()
        assert output_path.read_bytes().startswith(b"\x1f\x8b")
    finally:
        logger.set_log_path(previous_log_path)
        if Path.cwd() != cwd:
            import os
            os.chdir(cwd)
