from pathlib import Path

from snippy_ng.context import Context
from snippy_ng.stages import BaseOutput, BaseStage


class DirOutput(BaseOutput):
    current_directory: Path


class FailingStage(BaseStage):
    @property
    def output(self) -> BaseOutput:
        return DirOutput(current_directory=Path("."))

    def create_commands(self, ctx):
        return [
            self.python_cmd(
                func=self._fail,
                description="Fail intentionally",
            )
        ]

    @staticmethod
    def _fail():
        raise AssertionError("boom")


def test_stage_cleanup_skips_current_working_directory(tmp_path):
    stage = FailingStage(prefix="snippy")

    try:
        stage.run(Context(outdir=tmp_path))
    except AssertionError as exc:
        assert str(exc) == "boom"
    else:
        raise AssertionError("Expected stage to fail")

    assert tmp_path.exists()
