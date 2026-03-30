A stage is the smallest execution unit in Snippy-NG. A stage should do one job, declare its inputs and outputs explicitly, and be reusable across pipelines. Stages inherit from `BaseStage`, define a typed `output` property, and implement `create_commands()` to return shell commands, shell pipelines, or Python functions to run. For example:

```python
class ExampleStageOutput(BaseOutput):
    result: Path = Field(..., description="Output file")


class ExampleStage(BaseStage):
    input_file: Path = Field(..., description="Input file")

    @property
    def output(self) -> ExampleStageOutput:
        return ExampleStageOutput(
            result=Path(f"{self.prefix}.result.txt")
        )

    def create_commands(self, ctx) -> list[ShellCommand]:
        return [
            self.shell_cmd(
                ["cp", str(self.input_file), str(self.output.result)],
                description="Copy input to output",
            )
        ]
```

In practice, a stage usually contains:

- typed input fields
- a typed output model
- optional dependency declarations in `_dependencies`
- one or more commands in `create_commands()`
- optional `test_*` methods for post-run validation

## Context

The `ctx` argument passed to `create_commands()` is a runtime `Context` object. This contains execution settings that are not part of the stage definition itself, such as:

- `outdir`
- `tmpdir`
- `cpus`
- `ram`
- `quiet`
- `skip_check`

Stages should use `ctx` for runtime behavior such as thread counts, temporary file locations, or memory limits. User inputs that define what the stage does should remain explicit stage fields, while execution settings that affect how the stage runs should come from `Context`.

Stages should not contain top-level CLI parsing logic. If a task can be reused in more than one workflow, it should usually be a stage.
