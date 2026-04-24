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

## Environment variable-backed stage fields

Stages can declare fields whose default values come from environment variables by using `EnvVar` from [src/snippy_ng/envvars.py](src/snippy_ng/envvars.py). This is useful for runtime toggles that should remain part of the typed stage model rather than being hard-coded in `create_commands()`.

For example, [src/snippy_ng/stages/consequences.py](src/snippy_ng/stages/consequences.py) defines `use_local_csq` like this:

```python
from snippy_ng.envvars import EnvVar


class BcftoolsConsequencesCaller(BaseStage):
    reference: Path = Field(..., description="Reference file")
    variants: Path = Field(..., description="Input VCF file")
    features: Path = Field(..., description="Input features file")
    use_local_csq: bool = EnvVar(
        "LOCAL_BCFTOOLS_CSQ",
        default=True,
        description="Whether to use bcftools csq's --local-csq mode for consequence annotation",
    )
```

This field reads the prefixed environment variable:

- `SNIPPY_NG_LOCAL_BCFTOOLS_CSQ`

and parse the value based on the field default type. For booleans, accepted truthy values include `1`, `true`, `yes`, and `on`; falsy values include `0`, `false`, `no`, `off`, and the empty string.

Use `EnvVar` when:

- the option should appear as a normal typed field on the stage
- the value should be resolved when the stage model is instantiated
- you want the field to participate in validation, defaults, and model introspection like other Pydantic fields

Use the typed helpers in [src/snippy_ng/envvars.py](src/snippy_ng/envvars.py) (`bool_envvar()`, `int_envvar()`, `float_envvar()`, `str_envvar()`) when you need to read an environment variable outside a stage field declaration.