A pipeline is an ordered collection of stages that together perform a complete analysis. Pipelines are usually assembled by a `PipelineBuilder`, which decides which stages to include based on validated user input. Snippy-NG pipelines are linear: any branching happens while building the pipeline, not while executing it. For example:

```python
class ExamplePipelineBuilder(PipelineBuilder):
    input_file: Path = Field(..., description="Input file")
    prefix: str = Field(default="example", description="Output prefix")

    def build(self) -> SnippyPipeline:
        stages = []

        first = ExampleStage(
            input_file=self.input_file,
            prefix=self.prefix,
        )
        stages.append(first)

        second = AnotherStage(
            source=first.output.result,
            prefix=self.prefix,
        )
        stages.append(second)

        return SnippyPipeline(stages=stages)
```

In practice, a pipeline builder usually does three things:

- validates workflow-level inputs
- wires stage outputs into downstream stage inputs
- chooses which files should be kept after cleanup

Pipelines should contain workflow composition logic, not low-level implementation details. If a unit of work is independently meaningful, it should usually live in a stage and then be composed into one or more pipelines.
