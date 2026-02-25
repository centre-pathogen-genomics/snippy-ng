All CLI commands follow a common structure. They are defined as Click commands that use the `CommandWithGlobals` class and included the `add_snippy_global_options` decorator to add the global options. The command function then imports the appropriate pipeline builder, constructs the pipeline with the necessary arguments, and runs it with a typed runtime context. For example:
```python
@click.command(cls=CommandWithGlobals, context_settings={"show_default": True})
@add_snippy_global_options()
# ...pipeline-specific options...
def <name>(..., prefix, **context):
    from snippy_ng.context import Context
    from snippy_ng.pipelines.<name> import <Name>PipelineBuilder

    pipeline = <Name>PipelineBuilder(
        # pipeline-specific args only
        ...,
        prefix=prefix,
    ).build()

    run_ctx = Context(**context)
    pipeline.run(run_ctx)
```