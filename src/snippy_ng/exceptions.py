class SnippyError(Exception):
    """Base class for all Snippy exceptions."""
    pass

class DependencyError(SnippyError):
    pass

class MissingDependencyError(DependencyError):
    pass

class InvalidDependencyError(DependencyError):
    pass

class InvalidDependencyVersionError(DependencyError):
    pass

class MissingInputError(SnippyError):
    pass

class MissingOutputError(SnippyError):
    pass

class InvalidCommandTypeError(SnippyError):
    pass

class StageExecutionError(SnippyError):
    pass

class StageTestFailure(StageExecutionError):
    pass

class PipelineExecutionError(SnippyError):
    pass

class InvalidReferenceError(SnippyError):
    pass
