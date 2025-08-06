from .stats import SeqKitReadStats, SeqKitReadStatsBasic, SeqKitReadStatsDetailed
from .downsample_reads import RasusaDownsampleReads, RasusaDownsampleReadsByCoverage, RasusaDownsampleReadsByCount


class AtRunTime:
    """Wrapper for values that should be computed at runtime rather than initialization.
    
    This class allows deferring the execution of functions until the value is actually
    needed, typically when converting to string or accessing the value during command
    generation.
    
    Args:
        func: Callable that returns the actual value when called.
        
    Example:
        >>> def get_genome_length():
        ...     # This will only be called when the value is needed
        ...     return read_genome_length_from_file()
        >>> genome_length = at_run_time(get_genome_length)
        >>> # Function not called yet
        >>> str(genome_length)  # Now function is called
        '197394'
    """
    
    def __init__(self, func):
        self._func = func
        self._cached_value = None
        self._computed = False
    
    def __str__(self):
        """Convert to string, computing the value if needed."""
        if not self._computed:
            self._cached_value = self._func()
            self._computed = True
        return str(self._cached_value)
    
    def __int__(self):
        """Convert to int, computing the value if needed."""
        if not self._computed:
            self._cached_value = self._func()
            self._computed = True
        return int(self._cached_value)
    
    def __float__(self):
        """Convert to float, computing the value if needed.
        
        This method can be called by Pydantic during validation.
        We'll allow it to compute the value since it's needed for type conversion.
        """
        if not self._computed:
            self._cached_value = self._func()
            self._computed = True
        return float(self._cached_value)
    
    def __bool__(self):
        """Return True without triggering computation (AtRunTime objects are always truthy)."""
        return True

    @property
    def value(self):
        """Get the computed value."""
        if not self._computed:
            self._cached_value = self._func()
            self._computed = True
        return self._cached_value


def at_run_time(func):
    """Create an AtRunTime wrapper for deferred function execution.
    
    Args:
        func: Callable that returns the value when called.
        
    Returns:
        AtRunTime: Wrapper that will call func when the value is needed.
        
    Example:
        >>> def get_value():
        ...     return read_from_file()
        >>> deferred_value = at_run_time(get_value)
        >>> # Value not computed yet
        >>> print(f"Value is: {deferred_value}")  # Now computed
    """
    return AtRunTime(func)


__all__ = [
    "SeqKitReadStats",
    "SeqKitReadStatsBasic", 
    "SeqKitReadStatsDetailed",
    "RasusaDownsampleReads",
    "RasusaDownsampleReadsByCoverage",
    "RasusaDownsampleReadsByCount",
    "at_run_time",
    "AtRunTime"
]