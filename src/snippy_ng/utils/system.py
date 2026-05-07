import os


def available_cpu_count() -> int:
    process_cpu_count = getattr(os, "process_cpu_count", None)
    if callable(process_cpu_count):
        try:
            value = process_cpu_count()
            if value and value > 0:
                return int(value)
        except (OSError, ValueError):
            pass

    sched_getaffinity = getattr(os, "sched_getaffinity", None)
    if callable(sched_getaffinity):
        try:
            value = len(sched_getaffinity(0))
            if value > 0:
                return int(value)
        except (OSError, ValueError):
            pass

    value = os.cpu_count() or 1
    return max(1, int(value))
