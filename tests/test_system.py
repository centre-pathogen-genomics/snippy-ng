from snippy_ng.utils import system


def test_available_cpu_count_prefers_process_cpu_count(monkeypatch):
    monkeypatch.setattr(system.os, "process_cpu_count", lambda: 8, raising=False)
    monkeypatch.setattr(system.os, "sched_getaffinity", lambda _pid: {0, 1, 2, 3}, raising=False)
    monkeypatch.setattr(system.os, "cpu_count", lambda: 64, raising=False)

    assert system.available_cpu_count() == 8


def test_available_cpu_count_falls_back_to_sched_getaffinity(monkeypatch):
    monkeypatch.setattr(system.os, "process_cpu_count", lambda: None, raising=False)
    monkeypatch.setattr(system.os, "sched_getaffinity", lambda _pid: {0, 1, 2, 3, 4}, raising=False)
    monkeypatch.setattr(system.os, "cpu_count", lambda: 64, raising=False)

    assert system.available_cpu_count() == 5


def test_available_cpu_count_falls_back_to_cpu_count(monkeypatch):
    monkeypatch.setattr(system.os, "process_cpu_count", lambda: 0, raising=False)
    monkeypatch.setattr(system.os, "sched_getaffinity", lambda _pid: set(), raising=False)
    monkeypatch.setattr(system.os, "cpu_count", lambda: 12, raising=False)

    assert system.available_cpu_count() == 12
