from pathlib import Path

from snippy_ng.cli.utils import absolute_path, absolute_path_callback


class _Ctx:
    def __init__(self, resilient_parsing: bool = False):
        self.resilient_parsing = resilient_parsing


def test_absolute_path_handles_single_value(tmp_path):
    rel = tmp_path / "a.txt"
    result = absolute_path(str(rel))
    assert isinstance(result, Path)
    assert result.is_absolute()


def test_absolute_path_handles_tuple_values(tmp_path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    result = absolute_path((str(a), str(b)))
    assert isinstance(result, tuple)
    assert all(isinstance(x, Path) for x in result)
    assert all(x.is_absolute() for x in result)


def test_absolute_path_callback_resilient_parsing_returns_original_value(tmp_path):
    value = (str(tmp_path / "x"),)
    result = absolute_path_callback(_Ctx(resilient_parsing=True), None, value)
    assert result == value


def test_absolute_path_callback_converts_tuple_for_nargs(tmp_path):
    value = (str(tmp_path / "x"), str(tmp_path / "y"))
    result = absolute_path_callback(_Ctx(resilient_parsing=False), None, value)
    assert isinstance(result, tuple)
    assert all(isinstance(x, Path) for x in result)
    assert all(x.is_absolute() for x in result)
