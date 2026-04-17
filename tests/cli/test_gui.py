import builtins
import sys
import types

from tests.cli.helpers import run_cli_command


def test_gui_help_is_available_without_gradio():
    result = run_cli_command(["gui", "--help"])

    assert result.exit_code == 0
    assert "Launch the Snippy-NG graphical interface" in result.output
    assert "--no-browser" in result.output
    assert "--temp-output" in result.output
    assert "--no-temp-output" not in result.output
    assert "--max-cpus" in result.output
    assert "--no-server-paths" in result.output


def test_gui_disables_temp_output_by_default(monkeypatch):
    captured = {}

    class FakeApp:
        def launch(self, **kwargs):
            captured["launch"] = kwargs

    def fake_create_app(**kwargs):
        captured["create_app"] = kwargs
        return FakeApp()

    fake_gui = types.ModuleType("snippy_ng.gui")
    fake_gui.create_app = fake_create_app
    monkeypatch.setitem(sys.modules, "gradio", types.ModuleType("gradio"))
    monkeypatch.setitem(sys.modules, "snippy_ng.gui", fake_gui)

    result = run_cli_command(["gui", "--no-browser"])

    assert result.exit_code == 0, result.output
    assert captured["create_app"]["temp_output"] is False
    assert captured["create_app"]["max_cpus"] is None
    assert captured["launch"]["inbrowser"] is False


def test_gui_temp_output_enables_temp_output(monkeypatch):
    captured = {}

    class FakeApp:
        def launch(self, **kwargs):
            captured["launch"] = kwargs

    def fake_create_app(**kwargs):
        captured["create_app"] = kwargs
        return FakeApp()

    fake_gui = types.ModuleType("snippy_ng.gui")
    fake_gui.create_app = fake_create_app
    monkeypatch.setitem(sys.modules, "gradio", types.ModuleType("gradio"))
    monkeypatch.setitem(sys.modules, "snippy_ng.gui", fake_gui)

    result = run_cli_command(["gui", "--no-browser", "--temp-output"])

    assert result.exit_code == 0, result.output
    assert captured["create_app"]["temp_output"] is True


def test_gui_max_cpus_is_passed_to_app(monkeypatch):
    captured = {}

    class FakeApp:
        def launch(self, **kwargs):
            captured["launch"] = kwargs

    def fake_create_app(**kwargs):
        captured["create_app"] = kwargs
        return FakeApp()

    fake_gui = types.ModuleType("snippy_ng.gui")
    fake_gui.create_app = fake_create_app
    monkeypatch.setitem(sys.modules, "gradio", types.ModuleType("gradio"))
    monkeypatch.setitem(sys.modules, "snippy_ng.gui", fake_gui)

    result = run_cli_command(["gui", "--no-browser", "--max-cpus", "12"])

    assert result.exit_code == 0, result.output
    assert captured["create_app"]["max_cpus"] == 12


def test_gui_reports_missing_gradio(monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "gradio":
            raise ModuleNotFoundError("No module named 'gradio'", name="gradio")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    result = run_cli_command(["gui"])

    assert result.exit_code == 1
    assert "requires the optional Gradio dependency" in result.output
    assert "pip install 'snippy-nextgen[gui]'" in result.output
    assert "pip install gradio" in result.output


def test_gui_reports_missing_gradio_dependency(monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "gradio":
            raise ModuleNotFoundError("No module named 'pytz'", name="pytz")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    result = run_cli_command(["gui"])

    assert result.exit_code == 1
    assert "requires the optional Gradio dependency" in result.output
    assert "Missing Python module: pytz" in result.output
    assert "pip install 'snippy-nextgen[gui]'" in result.output
