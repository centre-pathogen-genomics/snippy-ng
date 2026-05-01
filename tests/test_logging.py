from pathlib import Path

from snippy_ng.logging import DEBUG, Logger, derive_log_path


def test_logger_resets_existing_log_file(tmp_path):
    log_path = tmp_path / "run.log"
    logger = Logger()
    logger.set_log_path(log_path)
    log_path.write_text("old contents\n")
    logger.reset_log_file()

    logger.info("hello")
    logger.warning("world")

    contents = log_path.read_text()
    assert "old contents" not in contents
    assert "hello" in contents
    assert "world" in contents


def test_derive_log_path_uses_outdir_and_basename():
    base = Path("/tmp/custom.log")
    outdir = Path("/work/core")

    assert derive_log_path(base, outdir) == (outdir / "custom.log").absolute()


def test_derive_log_path_without_outdir_is_absolute(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    assert derive_log_path(Path("LOG.txt"), None) == (tmp_path / "LOG.txt").absolute()


def test_debug_feature_has_expected_default(monkeypatch):
    monkeypatch.delenv("SNIPPY_NG_DEBUG", raising=False)

    assert DEBUG.enabled is False
