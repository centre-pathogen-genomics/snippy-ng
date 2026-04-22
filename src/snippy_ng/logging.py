import os
import click
import time
from pathlib import Path

from snippy_ng.envvars import bool_envvar


DEBUG = bool_envvar(
    "DEBUG",
    description="Enable debug-oriented runtime behavior and verbose diagnostics",
)


def derive_log_path(log_path: Path | None, outdir: Path | None) -> Path | None:
    if log_path is None:
        return None
    log_path = Path(log_path)
    if outdir is None:
        return log_path.absolute()
    return (Path(outdir) / log_path.name).absolute()


class Logger():
    levels = {
            'INFO': click.style('INFO', fg='green', bold=True),
            'WARNING': click.style('WARNING', fg='yellow', bold=True),
            'DEBUG': click.style('DEBUG', fg='blue', bold=True),
            'ERROR': click.style('ERROR', fg='red', bold=True)
        }

    def __init__(self, log_path: Path | None = None):
        self._log_path: Path | None = Path(log_path) if log_path is not None else None
        self._last_log_path: Path | None = self._log_path

    def _format(self, level, msg):
        """Format the log message."""
        asctime = time.strftime("%H:%M:%S")
        level_styled = self.levels.get(level, level)
        return f"[{asctime} - {level_styled}] {msg}"

    def get_log_path(self) -> Path | None:
        return self._log_path

    def get_last_log_path(self) -> Path | None:
        return self._last_log_path

    def set_log_path(self, log_path: Path | None) -> None:
        self._log_path = Path(log_path) if log_path is not None else None
        if self._log_path is not None:
            self._last_log_path = self._log_path

    def reset_log_file(self) -> None:
        log_path = self.get_log_path()
        if log_path is None:
            return
        log_path.parent.mkdir(parents=True, exist_ok=True)
        if log_path.exists():
            log_path.unlink()

    def _append_to_file(self, message: str) -> None:
        log_path = self.get_log_path()
        if log_path is None:
            return
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(click.unstyle(message))
            if not message.endswith("\n"):
                handle.write("\n")

    def echo(self, message='', err=False, console=True, **kwargs):
        """Echo a message to the console."""
        if console:
            click.echo(message, err=err, **kwargs)
        self._append_to_file(message)

    def info(self, msg):
       self.echo(self._format("INFO", msg), err=True)
    
    def warning(self, msg):
        self.echo(self._format("WARNING", msg), err=True)

    def is_debug(self):
        return DEBUG.enabled

    def debug(self, msg):
        if self.is_debug():
            self.echo(self._format("DEBUG", msg), err=True)
    
    def error(self, msg):
        self.echo(self._format("ERROR", msg), err=True)

    def horizontal_rule(self, msg = "", style: str = '=', color: str = None):
        """Create a horizontal rule with a message in the middle."""
        try:
            terminal_width = os.get_terminal_size().columns
        except OSError:
            terminal_width = 80
        msg_length = len(msg)
        if msg_length != 0:
            msg = f" {msg} "
            msg_length += 2  # account for spaces around the message
        left_padding = max((terminal_width - msg_length) // 2, 5)  # ensure at least 5 characters padding
        right_padding = max(terminal_width - left_padding - msg_length, 5)  # ensure at least 5 characters padding
        
        if color is not None:
            msg = click.style(msg, fg=color)
        
        line = f"{style * left_padding}{msg}{style * right_padding}"
        self.echo(line, err=True)


logger = Logger()
