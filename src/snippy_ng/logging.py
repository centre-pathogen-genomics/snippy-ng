import os
import click
import time


class Logger():
    def __init__(self):
        self.levels = {
            'INFO': click.style('INFO', fg='green', bold=True),
            'WARNING': click.style('WARNING', fg='yellow', bold=True),
            'DEBUG': click.style('DEBUG', fg='blue', bold=True),
            'ERROR': click.style('ERROR', fg='red', bold=True)
        }

    def _format(self, level, msg):
        """Format the log message."""
        asctime = time.strftime("%H:%M:%S")
        level_styled = self.levels.get(level, level)
        return f"[{asctime} - {level_styled}] {msg}"

    def echo(self, message='', err=False, **kwargs):
        """Echo a message to the console."""
        click.echo(message, err=err, **kwargs)

    def info(self, msg):
       self.echo(self._format("INFO", msg), err=True)
    
    def warning(self, msg):
        self.echo(self._format("WARNING", msg), err=True)
    
    def debug(self, msg):
        if os.getenv("SNIPPY_DEBUG", "0") == "1":
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