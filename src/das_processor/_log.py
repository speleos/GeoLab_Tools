"""Logging setup shared by the package: a log file, a console handler driven by
exec_params['verbose'], and a TRACE level below DEBUG."""
from __future__ import annotations

# --- Standard library --------------------------------------------------------
import datetime
import logging
import os
import sys
from typing import Any

_LOG_FILE = None

TRACE = 5   # finer than DEBUG: per-channel chatter, shown when exec_params["verbose"] == 2
logging.addLevelName(TRACE, "TRACE")

# Package-level logger: the console handler is attached here, so records from every module
# in the package (each uses logging.getLogger(__name__)) reach it.
logger = logging.getLogger(__name__.rpartition(".")[0] or __name__)
logger.setLevel(TRACE)       # the handlers, not the logger, decide what is shown
_CONSOLE_HANDLER = None


class _ConsoleFormatter(logging.Formatter):
    """Message only; warnings and errors get a level prefix unless the text already has one."""

    def format(self, record: logging.LogRecord) -> str:
        msg = super().format(record)
        if record.levelno >= logging.WARNING and not msg.upper().startswith(record.levelname):
            return f"{record.levelname}: {msg}"
        return msg


def _console_level(verbose: Any) -> int:
    """Console level for an exec_params ``verbose`` value: 0/False -> INFO, 1/True -> DEBUG, 2 -> TRACE."""
    try:
        v = int(verbose or 0)
    except (TypeError, ValueError):
        v = 1
    return TRACE if v >= 2 else logging.DEBUG if v == 1 else logging.INFO


def _ensure_console_handler(verbose: Any = None) -> None:
    """Show this module's log messages on stdout. Safe to call repeatedly.

    ``verbose`` (when given) sets how much is shown, see ``_console_level``.
    If the application already installed its own console handler on the root
    logger, that one is used instead and nothing is added here.
    """
    global _CONSOLE_HANDLER
    if _CONSOLE_HANDLER is None:
        if any(type(h) is logging.StreamHandler for h in logging.getLogger().handlers):
            return
        _CONSOLE_HANDLER = logging.StreamHandler(sys.stdout)
        _CONSOLE_HANDLER.setFormatter(_ConsoleFormatter())
        _CONSOLE_HANDLER.setLevel(logging.INFO)
        logger.addHandler(_CONSOLE_HANDLER)
    if verbose is not None:
        _CONSOLE_HANDLER.setLevel(_console_level(verbose))


def _setup_logging(log_dir: str | os.PathLike = ".", level: int = logging.INFO, verbose: Any = 0) -> str:
    """Configure logging: a log file (once per session) and the console.

    The console shows INFO and above, plus DEBUG when ``verbose`` is 1/True and
    TRACE as well when it is 2. The log file always receives DEBUG and above
    from this module. Safe to call repeatedly; returns the log file path.
    """
    global _LOG_FILE
    _ensure_console_handler(verbose)
    if _LOG_FILE is not None:
        return _LOG_FILE

    run_time = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    _LOG_FILE = os.path.join(log_dir, f'DAS_processor_{run_time}.log')

    handler = logging.FileHandler(_LOG_FILE)
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter(
        '%(asctime)s %(levelname)s: %(message)s', '%Y-%m-%d %H:%M:%S'))

    root = logging.getLogger()
    root.addHandler(handler)
    if root.level == logging.NOTSET or root.level > level:
        root.setLevel(level)
    return _LOG_FILE
