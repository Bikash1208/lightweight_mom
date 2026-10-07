"""Phase 13 Part 2: one logging setup for the whole project (solver core stays free of any
``logging`` calls of its own -- this is opt-in instrumentation the app and any script can use, not
something threaded through `sparameters.py`/`solver.py`).

- ``configure_logging()`` sets up a single ``"lwmom"`` logger: a plain console handler plus a
  rotating file handler under a per-user data directory (``~/.lwmom/logs/lwmom.log``, 5 MB x 3
  backups). Level is controlled by the ``LWMOM_LOG_LEVEL`` environment variable (default
  ``INFO``), or an explicit argument.
- ``JsonFormatter`` emits one JSON object per line (timestamp, level, logger, message, run_id, and
  any ``extra=`` fields) -- used for the file handler so diagnostic bundles are machine-parseable;
  the console handler stays human-readable.
- **Run ID propagation**: a solve invoked from the app (or any script) can call
  ``set_run_id(new_run_id())`` once; every subsequent log record on any logger obtained via
  ``get_logger`` automatically carries that ``run_id`` (a `contextvars.ContextVar`, safe across
  threads/async tasks), without the solver's own functions needing a ``run_id`` parameter -- exactly
  the "propagated through solver calls" requirement without touching the solver's signatures.
"""

import contextvars
import json
import logging
import logging.handlers
import os
import uuid
from pathlib import Path

_LOGGER_NAME = "lwmom"
_LOG_DIR = Path.home() / ".lwmom" / "logs"
_LOG_FILE = _LOG_DIR / "lwmom.log"
_MAX_BYTES = 5 * 1024 * 1024
_BACKUP_COUNT = 3

_run_id_var: contextvars.ContextVar[str | None] = contextvars.ContextVar("lwmom_run_id", default=None)

_configured = False


def new_run_id():
    """A fresh run ID (UUID4 hex, no dashes) -- call once per job/solve invocation."""
    return uuid.uuid4().hex


def set_run_id(run_id):
    """Binds `run_id` to every log record emitted from this context (thread/async task) until
    changed again. Returns the `contextvars.Token` so a caller can reset it when the job ends."""
    return _run_id_var.set(run_id)


def reset_run_id(token):
    _run_id_var.reset(token)


def current_run_id():
    return _run_id_var.get()


class _RunIdFilter(logging.Filter):
    def filter(self, record):
        record.run_id = current_run_id()
        return True


class JsonFormatter(logging.Formatter):
    def format(self, record):
        payload = {
            "timestamp": self.formatTime(record, "%Y-%m-%dT%H:%M:%S"),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
            "run_id": getattr(record, "run_id", None),
        }
        if record.exc_info:
            payload["exception"] = self.formatException(record.exc_info)
        # Any caller-supplied `extra={...}` fields (skip the standard LogRecord attributes).
        standard = set(logging.makeLogRecord({}).__dict__) | {"message", "asctime", "run_id"}
        for key, value in record.__dict__.items():
            if key not in standard:
                payload[key] = value
        return json.dumps(payload, default=str)


def configure_logging(level=None, log_dir=None):
    """Idempotent: safe to call more than once (e.g. once per app worker); re-configuring with a
    different `level` updates the existing handlers rather than duplicating them."""
    global _configured
    level_name = level or os.environ.get("LWMOM_LOG_LEVEL", "INFO")
    resolved_level = getattr(logging, str(level_name).upper(), logging.INFO)

    logger = logging.getLogger(_LOGGER_NAME)
    logger.setLevel(resolved_level)

    if _configured:
        for handler in logger.handlers:
            handler.setLevel(resolved_level)
        return logger

    log_dir = Path(log_dir) if log_dir else _LOG_DIR
    log_dir.mkdir(parents=True, exist_ok=True)

    console = logging.StreamHandler()
    console.setFormatter(logging.Formatter("%(asctime)s %(levelname)-8s %(name)s [%(run_id)s] %(message)s"))
    console.addFilter(_RunIdFilter())
    console.setLevel(resolved_level)

    file_handler = logging.handlers.RotatingFileHandler(
        log_dir / "lwmom.log",
        maxBytes=_MAX_BYTES,
        backupCount=_BACKUP_COUNT,
    )
    file_handler.setFormatter(JsonFormatter())
    file_handler.addFilter(_RunIdFilter())
    file_handler.setLevel(resolved_level)

    logger.addHandler(console)
    logger.addHandler(file_handler)
    logger.propagate = False
    _configured = True
    return logger


def get_logger(name=None):
    """A child of the `"lwmom"` logger (`"lwmom.<name>"`), or the root `"lwmom"` logger itself if
    `name` is omitted. Calls `configure_logging()` with defaults if nothing has configured it yet,
    so a bare `get_logger(__name__).info(...)` always works."""
    if not _configured:
        configure_logging()
    return logging.getLogger(f"{_LOGGER_NAME}.{name}" if name else _LOGGER_NAME)


def log_file_path(log_dir=None):
    return (Path(log_dir) if log_dir else _LOG_DIR) / "lwmom.log"
