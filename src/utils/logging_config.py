"""Central logging configuration.

Usage:
    from utils.logging_config import configure_logging
    configure_logging(json_logs=False, level="INFO")

Idempotent: safe to call multiple times.
"""
from __future__ import annotations
import logging
import sys
import os
from logging.handlers import RotatingFileHandler
from typing import Optional, Dict, Any

try:
    import json_log_formatter  # type: ignore
except Exception:  # pragma: no cover
    json_log_formatter = None  # type: ignore

_CONFIGURED = False

class _PlainFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:  # noqa: D401
        return True

def configure_logging(json_logs: bool = False, level: str = "INFO") -> None:
    global _CONFIGURED
    if _CONFIGURED:
        return
    root = logging.getLogger()
    root.setLevel(level.upper())
    for h in list(root.handlers):
        root.removeHandler(h)
    # Console handler configuration
    if json_logs and json_log_formatter:
        console_handler = logging.StreamHandler(sys.stdout)
        console_formatter = json_log_formatter.JSONFormatter()
        console_handler.setFormatter(console_formatter)
    else:
        console_handler = logging.StreamHandler(sys.stdout)
        console_formatter = logging.Formatter(
            fmt="%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        console_handler.setFormatter(console_formatter)
    console_handler.addFilter(_PlainFilter())
    root.addHandler(console_handler)

    # Optional structured JSON file (rotating) controlled by env MOLBRIDGE_LOG_FILE
    log_file = os.getenv('MOLBRIDGE_LOG_FILE')
    if log_file:
        try:
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
        except Exception:
            pass
        file_handler = RotatingFileHandler(log_file, maxBytes=5_000_000, backupCount=3)
        if json_log_formatter:
            file_formatter = json_log_formatter.JSONFormatter()
        else:
            file_formatter = logging.Formatter(
                fmt="%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )
        file_handler.setFormatter(file_formatter)
        file_handler.addFilter(_PlainFilter())
        root.addHandler(file_handler)
    _CONFIGURED = True

def emit_metrics(record: Dict[str, Any]) -> None:
    """Emit a metrics JSON line to a dedicated metrics log file if configured.

    Controlled by env MOLBRIDGE_METRICS_FILE. Falls back to standard logger
    if metrics file handler not present. Idempotent lightweight helper.
    """
    import json, time
    logger = logging.getLogger('metrics')
    payload = {
        'ts': time.time(),
        **record
    }
    try:
        logger.info(json.dumps(payload))
    except Exception:
        logging.getLogger(__name__).debug('Failed to serialize metrics record')
    # If metrics file was configured after initial import (env set at runtime), ensure handler exists
    if os.getenv('MOLBRIDGE_METRICS_FILE') and not any(isinstance(h, RotatingFileHandler) for h in logging.getLogger('metrics').handlers):  # pragma: no cover
        try:
            _configure_metrics_channel()
            logger.info(json.dumps(payload))  # re-emit to file after channel setup
        except Exception:
            pass

# Configure metrics file handler lazily on import if env present
def _configure_metrics_channel():  # pragma: no cover
    if os.getenv('MOLBRIDGE_METRICS_FILE'):
        path = os.getenv('MOLBRIDGE_METRICS_FILE')
        try:
            os.makedirs(os.path.dirname(path), exist_ok=True)
        except Exception:
            pass
        metrics_logger = logging.getLogger('metrics')
        if not metrics_logger.handlers:
            handler = RotatingFileHandler(path, maxBytes=5_000_000, backupCount=2)
            handler.setFormatter(logging.Formatter('%(message)s'))
            metrics_logger.addHandler(handler)
            metrics_logger.setLevel(logging.INFO)

_configure_metrics_channel()
