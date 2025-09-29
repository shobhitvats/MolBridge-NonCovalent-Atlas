"""Structured logging support (JSON) gated by environment variable.

Set MOLBRIDGE_LOG_FORMAT=json to emit JSON lines. Falls back to standard
logging if not set. Integrates with existing 'logging' module so external
handlers still function.
"""
from __future__ import annotations

import json
import logging
import os
import time
from typing import Any, Dict


class JsonFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:  # noqa: D401
        payload: Dict[str, Any] = {
            "ts": time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime(record.created)),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
        }
        if record.exc_info:
            payload["exc_info"] = self.formatException(record.exc_info)
        for k, v in getattr(record, 'extra', {}).items():  # type: ignore[attr-defined]
            if k not in payload:
                payload[k] = v
        return json.dumps(payload, ensure_ascii=False)


def enable_structured_logging():
    root = logging.getLogger()
    if not root.handlers:
        logging.basicConfig(level=logging.INFO)
    fmt = JsonFormatter()
    for h in root.handlers:
        h.setFormatter(fmt)


def configure_logging_if_requested():
    if os.getenv("MOLBRIDGE_LOG_FORMAT", "").lower() == "json":
        enable_structured_logging()


configure_logging_if_requested()
