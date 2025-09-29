"""Lightweight timing utilities for instrumentation.

Provides:
  - time_block context manager
  - time_function decorator
  - global TimingCollector (thread-safe)

Design goals:
  * Near-zero overhead when disabled.
  * Thread-safe accumulation of wall time, call count, and optional item counts.
  * Simple export for embedding into analysis result summaries.
Enable/disable globally via env var MOLBRIDGE_ENABLE_TIMING=0/1 (default 1).
"""
from __future__ import annotations
import os
import time
import threading
from contextlib import contextmanager
from typing import Callable, Dict, Any, Optional

_ENABLE = os.getenv("MOLBRIDGE_ENABLE_TIMING", "1") != "0"
_lock = threading.Lock()

class TimingCollector:
    def __init__(self):
        self._data: Dict[str, Dict[str, float]] = {}

    def add(self, key: str, duration: float, items: Optional[int] = None):
        if not _ENABLE:
            return
        with _lock:
            rec = self._data.setdefault(key, {"time": 0.0, "calls": 0.0, "items": 0.0})
            rec["time"] += float(duration)
            rec["calls"] += 1.0
            if items is not None:
                rec["items"] += float(items)

    def snapshot(self) -> Dict[str, Any]:
        # Return a shallow copy safe for serialization
        with _lock:
            out = {}
            for k, rec in self._data.items():
                avg = rec["time"] / rec["calls"] if rec["calls"] else 0.0
                rate = rec["items"] / rec["time"] if rec["time"] and rec["items"] else None
                out[k] = {
                    "total_time": round(rec["time"], 6),
                    "calls": int(rec["calls"]),
                    "avg_time": round(avg, 6),
                    **({"total_items": int(rec["items"]) } if rec["items"] else {}),
                    **({"items_per_sec": round(rate, 3)} if rate else {}),
                }
            return out

    def clear(self):
        with _lock:
            self._data.clear()

TIMINGS = TimingCollector()

@contextmanager
def time_block(name: str, items: Optional[int] = None):
    if not _ENABLE:
        yield
        return
    start = time.perf_counter()
    try:
        yield
    finally:
        TIMINGS.add(name, time.perf_counter() - start, items=items)

def time_function(name: Optional[str] = None, items_attr: Optional[str] = None):
    """Decorator to time a function.

    Args:
        name: Logical timing bucket (default: function.__name__).
        items_attr: If provided and the function result has this attribute or key, its
                    numeric value will be added as processed item count.
    """
    def deco(fn: Callable):
        bucket = name or fn.__name__
        def wrapper(*args, **kwargs):
            if not _ENABLE:
                return fn(*args, **kwargs)
            start = time.perf_counter()
            result = fn(*args, **kwargs)
            dur = time.perf_counter() - start
            items_val = None
            if items_attr and result is not None:
                try:
                    if isinstance(result, dict) and items_attr in result:
                        items_val = float(result[items_attr])
                    else:
                        items_val = float(getattr(result, items_attr))
                except Exception:
                    items_val = None
            TIMINGS.add(bucket, dur, items=items_val)
            return result
        wrapper.__name__ = fn.__name__
        wrapper.__doc__ = fn.__doc__
        return wrapper
    return deco

__all__ = ["time_block", "time_function", "TIMINGS", "TimingCollector"]
