"""Lightweight reusable workspace buffers for geometry operations.

Currently minimal; can be extended to implement shape-based pool reuse to avoid
repeated large allocations in hot vectorized detectors.
"""
from __future__ import annotations

from typing import Dict, Tuple
import numpy as _np

_POOL: Dict[Tuple[str, Tuple[int, ...], str], _np.ndarray] = {}

def get_workspace(name: str, shape: Tuple[int, ...], dtype: str = 'float32'):
    key = (name, shape, dtype)
    arr = _POOL.get(key)
    if arr is None:
        arr = _np.empty(shape, dtype=dtype)
        _POOL[key] = arr
    return arr

def clear_workspace(prefix: str | None = None):  # pragma: no cover - maintenance utility
    if prefix is None:
        _POOL.clear()
        return
    remove = [k for k in _POOL if k[0].startswith(prefix)]
    for k in remove:
        _POOL.pop(k, None)
