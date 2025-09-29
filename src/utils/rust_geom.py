"""Optional wrapper for Rust geometry kernels.

Exposes pairwise_sq_dists if the compiled extension (molbridge_geom) is present.
Falls back to a NumPy implementation otherwise.
"""
from __future__ import annotations
import importlib
from typing import Optional

def pairwise_sq_dists(a, b):  # a,b: (m,3)/(n,3) float32 arrays
    mod: Optional[object]
    try:
        mod = importlib.import_module('molbridge_geom')
    except Exception:
        mod = None
    if mod and hasattr(mod, 'pairwise_sq_dists'):
        # Flatten to list for FFI simplicity
        m, n = a.shape[0], b.shape[0]
        return mod.pairwise_sq_dists(a.ravel().tolist(), b.ravel().tolist(), m, n)  # type: ignore
    # Fallback
    import numpy as np
    diff = a[:, None, :] - b[None, :, :]
    return (diff * diff).sum(axis=-1).ravel().tolist()

__all__ = ["pairwise_sq_dists"]
