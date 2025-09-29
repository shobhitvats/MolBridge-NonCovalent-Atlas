"""Vectorized angle & geometric filtering utilities.

Provides bulk dot-product -> angle computation avoiding per-pair Python loops.
"""
from __future__ import annotations
import numpy as np
from typing import Tuple

def angles_between(v1: np.ndarray, v2: np.ndarray, normalize: bool = True) -> np.ndarray:
    """Compute angles (degrees) between corresponding vectors in v1 and v2.

    v1, v2: shape (N,3). If shapes broadcast, the smaller is broadcast.
    Returns: shape (N,) array of angles in degrees.
    """
    if normalize:
        v1 = v1 / (np.linalg.norm(v1, axis=-1, keepdims=True) + 1e-12)
        v2 = v2 / (np.linalg.norm(v2, axis=-1, keepdims=True) + 1e-12)
    dots = (v1 * v2).sum(-1)
    # Clip and then snap values extremely close to +/-1 to exactly those to avoid acos drift
    dots = np.clip(dots, -1.0, 1.0)
    # Snap tolerance
    close_to_one = np.abs(1.0 - dots) < 1e-6
    close_to_neg_one = np.abs((-1.0) - dots) < 1e-6
    dots[close_to_one] = 1.0
    dots[close_to_neg_one] = -1.0
    ang = np.degrees(np.arccos(dots))
    return ang

def filter_by_angle(v1: np.ndarray, v2: np.ndarray, max_angle: float) -> np.ndarray:
    """Return boolean mask where angle between v1[i], v2[i] <= max_angle."""
    ang = angles_between(v1, v2)
    return ang <= max_angle

__all__ = ["angles_between", "filter_by_angle"]