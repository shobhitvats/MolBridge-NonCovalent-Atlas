"""Geometry utility scaffold (non-breaking, not yet integrated).

Provides vectorized distance and angle helpers; optional KD-tree when scipy available.
"""
from __future__ import annotations
from typing import Tuple, Optional
import math
import numpy as np
import os

_USE_NUMBA = os.getenv('MOLBRIDGE_USE_NUMBA', '0') in {'1','true','True'}
_HAVE_NUMBA = False
if _USE_NUMBA:
    try:  # pragma: no cover - optional acceleration
        from numba import njit, prange  # type: ignore
        _HAVE_NUMBA = True
    except Exception:  # pragma: no cover
        _HAVE_NUMBA = False

try:  # optional dependency
    from scipy.spatial import cKDTree  # type: ignore
except Exception:  # pragma: no cover
    cKDTree = None  # type: ignore

Array = np.ndarray


if _HAVE_NUMBA:
    @njit(parallel=True, fastmath=True)
    def _pairwise_distances_numba(a: Array, b: Array) -> Array:  # type: ignore[misc]
        n = a.shape[0]
        m = b.shape[0]
        out = np.empty((n, m), dtype=np.float32)
        for i in prange(n):
            ax = a[i,0]; ay = a[i,1]; az = a[i,2]
            for j in range(m):
                dx = ax - b[j,0]
                dy = ay - b[j,1]
                dz = az - b[j,2]
                out[i,j] = (dx*dx + dy*dy + dz*dz) ** 0.5
        return out

def pairwise_distances(a: Array, b: Array) -> Array:
    """Compute full pairwise distances between two point sets.
    a: (N,3), b: (M,3) -> (N,M)
    Falls back to numpy broadcast if numba disabled or unavailable.
    """
    if _HAVE_NUMBA:
        return _pairwise_distances_numba(a, b)
    diff = a[:, None, :] - b[None, :, :]
    return np.linalg.norm(diff, axis=-1)


def cosine_angle(v1: Array, v2: Array) -> float:
    n1 = v1 / (np.linalg.norm(v1) + 1e-12)
    n2 = v2 / (np.linalg.norm(v2) + 1e-12)
    return float(np.clip(np.dot(n1, n2), -1.0, 1.0))


def angle_degrees(v1: Array, v2: Array) -> float:
    return math.degrees(math.acos(cosine_angle(v1, v2)))


def build_kdtree(coords: Array):  # pragma: no cover - tree building trivial
    if cKDTree is None:
        return None
    return cKDTree(coords)


def radius_query(tree, point: Array, radius: float):  # pragma: no cover
    if tree is None:
        return None
    return tree.query_ball_point(point, r=radius)
