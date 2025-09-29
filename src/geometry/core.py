"""High-performance geometry primitives (vectorized).

Centralized utilities to reduce duplicated O(N^2) distance / angle loops across
interaction detectors. Each function is NumPy-first with an optional SciPy KD-tree
fast path when available. All arrays are coerced to float32 to reduce memory
footprint and improve cache locality.

Design goals:
  * Reuse scratch buffers (via workspace module) where shapes are stable.
  * Provide unified API so detectors can migrate incrementally.
  * Avoid hard dependency on SciPy – silently fall back to pure NumPy.
"""
from __future__ import annotations

from typing import Tuple, Optional, List
import numpy as _np

try:  # optional SciPy
    from scipy.spatial import cKDTree as _KDTree  # type: ignore
    _HAVE_SCIPY = True
except Exception:  # pragma: no cover
    _KDTree = None  # type: ignore
    _HAVE_SCIPY = False

from .workspace import get_workspace

def as_f32(a):  # small helper
    return _np.asarray(a, dtype=_np.float32)

_KDTREE_CACHE: dict[tuple[int, int, float], tuple[object, int]] = {}

def _get_kdtree(arr: _np.ndarray, cutoff: float):  # simple reuse keyed by (id, len, cutoff_bucket)
    if not _HAVE_SCIPY:
        return None
    # Bucket cutoff to 0.1 Å to allow reuse across slight numeric variation
    bucket = round(float(cutoff) * 10.0) / 10.0
    key = (id(arr), arr.shape[0], bucket)
    entry = _KDTREE_CACHE.get(key)
    if entry is not None:
        tree, size = entry
        if size == arr.shape[0]:
            return tree
    tree = _KDTree(arr)
    _KDTREE_CACHE[key] = (tree, arr.shape[0])
    # Simple LRU-ish trim (bounded)
    if len(_KDTREE_CACHE) > 32:  # pragma: no cover - rarely hit in tests
        for _ in range(len(_KDTREE_CACHE) - 32):
            _KDTREE_CACHE.pop(next(iter(_KDTREE_CACHE)))
    return tree

def pairwise_within_cutoff(coords_a, coords_b, cutoff: float, use_kdtree: bool = True, return_distances: bool = False) -> Tuple[_np.ndarray, _np.ndarray, Optional[_np.ndarray]]:
    """Return index arrays (ia, ib) of pairs within distance <= cutoff.

    Heuristic: if NA * NB > 25_000 and SciPy available -> KD-tree per side.
    Otherwise use vectorized broadcasting (memory heavy but faster for small sets).
    """
    A = as_f32(coords_a); B = as_f32(coords_b)
    if A.size == 0 or B.size == 0:
        if return_distances:
            return _np.empty(0, dtype=_np.int32), _np.empty(0, dtype=_np.int32), _np.empty(0, dtype=_np.float32)
        return _np.empty(0, dtype=_np.int32), _np.empty(0, dtype=_np.int32), None
    na, nb = A.shape[0], B.shape[0]
    # Quick identical reference optimization
    same = A is B or (na == nb and A.base is B.base)
    threshold = 25_000
    if use_kdtree and _HAVE_SCIPY and na * nb > threshold:
        tree = _get_kdtree(B, cutoff)
        ia_buf = []
        ib_buf = []
        dist_buf = [] if return_distances else None
        for i, pt in enumerate(A):  # pragma: no cover (loop path not critical for correctness tests)
            idxs = tree.query_ball_point(pt, cutoff)
            if not idxs:
                continue
            if same:
                for j in idxs:
                    if j <= i:
                        continue
                    ia_buf.append(i); ib_buf.append(j)
                    if dist_buf is not None:
                        dist_buf.append(_np.linalg.norm(pt - B[j]))
            else:
                for j in idxs:
                    ia_buf.append(i); ib_buf.append(j)
                    if dist_buf is not None:
                        dist_buf.append(_np.linalg.norm(pt - B[j]))
        if not ia_buf:
            if return_distances:
                return (_np.empty(0, dtype=_np.int32), _np.empty(0, dtype=_np.int32), _np.empty(0, dtype=_np.float32))
            return _np.empty(0, dtype=_np.int32), _np.empty(0, dtype=_np.int32), None
        ia_arr = _np.asarray(ia_buf, dtype=_np.int32)
        ib_arr = _np.asarray(ib_buf, dtype=_np.int32)
        if return_distances:
            return ia_arr, ib_arr, _np.asarray(dist_buf, dtype=_np.float32) if dist_buf is not None else _np.empty(0, dtype=_np.float32)
        return ia_arr, ib_arr, None
    # Broadcasting path
    # Use workspace buffer if shapes stable (not reusing across variable dims yet)
    diff = A[:, None, :] - B[None, :, :]
    dist2 = _np.sum(diff * diff, axis=-1, dtype=_np.float32)
    mask = dist2 <= (cutoff * cutoff)
    if same:
        # keep upper triangle only
        iu, ju = _np.where(_np.triu(mask, k=1))
        if return_distances:
            d = _np.sqrt(dist2[iu, ju], dtype=_np.float32)
            return iu.astype(_np.int32), ju.astype(_np.int32), d
        return iu.astype(_np.int32), ju.astype(_np.int32), None
    i_idx, j_idx = _np.where(mask)
    if return_distances:
        d = _np.sqrt(dist2[i_idx, j_idx], dtype=_np.float32)
        return i_idx.astype(_np.int32), j_idx.astype(_np.int32), d
    return i_idx.astype(_np.int32), j_idx.astype(_np.int32), None

def norms(v: _np.ndarray) -> _np.ndarray:
    return _np.sqrt(_np.sum(v * v, axis=-1, dtype=_np.float32))

def cosine_angles(vecs_a: _np.ndarray, vecs_b: _np.ndarray, clamp: bool = True) -> _np.ndarray:
    """Compute cosine of angle between paired vectors (row-wise)."""
    a = as_f32(vecs_a); b = as_f32(vecs_b)
    num = _np.sum(a * b, axis=-1)
    den = norms(a) * norms(b)
    with _np.errstate(divide='ignore', invalid='ignore'):
        cos = _np.where(den > 0, num / den, 0.0)
    if clamp:
        cos = _np.clip(cos, -1.0, 1.0)
    return cos

def angle_degrees(cos_values: _np.ndarray) -> _np.ndarray:
    return _np.degrees(_np.arccos(cos_values))

def perpendicular_offset(vecs: _np.ndarray, normals: _np.ndarray) -> _np.ndarray:
    """Magnitude of component of vecs orthogonal to normals (paired rows)."""
    v = as_f32(vecs); n = as_f32(normals)
    proj = _np.sum(v * n, axis=-1, keepdims=True) * n
    perp = v - proj
    return norms(perp)

def vector_fields(indices_a: _np.ndarray, indices_b: _np.ndarray, coords_a: _np.ndarray, coords_b: _np.ndarray) -> Tuple[_np.ndarray, _np.ndarray, _np.ndarray]:
    """Return (vectors, distances, unit_vectors) for paired indices.

    Parameters
    ----------
    indices_a, indices_b : int32 arrays of same length
    coords_a, coords_b : (NA,3) and (NB,3) float arrays
    """
    if indices_a.size == 0:
        z = _np.zeros((0,3), dtype=_np.float32)
        return z, _np.zeros(0, dtype=_np.float32), z
    vecs = as_f32(coords_b[indices_b] - coords_a[indices_a])
    d = norms(vecs)
    with _np.errstate(divide='ignore', invalid='ignore'):
        uv = _np.where(d[:, None] > 0, vecs / d[:, None], 0.0)
    return vecs, d, uv

__all__ = [
    'pairwise_within_cutoff',
    'norms',
    'cosine_angles',
    'angle_degrees',
    'perpendicular_offset',
    'vector_fields'
]
