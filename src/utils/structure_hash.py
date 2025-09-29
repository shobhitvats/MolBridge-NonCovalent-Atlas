"""Structure coordinate hashing utilities.

Provides a lightweight, stable hash of all atomic coordinates of a Bio.PDB
`Structure` (or any object exposing `get_atoms()` iterator with atoms that
implement `get_coord()`).

Rationale:
  * Used to key persistent caching of already-computed interaction bundles or
    compact JSON blobs so that re-serialization can be skipped if neither the
    interaction parameter signature nor the underlying coordinates changed.
  * Float32 down‑casting is applied before hashing to guarantee stable byte
    representation across platforms and to reduce hashing cost / memory.

The produced hash is a truncated SHA256 hex string (default 16 chars) to keep
cache keys short while making collisions highly improbable for typical usage
sizes (16 hex chars ~ 64 bits of entropy).
"""
from __future__ import annotations

from typing import Iterable, Protocol
import hashlib
import numpy as np

class _HasCoord(Protocol):  # pragma: no cover - structural typing helper
    def get_coord(self): ...  # noqa: D401,E701 - simple protocol stub

def _iter_atom_coords(structure) -> Iterable["np.ndarray"]:
    """Yield coordinate arrays from a Structure-like object.

    Falls back gracefully returning empty iterator if object does not support
    `get_atoms` (in which case caller receives a deterministic hash of empty
    byte string).
    """
    get_atoms = getattr(structure, "get_atoms", None)
    if get_atoms is None:  # pragma: no cover - defensive
        return []
    try:
        for atom in get_atoms():  # type: ignore[call-arg]
            if atom is None:  # pragma: no cover
                continue
            get_coord = getattr(atom, "get_coord", None)
            if get_coord is None:  # pragma: no cover
                continue
            try:
                yield get_coord()
            except Exception:  # pragma: no cover
                continue
    except Exception:  # pragma: no cover
        return []

def compute_structure_coord_hash(structure, length: int = 16) -> str:
    """Return truncated SHA256 hash of float32-packed coordinate bytes.

    Parameters
    ----------
    structure : Any
        Bio.PDB Structure or compatible object.
    length : int, default 16
        Number of hex characters of the digest to return (<= 64).
    """
    # Reuse cached hash if present
    cached = getattr(structure, '_molbridge_coord_hash', None)
    if cached is not None and isinstance(cached, tuple):
        try:
            prev_len, prev_hash, prev_atom_count = cached
            # Fast path: if atom count stable, reuse existing (length may differ -> truncate)
            if isinstance(prev_atom_count, int):
                # We still confirm count quickly without materializing all coords again by attempting early stop
                count = 0
                for _ in _iter_atom_coords(structure):
                    count += 1
                    if count > prev_atom_count:
                        break
                if count == prev_atom_count:
                    return prev_hash[:length]
        except Exception:  # pragma: no cover
            pass
    coords = list(_iter_atom_coords(structure))
    if not coords:
        return "0" * min(length, 64)
    try:
        arr = np.asarray(coords, dtype="float32")  # shape (N,3) or (N,>=3)
    except Exception:  # pragma: no cover
        # Fallback: hash concatenated repr strings (rare path)
        h = hashlib.sha256("".join(map(repr, coords)).encode()).hexdigest()
        return h[:length]
    h = hashlib.sha256(arr.tobytes()).hexdigest()
    try:  # cache tuple (stored_length, full_hash, atom_count)
        setattr(structure, '_molbridge_coord_hash', (length, h, arr.shape[0]))
    except Exception:  # pragma: no cover
        pass
    return h[:length]

__all__ = ["compute_structure_coord_hash"]
