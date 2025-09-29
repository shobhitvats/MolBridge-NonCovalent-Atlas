"""Unified hashing utilities.

Provides a stable, reasonably fast hash for heterogeneous Python objects.
Primary function: stable_hash(obj, length=12) returning hex digest prefix.

Priority order:
  1. xxhash (if installed) for speed
  2. blake2b (stdlib hashlib) fallback
  3. sha256 fallback
Serialization strategy: compact JSON with sorted keys for dict-like objects; for
non-JSON-serializable objects falls back to repr().
"""
from __future__ import annotations
from typing import Any
import json, hashlib

try:  # optional fast hash
    import xxhash  # type: ignore
    _HAVE_XXHASH = True
except Exception:  # pragma: no cover
    xxhash = None  # type: ignore
    _HAVE_XXHASH = False

def _to_stable_bytes(obj: Any) -> bytes:
    if obj is None:
        return b"null"
    # Fast path for common primitives
    if isinstance(obj, (str, bytes)):
        return obj.encode() if isinstance(obj, str) else obj
    if isinstance(obj, (int, float, bool)):
        return str(obj).encode()
    # Attempt JSON
    try:
        return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()
    except Exception:
        return repr(obj).encode()

def stable_hash(obj: Any, length: int = 12) -> str:
    data = _to_stable_bytes(obj)
    if _HAVE_XXHASH:
        h = xxhash.xxh64(data).hexdigest()
    else:  # stdlib fallback
        try:
            h = hashlib.blake2b(data, digest_size=16).hexdigest()
        except Exception:  # pragma: no cover
            h = hashlib.sha256(data).hexdigest()
    return h[:length]

__all__ = ["stable_hash"]
