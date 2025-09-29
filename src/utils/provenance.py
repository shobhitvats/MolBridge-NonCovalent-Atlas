"""Provenance utilities for optional reproducibility metadata.

All functions are lightweight and only used when enable_provenance flag is set
in Settings (environment variable: MOLBRIDGE_ENABLE_PROVENANCE=1).

We avoid importing heavy libs; only stdlib hashlib used. Hashes are short
(to avoid bloating cache keys) but stable.
"""
from __future__ import annotations
from typing import Any, Dict
import hashlib
import json

SCHEMA_VERSION = "1.0.0"  # Increment if normalized interaction shape changes

_DEF_ENCODER_SEP = "|"  # simple separator to join parts before hashing


def _stable_json(data: Any) -> str:
    """Serialize to a stable JSON string (sorted keys, no whitespace)."""
    try:
        return json.dumps(data, sort_keys=True, separators=(",", ":"))
    except Exception:
        return repr(data)


def hash_params(params: Dict[str, Any]) -> str:
    """Return short hash (10 hex chars) of parameter dictionary."""
    h = hashlib.sha256(_stable_json(params).encode("utf-8")).hexdigest()
    return h[:10]


def hash_structure_identifier(pdb_id: str, assembly: str | int | None = None) -> str:
    """Hash structural identity (pdb id + assembly)."""
    base = f"{pdb_id.lower()}{_DEF_ENCODER_SEP}{assembly if assembly is not None else 'default'}"
    return hashlib.sha256(base.encode("utf-8")).hexdigest()[:10]


def build_provenance(pdb_id: str, params: Dict[str, Any], assembly: str | int | None = None) -> Dict[str, Any]:
    """Assemble provenance block.

    Returns dict with deterministic keys; safe to merge into result metadata.
    """
    return {
        "schema_version": SCHEMA_VERSION,
        "structure_hash": hash_structure_identifier(pdb_id, assembly),
        "params_hash": hash_params(params),
        "param_count": len(params),
    }
