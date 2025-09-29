"""Utilities for stable interaction parameter hashing & cache keys.

Centralizes construction so batch processors / detectors can import without
dragging in full config logic each time.
"""
from __future__ import annotations

from typing import Dict
from utils.config import InteractionConfig

def interaction_params_signature(cfg: InteractionConfig) -> str:
    """Return the current interaction parameter signature (version + hash).

    Format: v{_version}:{hash}
    Example: v1:3fa09c1b2d
    """
    return f"v{cfg._version}:{cfg.param_hash}"  # noqa: SLF001 (intentional internal access)

def compose_cache_key(base: str, cfg: InteractionConfig, extras: Dict[str, str] | None = None) -> str:
    """Compose a cache key using a base label, param signature, and optional extras.

    Extras dict values should be short tokens (already hashed if large).
    """
    parts = [base, interaction_params_signature(cfg)]
    if extras:
        for k in sorted(extras):
            parts.append(f"{k}={extras[k]}")
    return "|".join(parts)
