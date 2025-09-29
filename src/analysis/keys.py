"""Canonical interaction keys & alias resolution (non-breaking).

Existing code continues to use legacy keys; this module provides a resolver
for future integration while preserving current behavior.
"""
from __future__ import annotations
from typing import Dict

# Canonical snake_case keys
CANONICAL_KEYS = {
    "hydrogen_bond": "hydrogen_bond",
    "halogen_bond": "halogen_bond",
    "ch_pi": "ch_pi",
    "pi_pi": "pi_pi",
    "ionic_interaction": "ionic_interaction",
    "hydrophobic_contact": "hydrophobic_contact",
    "chalcogen_bond": "chalcogen_bond",
    "pnictogen_bond": "pnictogen_bond",
    "tetrel_bond": "tetrel_bond",
    "anion_pi": "anion_pi",
    "n_pi_star": "n_pi_star",
    "dispersion": "dispersion",
    "cation_pi": "cation_pi",
    "salt_bridge": "salt_bridge",
    "sulfur_pi": "sulfur_pi",
    "metal_coordination": "metal_coordination",
}

# Legacy / UI / API variants mapped to canonical
ALIASES: Dict[str, str] = {
    # legacy concatenated forms
    "hydrogenbond": "hydrogen_bond",
    "halogenbond": "halogen_bond",
    "ionicinteraction": "ionic_interaction",
    "hydrophobiccontact": "hydrophobic_contact",
    "chalcogenbond": "chalcogen_bond",
    "pnictogenbond": "pnictogen_bond",
    "tetrelbond": "tetrel_bond",
    "anionpi": "anion_pi",
    "npistar": "n_pi_star",
    # alternate underscores vs concatenation already canonical
    "pipi": "pi_pi",
    "chpi": "ch_pi",
    # already canonical but include identity mapping for clarity
    **{k: v for k, v in CANONICAL_KEYS.items()}
}

def resolve_key(key: str) -> str:
    """Resolve arbitrary key to canonical form (fallback: original)."""
    if not key:
        return key
    lowered = key.lower().strip()
    return ALIASES.get(lowered, lowered)
