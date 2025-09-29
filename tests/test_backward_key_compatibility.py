"""Tests ensuring legacy interaction keys resolve to canonical forms.

Guards against breaking existing API/UI clients that use older key variants.
"""
from analysis.keys import resolve_key

LEGACY_TO_CANONICAL = {
    "hydrogenbond": "hydrogen_bond",
    "halogenbond": "halogen_bond",
    "chpi": "ch_pi",
    "pipi": "pi_pi",
    "ionicinteraction": "ionic_interaction",
    "hydrophobiccontact": "hydrophobic_contact",
    "chalcogenbond": "chalcogen_bond",
    "pnictogenbond": "pnictogen_bond",
    "tetrelbond": "tetrel_bond",
    "anionpi": "anion_pi",
    "npistar": "n_pi_star",
}

def test_legacy_key_resolution():
    for legacy, canonical in LEGACY_TO_CANONICAL.items():
        assert resolve_key(legacy) == canonical

def test_identity_mapping_preserved():
    for canonical in set(LEGACY_TO_CANONICAL.values()):
        assert resolve_key(canonical) == canonical