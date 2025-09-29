"""Detector registry consolidating mapping of interaction keys to detector classes & method names.
Non-breaking: existing batch_processor still builds its own list; this layer can be adopted gradually.
"""
from __future__ import annotations
from typing import Dict, Type, Any

# Import detector classes (lazy-friendly but explicit for now)
from .hydrogen_bonds import HydrogenBondDetector
from .halogen_bonds import HalogenBondDetector
from .ch_pi_interactions import CHPiDetector
from .pi_pi_stacking import PiPiDetector
from .ionic_interactions import IonicInteractionDetector
from .hydrophobic_contacts import HydrophobicContactDetector
from .chalcogen_bonds import ChalcogenBondDetector
from .pnictogen_bonds import PnictogenBondDetector
from .tetrel_bonds import TetrelBondDetector
from .anion_pi_interactions import AnionPiDetector
from .n_pi_star_interactions import NPiStarDetector
from .london_dispersion import DispersionDetector
from .cation_pi_interactions import CationPiDetector
from .salt_bridges import SaltBridgeDetector
from .sulfur_pi_interactions import SulfurPiDetector
from .metal_coordination import MetalCoordinationDetector

# Key -> (Detector class, method name returning interactions)
DETECTOR_REGISTRY: Dict[str, tuple[type, str]] = {
    'hydrogenbond': (HydrogenBondDetector, 'detect_hydrogen_bonds'),
    'halogenbond': (HalogenBondDetector, 'detect_halogen_bonds'),
    'chpi': (CHPiDetector, 'detect_ch_pi_interactions'),
    'pipi': (PiPiDetector, 'detect_pi_pi_interactions'),
    'ionicinteraction': (IonicInteractionDetector, 'detect_ionic_interactions'),
    'hydrophobiccontact': (HydrophobicContactDetector, 'detect_hydrophobic_contacts'),
    'chalcogenbond': (ChalcogenBondDetector, 'detect_chalcogen_bonds'),
    'pnictogenbond': (PnictogenBondDetector, 'detect_pnictogen_bonds'),
    'tetrelbond': (TetrelBondDetector, 'detect_tetrel_bonds'),
    'anionpi': (AnionPiDetector, 'detect_anion_pi_interactions'),
    'npistar': (NPiStarDetector, 'detect_n_pi_star_interactions'),
    'dispersion': (DispersionDetector, 'detect_dispersion_interactions'),
    'cation_pi': (CationPiDetector, 'detect_cation_pi_interactions'),
    'salt_bridge': (SaltBridgeDetector, 'detect_salt_bridges'),
    'sulfur_pi': (SulfurPiDetector, 'detect_sulfur_pi_interactions'),
    'metal_coordination': (MetalCoordinationDetector, 'detect_metal_coordination')
}

def list_interaction_keys():
    return list(DETECTOR_REGISTRY.keys())

_def_fallback_config = None

def get_detector(key: str, config: Any | None = None):
    """Instantiate detector for key; returns (instance, method_name) or (None, None)."""
    entry = DETECTOR_REGISTRY.get(key)
    if not entry:
        return None, None
    cls, method_name = entry
    if config is None:
        global _def_fallback_config
        if _def_fallback_config is None:
            try:
                from utils.config import AppConfig
                _def_fallback_config = AppConfig()
            except Exception:
                _def_fallback_config = object()
        config = _def_fallback_config
    try:
        inst = cls(config)
    except Exception:
        inst = cls.__new__(cls)  # last resort, though detectors usually expect config
    return inst, method_name
