"""Central detector metadata registry.

Provides a single source of truth for mapping detector registry keys
(corresponding to @register_detector identifiers) to auxiliary metadata:
- human label
- adaptive threshold key (used by utils.kdtree_thresholds)
- environment variable name that can override KD-tree threshold (when applicable)
- short description

This indirection removes ad‑hoc string literals scattered across detector
implementations and enables higher‑level aggregation/visualization code
to reason about detectors uniformly.

Adding a detector:
1. Ensure it is registered via @register_detector(<key>, method=...)
2. Add an entry below with at least a label. Threshold metadata is optional.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Optional, Iterable

@dataclass(frozen=True)
class DetectorMeta:
    key: str                  # registry key
    label: str                # human readable label
    threshold_key: Optional[str] = None  # key used by kdtree_thresholds module
    env_threshold_var: Optional[str] = None  # environment variable override name
    description: Optional[str] = None

# Ordered (in roughly common usage frequency)
_DETECTOR_META: Dict[str, DetectorMeta] = {
    'hydrogenbond': DetectorMeta(
        key='hydrogenbond', label='Hydrogen Bonds', threshold_key='hbond',
        env_threshold_var='MOLBRIDGE_KDTREE_HBOND_THRESHOLD',
        description='Conventional hydrogen bonds including sulfur-mediated variants.'
    ),
    'ionic_interaction': DetectorMeta(
        key='ionic_interaction', label='Ionic Interactions', threshold_key='ionic',
        env_threshold_var='MOLBRIDGE_KDTREE_IONIC_THRESHOLD',
        description='Electrostatic contacts between oppositely charged side chains.'
    ),
    'hydrophobic_contact': DetectorMeta(
        key='hydrophobic_contact', label='Hydrophobic Contacts', threshold_key='hydro',
        env_threshold_var='MOLBRIDGE_KDTREE_HYDRO_THRESHOLD',
        description='Non‑polar side chain proximity interactions.'
    ),
    'pipi': DetectorMeta(
        key='pipi', label='π–π Stacking', threshold_key='pipi',
        description='Aromatic ring stacking interactions (parallel, T-shaped, offset).'
    ),
    'chalcogen_bond': DetectorMeta(
        key='chalcogen_bond', label='Chalcogen Bonds', threshold_key='chalcogen',
        description='σ‑hole interactions involving sulfur (and future Se/Te).' 
    ),
    'tetrel_bond': DetectorMeta(
        key='tetrel_bond', label='Tetrel Bonds', threshold_key='tetrel',
        description='Carbon σ‑hole contacts with electron-rich acceptors.'
    ),
    'halogenbond': DetectorMeta(
        key='halogenbond', label='Halogen Bonds', threshold_key='halogen',
        description='X···A σ‑hole interactions (Cl/Br/I donors).' 
    ),
    'pnictogen_bond': DetectorMeta(
        key='pnictogen_bond', label='Pnictogen Bonds', threshold_key='pnictogen',
        description='Group 15 (N/P/As) directional noncovalent interactions.'
    ),
    # Additional detectors (placeholders without thresholds for now)
    'anion_pi': DetectorMeta(
        key='anion_pi', label='Anion–π Interactions',
        description='Anionic center interacting with an aromatic π system.'
    ),
    'cation_pi': DetectorMeta(
        key='cation_pi', label='Cation–π Interactions',
        description='Cationic center engaging an aromatic π system.'
    ),
    'sulfur_pi': DetectorMeta(
        key='sulfur_pi', label='Sulfur–π Interactions',
        description='Sulfur lone pair / aromatic ring interactions.'
    ),
    'ch_pi': DetectorMeta(
        key='ch_pi', label='C–H···π Interactions',
        description='Weak hydrogen bond between aliphatic/aromatic C–H and π system.'
    ),
    'n_pi_star': DetectorMeta(
        key='n_pi_star', label='n→π* Interactions',
        description='Lone pair to antibonding π* orbital delocalization geometries.'
    ),
    'london_dispersion': DetectorMeta(
        key='london_dispersion', label='London Dispersion',
        description='Van der Waals dispersion contacts (broad, weak).' 
    ),
    'salt_bridge': DetectorMeta(
        key='salt_bridge', label='Salt Bridges',
        description='Paired cationic/anionic side-chain ion pair contacts.'
    ),
    'metal_coordination': DetectorMeta(
        key='metal_coordination', label='Metal Coordination',
        description='Protein side chain / ligand coordination to metal ion.'
    ),
}


def get_meta(key: str) -> Optional[DetectorMeta]:
    """Return metadata for detector key (None if unknown)."""
    return _DETECTOR_META.get(key)


def list_detector_keys() -> Iterable[str]:
    return _DETECTOR_META.keys()


def all_metadata() -> Dict[str, DetectorMeta]:
    return dict(_DETECTOR_META)

__all__ = [
    'DetectorMeta', 'get_meta', 'list_detector_keys', 'all_metadata'
]
