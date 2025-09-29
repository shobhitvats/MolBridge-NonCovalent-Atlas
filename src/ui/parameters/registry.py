"""Central interaction parameter registry.

Separated from the main interface for reuse and faster import.
Each entry maps an interaction type key (matching internal identifiers)
 to a list of slider/parameter metadata dicts consumed by the UI.
"""
from __future__ import annotations
from typing import Dict, List, Any

# Public: REGISTRY constant
REGISTRY: Dict[str, List[Dict[str, Any]]] = {
    'hydrogenbond': [
        {'field': 'hbond_distance_cutoff', 'label': 'Distance (Å)', 'min': 2.5, 'max': 5.0, 'step': 0.05, 'help': 'Max donor-acceptor heavy atom distance'},
        {'field': 'hbond_angle_cutoff', 'label': 'Angle (°)', 'min': 90.0, 'max': 180.0, 'step': 1.0, 'help': 'Minimum D–H···A angle'}
    ],
    'halogenbond': [
        {'field': 'halogen_distance_cutoff', 'label': 'Distance (Å)', 'min': 3.0, 'max': 6.0, 'step': 0.1, 'help': 'Halogen to acceptor distance'},
        {'field': 'halogen_angle_cutoff', 'label': 'Angle (°)', 'min': 100.0, 'max': 180.0, 'step': 1.0, 'help': 'C–X···A angle threshold'}
    ],
    'chalcogenbond': [
        {'field': 'chalcogen_distance_cutoff', 'label': 'Distance (Å)', 'min': 3.0, 'max': 6.0, 'step': 0.1, 'help': 'Chalcogen to acceptor distance'},
        {'fields': ('chalcogen_theta_min', 'chalcogen_theta_max'), 'label': 'Theta Range (°)', 'min': 90.0, 'max': 180.0, 'step': 1.0, 'help': 'Allowed θ window (centroid→S vs S→acceptor). Drag ends to set min & max.'},
        {'field': 'chalcogen_phi_max', 'label': 'Phi Max (°)', 'min': 10.0, 'max': 90.0, 'step': 1.0, 'help': 'Maximum |φ| (out-of-plane / delta deviation)'}
    ],
    'pnictogenbond': [
        {'field': 'pnictogen_distance_cutoff', 'label': 'Distance (Å)', 'min': 3.0, 'max': 6.0, 'step': 0.1, 'help': 'Pnictogen to acceptor distance'},
        {'field': 'pnictogen_angle_cutoff', 'label': 'Angle (°)', 'min': 100.0, 'max': 180.0, 'step': 1.0, 'help': 'σ-hole alignment angle'}
    ],
    'tetrelbond': [
        {'field': 'tetrel_distance_cutoff', 'label': 'Distance (Å)', 'min': 3.0, 'max': 6.0, 'step': 0.1, 'help': 'Tetrel to acceptor distance'},
        {'field': 'tetrel_angle_cutoff', 'label': 'Angle (°)', 'min': 100.0, 'max': 180.0, 'step': 1.0, 'help': 'σ-hole alignment angle'}
    ],
    'chpi': [
        {'field': 'ch_pi_min_distance', 'label': 'Min Dist (Å)', 'min': 1.5, 'max': 3.0, 'step': 0.05, 'help': 'Lower bound C···centroid'},
        {'field': 'ch_pi_max_distance', 'label': 'Max Dist (Å)', 'min': 3.5, 'max': 6.0, 'step': 0.05, 'help': 'Upper bound C···centroid'},
        {'field': 'ch_pi_angle_cutoff', 'label': 'Angle (°)', 'min': 60.0, 'max': 120.0, 'step': 1.0, 'help': 'Max C-H···centroid angle'},
        {'field': 'ch_pi_max_height', 'label': 'Max Height (Å)', 'min': 1.0, 'max': 4.0, 'step': 0.05, 'help': 'Perpendicular distance limit'}
    ],
    'pipi': [
        {'field': 'pi_pi_distance_cutoff', 'label': 'Distance (Å)', 'min': 3.5, 'max': 8.0, 'step': 0.1, 'help': 'Centroid–centroid distance'},
        {'field': 'pi_pi_angle_cutoff', 'label': 'Angle (°)', 'min': 0.0, 'max': 60.0, 'step': 1.0, 'help': 'Planar angle tolerance'}
    ],
    'anionpi': [
        {'field': 'anion_pi_distance_cutoff', 'label': 'Distance (Å)', 'min': 3.0, 'max': 7.0, 'step': 0.1, 'help': 'Anion–centroid distance'}
    ],
    'npistar': [
        {'field': 'n_pi_star_distance_cutoff', 'label': 'Distance (Å)', 'min': 3.0, 'max': 5.0, 'step': 0.1, 'help': 'Carbonyl O to carbonyl C distance'},
        {'field': 'n_pi_star_angle_cutoff', 'label': 'Angle (°)', 'min': 90.0, 'max': 150.0, 'step': 1.0, 'help': 'n→π* approach angle'}
    ],
    'dispersion': [
        {'field': 'dispersion_min_distance', 'label': 'Min Dist (Å)', 'min': 3.0, 'max': 4.0, 'step': 0.05, 'help': 'Lower van der Waals window'},
        {'field': 'dispersion_distance_cutoff', 'label': 'Max Dist (Å)', 'min': 4.0, 'max': 6.0, 'step': 0.05, 'help': 'Upper van der Waals window'}
    ],
    'ionicinteraction': [
        {'field': 'ionic_distance_cutoff', 'label': 'Distance (Å)', 'min': 4.0, 'max': 10.0, 'step': 0.1, 'help': 'Charged group centroid distance'}
    ],
    'hydrophobiccontact': [
        {'field': 'hydrophobic_distance_cutoff', 'label': 'Distance (Å)', 'min': 3.5, 'max': 7.0, 'step': 0.1, 'help': 'Carbon–carbon cutoff'}
    ],
    'cation_pi': [
        {'field': 'cation_pi_distance_cutoff', 'label': 'Distance (Å)', 'min': 4.0, 'max': 8.0, 'step': 0.1, 'help': 'Cation–centroid distance'}
    ],
    'salt_bridge': [
        {'field': 'salt_bridge_distance_cutoff', 'label': 'Atom Dist (Å)', 'min': 3.0, 'max': 6.0, 'step': 0.1, 'help': 'Atom-atom distance'},
        {'field': 'salt_bridge_centroid_cutoff', 'label': 'Centroid Dist (Å)', 'min': 4.0, 'max': 8.0, 'step': 0.1, 'help': 'Charged group centroid distance'}
    ],
    'sulfur_pi': [
        {'field': 'sulfur_pi_distance_cutoff', 'label': 'Distance (Å)', 'min': 4.0, 'max': 8.0, 'step': 0.1, 'help': 'Sulfur–centroid distance'},
        {'field': 'sulfur_pi_max_perp', 'label': 'Max Perp (Å)', 'min': 1.0, 'max': 5.0, 'step': 0.05, 'help': 'Max perpendicular offset'}
    ],
    'metal_coordination': [
        {'field': 'metal_coordination_primary_cutoff', 'label': 'Primary (Å)', 'min': 1.5, 'max': 3.5, 'step': 0.05, 'help': 'Primary coordination shell'},
        {'field': 'metal_coordination_extended_cutoff', 'label': 'Extended (Å)', 'min': 2.0, 'max': 4.0, 'step': 0.05, 'help': 'Extended shell'}
    ],
}

__all__ = ["REGISTRY", "get_parameters_for"]

def get_parameters_for(interaction_type: str):
    """Return parameter metadata list for an interaction type (empty list if none)."""
    return REGISTRY.get(interaction_type, [])
