"""Outlier / threshold-adjacent interaction detection.

Flags interactions whose geometric parameters are near configured cutoffs.
"""
from typing import Dict, Any, List

# Default relative windows
DIST_WINDOW = 0.2  # Å within cutoff
ANGLE_WINDOW = 5.0  # degrees near cutoff


def compute(result: Dict[str, Any], config) -> Dict[str, Any]:
    interactions = result.get('interactions', {})
    outlier_records: List[Dict[str, Any]] = []

    # Mapping from interaction type key fragments to config attribute patterns
    distance_map = {
        'hydrogen_bond': ('hbond_distance_cutoff', 'hbond_angle_cutoff'),
        'halogen_bond': ('halogen_distance_cutoff', 'halogen_angle_cutoff'),
        'chalcogen_bond': ('chalcogen_distance_cutoff', 'chalcogen_angle_cutoff'),
        'pnictogen_bond': ('pnictogen_distance_cutoff', 'pnictogen_angle_cutoff'),
        'tetrel_bond': ('tetrel_distance_cutoff', 'tetrel_angle_cutoff'),
        'ch_pi': ('ch_pi_max_distance', 'ch_pi_angle_cutoff'),
        'pi_pi': ('pi_pi_distance_cutoff', 'pi_pi_angle_cutoff'),
        'anion_pi': ('anion_pi_distance_cutoff', None),
        'n_pi_star': ('n_pi_star_distance_cutoff', 'n_pi_star_angle_cutoff'),
        'ionic': ('ionic_distance_cutoff', None),
        'hydrophobic': ('hydrophobic_distance_cutoff', None)
    }

    cfg = config.interactions

    for int_type, items in interactions.items():
        dist_attr, angle_attr = distance_map.get(int_type, (None, None))
        cutoff_d = getattr(cfg, dist_attr, None) if dist_attr else None
        cutoff_a = getattr(cfg, angle_attr, None) if angle_attr else None
        for it in items:
            dist = _prop(it, 'distance')
            ang = _prop(it, 'angle')
            near = False
            notes = []
            if cutoff_d and isinstance(dist, (int, float)) and dist > (cutoff_d - DIST_WINDOW):
                near = True
                notes.append(f"distance {dist:.2f} close to cutoff {cutoff_d}")
            if cutoff_a and isinstance(ang, (int, float)) and ang < (cutoff_a + ANGLE_WINDOW):
                near = True
                notes.append(f"angle {ang:.1f} near minimum {cutoff_a}")
            if near:
                outlier_records.append({
                    'interaction_type': int_type,
                    'residue1': _prop(it, 'residue1') or _prop(it, 'donor_residue') or _prop(it, 'ring1_residue'),
                    'residue2': _prop(it, 'residue2') or _prop(it, 'acceptor_residue') or _prop(it, 'ring2_residue'),
                    'distance': dist,
                    'angle': ang,
                    'notes': '; '.join(notes)
                })

    return {
        'count': len(outlier_records),
        'records': outlier_records
    }


def _prop(obj: Any, name: str):
    if hasattr(obj, name):
        return getattr(obj, name, None)
    if isinstance(obj, dict):
        return obj.get(name)
    return None
