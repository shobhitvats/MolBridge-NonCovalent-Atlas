"""π-π stacking refinement module.

Classifies existing π-π stacking interactions into subtypes based on ring plane angle if available:
- parallel: angle < 15°
- t_shaped: 60° <= angle <= 120°
- intermediate: otherwise
Falls back gracefully if angle not present.
Stores counts and per-interaction annotations without mutating original interaction list.
"""
from typing import Dict, Any, List


def compute(result: Dict[str, Any], config) -> Dict[str, Any]:
    interactions = result.get('interactions', {}) or {}
    # Collect possible pi-pi keys (support both UI name 'pipi' and detector name 'pi_pi')
    candidate_keys = ('pipi','pi_pi','pi-pi')
    pi_keys = [k for k in interactions.keys() if k in candidate_keys]
    if not pi_keys:
        return {
            'count': 0,
            'subtype_counts': {'parallel': 0, 't_shaped': 0, 'intermediate': 0},
            'interactions': [],
            'note': 'No pi-pi interactions found in interactions dict'
        }

    records: List[Dict[str, Any]] = []
    subtype_counts = {'parallel': 0, 't_shaped': 0, 'intermediate': 0}
    total = 0
    for key in pi_keys:
        for rec in interactions.get(key, []):
            total += 1
            angle = None
            if isinstance(rec, dict):
                angle = rec.get('angle') or rec.get('dihedral') or rec.get('tilt_angle')
            else:
                angle = getattr(rec, 'angle', None) or getattr(rec, 'dihedral', None) or getattr(rec, 'tilt_angle', None)
            subtype = 'intermediate'
            if angle is not None:
                try:
                    if angle < 15:
                        subtype = 'parallel'
                    elif 60 <= angle <= 120:
                        subtype = 't_shaped'
                except Exception:
                    pass
            subtype_counts[subtype] += 1
            records.append({'angle': angle, 'subtype': subtype})
    # Keep zero categories so UI can show explicit zeros instead of blank
    return {
        'count': total,
        'subtype_counts': subtype_counts,
        'interactions': records,
        'note': 'Heuristic pi-pi orientation classification'
    }
