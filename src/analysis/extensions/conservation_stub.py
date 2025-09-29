"""Conservation placeholder.

Without external MSA computation, we approximate conservation by:
- If per-residue interaction frequency (from residue_profiles extension) exists across batch runs, reuse those counts if provided in result['batch_summary'].
- Otherwise, mark all residues with conservation_score = None.
Future: integrate with external alignment service or provided MSA.
"""
from typing import Dict, Any


def compute(result: Dict[str, Any], config) -> Dict[str, Any]:
    batch = result.get('batch_summary') or {}
    residue_profiles = (result.get('extensions', {}) or {}).get('residue_profiles') or {}
    residues = residue_profiles.get('profiles') or {}
    conservation = {}
    # If batch has aggregated counts, use them to rank conservation
    agg = batch.get('residue_usage') or {}
    if agg:
        max_c = max(agg.values()) or 1
        for rid in residues.keys():
            c = agg.get(rid, 0)
            conservation[rid] = {'raw_count': c, 'normalized': c/max_c}
        note = 'Approx conservation from batch residue_usage'
    else:
        for rid in residues.keys():
            conservation[rid] = {'raw_count': None, 'normalized': None}
        note = 'No batch aggregate; placeholder only'
    return {'residues': conservation, 'note': note}
