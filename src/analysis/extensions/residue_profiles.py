"""Residue profile extraction for protein interaction analysis.

Aggregates per-residue interaction counts and strength distributions from the
standard `result['interactions']` structure.
"""
from typing import Dict, Any, List
from collections import defaultdict


def compute(result: Dict[str, Any], config) -> Dict[str, Any]:
    interactions = result.get('interactions', {})
    residue_data: Dict[str, Dict[str, Any]] = {}

    for int_type, items in interactions.items():
        for it in items:
            r1 = _prop(it, 'residue1') or _prop(it, 'donor_residue') or _prop(it, 'ring1_residue')
            r2 = _prop(it, 'residue2') or _prop(it, 'acceptor_residue') or _prop(it, 'ring2_residue')
            s = _prop(it, 'strength') or 'unknown'
            for r in (r1, r2):
                if not r:
                    continue
                rec = residue_data.setdefault(r, {'total': 0, 'by_type': defaultdict(int), 'strengths': defaultdict(int)})
                rec['total'] += 1
                rec['by_type'][int_type] += 1
                rec['strengths'][s] += 1

    # Convert defaultdicts to normal dicts
    for r, rec in residue_data.items():
        rec['by_type'] = dict(rec['by_type'])
        rec['strengths'] = dict(rec['strengths'])

    return {
        'residue_count': len(residue_data),
        'profiles': residue_data
    }


def _prop(obj: Any, name: str):
    if hasattr(obj, name):
        return getattr(obj, name, None)
    if isinstance(obj, dict):
        return obj.get(name)
    return None
