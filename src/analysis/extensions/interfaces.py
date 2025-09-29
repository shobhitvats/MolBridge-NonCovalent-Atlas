"""Interface analysis between protein chains.

Estimates chain-chain interfaces using interaction connectivity and optional
residue co-occurrence. This is a lightweight heuristic (not full SASA/BSA).
"""
from typing import Dict, Any, Tuple, Set
from collections import defaultdict


def compute(result: Dict[str, Any], config) -> Dict[str, Any]:
    interactions = result.get('interactions', {})
    chain_pair_counts: Dict[Tuple[str, str], int] = defaultdict(int)
    chain_residues: Dict[str, Set[str]] = defaultdict(set)

    for int_type, items in interactions.items():
        for it in items:
            r1 = _prop(it, 'residue1') or _prop(it, 'donor_residue') or _prop(it, 'ring1_residue')
            r2 = _prop(it, 'residue2') or _prop(it, 'acceptor_residue') or _prop(it, 'ring2_residue')
            c1 = _prop(it, 'chain1') or _prop(it, 'donor_chain') or _prop(it, 'ring1_chain')
            c2 = _prop(it, 'chain2') or _prop(it, 'acceptor_chain') or _prop(it, 'ring2_chain')
            if not (r1 and r2 and c1 and c2):
                continue
            chain_residues[c1].add(r1)
            chain_residues[c2].add(r2)
            if c1 != c2:
                key = tuple(sorted((c1, c2)))
                chain_pair_counts[key] += 1

    interface_list = []
    for (c1, c2), count in chain_pair_counts.items():
        interface_list.append({
            'chain_pair': f"{c1}-{c2}",
            'interaction_count': count,
            'approx_interface_residues': count  # placeholder heuristic
        })

    interface_list.sort(key=lambda x: x['interaction_count'], reverse=True)

    return {
        'interfaces_found': len(interface_list),
        'interfaces': interface_list
    }


def _prop(obj: Any, name: str):
    if hasattr(obj, name):
        return getattr(obj, name, None)
    if isinstance(obj, dict):
        return obj.get(name)
    return None
