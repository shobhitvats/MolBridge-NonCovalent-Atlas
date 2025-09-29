"""Disulfide bond detection module.

Detects S-S bonds between cysteine SG atoms within 2.2Å - 2.3Å (flexible tolerance up to 2.5Å).
Outputs list of disulfide pairs and simple statistics.
"""
from typing import Dict, Any
import math


def compute(result: Dict[str, Any], config) -> Dict[str, Any]:
    structures = result.get('structures') or []
    sg_atoms = []  # (res_id, coord)
    for s in structures:
        for res in s.get('residues', []):
            name = res.get('name', '').upper()
            if name not in ('CYS','CYX','CSS'):
                continue
            res_id = res.get('id') or f"{name}:{res.get('seq')}"
            for atom in res.get('atoms', []):
                if isinstance(atom, dict) and atom.get('name','').upper() in ('SG','S'):
                    c = atom.get('coord')
                    if c and len(c) == 3:
                        sg_atoms.append((res_id, c))
    pairs = []
    n = len(sg_atoms)
    for i in range(n):
        ri, ci = sg_atoms[i]
        for j in range(i+1, n):
            rj, cj = sg_atoms[j]
            dx = ci[0]-cj[0]; dy = ci[1]-cj[1]; dz = ci[2]-cj[2]
            d2 = dx*dx + dy*dy + dz*dz
            if d2 < 2.5*2.5:
                dist = math.sqrt(d2)
                if 1.8 <= dist <= 2.5:
                    pairs.append({'cys1': ri, 'cys2': rj, 'distance': dist})
    return {
        'disulfides': pairs,
        'count': len(pairs),
        'note': 'CYS SG-SG within 2.5A'
    }
