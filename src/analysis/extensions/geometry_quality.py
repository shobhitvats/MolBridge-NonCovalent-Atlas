"""Geometry quality metrics: Ramachandran classification, simple clash detection, torsion outliers.

Simplifications:
- Ramachandran region classification uses broad allowed ranges:
  * Favored: phi in (-160,-40), psi in (-80, 130)
  * Allowed: phi in (-180, 180), psi in (-180, 180) else -> none (should not happen)
- Clash detection: any non-bonded heavy atom pair < 2.4A flagged (cap at first 200 to limit size)
- Torsion outliers: reuse secondary_structure torsions and flag |omega| > 30 (cis/strained) or phi/psi outside broad region.

This is a heuristic placeholder; integrate full libraries (e.g., MolProbity) later if desired.
"""
from typing import Dict, Any, List, Tuple
import math


def compute(result: Dict[str, Any], config) -> Dict[str, Any]:
    ss_ext = (result.get('extensions', {}) or {}).get('secondary_structure')
    torsions = ss_ext.get('residues') if ss_ext else {}

    rama = {'favored': 0, 'allowed': 0, 'outlier': 0}
    torsion_outliers = []
    for rid, vals in torsions.items():
        phi = vals.get('phi'); psi = vals.get('psi'); omega = vals.get('omega')
        if phi is None or psi is None:
            continue
        if -160 < phi < -40 and -80 < psi < 130:
            rama['favored'] += 1
        elif -180 <= phi <= 180 and -180 <= psi <= 180:
            rama['allowed'] += 1
        else:
            rama['outlier'] += 1
        if omega is not None and abs(omega) > 30:
            torsion_outliers.append({'residue': rid, 'phi': phi, 'psi': psi, 'omega': omega, 'type': 'omega_strain'})
        if not (-180 <= phi <= 180 and -180 <= psi <= 180):
            torsion_outliers.append({'residue': rid, 'phi': phi, 'psi': psi, 'omega': omega, 'type': 'phi_psi_out_of_range'})

    # clash detection
    structures = result.get('structures') or []
    coords = []
    residues = []
    for s in structures:
        for res in s.get('residues', []):
            res_id = res.get('id') or f"{res.get('name')}:{res.get('seq')}"
            for atom in res.get('atoms', []):
                if isinstance(atom, dict):
                    name = atom.get('name','')
                    if name and name[0] != 'H':  # skip hydrogens
                        c = atom.get('coord')
                        if c and len(c) == 3:
                            coords.append((res_id, name, c))
                            residues.append(res_id)
    clashes = []
    n = len(coords)
    for i in range(n):
        ri, ai, ci = coords[i]
        for j in range(i+1, n):
            rj, aj, cj = coords[j]
            dx = ci[0]-cj[0]; dy = ci[1]-cj[1]; dz = ci[2]-cj[2]
            d2 = dx*dx + dy*dy + dz*dz
            if d2 < 2.4*2.4:
                clashes.append({'res1': ri, 'atom1': ai, 'res2': rj, 'atom2': aj, 'distance': math.sqrt(d2)})
                if len(clashes) >= 200:
                    break
        if len(clashes) >= 200:
            break

    total_eval = sum(rama.values()) or 1
    rama_fractions = {k: v/total_eval for k,v in rama.items()}
    return {
        'ramachandran': rama,
        'ramachandran_fractions': rama_fractions,
        'torsion_outliers': torsion_outliers,
        'clashes': clashes,
        'clash_count': len(clashes),
        'note': 'Heuristic geometry metrics'
    }
