"""Secondary structure & torsion heuristic module.

Attempts to compute phi/psi/omega from backbone atom coordinates if present in result['structures']
(assuming a normalized structure dict). If unavailable, returns a placeholder.
A very lightweight H/E/C assignment:
- Helix (H): phi in (-100,-30) and psi in (-80, -5)
- Sheet (E): phi in (-160,-80) and psi in (90,180)
- Coil (C): else
"""
from typing import Dict, Any, List, Tuple
import math

# Helper to compute torsion angle between four 3D points

def _torsion(p1, p2, p3, p4):
    import numpy as np
    b1 = np.array(p2) - np.array(p1)
    b2 = np.array(p3) - np.array(p2)
    b3 = np.array(p4) - np.array(p3)
    # normalize b2 for numerical stability
    b2n = b2 / (np.linalg.norm(b2) + 1e-9)
    v = b1 - (b1 @ b2n) * b2n
    w = b3 - (b3 @ b2n) * b2n
    x = v @ w
    y = (np.cross(b2n, v) @ w)
    ang = math.degrees(math.atan2(y, x))
    return ang


def _assign_ss(phi: float, psi: float) -> str:
    if phi is None or psi is None:
        return 'C'
    if -100 < phi < -30 and -80 < psi < -5:
        return 'H'
    if -160 < phi < -80 and 90 < psi < 180:
        return 'E'
    return 'C'


def compute(result: Dict[str, Any], config) -> Dict[str, Any]:
    structures = result.get('structures') or []
    if not structures:
        return {'residues': {}, 'note': 'No structure coordinate data present'}

    # Expect each structure has residues list with atoms [{'name': 'N','coord':(x,y,z)}, ...]
    residues_out = {}
    for s in structures:
        residues = s.get('residues', [])
        # Build index by (res_id)
        for i, res in enumerate(residues):
            atoms = {a['name']: a.get('coord') for a in res.get('atoms', []) if isinstance(a, dict)}
            # Need previous C, current N, CA, C, next N for torsions
            if i == 0 or i == len(residues) - 1:
                phi = psi = omega = None
            else:
                prev_atoms = {a['name']: a.get('coord') for a in residues[i-1].get('atoms', []) if isinstance(a, dict)}
                next_atoms = {a['name']: a.get('coord') for a in residues[i+1].get('atoms', []) if isinstance(a, dict)}
                try:
                    if all(k in prev_atoms for k in ['C']) and all(k in atoms for k in ['N','CA','C']) and 'N' in next_atoms:
                        phi = _torsion(prev_atoms['C'], atoms['N'], atoms['CA'], atoms['C'])
                        psi = _torsion(atoms['N'], atoms['CA'], atoms['C'], next_atoms['N']) if next_atoms.get('N') else None
                        # omega: preceding peptide bond torsion (C- N - CA - C)
                        omega = phi - psi if (phi is not None and psi is not None) else None
                    else:
                        phi = psi = omega = None
                except Exception:
                    phi = psi = omega = None
            ss = _assign_ss(phi, psi)
            residues_out[res.get('id') or f"{res.get('name')}:{res.get('seq', i)}"] = {
                'phi': phi,
                'psi': psi,
                'omega': omega,
                'ss': ss
            }
    # Summaries
    counts = {'H':0,'E':0,'C':0}
    for r in residues_out.values():
        counts[r['ss']] = counts.get(r['ss'],0)+1
    total = sum(counts.values()) or 1
    fractions = {k: v/total for k,v in counts.items()}
    # Provide both simple keys (H,E,C) and suffixed keys (H_frac, etc.) expected by UI/export
    fractions_augmented = {**fractions, **{f"{k}_frac": fractions[k] for k in counts}}
    return {
        'residues': residues_out,
        'counts': counts,
        'fractions': fractions_augmented,
        'note': 'Heuristic secondary structure (no DSSP)'
    }
