"""SASA & BSA heuristic module.

Approach:
- Attempt to import freesasa if available for accurate SASA.
- Fallback: very rough approximation using count of exposed heavy atoms (atoms with < X neighbors within 5A) * constant.
- BSA: If multiple chains, estimate interface buried surface by: BSA = sum(SASA(chain)) - SASA(all)
  (Only if freesasa present; fallback leaves BSA None.)
"""
from typing import Dict, Any


def _try_freesasa(structures):
    try:
        import freesasa  # type: ignore
    except Exception:
        return None
    # Build PDB string from structures minimal
    pdb_lines = []
    atom_serial = 1
    for s in structures:
        chain = s.get('chain_id', 'A')
        for res in s.get('residues', []):
            resn = res.get('name', 'UNK')[:3]
            resid = res.get('seq') or res.get('id') or 0
            for atom in res.get('atoms', []):
                if not isinstance(atom, dict):
                    continue
                name = atom.get('name','?')
                if len(name) < 4:
                    name_fmt = f" {name:<3}"  # right align in 4 cols
                else:
                    name_fmt = name[:4]
                x,y,z = atom.get('coord',(0.0,0.0,0.0))
                pdb_lines.append(f"ATOM  {atom_serial:5d} {name_fmt} {resn:>3} {chain:1}{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C")
                atom_serial += 1
    pdb_str = "\n".join(pdb_lines)
    try:
        structure = freesasa.Structure(pdb_str)
        result = freesasa.calc(structure)
    except Exception:
        return None
    # Per-chain and total
    chain_sasa = {}
    for ch in structure.chainLabels():
        sel = freesasa.selectArea({'sel': f"chain {ch}"}, structure, result)
        chain_sasa[ch] = sel['sel']
    total = result.totalArea()

    bsa = None
    if len(chain_sasa) > 1:
        # naive interface buried surface = sum(chain SASA) - complex SASA
        sum_ch = sum(chain_sasa.values())
        bsa = sum_ch - total
    return {'total_sasa': total, 'chain_sasa': chain_sasa, 'bsa_estimate': bsa, 'method': 'freesasa'}


def _rough_fallback(structures):
    # Count atoms with few neighbors as 'exposed'. Very crude.
    import math
    coords = []
    for s in structures:
        for res in s.get('residues', []):
            for atom in res.get('atoms', []):
                if isinstance(atom, dict):
                    coords.append(atom.get('coord'))
    coords = [c for c in coords if c and len(c) == 3]
    exposed = 0
    for i, c in enumerate(coords):
        neighbors = 0
        for j, d in enumerate(coords):
            if i == j:
                continue
            dx = c[0]-d[0]; dy = c[1]-d[1]; dz = c[2]-d[2]
            if dx*dx + dy*dy + dz*dz < 25:  # 5A^2
                neighbors += 1
                if neighbors > 25:  # crowded
                    break
        if neighbors < 15:
            exposed += 1
    total_sasa = exposed * 25.0  # arbitrary scale factor
    return {'total_sasa': total_sasa, 'chain_sasa': {}, 'bsa_estimate': None, 'method': 'heuristic_exposed_atoms'}


def compute(result: Dict[str, Any], config) -> Dict[str, Any]:
    structures = result.get('structures') or []
    if not structures:
        return {'total_sasa': None, 'chain_sasa': {}, 'bsa_estimate': None, 'method': 'no_structure'}
    data = _try_freesasa(structures)
    if data:
        return data
    return _rough_fallback(structures)
