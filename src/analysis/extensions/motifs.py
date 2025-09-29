"""Experimental structural motifs detection.
Currently implements a simple heuristic for 'aromatic cage' motifs: 
Three or more aromatic residues (PHE, TYR, TRP, HIS) with ring centroids within 6 Å of a central cationic residue (ARG or LYS) side-chain nitrogen centroid.
This is a coarse heuristic; future work can refine geometric criteria.
"""
from typing import Dict, Any, List, Tuple
import math

AROMATICS = {"PHE", "TYR", "TRP", "HIS"}
CATIONICS = {"ARG", "LYS"}
MAX_RING_DISTANCE = 6.0
MIN_AROMATICS = 3

# Utility centroid placeholder; real implementation would use actual atom coordinates.

def _centroid(coord_list: List[Tuple[float, float, float]]):
    if not coord_list:
        return (0.0, 0.0, 0.0)
    x = sum(c[0] for c in coord_list) / len(coord_list)
    y = sum(c[1] for c in coord_list) / len(coord_list)
    z = sum(c[2] for c in coord_list) / len(coord_list)
    return (x, y, z)

def _distance(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def compute(result: Dict[str, Any], config) -> Dict[str, Any]:
    structures = result.get('structures', [])
    motifs = []
    # Expect that each structure contains residues with minimal coordinate info; if absent, abort gracefully.
    for s in structures:
        residues = s.get('residues', [])
        aromatic_res = [r for r in residues if r.get('name') in AROMATICS]
        cationic_res = [r for r in residues if r.get('name') in CATIONICS]
        if not aromatic_res or not cationic_res:
            continue
        # Placeholder coordinate retrieval; using 'centroid' key if present
        for cat in cationic_res:
            cat_centroid = cat.get('centroid')
            if not cat_centroid:
                continue
            close_aromatics = []
            for ar in aromatic_res:
                ar_centroid = ar.get('centroid')
                if not ar_centroid:
                    continue
                if _distance(cat_centroid, ar_centroid) <= MAX_RING_DISTANCE:
                    close_aromatics.append(ar)
            if len(close_aromatics) >= MIN_AROMATICS:
                motifs.append({
                    'type': 'aromatic_cage',
                    'structure_id': s.get('id'),
                    'central_residue': cat.get('id'),
                    'aromatic_residues': [r.get('id') for r in close_aromatics],
                    'count': len(close_aromatics)
                })
    return {
        'motif_count': len(motifs),
        'motifs': motifs,
        'parameters': {
            'max_ring_distance': MAX_RING_DISTANCE,
            'min_aromatics': MIN_AROMATICS
        }
    }
