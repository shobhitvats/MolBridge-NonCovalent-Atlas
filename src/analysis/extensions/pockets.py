"""Pocket / cavity heuristic detection.

Extremely lightweight centroid grid approach:
- Collect all heavy atom coordinates.
- Compute bounding box & coarse 3A grid.
- For each grid point inside box, evaluate distance to nearest atom.
- Points farther than threshold (e.g., >3.5A) from any atom but within inflated bounding box treated as cavity points.
- Cluster cavity points by simple BFS using grid adjacency.
- Report top clusters by size with approximate volume (#points * 27 Å^3 for 3A spacing).

This is a placeholder; replace with proper pocket detection (fpocket, castp, etc.) later.
"""
from typing import Dict, Any, List, Tuple
import math

GRID = 3.0
NEAR_THRESH2 = 3.5 * 3.5
MAX_POINTS = 40000


def _cluster(points):
    # points: list of (x,y,z)
    clusters = []
    visited = set()
    idx_map = {p: i for i,p in enumerate(points)}
    # Map grid coords for adjacency (round to nearest GRID integer index)
    grid_index = {}
    for p in points:
        gx = round(p[0]/GRID); gy = round(p[1]/GRID); gz = round(p[2]/GRID)
        grid_index.setdefault((gx,gy,gz), []).append(p)

    def neighbors(p):
        gx = round(p[0]/GRID); gy = round(p[1]/GRID); gz = round(p[2]/GRID)
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for dz in (-1,0,1):
                    for q in grid_index.get((gx+dx,gy+dy,gz+dz), []):
                        yield q

    for p in points:
        if p in visited:
            continue
        stack = [p]
        cluster = []
        visited.add(p)
        while stack:
            cur = stack.pop()
            cluster.append(cur)
            for q in neighbors(cur):
                if q not in visited:
                    visited.add(q)
                    stack.append(q)
        clusters.append(cluster)
    return clusters


def compute(result: Dict[str, Any], config) -> Dict[str, Any]:
    structures = result.get('structures') or []
    coords = []
    for s in structures:
        for res in s.get('residues', []):
            for atom in res.get('atoms', []):
                if isinstance(atom, dict):
                    name = atom.get('name','')
                    if name and name[0] != 'H':
                        c = atom.get('coord')
                        if c and len(c) == 3:
                            coords.append(c)
    if not coords:
        return {'pockets': [], 'note': 'No coordinates'}
    if len(coords) > 20000:
        return {'pockets': [], 'note': 'Structure too large for heuristic'}

    xs = [c[0] for c in coords]; ys = [c[1]]; zs = [c[2] for c in coords]
    minx = min(x[0] for x in coords); maxx = max(x[0] for x in coords)
    miny = min(x[1] for x in coords); maxy = max(x[1] for x in coords)
    minz = min(x[2] for x in coords); maxz = max(x[2] for x in coords)
    # Inflate bounding box by 2A
    minx -= 2; miny -= 2; minz -= 2
    maxx += 2; maxy += 2; maxz += 2

    cavity_pts = []
    # Precompute for speed
    for gx in range(int(minx//GRID), int(maxx//GRID)+1):
        x = gx * GRID
        for gy in range(int(miny//GRID), int(maxy//GRID)+1):
            y = gy * GRID
            for gz in range(int(minz//GRID), int(maxz//GRID)+1):
                z = gz * GRID
                # find nearest atom squared distance
                mind2 = 1e9
                for c in coords:
                    dx = x-c[0]; dy=y-c[1]; dz=z-c[2]
                    d2 = dx*dx + dy*dy + dz*dz
                    if d2 < mind2:
                        mind2 = d2
                        if mind2 < NEAR_THRESH2:
                            break
                if mind2 > NEAR_THRESH2:
                    cavity_pts.append((x,y,z))
                    if len(cavity_pts) > MAX_POINTS:
                        break
            if len(cavity_pts) > MAX_POINTS:
                break
        if len(cavity_pts) > MAX_POINTS:
            break

    if not cavity_pts:
        return {'pockets': [], 'note': 'No cavity points detected'}

    clusters = _cluster(cavity_pts)
    pockets = []
    for cl in clusters:
        if len(cl) < 4:
            continue
        xs = [p[0] for p in cl]; ys = [p[1] for p in cl]; zs = [p[2] for p in cl]
        centroid = (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))
        volume = len(cl) * (GRID**3)
        pockets.append({'size_points': len(cl), 'centroid': centroid, 'approx_volume': volume})
    pockets.sort(key=lambda x: x['size_points'], reverse=True)
    return {'pockets': pockets[:10], 'note': 'Heuristic grid-based pockets'}
