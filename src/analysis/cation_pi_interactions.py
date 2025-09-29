"""Cation–π interaction detection.

Heuristic detection:
- Cation residues: LYS (NZ), ARG (guanidinium CZ center), HIS (if protonated form unknown: use ring centroid), occasionally protonated termini (not handled here), plus ligand atoms with +1 formal charge (not implemented yet).
- Aromatic rings: PHE, TYR, TRP, HIS (same ring logic as pi-pi detector); reuse simplified centroid computation.
- Distance criterion: centroid(cation group) to aromatic ring centroid <= 6.0 Å (configurable via interactions.cation_pi_distance_cutoff, fallback 6.0).
- Optional vertical geometry (angle between ring normal and vector to cation) is recorded but not filtered strictly; we supply it for future refinement.

Returns list of dict records consistent with other interaction detectors: keys include
residue1, residue2, chain1, chain2, distance, angle(optional), interaction_type, strength.

Strength heuristic:
  distance < 4.5 Å => strong
  4.5–5.5 Å => moderate
  else weak
"""
from typing import List, Dict, Any, Tuple
import numpy as np
from Bio.PDB import Structure, Residue, Atom
from utils.settings import get_settings
from utils.instrumentation import init_funnel, update_counts, finalize_funnel
from analysis.feature_store import get_feature_store
from .base_detector import register_detector
try:  # optional vector index
    from analysis.structure_index import build_aromatic_ring_index
except Exception:  # pragma: no cover
    build_aromatic_ring_index = None  # type: ignore

AROMATIC_RESIDUES = {"PHE","TYR","TRP","HIS"}
CATION_RESIDUES = {"LYS","ARG","HIS"}  # HIS ambiguous; still include heuristically

@register_detector("cation_pi", method="detect_cation_pi")
class CationPiDetector:
    def __init__(self, config):
        self.config = config
        self.distance_cutoff = getattr(self.config.interactions, 'cation_pi_distance_cutoff', 6.0)

    # ---------------- Utility geometry helpers -----------------
    def _residue_centroid(self, residue: Residue.Residue, atom_names: List[str]) -> np.ndarray:
        pts = []
        for name in atom_names:
            if name in residue:
                pts.append(residue[name].get_coord())
        if not pts:
            # fallback all heavy atoms
            pts = [atom.get_coord() for atom in residue if atom.element != 'H']
        if not pts:
            return np.array([0.0,0.0,0.0])
        return np.mean(np.vstack(pts), axis=0)

    def _aromatic_ring_centroid(self, residue: Residue.Residue) -> Tuple[np.ndarray, np.ndarray]:
        # approximate ring atoms similar to pi-pi module
        mapping = {
            'PHE': ['CG','CD1','CD2','CE1','CE2','CZ'],
            'TYR': ['CG','CD1','CD2','CE1','CE2','CZ'],
            'TRP': ['CG','CD1','NE1','CE2','CD2','CE3','CZ2','CZ3','CH2'],
            'HIS': ['CG','ND1','CD2','CE1','NE2']
        }
        atom_names = mapping.get(residue.get_resname(), [])
        coords = []
        for n in atom_names:
            if n in residue:
                coords.append(residue[n].get_coord())
        if len(coords) < 3:  # insufficient atoms
            return np.array([0.0,0.0,0.0]), np.array([0,0,1])
        coords_arr = np.vstack(coords)
        center = np.mean(coords_arr, axis=0)
        # Normal via SVD
        centered = coords_arr - center
        try:
            _, _, vh = np.linalg.svd(centered)
            normal = vh[-1]
            if normal[2] < 0:
                normal = -normal
        except Exception:
            normal = np.array([0,0,1])
        return center, normal

    def _cation_group_centroid(self, residue: Residue.Residue) -> np.ndarray:
        resname = residue.get_resname()
        if resname == 'LYS':
            return self._residue_centroid(residue, ['NZ','CE','CD'])
        if resname == 'ARG':
            return self._residue_centroid(residue, ['CZ','NH1','NH2'])
        if resname == 'HIS':  # if protonated; approximate using ring
            ctr, _ = self._aromatic_ring_centroid(residue)
            return ctr
        # fallback heavy centroid
        return self._residue_centroid(residue, [])

    def detect_cation_pi(self, structure: Structure.Structure) -> List[Any]:
        """Public entrypoint with feature-store backed vector path.

        Vector path triggers when MOLBRIDGE_ENABLE_VECTOR_GEOM=1 (settings.enable_vector_geom)
        and we can obtain at least one cation & one aromatic ring. Uses either the
        optional structure_index (if available) or falls back to FeatureStore cached
        ring geometry (centroid + normal arrays).
        """
        settings = get_settings()
        if getattr(settings, 'enable_vector_geom', False):
            try:
                return self._vector_detect(structure)
            except Exception:  # pragma: no cover - fail safe to legacy
                return self._legacy_detect(structure)
        return self._legacy_detect(structure)

    def _legacy_detect(self, structure: Structure.Structure) -> List[Any]:
        interactions: List[Dict[str, Any]] = []
        try:
            model = structure[0]
        except Exception:
            return interactions
        aromatics = []
        cations = []
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname in AROMATIC_RESIDUES:
                    c, n = self._aromatic_ring_centroid(residue)
                    aromatics.append((residue, c, n))
                if resname in CATION_RESIDUES:
                    cg = self._cation_group_centroid(residue)
                    cations.append((residue, cg))
        # Collect counts for instrumentation
        raw_pairs = len(cations) * len(aromatics)
        candidate_pairs = 0
        for res_c, c_cent in cations:
            for res_a, a_cent, a_norm in aromatics:
                if res_c is res_a:
                    continue
                dist = float(np.linalg.norm(c_cent - a_cent))
                if dist > self.distance_cutoff:
                    continue
                candidate_pairs += 1
                vec = c_cent - a_cent
                if np.linalg.norm(vec) > 1e-6:
                    vec_n = vec / np.linalg.norm(vec)
                    angle = np.degrees(np.arccos(np.clip(abs(np.dot(vec_n, a_norm)),0,1)))
                else:
                    angle = None
                if dist < 4.5:
                    strength = 'strong'
                elif dist < 5.5:
                    strength = 'moderate'
                else:
                    strength = 'weak'
                interactions.append({
                    'residue1': f"{res_c.get_resname()}{res_c.get_id()[1]}",
                    'residue2': f"{res_a.get_resname()}{res_a.get_id()[1]}",
                    'chain1': res_c.get_parent().get_id(),
                    'chain2': res_a.get_parent().get_id(),
                    'distance': round(dist,3),
                    'angle': round(angle,1) if angle is not None else None,
                    'interaction_type': 'cation_pi',
                    'strength': strength
                })
        # Legacy instrumentation (single phase)
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=len(interactions),
            core_pair_generation=False,
            extra={
                'cations': len(cations),
                'rings': len(aromatics)
            }
        )
        finalize_funnel(self.instrumentation, pair_gen_seconds=0.0, eval_seconds=0.0, build_seconds=0.0)
        return interactions

    def _vector_detect(self, structure: Structure.Structure) -> List[Any]:
        """Vector candidate pruning using either structure_index or FeatureStore caches.

        Fallback chain is:
          1. If build_aromatic_ring_index available -> use it (may include richer metadata)
          2. Else use FeatureStore.ensure_rings(structure) to get ring centroids/normals
          3. If no rings or cations -> return []
        """
        try:
            model = structure[0]
        except Exception:
            return []
        fs = get_feature_store()
        ring_centers: np.ndarray
        ring_normals: np.ndarray
        ring_records: List[Any]
        if build_aromatic_ring_index is not None:
            ring_index = build_aromatic_ring_index(model)
            if len(ring_index) == 0:
                return []
            ring_centers = ring_index.centers.astype('float32')
            ring_normals = ring_index.normals.astype('float32')
            ring_records = ring_index.rings
        else:
            rings_cached = fs.ensure_rings(structure)
            if not rings_cached:
                return []
            ring_centers = np.vstack([r['centroid'] for r in rings_cached]).astype('float32')
            ring_normals = np.vstack([r['normal'] for r in rings_cached]).astype('float32')
            ring_records = rings_cached  # dicts
        # Collect cation groups (reuse FeatureStore charged centers if possible)
        charged_centers = fs.ensure_charged_centers(structure)
        # Filter only residues we consider cationic
        charged_map = {}
        for entry in charged_centers:
            if entry['resname'] in CATION_RESIDUES:
                charged_map[id(entry['residue'])] = entry
        cation_residues: List[Residue.Residue] = []
        c_centers: List[np.ndarray] = []
        reused = 0
        recomputed = 0
        for chain in model:
            for residue in chain:
                if residue.get_resname() in CATION_RESIDUES:
                    cation_residues.append(residue)
                    key = id(residue)
                    if key in charged_map:
                        c_centers.append(charged_map[key]['centroid'])
                        reused += 1
                    else:
                        c_centers.append(self._cation_group_centroid(residue))
                        recomputed += 1
        if not c_centers:
            return []
        c_centers_arr = np.vstack(c_centers).astype('float32')
        ring_centers = ring_centers.astype('float32')
        # Use geometry core for pair pruning (KD-tree aware)
        import time
        t_pair_start = time.time()
        try:
            pairs = fs.neighbor_between(structure, float(self.distance_cutoff), c_centers_arr, ring_centers)
            if not pairs:
                return []
            k_idx = np.asarray([p[0] for p in pairs], dtype='int32')
            r_idx = np.asarray([p[1] for p in pairs], dtype='int32')
            used_core = True
        except Exception:
            diff = c_centers_arr[:, None, :] - ring_centers[None, :, :]
            dist = np.linalg.norm(diff, axis=-1)
            ii, jj = np.where((dist <= self.distance_cutoff))
            k_idx, r_idx = ii.astype(int), jj.astype(int)
            used_core = False
        t_pair_end = time.time()
        # On-demand distance compute if not precomputed
        have_dist_matrix = 'dist' in locals()
        interactions: List[Dict[str, Any]] = []
        # Instrumentation via helper
        raw_pairs = int(c_centers_arr.shape[0] * ring_centers.shape[0])
        candidate_pairs = int(len(k_idx))
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=0,
            core_pair_generation=used_core,
            extra={
                'cations': len(c_centers_arr),
                'rings': ring_centers.shape[0],
                'pairs_considered': raw_pairs,
                'pairs_within_cutoff': candidate_pairs,
                'charged_centers_reused': reused,
                'charged_centers_recomputed': recomputed
            }
        )
        t_eval_start = time.time()
        for ci, ri in zip(k_idx.tolist(), r_idx.tolist()):
            res_c = cation_residues[ci]
            ring_obj = ring_records[ri]
            # unify access for ring residue, normal
            if build_aromatic_ring_index is not None:
                ring_res = ring_obj.residue  # type: ignore[attr-defined]
                ring_normal = ring_normals[ri]
            else:  # dict
                ring_res = ring_obj['residue']
                ring_normal = ring_normals[ri]
            if res_c is ring_res:
                continue
            if have_dist_matrix:
                d = float(dist[ci, ri])  # type: ignore[index-defined]
            else:
                vec_tmp = c_centers_arr[ci] - ring_centers[ri]
                d = float(np.linalg.norm(vec_tmp))
            vec = c_centers_arr[ci] - ring_centers[ri]
            if np.linalg.norm(vec) > 1e-6:
                vec_n = vec / np.linalg.norm(vec)
                angle = float(np.degrees(np.arccos(np.clip(abs(np.dot(vec_n, ring_normal)), 0, 1))))
            else:
                angle = None
            if d < 4.5:
                strength = 'strong'
            elif d < 5.5:
                strength = 'moderate'
            else:
                strength = 'weak'
            interactions.append({
                'residue1': f"{res_c.get_resname()}{res_c.get_id()[1]}",
                'residue2': f"{ring_res.get_resname()}{ring_res.get_id()[1]}",
                'chain1': res_c.get_parent().get_id(),
                'chain2': ring_res.get_parent().get_id(),
                'distance': round(d, 3),
                'angle': round(angle, 1) if angle is not None else None,
                'interaction_type': 'cation_pi',
                'strength': strength
            })
        t_eval_end = time.time()
        update_counts(self.instrumentation, accepted=len(interactions))
        finalize_funnel(
            self.instrumentation,
            pair_gen_seconds=(t_pair_end - t_pair_start),
            eval_seconds=(t_eval_end - t_eval_start),
            build_seconds=0.0
        )
        return interactions

    def to_dict_list(self, interactions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        return interactions
