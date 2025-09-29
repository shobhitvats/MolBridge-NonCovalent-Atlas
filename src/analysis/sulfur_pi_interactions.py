"""Sulfur–π interaction detection.

Detects interactions between sulfur-containing side chains (MET, CYS disulfide or free) and aromatic rings (PHE, TYR, TRP, HIS).

Criteria (heuristic):
  - Sulfur atom (SG in CYS, SD in MET) to aromatic ring centroid distance <= 6.0 Å (configurable: sulfur_pi_distance_cutoff)
  - Vertical offset (projection of vector onto ring normal) absolute value <= 3.0 Å (sulfur_pi_max_perp) to ensure approach is near face region, but not compulsory for initial inclusion; we record it and filter softly.
  - Lateral (in-plane) distance <= 5.5 Å (derived from sqrt(total^2 - perp^2)).

Returned dict keys:
  residue_s, residue_ring, chain_s, chain_ring, distance, perp_offset, lateral_distance,
  interaction_type='sulfur_pi', strength (distance-based classification)

Strength heuristic:
  distance < 4.0 Å strong; 4.0–5.0 moderate; else weak.

Future refinements:
  - Distinguish disulfide-bonded cysteines by checking SG-SG connectivity.
  - Angle criteria relative to ring plane orientation stringency.
"""
from typing import List, Dict, Any, Tuple
import time
import numpy as np
from Bio.PDB import Structure, Residue
from utils.settings import get_settings
from utils.instrumentation import init_funnel, update_counts, finalize_funnel
from analysis.feature_store import get_feature_store
from geometry.core import pairwise_within_cutoff, norms, vector_fields
from .base_detector import register_detector

AROMATIC_RESIDUES = {"PHE","TYR","TRP","HIS"}
SULFUR_RESIDUES = {"CYS","MET"}

RING_ATOMS = {
    'PHE': ['CG','CD1','CD2','CE1','CE2','CZ'],
    'TYR': ['CG','CD1','CD2','CE1','CE2','CZ'],
    'TRP': ['CG','CD1','NE1','CE2','CD2','CE3','CZ2','CZ3','CH2'],
    'HIS': ['CG','ND1','CD2','CE1','NE2']
}

@register_detector("sulfur_pi", method="detect_sulfur_pi")
class SulfurPiDetector:
    def __init__(self, config):
        self.config = config
        self.distance_cutoff = getattr(self.config.interactions, 'sulfur_pi_distance_cutoff', 6.0)
        self.max_perp = getattr(self.config.interactions, 'sulfur_pi_max_perp', 3.0)
    
    def _ring_centroid_normal(self, residue: Residue.Residue) -> Tuple[np.ndarray, np.ndarray]:
        atoms = RING_ATOMS.get(residue.get_resname(), [])
        coords = []
        for n in atoms:
            if n in residue:
                coords.append(residue[n].get_coord())
        if len(coords) < 3:
            return np.array([0.0,0.0,0.0]), np.array([0,0,1])
        arr = np.vstack(coords)
        center = arr.mean(axis=0)
        centered = arr - center
        try:
            _,_,vh = np.linalg.svd(centered)
            normal = vh[-1]
            if normal[2] < 0:
                normal = -normal
        except Exception:
            normal = np.array([0,0,1])
        return center, normal
    
    def _sulfur_atom(self, residue: Residue.Residue):
        resname = residue.get_resname()
        if resname == 'CYS' and 'SG' in residue:
            return residue['SG']
        if resname == 'MET' and 'SD' in residue:
            return residue['SD']
        return None
    
    def detect_sulfur_pi(self, structure: Structure.Structure) -> List[Dict[str, Any]]:
        """Detect sulfur–π interactions (vector fast path with FeatureStore reuse).

        Fast path enabled when settings.enable_vector_geom is True; otherwise we
        fallback to legacy nested loops (still using FeatureStore ring cache for reuse).
        """
        settings = get_settings()
        if getattr(settings, 'enable_vector_geom', False):
            try:
                return self._vector_detect(structure)
            except Exception:  # pragma: no cover
                return self._legacy_detect(structure)
        return self._legacy_detect(structure)

    def _legacy_detect(self, structure: Structure.Structure) -> List[Dict[str, Any]]:
        interactions: List[Dict[str, Any]] = []
        t_pair_start = time.time()
        try:
            model = structure[0]
        except Exception:
            return interactions
        fs = get_feature_store()
        ring_cache = fs.ensure_rings(structure)
        rings = [(r['residue'].get_parent(), r['residue'], r['centroid'], r['normal']) for r in ring_cache]
        sulfurs = []
        for chain in model:
            for residue in chain:
                if residue.get_resname() in SULFUR_RESIDUES:
                    atom = self._sulfur_atom(residue)
                    if atom is not None:
                        sulfurs.append((chain,residue,atom))
        t_pair_end = time.time()
        t_eval_start = t_pair_end
        for ch_s, res_s, s_atom in sulfurs:
            s_coord = s_atom.get_coord()
            for ch_r, res_r, ring_center, ring_normal in rings:
                if res_s is res_r:
                    continue
                vec = s_coord - ring_center
                dist = float(np.linalg.norm(vec))
                if dist > self.distance_cutoff:
                    continue
                perp = abs(np.dot(vec, ring_normal))
                lateral_sq = max(dist**2 - perp**2, 0.0)
                lateral = lateral_sq**0.5
                if dist < 4.0:
                    strength = 'strong'
                elif dist < 5.0:
                    strength = 'moderate'
                else:
                    strength = 'weak'
                interactions.append({
                    'residue_s': f"{res_s.get_resname()}{res_s.get_id()[1]}",
                    'residue_ring': f"{res_r.get_resname()}{res_r.get_id()[1]}",
                    'chain_s': ch_s.get_id(),
                    'chain_ring': ch_r.get_id(),
                    'distance': round(dist,3),
                    'perp_offset': round(perp,3),
                    'lateral_distance': round(lateral,3),
                    'interaction_type': 'sulfur_pi',
                    'strength': strength
                })
        t_eval_end = time.time()
        raw_pairs = len(sulfurs) * len(rings)
        candidate_pairs = raw_pairs  # no pruning in legacy implementation
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=len(interactions),
            core_pair_generation=False,
            extra={'sulfurs': len(sulfurs), 'rings': len(rings), 'candidate_density': 1.0}
        )
        finalize_funnel(
            self.instrumentation,
            pair_gen_seconds=(t_pair_end - t_pair_start),
            eval_seconds=(t_eval_end - t_eval_start),
            build_seconds=0.0
        )
        return interactions

    def _vector_detect(self, structure: Structure.Structure) -> List[Dict[str, Any]]:
        t_pair_start = time.time()
        try:
            model = structure[0]
        except Exception:
            return []
        fs = get_feature_store()
        rings_cached = fs.ensure_rings(structure)
        if not rings_cached:
            return []
        ring_centers = np.vstack([r['centroid'] for r in rings_cached]).astype('float32')
        ring_normals = np.vstack([r['normal'] for r in rings_cached]).astype('float32')
        sulfur_atoms = []
        sulfur_meta = []
        for chain in model:
            for residue in chain:
                if residue.get_resname() in SULFUR_RESIDUES:
                    satom = self._sulfur_atom(residue)
                    if satom is not None:
                        sulfur_atoms.append(satom.get_coord().astype('float32') if hasattr(satom.get_coord(), 'astype') else satom.get_coord())
                        sulfur_meta.append((chain, residue, satom))
        if not sulfur_atoms:
            return []
        s_arr = np.vstack(sulfur_atoms).astype('float32')
        try:
            si, ri, dists = pairwise_within_cutoff(s_arr, ring_centers, float(self.distance_cutoff), use_kdtree=True, return_distances=True)
            core = True
        except Exception:  # pragma: no cover
            diff = s_arr[:, None, :] - ring_centers[None, :, :]
            dist2 = np.sum(diff*diff, axis=-1)
            mask = dist2 <= (self.distance_cutoff ** 2)
            si, ri = np.where(mask)
            si = si.astype('int32'); ri = ri.astype('int32')
            dists = np.sqrt(dist2[si, ri]).astype('float32')
            core = False
        t_pair_end = time.time()
        t_eval_start = t_pair_end
        interactions: List[Dict[str, Any]] = []
        raw_pairs = int(s_arr.shape[0] * ring_centers.shape[0])
        candidate_pairs = int(len(si))
        accepted = 0
        for idx in range(candidate_pairs):
            s_idx = int(si[idx]); r_idx = int(ri[idx])
            chain_s, res_s, s_atom = sulfur_meta[s_idx]
            ring_record = rings_cached[r_idx]
            res_r = ring_record['residue']
            if res_s is res_r:
                continue
            d = float(dists[idx])
            vec = s_arr[s_idx] - ring_centers[r_idx]
            ring_normal = ring_normals[r_idx]
            perp = abs(float(np.dot(vec, ring_normal)))
            lateral_sq = max(d**2 - perp**2, 0.0)
            lateral = lateral_sq**0.5
            if d < 4.0:
                strength = 'strong'
            elif d < 5.0:
                strength = 'moderate'
            else:
                strength = 'weak'
            interactions.append({
                'residue_s': f"{res_s.get_resname()}{res_s.get_id()[1]}",
                'residue_ring': f"{res_r.get_resname()}{res_r.get_id()[1]}",
                'chain_s': chain_s.get_id(),
                'chain_ring': res_r.get_parent().get_id(),
                'distance': round(d,3),
                'perp_offset': round(perp,3),
                'lateral_distance': round(lateral,3),
                'interaction_type': 'sulfur_pi',
                'strength': strength
            })
            accepted += 1
        t_eval_end = time.time()
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=accepted,
            core_pair_generation=core,
            extra={
                'sulfurs': s_arr.shape[0],
                'rings': ring_centers.shape[0],
                'candidate_density': (candidate_pairs / raw_pairs) if raw_pairs else 0.0
            }
        )
        finalize_funnel(
            self.instrumentation,
            pair_gen_seconds=(t_pair_end - t_pair_start),
            eval_seconds=(t_eval_end - t_eval_start),
            build_seconds=0.0
        )
        return interactions
    
    def to_dict_list(self, interactions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        return interactions
