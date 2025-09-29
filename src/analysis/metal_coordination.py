"""Metal coordination detection.

Identifies coordination interactions between metal ions (common biologically: Zn, Fe, Mg, Mn, Ca, Cu, Co, Ni, Na, K) and nearby donor atoms from amino acid side chains or ligands.

Heuristic approach:
  - Metals recognized by element symbol (case-insensitive) in atom.element or residue name (e.g., 'ZN', 'FE').
  - Donor atoms from residues:
      * CYS: SG
      * HIS: ND1, NE2
      * ASP: OD1, OD2
      * GLU: OE1, OE2
      * TYR: OH
      * SER: OG
      * THR: OG1
      * MET: SD
      * ASN: OD1
      * GLN: OE1
      * LYS: NZ (rare but include heuristically)
      * backbone carbonyl O (O) optionally (enabled via flag)
  - Distance cutoff: metal_coordination_distance_cutoff (default 2.6 Å for strong inner-sphere; we extend to 3.0 Å upper bound to include borderline). We'll implement two tiers: primary <=2.6, extended <=3.0.

Returned dict fields:
  metal_id, metal_element, residue_donor, chain_donor, donor_atom, distance, tier (primary/extended), interaction_type='metal_coordination'

Future improvements (not implemented):
  - Geometry validation (tetrahedral, octahedral scoring)
  - Water-mediated links
  - Charge / oxidation state inference
"""
from typing import List, Dict, Any, Tuple
import numpy as np
from Bio.PDB import Structure, Residue, Atom
from utils.settings import get_settings
from utils.instrumentation import init_funnel, update_counts, finalize_funnel
from geometry.core import pairwise_within_cutoff, norms

METAL_ELEMENTS = {"ZN","FE","MG","MN","CA","CU","CO","NI","NA","K","CD"}
DONOR_MAP = {
    'CYS': ['SG'],
    'HIS': ['ND1','NE2'],
    'ASP': ['OD1','OD2'],
    'GLU': ['OE1','OE2'],
    'TYR': ['OH'],
    'SER': ['OG'],
    'THR': ['OG1'],
    'MET': ['SD'],
    'ASN': ['OD1'],
    'GLN': ['OE1'],
    'LYS': ['NZ']
}

class MetalCoordinationDetector:
    def __init__(self, config):
        self.config = config
        self.primary_cutoff = getattr(self.config.interactions, 'metal_coordination_primary_cutoff', 2.6)
        self.extended_cutoff = getattr(self.config.interactions, 'metal_coordination_extended_cutoff', 3.0)
        self.include_backbone = getattr(self.config.interactions, 'metal_include_backbone_carbonyl', True)
    
    def _is_metal_atom(self, atom: Atom.Atom, residue: Residue.Residue) -> bool:
        element = (atom.element or '').upper()
        resname = residue.get_resname().upper()
        return element in METAL_ELEMENTS or resname in METAL_ELEMENTS
    
    def _collect_metals(self, structure: Structure.Structure):
        metals = []
        try:
            model = structure[0]
        except Exception:
            return metals
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if self._is_metal_atom(atom, residue):
                        metals.append((chain, residue, atom))
                        break
        return metals
    
    def _collect_donors(self, structure: Structure.Structure):
        donors = []
        try:
            model = structure[0]
        except Exception:
            return donors
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                # side chain donors
                if resname in DONOR_MAP:
                    for name in DONOR_MAP[resname]:
                        if name in residue:
                            donors.append((chain, residue, residue[name]))
                # backbone carbonyl O
                if self.include_backbone and 'O' in residue:
                    donors.append((chain, residue, residue['O']))
        return donors
    
    def detect_metal_coordination(self, structure: Structure.Structure) -> List[Dict[str, Any]]:
        settings = get_settings()
        if getattr(settings, 'enable_vector_geom', False):
            try:
                return self._vector_detect(structure)
            except Exception:  # pragma: no cover
                pass
        return self._legacy_detect(structure)

    def _legacy_detect(self, structure: Structure.Structure) -> List[Dict[str, Any]]:
        import time as _t
        t_pair_start = _t.time()
        interactions: List[Dict[str, Any]] = []
        metals = self._collect_metals(structure)
        donors = self._collect_donors(structure)
        if not metals or not donors:
            self.instrumentation = init_funnel(
                raw_pairs=len(metals)*len(donors),
                candidate_pairs=0,
                accepted_pairs=0,
                core_pair_generation=False,
                extra={
                    'metals': len(metals),
                    'donors': len(donors),
                    'candidate_density': 0.0
                }
            )
            finalize_funnel(self.instrumentation, pair_gen_seconds=0.0, eval_seconds=0.0, build_seconds=0.0)
            return interactions
        candidate_pairs = 0
        t_pair_end = _t.time()
        t_eval_start = t_pair_end
        for ch_m, res_m, atom_m in metals:
            m_coord = atom_m.get_coord()
            for ch_d, res_d, atom_d in donors:
                if res_m is res_d:
                    continue
                dist = float(np.linalg.norm(m_coord - atom_d.get_coord()))
                if dist <= self.extended_cutoff:
                    candidate_pairs += 1
                    tier = 'primary' if dist <= self.primary_cutoff else 'extended'
                    interactions.append({
                        'metal_id': f"{res_m.get_resname()}{res_m.get_id()[1]}",
                        'metal_element': (atom_m.element or res_m.get_resname()).upper(),
                        'residue_donor': f"{res_d.get_resname()}{res_d.get_id()[1]}",
                        'chain_donor': ch_d.get_id(),
                        'donor_atom': atom_d.get_name(),
                        'distance': round(dist,3),
                        'tier': tier,
                        'interaction_type': 'metal_coordination'
                    })
        raw_pairs = len(metals)*len(donors)
        t_eval_end = _t.time()
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=len(interactions),
            core_pair_generation=False,
            extra={
                'metals': len(metals),
                'donors': len(donors),
                'candidate_density': (candidate_pairs/raw_pairs) if raw_pairs else 0.0
            }
        )
        finalize_funnel(
            self.instrumentation,
            pair_gen_seconds=(t_pair_end - t_pair_start),
            eval_seconds=(t_eval_end - t_eval_start),
            build_seconds=0.0
        )
        return interactions

    def _vector_detect(self, structure: Structure.Structure) -> List[Dict[str, Any]]:
        import time as _t
        t_pair_start = _t.time()
        interactions: List[Dict[str, Any]] = []
        metals = self._collect_metals(structure)
        donors = self._collect_donors(structure)
        if not metals or not donors:
            self.instrumentation = init_funnel(
                raw_pairs=len(metals)*len(donors),
                candidate_pairs=0,
                accepted_pairs=0,
                core_pair_generation=False,
                extra={
                    'metals': len(metals),
                    'donors': len(donors),
                    'candidate_density': 0.0
                }
            )
            finalize_funnel(self.instrumentation, pair_gen_seconds=0.0, eval_seconds=0.0, build_seconds=0.0)
            return interactions
        m_coords = np.vstack([atm.get_coord() for _,_,atm in metals]).astype('float32')
        d_coords = np.vstack([atm.get_coord() for _,_,atm in donors]).astype('float32')
        # Pair pruning with extended cutoff
        try:
            mi, di = pairwise_within_cutoff(m_coords, d_coords, self.extended_cutoff, use_kdtree=True)
            core = True
        except Exception:  # pragma: no cover
            diff = m_coords[:, None, :] - d_coords[None, :, :]
            dist2 = np.sum(diff*diff, axis=-1)
            mask = dist2 <= (self.extended_cutoff ** 2)
            mi, di = np.where(mask)
            mi = mi.astype('int32'); di = di.astype('int32')
            core = False
        raw_pairs = int(m_coords.shape[0] * d_coords.shape[0])
        t_pair_end = _t.time()
        t_eval_start = t_pair_end
        vecs = d_coords[di] - m_coords[mi]
        distances = norms(vecs)
        # Filter distances (<= extended) already guaranteed, just iterate
        candidate_pairs = int(len(mi))
        accepted = 0
        for idx in range(candidate_pairs):
            m_idx = int(mi[idx]); d_idx = int(di[idx])
            ch_m, res_m, atom_m = metals[m_idx]
            ch_d, res_d, atom_d = donors[d_idx]
            if res_m is res_d:
                continue
            dist = float(distances[idx])
            tier = 'primary' if dist <= self.primary_cutoff else 'extended'
            interactions.append({
                'metal_id': f"{res_m.get_resname()}{res_m.get_id()[1]}",
                'metal_element': (atom_m.element or res_m.get_resname()).upper(),
                'residue_donor': f"{res_d.get_resname()}{res_d.get_id()[1]}",
                'chain_donor': ch_d.get_id(),
                'donor_atom': atom_d.get_name(),
                'distance': round(dist,3),
                'tier': tier,
                'interaction_type': 'metal_coordination'
            })
            accepted += 1
        t_eval_end = _t.time()
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=accepted,
            core_pair_generation=core,
            extra={
                'metals': m_coords.shape[0],
                'donors': d_coords.shape[0],
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
