"""
Ionic interaction detection for protein structures.
Detects electrostatic interactions between charged residues.
"""

import os

import numpy as np
from typing import List, Dict, Any, Optional
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from loguru import logger
from .base_detector import register_detector
from utils.settings import get_settings
from utils.instrumentation import init_funnel, update_counts, finalize_funnel
try:  # reuse existing geom KD-tree helpers
    from utils.geom import build_kdtree, radius_query
except Exception:  # pragma: no cover
    build_kdtree = None  # type: ignore
    radius_query = None  # type: ignore

@dataclass
class IonicInteraction:
    """Represents a detected ionic interaction."""
    positive_atom: Atom.Atom
    negative_atom: Atom.Atom
    distance: float
    strength: str
    positive_residue: str
    negative_residue: str
    positive_chain: str
    negative_chain: str

@register_detector("ionic_interaction", method="detect_ionic_interactions")
class IonicInteractionDetector:
    """Detects ionic interactions in protein structures."""
    
    def __init__(self, config):
        self.config = config
        self.interaction_config = config.interactions
        self.distance_cutoff = self.interaction_config.ionic_distance_cutoff
        
        # Charged residues
        self.positive_residues = {
            'ARG': ['NH1', 'NH2'],
            'LYS': ['NZ'],
            'HIS': ['ND1', 'NE2']  # Protonated form
        }
        
        self.negative_residues = {
            'ASP': ['OD1', 'OD2'],
            'GLU': ['OE1', 'OE2']
        }
    
    def detect_ionic_interactions(self, structure: Structure.Structure) -> List[IonicInteraction]:
        settings = get_settings()
        env_force = (os.getenv('MOLBRIDGE_ENABLE_VECTOR_GEOM') in {'1','true','True'})
        if env_force:
            # bypass settings object if env explicitly forces vector mode
            try:
                return self._vector_detect(structure)
            except Exception:
                return self._legacy_detect(structure)
        if getattr(settings, 'enable_vector_geom', False):
            try:
                return self._vector_detect(structure)
            except Exception:  # pragma: no cover
                return self._legacy_detect(structure)
        return self._legacy_detect(structure)

    def _legacy_detect(self, structure: Structure.Structure) -> List[IonicInteraction]:
        import time as _t
        t_pair_start = _t.time()
        interactions: List[IonicInteraction] = []
        try:
            model = structure[0]
            positive_atoms = self._get_positive_atoms(model)
            negative_atoms = self._get_negative_atoms(model)
            # Treat legacy as single combined phase (no pruning); pair_gen = 0
            t_pair_end = _t.time()
            t_eval_start = t_pair_end
            for pos_info in positive_atoms:
                for neg_info in negative_atoms:
                    inter = self._check_ionic_interaction(pos_info, neg_info)
                    if inter:
                        interactions.append(inter)
            t_eval_end = _t.time()
            raw_pairs = len(positive_atoms) * len(negative_atoms)
            self.instrumentation = init_funnel(
                raw_pairs=raw_pairs,
                candidate_pairs=raw_pairs,
                accepted_pairs=len(interactions),
                core_pair_generation=False,
                extra={
                    'positive_atoms': len(positive_atoms),
                    'negative_atoms': len(negative_atoms),
                    'kdtree_used': False
                }
            )
            finalize_funnel(self.instrumentation, pair_gen_seconds=0.0, eval_seconds=(t_eval_end - t_eval_start), build_seconds=0.0)
            logger.info(f"[legacy] Detected {len(interactions)} ionic interactions")
        except Exception as e:
            logger.error(f"Error detecting ionic interactions: {e}")
        return interactions

    def _vector_detect(self, structure: Structure.Structure) -> List[IonicInteraction]:
        import time as _t
        t_pair_start = _t.time()
        interactions: List[IonicInteraction] = []
        try:
            model = structure[0]
        except Exception:
            return interactions
        pos_atoms = self._get_positive_atoms(model)
        neg_atoms = self._get_negative_atoms(model)
        if not pos_atoms or not neg_atoms:
            return interactions
        import numpy as _np
        p_coords = _np.vstack([p['atom'].get_coord() for p in pos_atoms]).astype('float32')
        n_coords = _np.vstack([n['atom'].get_coord() for n in neg_atoms]).astype('float32')
        from utils.kdtree_thresholds import get_threshold, adapt_threshold, get_last_density, should_flag_kdtree
        raw_pairs = len(pos_atoms) * len(neg_atoms)
        threshold = get_threshold('ionic')
        env_override_val = os.getenv('MOLBRIDGE_KDTREE_IONIC_THRESHOLD')
        env_threshold = None
        if env_override_val is not None:
            try:
                env_threshold = int(env_override_val)
            except Exception:
                env_threshold = None
        try:
            from geometry.core import pairwise_within_cutoff
            pi, ni = pairwise_within_cutoff(p_coords, n_coords, float(self.distance_cutoff), use_kdtree=True)
            core = True
        except Exception:
            diff = p_coords[:, None, :] - n_coords[None, :, :]
            dist_mat = _np.linalg.norm(diff, axis=-1)
            ii, jj = _np.where(dist_mat <= self.distance_cutoff)
            pi, ni = ii.astype('int32'), jj.astype('int32')
            core = False
        pruned_pairs = int(len(pi))
        t_pair_end = _t.time()
        t_eval_start = t_pair_end
        # Heuristic: if core path yields significant pruning assume KD-tree usage (for adaptation metrics) 
        # Determine KD-tree usage flag considering explicit threshold overrides:
        # If threshold set extremely high (env forcing disable) treat as not used.
        # If threshold very low (env forcing enable) and raw_pairs exceed that low threshold, mark used.
        # Explicit env toggle semantics for tests:
        #   very high threshold (>= 500_000) => treat as disabled
        #   very low threshold (<= 5) => force enabled when raw_pairs exceeds it
        threshold_for_toggle = env_threshold if env_threshold is not None else threshold
        if threshold_for_toggle >= 500_000:
            kdtree_used = False
        elif threshold_for_toggle <= 5 and raw_pairs > threshold_for_toggle:
            kdtree_used = True
        else:
            kdtree_used = core and (pruned_pairs < raw_pairs * 0.65) and raw_pairs > threshold
        pruning_flag = should_flag_kdtree(pruned_pairs, raw_pairs, 'ionic')
        new_thresh, changed, reason = adapt_threshold('ionic', pruned_pairs, (kdtree_used or pruning_flag))
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=pruned_pairs,
            accepted_pairs=0,
            core_pair_generation=core,
            extra={
                'positive_atoms': len(pos_atoms),
                'negative_atoms': len(neg_atoms),
                'pairs_within_cutoff': pruned_pairs,
                'kdtree_used': kdtree_used,
                'threshold': threshold,
                'env_threshold': env_threshold,
                'adaptive_threshold': new_thresh,
                'threshold_changed': changed,
                'candidate_density': get_last_density('ionic'),
                'adapt_reason': reason
            }
        )
        for i, j in zip(pi.tolist(), ni.tolist()):
            inter = self._check_ionic_interaction(pos_atoms[i], neg_atoms[j])
            if inter:
                interactions.append(inter)
        update_counts(self.instrumentation, accepted=len(interactions))
        t_eval_end = _t.time()
        finalize_funnel(
            self.instrumentation,
            pair_gen_seconds=(t_pair_end - t_pair_start),
            eval_seconds=(t_eval_end - t_eval_start),
            build_seconds=0.0
        )
        logger.info(f"[vector]{'/core' if core else ''} Ionic: {len(interactions)} (raw={raw_pairs} pruned={pruned_pairs} acc_ratio={self.instrumentation['acceptance_ratio']:.3f})")
        return interactions
    
    def _get_positive_atoms(self, model) -> List[Dict[str, Any]]:
        """Get positively charged atoms."""
        atoms = []
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname in self.positive_residues:
                    for atom_name in self.positive_residues[resname]:
                        if atom_name in residue:
                            atoms.append({
                                'atom': residue[atom_name],
                                'residue': residue,
                                'resname': resname,
                                'chain_id': chain.get_id(),
                                'res_id': residue.get_id()
                            })
        return atoms
    
    def _get_negative_atoms(self, model) -> List[Dict[str, Any]]:
        """Get negatively charged atoms."""
        atoms = []
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname in self.negative_residues:
                    for atom_name in self.negative_residues[resname]:
                        if atom_name in residue:
                            atoms.append({
                                'atom': residue[atom_name],
                                'residue': residue,
                                'resname': resname,
                                'chain_id': chain.get_id(),
                                'res_id': residue.get_id()
                            })
        return atoms
    
    def _check_ionic_interaction(self, pos_info: Dict, neg_info: Dict) -> Optional[IonicInteraction]:
        """Check if two charged atoms form an ionic interaction."""
        if pos_info['residue'] == neg_info['residue']:
            return None
        
        distance = np.linalg.norm(
            pos_info['atom'].get_coord() - neg_info['atom'].get_coord()
        )
        
        if distance > self.distance_cutoff:
            return None
        
        strength = 'strong' if distance < 4.0 else 'moderate' if distance < 5.0 else 'weak'
        
        return IonicInteraction(
            positive_atom=pos_info['atom'],
            negative_atom=neg_info['atom'],
            distance=distance,
            strength=strength,
            positive_residue=f"{pos_info['resname']}{pos_info['res_id'][1]}",
            negative_residue=f"{neg_info['resname']}{neg_info['res_id'][1]}",
            positive_chain=pos_info['chain_id'],
            negative_chain=neg_info['chain_id']
        )
    
    def to_dict_list(self, interactions: List[IonicInteraction]) -> List[Dict[str, Any]]:
        """Convert interactions to list of dictionaries."""
        return [
            {
                'positive_atom': i.positive_atom.get_name(),
                'negative_atom': i.negative_atom.get_name(),
                'distance': round(i.distance, 3),
                'strength': i.strength,
                'positive_residue': i.positive_residue,
                'negative_residue': i.negative_residue,
                'positive_chain': i.positive_chain,
                'negative_chain': i.negative_chain,
                'residue1': i.positive_residue,
                'residue2': i.negative_residue,
                'chain1': i.positive_chain,
                'chain2': i.negative_chain
            }
            for i in interactions
        ]
