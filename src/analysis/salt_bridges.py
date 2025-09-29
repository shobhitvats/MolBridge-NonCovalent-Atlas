"""Salt bridge detection leveraging FeatureStore centroids.

This implementation reuses charged (positive) and acidic (negative) group centroids
produced for ionic / cation-π workflows to avoid per-atom nested loops.
"""
from typing import List, Dict, Any
import numpy as np
from .base_detector import register_detector
from utils.instrumentation import init_funnel, update_counts, finalize_funnel
from loguru import logger

@register_detector("salt_bridge", method="detect_salt_bridges")
class SaltBridgeDetector:
    def __init__(self, config):
        self.config = config
        inter = config.interactions
        self.centroid_cutoff = getattr(inter, 'salt_bridge_centroid_cutoff', 5.0)
        self.strong_cut = 3.5
        self.moderate_cut = 4.5

    def detect_salt_bridges(self, structure) -> List[Dict[str, Any]]:
        import time as _t
        t_pair_start = _t.time()
        try:
            from analysis.feature_store import get_feature_store
            fs = get_feature_store()
            positives = fs.ensure_charged_centers(structure) or []
            negatives = fs.ensure_acidic_centers(structure) or []
        except Exception as e:
            logger.error(f"Salt bridge centroid caching failed: {e}")
            return []
        if not positives or not negatives:
            return []
        t_pair_end = _t.time()
        t_eval_start = t_pair_end
        p = np.vstack([c['centroid'] for c in positives]).astype('float32')
        n = np.vstack([c['centroid'] for c in negatives]).astype('float32')
        try:
            from geometry.core import pairwise_within_cutoff
            pi, ni = pairwise_within_cutoff(p, n, float(self.centroid_cutoff), use_kdtree=True)
            dist_matrix = None
            core = True
        except Exception:
            diff = p[:, None, :] - n[None, :, :]
            dist_matrix = np.linalg.norm(diff, axis=-1)
            ii, jj = np.where(dist_matrix <= self.centroid_cutoff)
            pi, ni = ii.astype(int), jj.astype(int)
            core = False
        bridges: List[Dict[str, Any]] = []
        have_matrix = dist_matrix is not None
        for idx, (i, j) in enumerate(zip(pi.tolist(), ni.tolist())):
            if have_matrix:
                d = float(dist_matrix[i, j])
            else:
                d = float(np.linalg.norm(p[i] - n[j]))
            if d < self.strong_cut:
                strength = 'strong'
            elif d < self.moderate_cut:
                strength = 'moderate'
            else:
                strength = 'weak'
            pos = positives[i]
            neg = negatives[j]
            bridges.append({
                'positive_residue': f"{pos['resname']}{pos['residue'].get_id()[1]}",
                'negative_residue': f"{neg['resname']}{neg['residue'].get_id()[1]}",
                'positive_chain': pos['chain_id'],
                'negative_chain': neg['chain_id'],
                'distance': round(d, 3),
                'strength': strength,
                'type': 'salt_bridge'
            })
        t_eval_end = _t.time()
        raw_pairs = len(positives) * len(negatives)
        candidate_pairs = len(pi)
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=len(bridges),
            core_pair_generation=core,
            extra={
                'positive_centroids': len(positives),
                'negative_centroids': len(negatives),
                'pairs_considered': raw_pairs,
                'within_cutoff': len(bridges),
                'cutoff': self.centroid_cutoff
            }
        )
        finalize_funnel(
            self.instrumentation,
            pair_gen_seconds=(t_pair_end - t_pair_start),
            eval_seconds=(t_eval_end - t_eval_start),
            build_seconds=0.0
        )
        return bridges

    def to_dict_list(self, interactions: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        return interactions
