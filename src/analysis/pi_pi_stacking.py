"""
π-π stacking interaction detection for protein structures.
Detects aromatic stacking interactions between aromatic rings.
"""

import numpy as np
import time
from typing import List, Dict, Any, Tuple, Optional
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from loguru import logger
from utils.settings import get_settings  # feature flag gating
from analysis.feature_store import get_feature_store
from utils.instrumentation import init_funnel, update_counts, finalize_funnel
try:  # optional index (new vector path)
    from analysis.structure_index import build_aromatic_ring_index, AromaticRingIndex
except Exception:  # pragma: no cover - if import fails we just skip vector path
    build_aromatic_ring_index = None  # type: ignore
    AromaticRingIndex = None  # type: ignore

@dataclass
class PiPiInteraction:
    """Represents a detected π-π stacking interaction."""
    ring1_residue: str
    ring2_residue: str
    ring1_chain: str
    ring2_chain: str
    ring1_center: np.ndarray
    ring2_center: np.ndarray
    distance: float
    angle: float
    offset: float
    interaction_type: str  # 'face_to_face', 'edge_to_face', 'offset'
    strength: str

@dataclass
class AromaticRing:
    """Represents an aromatic ring."""
    residue: Residue.Residue
    ring_type: str  # 'PHE', 'TYR', 'TRP', 'HIS'
    atoms: List[Atom.Atom]
    center: np.ndarray
    normal: np.ndarray
    plane_coeffs: np.ndarray

from .base_detector import register_detector

@register_detector("pipi", method="detect_pi_pi_interactions")
class PiPiDetector:
    """Detects π-π stacking interactions in protein structures.

    Vector fast-path: Enabled when MOLBRIDGE_ENABLE_VECTOR_GEOM=1 (via utils.settings)
    AND the structure_index module is importable. Falls back to legacy nested loops
    otherwise. This preserves prior behavior while providing a performance lift
    for large structures with many aromatics.
    """

    def __init__(self, config):
        self.config = config
        self.interaction_config = config.interactions
        self.distance_cutoff = self.interaction_config.pi_pi_distance_cutoff
        self.angle_cutoff = self.interaction_config.pi_pi_angle_cutoff
        self.offset_cutoff = getattr(self.interaction_config, 'pi_pi_offset_cutoff', 3.0)

        # Retain legacy ring atom definitions for fallback path
        self.aromatic_atoms = {
            'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TRP': {
                'indole': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
                'pyrrole': ['CG', 'CD1', 'NE1', 'CE2', 'CD2'],
                'benzene': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']
            },
            'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2']
        }
    
    def detect_pi_pi_interactions(self, structure: Structure.Structure) -> List[PiPiInteraction]:
        """Public entrypoint resolving to vector or legacy path."""
        settings = get_settings()
        use_vector = bool(getattr(settings, 'enable_vector_geom', False)) and build_aromatic_ring_index is not None
        if not use_vector:
            return self._legacy_detect(structure)
        try:
            return self._vector_detect(structure)
        except Exception as e:  # safety net
            logger.error(f"Vector π-π path failed ({e}); falling back to legacy")
            return self._legacy_detect(structure)

    # ---------------- Legacy Path (unchanged logic) -----------------
    def _legacy_detect(self, structure: Structure.Structure) -> List[PiPiInteraction]:
        interactions: List[PiPiInteraction] = []
        try:
            model = structure[0]
            rings = self._get_aromatic_rings(model)
            logger.info(f"[legacy] Found {len(rings)} aromatic rings")
            for i, r1 in enumerate(rings):
                for r2 in rings[i+1:]:
                    inter = self._check_pi_pi_interaction(r1, r2)
                    if inter:
                        interactions.append(inter)
            logger.info(f"[legacy] Detected {len(interactions)} π-π interactions")
        except Exception as e:
            logger.error(f"Error (legacy) detecting π-π interactions: {e}")
        return interactions

    # ---------------- Vector Path -----------------
    def _vector_detect(self, structure: Structure.Structure) -> List[PiPiInteraction]:
        interactions: List[PiPiInteraction] = []
        t_pair_start = time.time()
        model = structure[0]
        fs = get_feature_store()
        # Attempt to reuse cached rings (ensure_rings returns list of dicts)
        rings_cached = fs.ensure_rings(structure)
        if rings_cached and build_aromatic_ring_index is None:
            # Build a lightweight centers/normals array from cached rings
            centers = np.vstack([r['centroid'] for r in rings_cached]) if len(rings_cached) else np.zeros((0,3))
            normals = np.vstack([r['normal'] for r in rings_cached]) if len(rings_cached) else np.zeros((0,3))
            class _RingWrapper:
                def __init__(self, rings, centers, normals):
                    self.rings = rings
                    self.centers = centers
                    self.normals = normals
                def __len__(self):
                    return len(self.rings)
            ring_index = _RingWrapper(rings_cached, centers, normals)
        else:
            if build_aromatic_ring_index is None:  # fallback to legacy if no index and no cached rings
                return self._legacy_detect(structure)
            ring_index = build_aromatic_ring_index(model)
        n = len(ring_index)
        if n < 2:
            return []
        centers = ring_index.centers
        normals = ring_index.normals
        # Candidate pruning via shared FeatureStore subset neighbor API
        from analysis.feature_store import get_feature_store as _gfs
        fs2 = _gfs()
        cutoff = float(self.distance_cutoff)
        pairs = fs2.neighbor_within(structure, cutoff, points=centers.astype('float32'))
        if not pairs:
            return []
        ia = np.asarray([p[0] for p in pairs], dtype='int32')
        ib = np.asarray([p[1] for p in pairs], dtype='int32')
        # distance compute on demand
        dists = np.linalg.norm(centers[ib] - centers[ia], axis=1)
        core_pairs = True  # treat as pruning path for instrumentation
        # Adaptive KD-tree threshold style decision (even if we don't use a KD-tree yet, we track complexity)
        from utils.kdtree_thresholds import get_threshold, adapt_threshold, get_last_density, should_flag_kdtree
        raw_pairs = int(n * (n - 1) / 2)
        threshold = get_threshold('pipi')
        # Placeholder for future KD-tree (not implemented): treat large pairing as 'kdtree_used' if above threshold to drive adaptation
        kdtree_used = raw_pairs > threshold * 1.3  # synthetic heuristic
        pruned_pairs = len(ia)
        pruning_flag = should_flag_kdtree(pruned_pairs, raw_pairs, 'pipi')
        new_thresh, changed, reason = adapt_threshold('pipi', pruned_pairs, (kdtree_used or pruning_flag))
        t_pair_end = time.time()
        t_eval_start = time.time()
        for idx, (i_idx, j_idx) in enumerate(zip(ia.tolist(), ib.tolist())):
            r1 = ring_index.rings[i_idx]
            r2 = ring_index.rings[j_idx]
            # Skip same residue (shouldn't happen but keep parity with legacy)
            if r1.residue == r2.residue:
                continue
            distance = float(dists[idx])
            cos_angle = abs(float(np.clip(np.dot(normals[i_idx], normals[j_idx]), -1.0, 1.0)))
            ang = float(np.degrees(np.arccos(cos_angle)))
            angle = min(ang, 180 - ang)
            center_vec = centers[j_idx] - centers[i_idx]
            parallel = np.dot(center_vec, normals[i_idx])
            perp_vec = center_vec - parallel * normals[i_idx]
            offset = float(np.linalg.norm(perp_vec))
            interaction_type = self._classify_pi_pi_type(angle, offset, distance)
            if not self._meets_criteria(interaction_type, angle, offset):
                continue
            strength = self._calculate_pi_pi_strength(distance, angle, offset, interaction_type)
            interactions.append(PiPiInteraction(
                ring1_residue=f"{r1.residue.get_resname()}{r1.residue.get_id()[1]}",
                ring2_residue=f"{r2.residue.get_resname()}{r2.residue.get_id()[1]}",
                ring1_chain=r1.residue.get_parent().get_id(),
                ring2_chain=r2.residue.get_parent().get_id(),
                ring1_center=r1.center,
                ring2_center=r2.center,
                distance=float(distance),
                angle=float(angle),
                offset=float(offset),
                interaction_type=interaction_type,
                strength=strength
            ))
        t_eval_end = time.time()
        # Funnel instrumentation using shared helpers
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=pruned_pairs,
            accepted_pairs=len(interactions),
            core_pair_generation=core_pairs,
            extra={
                'rings': n,
                'kdtree_used': int(kdtree_used),
                'threshold': threshold,
                'adaptive_threshold': new_thresh,
                'threshold_changed': changed,
                'candidate_density': get_last_density('pipi'),
                'adapt_reason': reason
            }
        )
        finalize_funnel(
            self.instrumentation,
            pair_gen_seconds=(t_pair_end - t_pair_start),
            eval_seconds=(t_eval_end - t_eval_start),
            build_seconds=0.0
        )
        logger.info(
            f"[vector]{'/adaptive' if changed else ''} π-π: {len(interactions)} (rings={n} raw={raw_pairs} pruned={pruned_pairs} acc_ratio={self.instrumentation['acceptance_ratio']:.3f})"
        )
        return interactions
    
    def _get_aromatic_rings(self, model) -> List[AromaticRing]:
        """Get all aromatic rings in the structure."""
        rings = []
        
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                
                if resname in ['PHE', 'TYR', 'HIS']:
                    ring = self._create_simple_ring(residue, resname)
                    if ring:
                        rings.append(ring)
                
                elif resname == 'TRP':
                    # Tryptophan has multiple rings
                    trp_rings = self._create_tryptophan_rings(residue)
                    rings.extend(trp_rings)
        
        return rings
    
    def _create_simple_ring(self, residue: Residue.Residue, resname: str) -> Optional[AromaticRing]:
        """Create aromatic ring for PHE, TYR, or HIS."""
        if resname not in self.aromatic_atoms:
            return None
        
        ring_atoms = []
        for atom_name in self.aromatic_atoms[resname]:
            if atom_name in residue:
                ring_atoms.append(residue[atom_name])
        
        if len(ring_atoms) < 3:  # Need at least 3 atoms to define a plane
            return None
        
        # Calculate ring center and normal
        center, normal, plane_coeffs = self._calculate_ring_geometry(ring_atoms)
        
        return AromaticRing(
            residue=residue,
            ring_type=resname,
            atoms=ring_atoms,
            center=center,
            normal=normal,
            plane_coeffs=plane_coeffs
        )
    
    def _create_tryptophan_rings(self, residue: Residue.Residue) -> List[AromaticRing]:
        """Create aromatic rings for tryptophan (indole system)."""
        rings = []
        
        # For simplicity, we'll treat the entire indole as one ring
        # In a more sophisticated implementation, you might separate pyrrole and benzene rings
        
        indole_atoms = []
        for atom_name in self.aromatic_atoms['TRP']['indole']:
            if atom_name in residue:
                indole_atoms.append(residue[atom_name])
        
        if len(indole_atoms) >= 6:  # Minimum for indole ring
            center, normal, plane_coeffs = self._calculate_ring_geometry(indole_atoms)
            
            rings.append(AromaticRing(
                residue=residue,
                ring_type='TRP_indole',
                atoms=indole_atoms,
                center=center,
                normal=normal,
                plane_coeffs=plane_coeffs
            ))
        
        return rings
    
    def _calculate_ring_geometry(self, atoms: List[Atom.Atom]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Calculate ring center, normal vector, and plane coefficients."""
        # Get atom coordinates
        coords = np.array([atom.get_coord() for atom in atoms])
        
        # Calculate center as centroid
        center = np.mean(coords, axis=0)
        
        # Calculate normal vector using SVD
        centered_coords = coords - center
        _, _, vh = np.linalg.svd(centered_coords)
        normal = vh[-1]  # Last row of V^T is normal to best-fit plane
        
        # Ensure normal points in consistent direction
        if normal[2] < 0:
            normal = -normal
        
        # Calculate plane coefficients (ax + by + cz + d = 0)
        d = -np.dot(normal, center)
        plane_coeffs = np.append(normal, d)
        
        return center, normal, plane_coeffs
    
    def _check_pi_pi_interaction(self, ring1: AromaticRing, ring2: AromaticRing) -> Optional[PiPiInteraction]:
        """Check if two aromatic rings form a π-π interaction."""
        # Skip if same residue
        if ring1.residue == ring2.residue:
            return None
        
        # Calculate center-to-center distance
        distance = np.linalg.norm(ring1.center - ring2.center)
        
        if distance > self.distance_cutoff:
            return None
        
        # Calculate angle between ring planes
        angle = self._calculate_ring_angle(ring1.normal, ring2.normal)
        
        # Calculate offset (perpendicular distance)
        offset = self._calculate_offset(ring1, ring2)
        
        # Classify interaction type
        interaction_type = self._classify_pi_pi_type(angle, offset, distance)
        
        # Check if interaction meets criteria
        if not self._meets_criteria(interaction_type, angle, offset):
            return None
        
        # Calculate strength
        strength = self._calculate_pi_pi_strength(distance, angle, offset, interaction_type)
        
        return PiPiInteraction(
            ring1_residue=f"{ring1.residue.get_resname()}{ring1.residue.get_id()[1]}",
            ring2_residue=f"{ring2.residue.get_resname()}{ring2.residue.get_id()[1]}",
            ring1_chain=ring1.residue.get_parent().get_id(),
            ring2_chain=ring2.residue.get_parent().get_id(),
            ring1_center=ring1.center,
            ring2_center=ring2.center,
            distance=distance,
            angle=angle,
            offset=offset,
            interaction_type=interaction_type,
            strength=strength
        )
    
    def _calculate_ring_angle(self, normal1: np.ndarray, normal2: np.ndarray) -> float:
        """Calculate angle between two ring planes."""
        cos_angle = np.abs(np.dot(normal1, normal2))
        cos_angle = np.clip(cos_angle, 0.0, 1.0)
        
        # Return acute angle
        angle = np.degrees(np.arccos(cos_angle))
        return min(angle, 180 - angle)
    
    def _calculate_offset(self, ring1: AromaticRing, ring2: AromaticRing) -> float:
        """Calculate offset between ring centers."""
        # Vector from ring1 center to ring2 center
        center_vector = ring2.center - ring1.center
        
        # Project onto ring1 normal to get perpendicular component
        parallel_component = np.dot(center_vector, ring1.normal)
        perpendicular_vector = center_vector - parallel_component * ring1.normal
        
        return np.linalg.norm(perpendicular_vector)
    
    def _classify_pi_pi_type(self, angle: float, offset: float, distance: float) -> str:
        """Classify the type of π-π interaction."""
        if angle < 30:  # Nearly parallel
            if offset < 2.0:
                return 'face_to_face'
            else:
                return 'offset'
        elif angle > 60:  # Nearly perpendicular
            return 'edge_to_face'
        else:
            return 'offset'
    
    def _meets_criteria(self, interaction_type: str, angle: float, offset: float) -> bool:
        """Check if interaction meets geometric criteria."""
        if interaction_type == 'face_to_face':
            return angle <= self.angle_cutoff and offset <= self.offset_cutoff
        
        elif interaction_type == 'edge_to_face':
            return angle >= (90 - self.angle_cutoff)
        
        elif interaction_type == 'offset':
            return angle <= 45  # Allow more flexibility for offset interactions
        
        return False
    
    def _calculate_pi_pi_strength(self, distance: float, angle: float, offset: float, interaction_type: str) -> str:
        """Calculate π-π interaction strength."""
        # Base strength on geometry
        distance_score = max(0, (6.0 - distance) / 6.0)
        
        if interaction_type == 'face_to_face':
            angle_score = max(0, (30 - angle) / 30)
            offset_score = max(0, (3.0 - offset) / 3.0)
            overall_score = (distance_score + angle_score + offset_score) / 3
        
        elif interaction_type == 'edge_to_face':
            angle_score = max(0, (angle - 60) / 30)
            overall_score = (distance_score + angle_score) / 2
        
        else:  # offset
            overall_score = distance_score * 0.7  # Generally weaker
        
        if overall_score > 0.7:
            return 'strong'
        elif overall_score > 0.4:
            return 'moderate'
        else:
            return 'weak'
    
    def _point_to_plane_distance(self, point: np.ndarray, plane_coeffs: np.ndarray) -> float:
        """Calculate distance from point to plane."""
        a, b, c, d = plane_coeffs
        return abs(a * point[0] + b * point[1] + c * point[2] + d) / np.sqrt(a**2 + b**2 + c**2)
    
    def get_statistics(self, interactions: List[PiPiInteraction]) -> Dict[str, Any]:
        """Get statistics about detected π-π interactions."""
        if not interactions:
            return {}
        
        stats = {
            'total_interactions': len(interactions),
            'interaction_types': {},
            'interaction_strengths': {},
            'residue_pairs': {},
            'distance_stats': {
                'mean': np.mean([i.distance for i in interactions]),
                'std': np.std([i.distance for i in interactions]),
                'min': min([i.distance for i in interactions]),
                'max': max([i.distance for i in interactions])
            },
            'angle_stats': {
                'mean': np.mean([i.angle for i in interactions]),
                'std': np.std([i.angle for i in interactions]),
                'min': min([i.angle for i in interactions]),
                'max': max([i.angle for i in interactions])
            },
            'offset_stats': {
                'mean': np.mean([i.offset for i in interactions]),
                'std': np.std([i.offset for i in interactions]),
                'min': min([i.offset for i in interactions]),
                'max': max([i.offset for i in interactions])
            }
        }
        
        # Count by type and strength
        for interaction in interactions:
            stats['interaction_types'][interaction.interaction_type] = \
                stats['interaction_types'].get(interaction.interaction_type, 0) + 1
            
            stats['interaction_strengths'][interaction.strength] = \
                stats['interaction_strengths'].get(interaction.strength, 0) + 1
            
            # Count residue pair types
            res1_type = interaction.ring1_residue[:3]
            res2_type = interaction.ring2_residue[:3]
            pair_key = f"{res1_type}-{res2_type}"
            stats['residue_pairs'][pair_key] = stats['residue_pairs'].get(pair_key, 0) + 1
        
        return stats
    
    def to_dict_list(self, interactions: List[PiPiInteraction]) -> List[Dict[str, Any]]:
        """Convert π-π interactions to list of dictionaries."""
        return [
            {
                'ring1_residue': interaction.ring1_residue,
                'ring2_residue': interaction.ring2_residue,
                'ring1_chain': interaction.ring1_chain,
                'ring2_chain': interaction.ring2_chain,
                'distance': round(interaction.distance, 3),
                'angle': round(interaction.angle, 1),
                'offset': round(interaction.offset, 3),
                'interaction_type': interaction.interaction_type,
                'strength': interaction.strength,
                'residue1': interaction.ring1_residue,
                'residue2': interaction.ring2_residue,
                'chain1': interaction.ring1_chain,
                'chain2': interaction.ring2_chain
            }
            for interaction in interactions
        ]
