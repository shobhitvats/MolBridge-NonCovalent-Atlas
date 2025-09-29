"""
C-H···π interaction detection for protein structures.
Detects interactions between C-H bonds and π-systems.
"""

from typing import List, Dict, Any, Optional, Set, Tuple
import time
import numpy as np
from Bio.PDB import Structure, Residue, Atom
from utils.settings import get_settings
from utils.instrumentation import init_funnel, update_counts, finalize_funnel
from analysis.feature_store import get_feature_store
from geometry.core import pairwise_within_cutoff, norms, vector_fields

from utils.config import AppConfig

from .base_detector import register_detector

@register_detector("chpi", method="detect_ch_pi_interactions")
class CHPiDetector:
    """Detects C-H···π interactions in protein structures."""
    
    def __init__(self, config: AppConfig):
        self.config = config
        self.interaction_config = config.interactions
        
        # Aromatic residues with ring atoms
        self.aromatic_rings = {
            'PHE': [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
            'TYR': [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
            'TRP': [
                ['CG', 'CD1', 'NE1', 'CE2', 'CD2'],  # pyrrole ring
                ['CD2', 'CE2', 'CE3', 'CZ3', 'CH2', 'CZ2']  # benzene ring
            ],
            'HIS': [['CG', 'ND1', 'CE1', 'NE2', 'CD2']],
            'HIP': [['CG', 'ND1', 'CE1', 'NE2', 'CD2']],
            'HIE': [['CG', 'ND1', 'CE1', 'NE2', 'CD2']],
            'HID': [['CG', 'ND1', 'CE1', 'NE2', 'CD2']]
        }
        
        # C-H donor atoms and their hydrogen atoms
        self.ch_donors = {
            'ALA': [('CB', ['HB1', 'HB2', 'HB3'])],
            'VAL': [('CB', ['HB']), ('CG1', ['HG11', 'HG12', 'HG13']), ('CG2', ['HG21', 'HG22', 'HG23'])],
            'LEU': [('CB', ['HB2', 'HB3']), ('CG', ['HG']), ('CD1', ['HD11', 'HD12', 'HD13']), ('CD2', ['HD21', 'HD22', 'HD23'])],
            'ILE': [('CB', ['HB']), ('CG1', ['HG12', 'HG13']), ('CG2', ['HG21', 'HG22', 'HG23']), ('CD1', ['HD11', 'HD12', 'HD13'])],
            'THR': [('CB', ['HB']), ('CG2', ['HG21', 'HG22', 'HG23'])],
            'SER': [('CB', ['HB2', 'HB3'])],
            'CYS': [('CB', ['HB2', 'HB3'])],
            'MET': [('CB', ['HB2', 'HB3']), ('CG', ['HG2', 'HG3']), ('CE', ['HE1', 'HE2', 'HE3'])],
            'ASN': [('CB', ['HB2', 'HB3'])],
            'GLN': [('CB', ['HB2', 'HB3']), ('CG', ['HG2', 'HG3'])],
            'ASP': [('CB', ['HB2', 'HB3'])],
            'GLU': [('CB', ['HB2', 'HB3']), ('CG', ['HG2', 'HG3'])],
            'LYS': [('CB', ['HB2', 'HB3']), ('CG', ['HG2', 'HG3']), ('CD', ['HD2', 'HD3']), ('CE', ['HE2', 'HE3'])],
            'ARG': [('CB', ['HB2', 'HB3']), ('CG', ['HG2', 'HG3']), ('CD', ['HD2', 'HD3'])],
            'PRO': [('CB', ['HB2', 'HB3']), ('CG', ['HG2', 'HG3']), ('CD', ['HD2', 'HD3'])],
            'GLY': [('CA', ['HA2', 'HA3'])],
            # Aromatic C-H donors
            'PHE': [('CD1', ['HD1']), ('CE1', ['HE1']), ('CZ', ['HZ']), ('CE2', ['HE2']), ('CD2', ['HD2'])],
            'TYR': [('CD1', ['HD1']), ('CE1', ['HE1']), ('CE2', ['HE2']), ('CD2', ['HD2'])],
            'TRP': [('CD1', ['HD1']), ('CE1', ['HE1']), ('CE3', ['HE3']), ('CZ3', ['HZ3']), ('CH2', ['HH2']), ('CZ2', ['HZ2'])],
            'HIS': [('CD2', ['HD2']), ('CE1', ['HE1'])],
            'HIP': [('CD2', ['HD2']), ('CE1', ['HE1'])],
            'HIE': [('CD2', ['HD2']), ('CE1', ['HE1'])],
            'HID': [('CE1', ['HE1'])]
        }
    
    def detect_ch_pi_interactions(self, structure: Structure) -> List[Dict[str, Any]]:
        """
        Detect C-H···π interactions in the structure.
        
        Args:
            structure: BioPython Structure object
            
        Returns:
            List of detected interactions with details
        """
        settings = get_settings()
        if getattr(settings, 'enable_vector_geom', False):
            try:
                return self._vector_detect(structure)
            except Exception:  # pragma: no cover
                pass
        interactions: List[Dict[str, Any]] = []
        raw_pairs = 0
        candidate_pairs = 0
        t_pair_start = time.time()

        # Get all residues
        residues = list(structure.get_residues())

        # Find aromatic residues (π-acceptors)
        aromatic_residues = []
        for residue in residues:
            if residue.resname in self.aromatic_rings:
                aromatic_residues.append(residue)
        
        # Find C-H donor residues
        ch_donor_residues = []
        for residue in residues:
            if residue.resname in self.ch_donors:
                ch_donor_residues.append(residue)
        
        # Check interactions between C-H donors and π-acceptors
        raw_pairs = len(ch_donor_residues) * len(aromatic_residues)
        for donor_res in ch_donor_residues:
            for acceptor_res in aromatic_residues:
                # Skip same residue
                if donor_res == acceptor_res:
                    continue
                
                # Skip very close residues in sequence (less than 3 residues apart)
                if (donor_res.parent == acceptor_res.parent and
                    abs(donor_res.id[1] - acceptor_res.id[1]) < 3):
                    continue
                
                # Check each aromatic ring in acceptor
                for ring_atoms in self.aromatic_rings[acceptor_res.resname]:
                    ring_center, ring_normal = self._get_ring_center_and_normal(acceptor_res, ring_atoms)
                    if ring_center is None:
                        continue
                    
                    # Check each C-H donor in donor residue
                    for carbon_name, hydrogen_names in self.ch_donors[donor_res.resname]:
                        ch_interactions = self._check_ch_ring_interaction(
                            donor_res, carbon_name, hydrogen_names,
                            acceptor_res, ring_center, ring_normal, ring_atoms
                        )
                        prev_len = len(interactions)
                        interactions.extend(ch_interactions)
                        if ch_interactions and len(interactions) > prev_len:
                            candidate_pairs += 1  # each productive ring counts once

        t_eval_end = time.time()
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=len(interactions),
            core_pair_generation=False,
            extra={
                'donor_residues': len(ch_donor_residues),
                'aromatic_residues': len(aromatic_residues),
                'candidate_density': (candidate_pairs/raw_pairs) if raw_pairs else 0.0
            }
        )
        finalize_funnel(
            self.instrumentation,
            pair_gen_seconds=0.0,
            eval_seconds=(t_eval_end - t_pair_start),
            build_seconds=0.0
        )
        from loguru import logger as _l
        _l.info(f"CH-π: {len(interactions)} (raw={raw_pairs} cand_units={candidate_pairs} acc_ratio={self.instrumentation['acceptance_ratio']:.3f})")
        return interactions

    # ---- Vector fast path leveraging ring FeatureStore and pair pruning ----
    def _vector_detect(self, structure: Structure) -> List[Dict[str, Any]]:
        interactions: List[Dict[str, Any]] = []
        t_pair_start = time.time()
        try:
            model = structure[0]
        except Exception:
            return interactions
        fs = get_feature_store()
        rings = fs.ensure_rings(structure)
        if not rings:
            return interactions
        ring_centers = np.vstack([r['centroid'] for r in rings]).astype('float32')
        ring_normals = np.vstack([r['normal'] for r in rings]).astype('float32')
        ring_meta = [(r['residue'].parent.id, r['residue']) for r in rings]
        donor_entries = []  # (residue, carbon_atom, hydrogen_atom or inferred_coord, hydrogen_name, is_inferred)
        # Collect all donor heavy atoms & hydrogens (explicit and infer where needed)
        for chain in model:
            for residue in chain:
                resname = residue.resname
                if resname not in self.ch_donors:
                    continue
                for carbon_name, hydrogen_names in self.ch_donors[resname]:
                    if carbon_name not in residue:
                        continue
                    carbon_atom = residue[carbon_name]
                    for h_name in hydrogen_names:
                        if h_name in residue:
                            donor_entries.append((residue, carbon_atom, residue[h_name], h_name, False))
                        else:
                            h_coord = self._infer_hydrogen_position(residue, carbon_atom, h_name)
                            if h_coord is not None:
                                # lightweight proxy storing coord
                                donor_entries.append((residue, carbon_atom, type('HProxy', (), {'name': h_name, 'coord': h_coord, 'element': 'H'})(), h_name, True))
        if not donor_entries:
            return interactions
        # Build coordinate arrays for hydrogen positions to prune vs ring centroids
        h_coords = np.vstack([d[2].coord for d in donor_entries]).astype('float32')
        # Broad pruning with max distance
        max_d = float(self.interaction_config.ch_pi_max_distance)
        try:
            hi, ri, dists = pairwise_within_cutoff(h_coords, ring_centers, max_d, use_kdtree=True, return_distances=True)
            core = True
        except Exception:  # pragma: no cover
            diff = h_coords[:, None, :] - ring_centers[None, :, :]
            dist2 = np.sum(diff*diff, axis=-1)
            mask = dist2 <= (max_d**2)
            hi, ri = np.where(mask)
            hi = hi.astype('int32'); ri = ri.astype('int32')
            dists = np.sqrt(dist2[hi, ri]).astype('float32')
            core = False
        raw_pairs = int(h_coords.shape[0] * ring_centers.shape[0])
        candidate_pairs = int(len(hi))
        t_pair_end = time.time()
        t_eval_start = t_pair_end
        # Precompute distances for candidates for early rejects on min distance and height later
    # dists already computed when return_distances True
        accepted_units = 0
        # Map hydrogen index to (residue, carbon_atom)
        for idx in range(candidate_pairs):
            h_index = int(hi[idx]); r_index = int(ri[idx])
            dist = float(dists[idx])
            if dist < self.interaction_config.ch_pi_min_distance or dist > self.interaction_config.ch_pi_max_distance:
                continue
            residue, carbon_atom, hydrogen_atom, h_name, inferred = donor_entries[h_index]
            ring_chain, ring_residue = ring_meta[r_index]
            # Skip same residue or near sequence proximity (<3 apart same chain)
            if residue is ring_residue:
                continue
            if (residue.parent == ring_residue.parent and abs(residue.id[1] - ring_residue.id[1]) < 3):
                continue
            ring_center = ring_centers[r_index]
            ring_normal = ring_normals[r_index]
            # Geometry evaluation
            ch_vector = hydrogen_atom.coord - carbon_atom.coord
            h_ring_vector = ring_center - hydrogen_atom.coord
            ch_norm = np.linalg.norm(ch_vector)
            hr_norm = np.linalg.norm(h_ring_vector)
            if ch_norm == 0 or hr_norm == 0:
                continue
            ch_unit = ch_vector / ch_norm
            hr_unit = h_ring_vector / hr_norm
            cos_angle = np.clip(np.dot(ch_unit, hr_unit), -1.0, 1.0)
            ch_ring_angle = float(np.degrees(np.arccos(cos_angle)))
            if ch_ring_angle > self.interaction_config.ch_pi_max_angle:
                continue
            # Height above plane
            height = abs(float(np.dot(hydrogen_atom.coord - ring_center, ring_normal)))
            if height > self.interaction_config.ch_pi_max_height:
                continue
            strength = self._calculate_ch_pi_strength(dist, ch_ring_angle, height)
            interaction_type = self._classify_ch_pi_type(residue, carbon_atom, ring_residue)
            interactions.append({
                'type': 'ch_pi',
                'subtype': interaction_type,
                'residue1': f"{residue.resname}{residue.id[1]}",
                'chain1': residue.parent.id,
                'residue2': f"{ring_residue.resname}{ring_residue.id[1]}",
                'chain2': ring_residue.parent.id,
                'atom1': carbon_atom.name,
                'atom2': 'RING_CENTER',
                'hydrogen': hydrogen_atom.name,
                'distance': dist,
                'angle': ch_ring_angle,
                'height': height,
                'strength': strength,
                'donor_coord': carbon_atom.coord.tolist(),
                'acceptor_coord': ring_center.tolist(),
                'hydrogen_coord': hydrogen_atom.coord.tolist(),
                'ring_atoms': [a for a in []],  # omit per-ring atoms list for speed (can be added)
            })
            accepted_units += 1
        t_eval_end = time.time()
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=accepted_units,
            core_pair_generation=core,
            extra={
                'hydrogens': h_coords.shape[0],
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
        from loguru import logger as _l
        _l.info(f"CH-π[vector]{'/core' if core else ''}: {accepted_units} raw={raw_pairs} pruned={candidate_pairs} acc={self.instrumentation['acceptance_ratio']:.3f}")
        return interactions
    
    def _get_ring_center_and_normal(self, residue: Residue, ring_atoms: List[str]) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """Calculate ring center and normal vector."""
        try:
            coords = []
            for atom_name in ring_atoms:
                if atom_name in residue:
                    coords.append(residue[atom_name].coord)
                else:
                    return None, None
            
            coords = np.array(coords)
            center = np.mean(coords, axis=0)
            
            # Calculate normal using first three atoms
            if len(coords) >= 3:
                v1 = coords[1] - coords[0]
                v2 = coords[2] - coords[0]
                normal = np.cross(v1, v2)
                normal = normal / np.linalg.norm(normal)
                return center, normal
            
            return center, None
            
        except Exception:
            return None, None
    
    def _check_ch_ring_interaction(self,
                                 donor_res: Residue,
                                 carbon_name: str,
                                 hydrogen_names: List[str],
                                 acceptor_res: Residue,
                                 ring_center: np.ndarray,
                                 ring_normal: Optional[np.ndarray],
                                 ring_atoms: List[str]) -> List[Dict[str, Any]]:
        """Check C-H···π interaction between specific C-H and ring."""
        interactions = []
        
        # Get carbon atom
        if carbon_name not in donor_res:
            return interactions
        
        carbon_atom = donor_res[carbon_name]
        carbon_coord = carbon_atom.coord
        
        # Distance check first (quick filter)
        ch_ring_distance = np.linalg.norm(carbon_coord - ring_center)
        if ch_ring_distance > self.interaction_config.ch_pi_max_distance:
            return interactions
        
        # Check with each hydrogen (if available) or infer H position
        for h_name in hydrogen_names:
            if h_name in donor_res:
                # Explicit hydrogen
                h_coord = donor_res[h_name].coord
                interaction = self._evaluate_ch_pi_geometry(
                    donor_res, carbon_atom, donor_res[h_name],
                    acceptor_res, ring_center, ring_normal, ring_atoms
                )
                if interaction:
                    interactions.append(interaction)
            else:
                # Infer hydrogen position
                h_coord = self._infer_hydrogen_position(donor_res, carbon_atom, h_name)
                if h_coord is not None:
                    # Create temporary hydrogen atom for evaluation
                    temp_h_atom = type('TempAtom', (), {
                        'name': h_name,
                        'coord': h_coord,
                        'element': 'H'
                    })()
                    
                    interaction = self._evaluate_ch_pi_geometry(
                        donor_res, carbon_atom, temp_h_atom,
                        acceptor_res, ring_center, ring_normal, ring_atoms
                    )
                    if interaction:
                        interactions.append(interaction)
        
        return interactions
    
    def _evaluate_ch_pi_geometry(self,
                               donor_res: Residue,
                               carbon_atom: Atom,
                               hydrogen_atom: Atom,
                               acceptor_res: Residue,
                               ring_center: np.ndarray,
                               ring_normal: Optional[np.ndarray],
                               ring_atoms: List[str]) -> Optional[Dict[str, Any]]:
        """Evaluate C-H···π interaction geometry."""
        
        c_coord = carbon_atom.coord
        h_coord = hydrogen_atom.coord
        
        # Calculate distances
        h_ring_distance = np.linalg.norm(h_coord - ring_center)
        c_ring_distance = np.linalg.norm(c_coord - ring_center)
        
        # Distance criteria
        if (h_ring_distance > self.interaction_config.ch_pi_max_distance or
            h_ring_distance < self.interaction_config.ch_pi_min_distance):
            return None
        
        # Angle criteria: C-H···π
        ch_vector = h_coord - c_coord
        h_ring_vector = ring_center - h_coord
        
        # C-H···ring angle
        try:
            ch_norm = ch_vector / np.linalg.norm(ch_vector)
            hr_norm = h_ring_vector / np.linalg.norm(h_ring_vector)
            
            cos_angle = np.dot(ch_norm, hr_norm)
            cos_angle = np.clip(cos_angle, -1.0, 1.0)
            ch_ring_angle = np.degrees(np.arccos(cos_angle))
            
            if ch_ring_angle > self.interaction_config.ch_pi_max_angle:
                return None
            
        except (ZeroDivisionError, ValueError):
            return None
        
        # Height above ring plane (if normal available)
        height = None
        if ring_normal is not None:
            ring_to_h = h_coord - ring_center
            height = abs(np.dot(ring_to_h, ring_normal))
            
            if height > self.interaction_config.ch_pi_max_height:
                return None
        
        # Determine interaction strength
        strength = self._calculate_ch_pi_strength(h_ring_distance, ch_ring_angle, height)
        
        # Classify interaction type
        interaction_type = self._classify_ch_pi_type(donor_res, carbon_atom, acceptor_res)
        
        return {
            'type': 'ch_pi',
            'subtype': interaction_type,
            'residue1': f"{donor_res.resname}{donor_res.id[1]}",
            'chain1': donor_res.parent.id,
            'residue2': f"{acceptor_res.resname}{acceptor_res.id[1]}",
            'chain2': acceptor_res.parent.id,
            'atom1': carbon_atom.name,
            'atom2': 'RING_CENTER',
            'hydrogen': hydrogen_atom.name,
            'distance': h_ring_distance,
            'angle': ch_ring_angle,
            'height': height,
            'strength': strength,
            'donor_coord': c_coord.tolist(),
            'acceptor_coord': ring_center.tolist(),
            'hydrogen_coord': h_coord.tolist(),
            'ring_atoms': ring_atoms
        }
    
    def _infer_hydrogen_position(self, residue: Residue, carbon_atom: Atom, h_name: str) -> Optional[np.ndarray]:
        """Infer hydrogen position based on carbon geometry."""
        try:
            c_coord = carbon_atom.coord
            
            # Get neighboring atoms to carbon
            neighbors = []
            for atom in residue:
                if atom != carbon_atom and atom.element != 'H':
                    distance = np.linalg.norm(atom.coord - c_coord)
                    if distance < 2.0:  # Typical C-X bond length
                        neighbors.append(atom.coord)
            
            if len(neighbors) == 0:
                return None
            
            # Simple geometric hydrogen placement
            if len(neighbors) == 1:
                # Linear case: place H opposite to neighbor
                direction = c_coord - neighbors[0]
                direction = direction / np.linalg.norm(direction)
                h_coord = c_coord + direction * 1.09  # C-H bond length
                
            elif len(neighbors) == 2:
                # Trigonal case: place H in plane opposite to neighbors
                v1 = neighbors[0] - c_coord
                v2 = neighbors[1] - c_coord
                
                # Average direction opposite to neighbors
                avg_direction = -(v1 + v2)
                avg_direction = avg_direction / np.linalg.norm(avg_direction)
                h_coord = c_coord + avg_direction * 1.09
                
            elif len(neighbors) >= 3:
                # Tetrahedral case: place H opposite to average of neighbors
                avg_neighbor = np.mean(neighbors, axis=0)
                direction = c_coord - avg_neighbor
                direction = direction / np.linalg.norm(direction)
                h_coord = c_coord + direction * 1.09
                
            else:
                return None
            
            return h_coord
            
        except Exception:
            return None
    
    def _calculate_ch_pi_strength(self, distance: float, angle: float, height: Optional[float]) -> str:
        """Calculate C-H···π interaction strength."""
        
        # Distance component
        if distance <= 3.0:
            distance_score = 3
        elif distance <= 3.5:
            distance_score = 2
        else:
            distance_score = 1
        
        # Angle component
        if angle <= 30:
            angle_score = 3
        elif angle <= 45:
            angle_score = 2
        else:
            angle_score = 1
        
        # Height component
        height_score = 2
        if height is not None:
            if height <= 1.0:
                height_score = 3
            elif height <= 2.0:
                height_score = 2
            else:
                height_score = 1
        
        total_score = distance_score + angle_score + height_score
        
        if total_score >= 8:
            return "strong"
        elif total_score >= 6:
            return "moderate"
        else:
            return "weak"
    
    def _classify_ch_pi_type(self, donor_res: Residue, carbon_atom: Atom, acceptor_res: Residue) -> str:
        """Classify C-H···π interaction type."""
        
        # Check if donor is aromatic
        if donor_res.resname in self.aromatic_rings:
            # Check if carbon is part of aromatic ring
            for ring_atoms in self.aromatic_rings[donor_res.resname]:
                if carbon_atom.name in ring_atoms:
                    return "aromatic_ch_pi"
            return "aromatic_aliphatic_ch_pi"
        else:
            return "aliphatic_ch_pi"
