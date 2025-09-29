"""
Chalcogen bond detection for protein structures.
Detects S/Se/Te σ-hole donor interactions with O/N/S acceptors using directionality angles θ (115–155°) and δ (±50°), distance ≤ 3.6 Å.
"""

import numpy as np
import time
import math
from typing import List, Dict, Any, Optional, Tuple
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from loguru import logger
from .base_detector import register_detector
from utils.settings import get_settings
from utils.instrumentation import init_funnel, update_counts, finalize_funnel
from geometry.core import pairwise_within_cutoff, norms

# Set debug level for detailed angle calculations
logger.add("chalcogen_debug.log", level="DEBUG", filter=lambda record: "chalcogen" in record["name"].lower())

@dataclass
class ChalcogenBond:
    """Represents a detected chalcogen bond."""
    chalcogen_atom: Atom.Atom
    acceptor_atom: Atom.Atom
    distance: float
    theta_angle: float  # θ angle (115–155°)
    delta_angle: float  # δ angle (±50°) 
    strength: str
    chalcogen_residue: str
    acceptor_residue: str
    chalcogen_chain: str
    acceptor_chain: str

@register_detector("chalcogen_bond", method="detect_chalcogen_bonds")
class ChalcogenBondDetector:
    """Detects chalcogen bonds in protein structures."""
    
    def __init__(self, config):
        """
        Initialize chalcogen bond detector.
        
        Args:
            config: AppConfig object containing interaction parameters
        """
        self.config = config
        self.interaction_config = config.interactions
        # Pull thresholds from unified interaction configuration (with backwards compatibility)
        self.distance_cutoff = getattr(self.interaction_config, 'chalcogen_distance_cutoff', 3.6)
        # Theta range now user-configurable
        self.theta_min = getattr(self.interaction_config, 'chalcogen_theta_min', 115.0)
        self.theta_max = getattr(self.interaction_config, 'chalcogen_theta_max', 155.0)
        # Phi (delta) absolute max deviation
        self.delta_max = getattr(self.interaction_config, 'chalcogen_phi_max', 50.0)
        
        # Chalcogen atoms with σ-holes (group 16 elements)
        self.chalcogen_donors = {
            'MET': ['SD'],   # Sulfur in methionine
            'CYS': ['SG']    # Sulfur in cysteine
            # Se and Te could be added for selenomethionine, etc.
        }
        
        # Nucleophilic acceptor atoms (electron-rich)
        self.acceptor_atoms = {
            'ASP': ['OD1', 'OD2'],    # Aspartate oxygens
            'GLU': ['OE1', 'OE2'],    # Glutamate oxygens
            'SER': ['OG'],            # Serine oxygen
            'THR': ['OG1'],           # Threonine oxygen
            'TYR': ['OH'],            # Tyrosine oxygen
            'ASN': ['OD1', 'ND2'],    # Asparagine O/N
            'GLN': ['OE1', 'NE2'],    # Glutamine O/N
            'HIS': ['ND1', 'NE2'],    # Histidine nitrogens
            'ARG': ['NE', 'NH1', 'NH2'],  # Arginine nitrogens
            'LYS': ['NZ'],            # Lysine nitrogen
            'TRP': ['NE1'],           # Tryptophan nitrogen
            'MET': ['SD'],            # Methionine sulfur
            'CYS': ['SG'],            # Cysteine sulfur
            'backbone': ['O', 'N']    # Backbone atoms
        }
    
    def detect_chalcogen_bonds(self, structure: Structure.Structure) -> List[ChalcogenBond]:
        """
        Detect chalcogen bonds in structure using S/Se/Te σ-hole donors.
        Based on the reference implementation with proper theta and delta calculations.
        """
        settings = get_settings()
        if getattr(settings, 'enable_vector_geom', False):
            try:
                return self._vector_detect(structure)
            except Exception:  # pragma: no cover
                pass  # fall through to legacy
        t_pair_start = time.time()
        bonds: List[ChalcogenBond] = []
        try:
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            if atom.element != 'S':
                                continue
                            sulfur = atom
                            bonded_atoms = [nbr for nbr in residue if nbr.element == 'C' and self._calculate_distance(sulfur, nbr) < 1.9]
                            if len(bonded_atoms) < 2:
                                continue
                            atoms_for_centroid = bonded_atoms + [sulfur]
                            x_atom = bonded_atoms[1]
                            for other_model in structure:
                                for other_chain in other_model:
                                    for other_residue in other_chain:
                                        if residue == other_residue:
                                            continue
                                        for other_atom in other_residue:
                                            if other_atom.element not in ['O','N']:
                                                continue
                                            distance = self._calculate_distance(sulfur, other_atom)
                                            if distance > self.distance_cutoff:
                                                continue
                                            try:
                                                theta = self._calculate_theta_ref(sulfur, atoms_for_centroid, other_atom)
                                                delta = self._calculate_delta_ref(x_atom, atoms_for_centroid, sulfur, other_atom)
                                            except Exception as e:
                                                logger.warning(f"Error calculating chalcogen angles: {e}; fallback defaults")
                                                theta, delta = 135.0, 0.0
                                            if (self.theta_min <= theta <= self.theta_max) and (abs(delta) <= self.delta_max):
                                                bonds.append(ChalcogenBond(
                                                    chalcogen_atom=sulfur,
                                                    acceptor_atom=other_atom,
                                                    distance=distance,
                                                    theta_angle=theta,
                                                    delta_angle=delta,
                                                    strength=self._classify_bond_strength(distance, theta),
                                                    chalcogen_residue=f"{residue.get_resname()}{residue.get_id()[1]}",
                                                    acceptor_residue=f"{other_residue.get_resname()}{other_residue.get_id()[1]}",
                                                    chalcogen_chain=chain.get_id(),
                                                    acceptor_chain=other_chain.get_id()
                                                ))
            t_eval_end = time.time()
            logger.info(f"Detected {len(bonds)} chalcogen bonds")
            # Legacy instrumentation via helper (single combined phase)
            donors_ct = 0
            acceptors_ct = 0
            try:
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                if atom.element == 'S':
                                    bonded = [nbr for nbr in residue if nbr.element == 'C' and self._calculate_distance(atom, nbr) < 1.9]
                                    if len(bonded) >= 2:
                                        donors_ct += 1
                                if atom.element in ['O', 'N']:
                                    acceptors_ct += 1
            except Exception:
                pass
            raw_pairs = donors_ct * acceptors_ct
            self.instrumentation = init_funnel(
                raw_pairs=raw_pairs,
                candidate_pairs=raw_pairs,
                accepted_pairs=len(bonds),
                core_pair_generation=False,
                extra={
                    'donors': donors_ct,
                    'acceptors': acceptors_ct,
                    'kdtree_used': False
                }
            )
            finalize_funnel(
                self.instrumentation,
                pair_gen_seconds=0.0,
                eval_seconds=(t_eval_end - t_pair_start),
                build_seconds=0.0
            )
        except Exception as e:
            logger.error(f"Error detecting chalcogen bonds: {e}")
        return bonds
    
    # ---- Vector candidate pruning (distance matrix) then angle refinement ----
    def _vector_detect(self, structure: Structure.Structure) -> List[ChalcogenBond]:
        bonds: List[ChalcogenBond] = []
        t_pair_start = time.time()
        # Gracefully handle empty structure
        try:
            model = structure[0]
        except Exception:
            return bonds

        # 1. Collect sulfur donors (S with two carbon neighbors <1.9 Å)
        donor_atoms = []
        donor_bonded_sets = []  # per-donor list of [C1, C2, S]
        donor_centroids = []
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == 'S':
                        bonded = [nbr for nbr in residue if nbr.element == 'C' and self._calculate_distance(atom, nbr) < 1.9]
                        if len(bonded) < 2:
                            continue
                        donor_atoms.append(atom)
                        donor_bonded_sets.append(bonded + [atom])
                        donor_centroids.append(self._calculate_centroid_ref(bonded + [atom]))
        if not donor_atoms:
            return bonds

        # 2. Collect acceptor atoms (O/N) from whole model
        acceptor_atoms = []
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element in ['O', 'N']:
                        acceptor_atoms.append(atom)
        if not acceptor_atoms:
            return bonds

        donor_coords = np.vstack([a.get_coord() for a in donor_atoms]).astype('float32')
        donor_centroids_arr = np.vstack(donor_centroids).astype('float32')
        # centroid -> sulfur vector (as in reference: centroid - S)
        vec_sc = donor_centroids_arr - donor_coords
        norms_sc = np.linalg.norm(vec_sc, axis=1)
        norms_sc[norms_sc == 0] = 1.0
        unit_sc = vec_sc / norms_sc[:, None]
        acceptor_coords = np.vstack([a.get_coord() for a in acceptor_atoms]).astype('float32')

        # 3. Pair pruning (KD-tree if available)
        try:
            donor_idx, acc_idx = pairwise_within_cutoff(donor_coords, acceptor_coords, float(self.distance_cutoff), use_kdtree=True)
            core = True
        except Exception:  # pragma: no cover
            diff = donor_coords[:, None, :] - acceptor_coords[None, :, :]
            dist2 = np.sum(diff * diff, axis=-1)
            mask = dist2 <= (self.distance_cutoff ** 2)
            donor_idx, acc_idx = np.where(mask)
            donor_idx = donor_idx.astype('int32'); acc_idx = acc_idx.astype('int32')
            core = False

        # 4. Vector theta calculation
        s_to_a_vec = acceptor_coords[acc_idx] - donor_coords[donor_idx]
        sa_norms = np.linalg.norm(s_to_a_vec, axis=1)
        sa_norms[sa_norms == 0] = 1.0
        unit_sa = s_to_a_vec / sa_norms[:, None]
        theta_vals = np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i', unit_sc[donor_idx], unit_sa), -1.0, 1.0)))

        # 5. Vector delta (φ) calculation
        ref_vecs = []
        plane_norms = []
        for di in donor_idx.tolist():
            bonded_set = donor_bonded_sets[di]  # [C1, C2, S]
            if len(bonded_set) < 3:
                c1 = bonded_set[0].get_coord()
                s_coord = donor_atoms[di].get_coord()
                ref_vecs.append(c1 - s_coord)
                plane_norms.append(np.array([0.0, 0.0, 1.0]))
                continue
            c1, c2, s_atom_local = bonded_set[0], bonded_set[1], bonded_set[-1]
            s_coord = s_atom_local.get_coord()
            v1 = c1.get_coord() - s_coord
            v2 = c2.get_coord() - s_coord
            ref_vecs.append(v1)
            plane_norms.append(np.cross(v1, v2))
        ref_vecs_arr = np.vstack(ref_vecs).astype('float32') if ref_vecs else np.zeros((0, 3), dtype='float32')
        plane_norms_arr = np.vstack(plane_norms).astype('float32') if plane_norms else np.zeros((0, 3), dtype='float32')
        ref_norm = np.linalg.norm(ref_vecs_arr, axis=1); ref_norm[ref_norm == 0] = 1.0
        ref_unit = ref_vecs_arr / ref_norm[:, None]
        plane_norm_mag = np.linalg.norm(plane_norms_arr, axis=1); plane_norm_mag[plane_norm_mag == 0] = 1.0
        plane_unit = plane_norms_arr / plane_norm_mag[:, None]
        sa_vec = s_to_a_vec
        ref_unit_exp = ref_unit[donor_idx]
        plane_unit_exp = plane_unit[donor_idx]
        dot_vn = np.einsum('ij,ij->i', sa_vec, plane_unit_exp)
        v_proj = sa_vec - dot_vn[:, None] * plane_unit_exp
        v_proj_norm = np.linalg.norm(v_proj, axis=1); v_proj_norm[v_proj_norm == 0] = 1.0
        v_proj_unit = v_proj / v_proj_norm[:, None]
        cross_r_v = np.cross(ref_unit_exp, v_proj_unit)
        dot_r_v = np.einsum('ij,ij->i', ref_unit_exp, v_proj_unit)
        dot_r_v = np.clip(dot_r_v, -1.0, 1.0)
        unsigned = np.degrees(np.arccos(dot_r_v))
        sign = np.sign(np.einsum('ij,ij->i', plane_unit_exp, cross_r_v))
        delta_vals = unsigned * sign

        # 6. Instrumentation (before evaluation loop) using shared helpers
        t_pair_end = time.time()
        t_eval_start = t_pair_end
        raw_pairs = int(donor_coords.shape[0] * acceptor_coords.shape[0])
        candidate_pairs = int(len(donor_idx))
        base_threshold = None
        try:
            from utils.kdtree_thresholds import get_threshold, adapt_threshold, get_last_density, should_flag_kdtree
            base_threshold = get_threshold('chalcogen')
            kdtree_used = should_flag_kdtree(candidate_pairs, raw_pairs, 'chalcogen')
            adaptive_threshold, changed, reason = adapt_threshold('chalcogen', candidate_pairs, kdtree_used)
            density_val = get_last_density('chalcogen')
        except Exception:  # pragma: no cover
            adaptive_threshold = None; changed = False; reason = 'unavailable'; density_val = None; kdtree_used = core
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=0,
            core_pair_generation=core,
            extra={
                'donors': donor_coords.shape[0],
                'acceptors': acceptor_coords.shape[0],
                'candidate_density': (candidate_pairs / raw_pairs) if raw_pairs else 0.0,
                'theta_vectorized_pairs': candidate_pairs,
                'delta_vectorized_pairs': candidate_pairs,
                'threshold': base_threshold,
                'adaptive_threshold': adaptive_threshold,
                'threshold_changed': changed,
                'adapt_reason': reason,
                'last_density': density_val,
                'kdtree_used': kdtree_used
            }
        )

        # 7. Filter & build objects
        dist_vec = norms(acceptor_coords[acc_idx] - donor_coords[donor_idx])
        for i in range(candidate_pairs):
            di = int(donor_idx[i]); ai = int(acc_idx[i])
            s_atom = donor_atoms[di]; a_atom = acceptor_atoms[ai]
            if s_atom.get_parent() is a_atom.get_parent():
                continue
            distance = float(dist_vec[i])
            theta = float(theta_vals[i])
            delta = float(delta_vals[i])
            if (self.theta_min <= theta <= self.theta_max) and (abs(delta) <= self.delta_max):
                bonds.append(ChalcogenBond(
                    chalcogen_atom=s_atom,
                    acceptor_atom=a_atom,
                    distance=distance,
                    theta_angle=theta,
                    delta_angle=delta,
                    strength=self._classify_bond_strength(distance, theta),
                    chalcogen_residue=f"{s_atom.get_parent().get_resname()}{s_atom.get_parent().id[1]}",
                    acceptor_residue=f"{a_atom.get_parent().get_resname()}{a_atom.get_parent().id[1]}",
                    chalcogen_chain=s_atom.get_parent().get_parent().id,
                    acceptor_chain=a_atom.get_parent().get_parent().id
                ))

        # 8. Finalize instrumentation
        update_counts(self.instrumentation, accepted=len(bonds))
        t_eval_end = time.time()
        finalize_funnel(
            self.instrumentation,
            pair_gen_seconds=(t_pair_end - t_pair_start),
            eval_seconds=(t_eval_end - t_eval_start),
            build_seconds=0.0
        )
        return bonds
    
    def _calculate_centroid_ref(self, bonded_atoms):
        """Calculate centroid of a list of atoms (reference implementation)."""
        coords = [atom.get_coord() for atom in bonded_atoms]
        return np.mean(coords, axis=0)
    
    def _calculate_theta_ref(self, sulfur, bonded_atoms, acceptor):
        """
        Calculate theta angle using reference implementation.
        Angle between centroid→sulfur vector and sulfur→acceptor vector.
        """
        try:
            centroid = self._calculate_centroid_ref(bonded_atoms)
            s_coord = np.array(sulfur.get_coord())
            centroid = np.array(centroid)
            acceptor_coord = np.array(acceptor.get_coord())
            
            # Vector from centroid to sulfur (exactly as in reference)
            vec_sc = centroid - s_coord
            # Vector from sulfur to acceptor (exactly as in reference)
            vec_sa = acceptor_coord - s_coord
            
            # Check for zero vectors
            if np.linalg.norm(vec_sc) == 0 or np.linalg.norm(vec_sa) == 0:
                logger.warning("Zero vector in theta calculation, returning default")
                return 135.0
            
            # Normalize vectors
            vec_sc_norm = vec_sc / np.linalg.norm(vec_sc)
            vec_sa_norm = vec_sa / np.linalg.norm(vec_sa)
            
            # Calculate angle
            dot_product = np.dot(vec_sc_norm, vec_sa_norm)
            dot_product = np.clip(dot_product, -1.0, 1.0)
            
            angle = np.degrees(np.arccos(dot_product))
            logger.debug(f"Theta calculation: centroid={centroid}, s_coord={s_coord}, acceptor_coord={acceptor_coord}, angle={angle}")
            return angle
            
        except Exception as e:
            logger.error(f"Error in theta calculation: {e}")
            return 135.0
    
    def _calculate_delta_ref(self, x_atom, bonded_atoms, sulfur, acceptor):
        """
        Calculate delta angle using reference implementation.
        Dihedral angle for out-of-plane deviation.
        """
        try:
            centroid = self._calculate_centroid_ref(bonded_atoms)
            
            # Define bond vectors exactly as in reference
            b1 = centroid - x_atom.get_coord()
            b2 = sulfur.get_coord() - centroid
            b3 = acceptor.get_coord() - sulfur.get_coord()
            
            # Check for zero vectors
            if (np.linalg.norm(b1) == 0 or np.linalg.norm(b2) == 0 or np.linalg.norm(b3) == 0):
                logger.warning("Zero vector in delta calculation, returning default")
                return 0.0
            
            # Normalize the bond vectors
            b1 = b1 / np.linalg.norm(b1)
            b2 = b2 / np.linalg.norm(b2)
            b3 = b3 / np.linalg.norm(b3)
            
            # Calculate normal vectors to the two planes
            n1 = np.cross(b1, b2)
            n2 = np.cross(b2, b3)
            
            # Check for zero normal vectors
            if np.linalg.norm(n1) == 0 or np.linalg.norm(n2) == 0:
                logger.warning("Zero normal vector in delta calculation, returning default")
                return 0.0
            
            # Normalize the normal vectors
            n1 = n1 / np.linalg.norm(n1)
            n2 = n2 / np.linalg.norm(n2)
            
            # Calculate the dihedral angle
            cos_delta = np.dot(n1, n2)
            cos_delta = np.clip(cos_delta, -1.0, 1.0)
            delta = np.degrees(np.arccos(cos_delta))
            
            # Determine the sign of the dihedral angle
            if np.dot(np.cross(n1, n2), b2) < 0:
                delta = -delta
                
            logger.debug(f"Delta calculation: centroid={centroid}, x_atom={x_atom.get_coord()}, sulfur={sulfur.get_coord()}, acceptor={acceptor.get_coord()}, delta={delta}")
            return delta
            
        except Exception as e:
            logger.error(f"Error in delta calculation: {e}")
            return 0.0
    
    def _get_chalcogen_donors(self, model) -> List[Dict[str, Any]]:
        """Get chalcogen atoms that can act as σ-hole donors."""
        donors = []
        
        for chain in model:
            chain_id = chain.get_id()
            
            for residue in chain:
                resname = residue.get_resname()
                res_id = residue.get_id()
                
                if resname in self.chalcogen_donors:
                    for atom_name in self.chalcogen_donors[resname]:
                        if atom_name in residue:
                            donors.append({
                                'atom': residue[atom_name],
                                'residue': residue,
                                'resname': resname,
                                'chain_id': chain_id,
                                'res_id': res_id
                            })
        
        return donors
    
    def _get_acceptor_atoms(self, model) -> List[Dict[str, Any]]:
        """Get nucleophilic acceptor atoms."""
        acceptors = []
        
        for chain in model:
            chain_id = chain.get_id()
            
            for residue in chain:
                resname = residue.get_resname()
                res_id = residue.get_id()
                
                # Check backbone acceptors
                if 'O' in residue:
                    acceptors.append({
                        'atom': residue['O'],
                        'residue': residue,
                        'resname': resname,
                        'chain_id': chain_id,
                        'res_id': res_id
                    })
                
                if 'N' in residue:
                    acceptors.append({
                        'atom': residue['N'],
                        'residue': residue,
                        'resname': resname,
                        'chain_id': chain_id,
                        'res_id': res_id
                    })
                
                # Check side chain acceptors
                if resname in self.acceptor_atoms:
                    for atom_name in self.acceptor_atoms[resname]:
                        if atom_name != 'backbone' and atom_name in residue:
                            acceptors.append({
                                'atom': residue[atom_name],
                                'residue': residue,
                                'resname': resname,
                                'chain_id': chain_id,
                                'res_id': res_id
                            })
        
        return acceptors
    
    def _calculate_distance(self, atom1: Atom.Atom, atom2: Atom.Atom) -> float:
        """Calculate distance between two atoms using reference implementation."""
        coord1 = atom1.get_coord()
        coord2 = atom2.get_coord()
        return math.sqrt((coord1[0] - coord2[0])**2 +
                        (coord1[1] - coord2[1])**2 +
                        (coord1[2] - coord2[2])**2)
    
    def _classify_bond_strength(self, distance: float, theta_angle: float) -> str:
        """Classify chalcogen bond strength."""
        # Strong: short distance and optimal angle
        if distance < 3.2 and 125 <= theta_angle <= 145:
            return 'strong'
        # Weak: long distance or poor angle
        elif distance > 3.4 or theta_angle < 120 or theta_angle > 150:
            return 'weak'
        # Moderate: in between
        else:
            return 'moderate'
    
    def to_dict_list(self, bonds: List[ChalcogenBond]) -> List[Dict[str, Any]]:
        """Convert bonds to list of dictionaries."""
        return [
            {
                'type': 'chalcogen_bond',
                'chalcogen_atom': f"{bond.chalcogen_atom.get_parent().get_resname()}{bond.chalcogen_atom.get_parent().id[1]}:{bond.chalcogen_atom.name}",
                'acceptor_atom': f"{bond.acceptor_atom.get_parent().get_resname()}{bond.acceptor_atom.get_parent().id[1]}:{bond.acceptor_atom.name}",
                'distance': round(bond.distance, 2),
                'theta_angle': round(bond.theta_angle, 1),
                'delta_angle': round(bond.delta_angle, 1),
                'strength': bond.strength,
                'chalcogen_residue': bond.chalcogen_residue,
                'acceptor_residue': bond.acceptor_residue,
                'chalcogen_chain': bond.chalcogen_chain,
                'acceptor_chain': bond.acceptor_chain,
                # Standard keys for UI compatibility
                'residue1': bond.chalcogen_residue,
                'residue2': bond.acceptor_residue,
                'chain1': bond.chalcogen_chain,
                'chain2': bond.acceptor_chain
            }
            for bond in bonds
        ]
