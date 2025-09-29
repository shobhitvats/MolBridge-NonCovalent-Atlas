"""
Tetrel (Carbon σ-hole) bond detection for protein structures.
Detects C atoms bound to electron-withdrawing groups (N, O, Cl, F, C=C, C=O, CH3, CF3) 
with O/N/S/F/Cl/P acceptors, C...acceptor 2.5–3.6 Å, angles θ ≥160° (σ-hole directionality).
"""

import numpy as np
import math
from typing import List, Dict, Any, Optional, Tuple
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from loguru import logger
from utils.settings import get_settings
import numpy as np
from .base_detector import register_detector
from geometry.core import pairwise_within_cutoff, norms
from utils.instrumentation import init_funnel, update_counts, finalize_funnel

@dataclass
class TetrelBond:
    """Represents a detected tetrel bond."""
    tetrel_atom: Atom.Atom
    acceptor_atom: Atom.Atom
    distance: float
    theta1_angle: float  # First C-Z...A angle (≥160°)
    theta2_angle: float  # Second C-Z...A angle (≥160°) if applicable
    nc_ratio: float     # Distance ratio to vdW sum
    strength: str
    tetrel_residue: str
    acceptor_residue: str
    tetrel_chain: str
    acceptor_chain: str

@register_detector("tetrel_bond", method="detect_tetrel_bonds")
class TetrelBondDetector:
    """Detects tetrel (carbon σ-hole) bonds in protein structures."""
    
    def __init__(self, config):
        """
        Initialize tetrel bond detector.
        
        Args:
            config: AppConfig object containing interaction parameters
        """
        self.config = config
        self.interaction_config = config.interactions
        self.distance_min = 2.5  # Minimum C...acceptor distance
        self.distance_max = 3.6  # Maximum C...acceptor distance
        self.angle_min = 160.0   # Minimum θ angle (can be relaxed to 140°)
        self.angle_relaxed = 140.0  # Relaxed θ angle criterion
        self.nc_ratio_max = 1.1  # Maximum Nc ratio (distance/vdW_sum)
        
        # Van der Waals radii (Å)
        self.vdw_radii = {
            'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
            'F': 1.47, 'Cl': 1.75, 'P': 1.80, 'H': 1.20
        }
        
        # Electron-withdrawing groups and residues that can have tetrel donors
        self.tetrel_residues = {
            'LYS': ['C', 'CA', 'CB', 'CG', 'CD', 'CE'],  # Lysine carbons (N-bound)
            'ALA': ['CA', 'CB'],                         # Alanine (methyl)
            'ARG': ['C', 'CA', 'CB', 'CG', 'CD', 'CZ'], # Arginine carbons (N-bound)
            'ASP': ['C', 'CA', 'CB', 'CG'],             # Aspartate carbons (O-bound)
            'GLU': ['C', 'CA', 'CB', 'CG', 'CD'],       # Glutamate carbons (O-bound)
            'MET': ['CA', 'CB', 'CG', 'CE'],            # Methionine carbons (S-bound)
            'CYS': ['CA', 'CB'],                        # Cysteine carbons (S-bound)
            'VAL': ['CA', 'CB', 'CG1', 'CG2'],         # Valine (methyl groups)
            'LEU': ['CA', 'CB', 'CG', 'CD1', 'CD2'],   # Leucine (methyl groups)
            'ILE': ['CA', 'CB', 'CG1', 'CG2', 'CD1'],  # Isoleucine (methyl groups)
            'THR': ['CA', 'CB', 'CG2'],                # Threonine (O-bound)
            'SER': ['CA', 'CB'],                        # Serine (O-bound)
            'TYR': ['CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],  # Tyrosine aromatics
            'PHE': ['CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],  # Phenylalanine aromatics
            'TRP': ['CA', 'CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],  # Tryptophan aromatics
            'HIS': ['CA', 'CB', 'CG', 'CD2', 'CE1'],   # Histidine aromatics (N-bound)
            'PRO': ['CA', 'CB', 'CG', 'CD'],           # Proline carbons
            'ASN': ['CA', 'CB', 'CG'],                 # Asparagine carbons (N/O-bound)
            'GLN': ['CA', 'CB', 'CG', 'CD'],          # Glutamine carbons (N/O-bound)
        }
        
        # Acceptor atoms with lone pairs
        self.acceptor_elements = ['O', 'N', 'S', 'F', 'Cl', 'P']
    
    def detect_tetrel_bonds(self, structure: Structure.Structure) -> List[TetrelBond]:
        """
        Detect tetrel bonds in structure using the new logic.
        
        Steps:
        1. Select candidate pairs: C atoms with nearby O/N/S/F/Cl/P within 3.6 Å
        2. Filter by chemical environment: C bound to electron-withdrawing groups
        3. Measure angle: C-Z...A angle ≥160° (relaxed to 140°)
        4. Exclude H-bonding: No H within vdW contact with acceptor
        5. Check Nc ratio: distance/vdW_sum ≤ 1.1
        """
        settings = get_settings()
        if getattr(settings, 'enable_vector_geom', False):
            try:
                return self._vector_detect(structure)
            except Exception:  # pragma: no cover
                pass
        import time as _t
        t_pair_start = _t.time()
        bonds = []
        
        try:
            model = structure[0]  # Use first model
            
            # Step 1: Find all candidate carbon-acceptor pairs within distance
            candidate_pairs = self._find_candidate_pairs(model)
            logger.info(f"Found {len(candidate_pairs)} candidate tetrel pairs")
            
            for carbon_atom, acceptor_atom in candidate_pairs:
                # Step 2: Check if carbon has electron-withdrawing environment
                if not self._has_electron_withdrawing_environment(carbon_atom):
                    continue
                
                distance = self._calculate_distance(carbon_atom, acceptor_atom)
                
                # Step 5: Check Nc ratio
                nc_ratio = self._calculate_nc_ratio(carbon_atom, acceptor_atom, distance)
                if nc_ratio > self.nc_ratio_max:
                    continue
                
                # Step 3: Measure angles for σ-hole directionality
                angles = self._calculate_tetrel_angles(carbon_atom, acceptor_atom)
                if not angles:
                    continue
                
                theta1, theta2 = angles
                
                # Check angle criteria (relaxed to 140°)
                if theta1 < self.angle_relaxed and (theta2 is None or theta2 < self.angle_relaxed):
                    continue
                
                # Step 4: Exclude conventional H-bonding
                if self._has_competing_hbond(carbon_atom, acceptor_atom):
                    continue
                
                # Calculate interaction strength
                strength = self._classify_bond_strength(distance, theta1, nc_ratio)
                
                bond = TetrelBond(
                    tetrel_atom=carbon_atom,
                    acceptor_atom=acceptor_atom,
                    distance=distance,
                    theta1_angle=theta1,
                    theta2_angle=theta2 if theta2 is not None else 0.0,
                    nc_ratio=nc_ratio,
                    strength=strength,
                    tetrel_residue=f"{carbon_atom.get_parent().get_resname()}{carbon_atom.get_parent().id[1]}",
                    acceptor_residue=f"{acceptor_atom.get_parent().get_resname()}{acceptor_atom.get_parent().id[1]}",
                    tetrel_chain=carbon_atom.get_parent().get_parent().id,
                    acceptor_chain=acceptor_atom.get_parent().get_parent().id
                )
                bonds.append(bond)
            
            t_eval_end = _t.time()
            logger.info(f"Detected {len(bonds)} tetrel bonds")
            self.instrumentation = init_funnel(
                raw_pairs=len(candidate_pairs),
                candidate_pairs=len(candidate_pairs),
                accepted_pairs=len(bonds),
                core_pair_generation=False
            )
            finalize_funnel(
                self.instrumentation,
                pair_gen_seconds=0.0,
                eval_seconds=(t_eval_end - t_pair_start),
                build_seconds=0.0
            )
            
        except Exception as e:
            logger.error(f"Error detecting tetrel bonds: {e}")
        
        return bonds

    # ---- Vector candidate pruning (carbon vs acceptor distance matrix) ----
    def _vector_detect(self, structure: Structure.Structure) -> List[TetrelBond]:
        import time as _t
        t_pair_start = _t.time()
        bonds: List[TetrelBond] = []
        try:
            model = structure[0]
        except Exception:
            return bonds
        carbon_atoms: List[Atom.Atom] = []
        carbon_bond_vectors: List[np.ndarray] = []
        for chain in model:
            for residue in chain:
                rname = residue.get_resname()
                if rname in self.tetrel_residues:
                    for an in self.tetrel_residues[rname]:
                        if an in residue and residue[an].element == 'C':
                            c_atom = residue[an]
                            carbon_atoms.append(c_atom)
                            bonded_vecs = []
                            for at in residue:
                                if at is not c_atom and at.element != 'H' and self._calculate_distance(c_atom, at) < 1.8:
                                    bonded_vecs.append(at.get_coord() - c_atom.get_coord())
                            carbon_bond_vectors.append(np.vstack(bonded_vecs) if bonded_vecs else np.zeros((0,3), dtype='float32'))
        acceptor_atoms: List[Atom.Atom] = []
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element in self.acceptor_elements:
                        acceptor_atoms.append(atom)
        if not carbon_atoms or not acceptor_atoms:
            return bonds
        c_coords = np.vstack([a.get_coord() for a in carbon_atoms]).astype('float32')
        a_coords = np.vstack([a.get_coord() for a in acceptor_atoms]).astype('float32')
        try:
            ci_all, ai_all = pairwise_within_cutoff(c_coords, a_coords, self.distance_max, use_kdtree=True)
            core = True
        except Exception:  # pragma: no cover
            diff = c_coords[:, None, :] - a_coords[None, :, :]
            dist2 = np.sum(diff*diff, axis=-1)
            mask = dist2 <= (self.distance_max ** 2)
            ci_all, ai_all = np.where(mask)
            ci_all = ci_all.astype('int32'); ai_all = ai_all.astype('int32')
            core = False
        t_pair_end = _t.time()
        raw_pairs = int(c_coords.shape[0] * a_coords.shape[0])
        vecs = a_coords[ai_all] - c_coords[ci_all]
        dists = norms(vecs)
        keep_mask = (dists >= self.distance_min) & (dists <= self.distance_max)
        ci = ci_all[keep_mask]; ai = ai_all[keep_mask]; dists_kept = dists[keep_mask]
        candidate_pairs = int(len(ci))
        from utils.kdtree_thresholds import adapt_threshold, get_last_density
        kdtree_used = (candidate_pairs / raw_pairs) < 0.55 if raw_pairs else False
        new_thresh, changed, reason = adapt_threshold('tetrel', candidate_pairs, kdtree_used)
        accepted = 0
        t_eval_start = t_pair_end
        for idx in range(candidate_pairs):
            c_idx = int(ci[idx]); a_idx = int(ai[idx])
            c_atom = carbon_atoms[c_idx]; a_atom = acceptor_atoms[a_idx]
            if c_atom.get_parent() is a_atom.get_parent():
                continue
            distance = float(dists_kept[idx])
            if not self._has_electron_withdrawing_environment(c_atom):
                continue
            nc_ratio = self._calculate_nc_ratio(c_atom, a_atom, distance)
            if nc_ratio > self.nc_ratio_max:
                continue
            bonded_mat = carbon_bond_vectors[c_idx]
            if bonded_mat.shape[0] == 0:
                continue
            vec_ca = (a_coords[a_idx] - c_coords[c_idx])
            norm_vec = np.linalg.norm(vec_ca) or 1.0
            unit_ca = vec_ca / norm_vec
            bnorms = np.linalg.norm(bonded_mat, axis=1); bnorms[bnorms == 0] = 1.0
            unit_cz = bonded_mat / bnorms[:, None]
            dots = np.clip(np.dot(unit_cz, unit_ca), -1.0, 1.0)
            acute = np.degrees(np.arccos(np.abs(dots)))
            sigma_angles = 180.0 - acute
            sigma_angles_sorted = np.sort(sigma_angles)[::-1]
            theta1 = float(sigma_angles_sorted[0]) if sigma_angles_sorted.size else 0.0
            theta2_val = float(sigma_angles_sorted[1]) if sigma_angles_sorted.size > 1 else None
            if theta1 < self.angle_relaxed and (theta2_val is None or theta2_val < self.angle_relaxed):
                continue
            if self._has_competing_hbond(c_atom, a_atom):
                continue
            strength = self._classify_bond_strength(distance, theta1, nc_ratio)
            bonds.append(TetrelBond(
                tetrel_atom=c_atom,
                acceptor_atom=a_atom,
                distance=distance,
                theta1_angle=theta1,
                theta2_angle=theta2_val if theta2_val is not None else 0.0,
                nc_ratio=nc_ratio,
                strength=strength,
                tetrel_residue=f"{c_atom.get_parent().get_resname()}{c_atom.get_parent().id[1]}",
                acceptor_residue=f"{a_atom.get_parent().get_resname()}{a_atom.get_parent().id[1]}",
                tetrel_chain=c_atom.get_parent().get_parent().id,
                acceptor_chain=a_atom.get_parent().get_parent().id
            ))
            accepted += 1
        t_eval_end = _t.time()
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=accepted,
            core_pair_generation=core,
            extra={
                'carbons': c_coords.shape[0],
                'acceptors': a_coords.shape[0],
                'candidate_density': (candidate_pairs / raw_pairs) if raw_pairs else 0.0,
                'threshold': new_thresh,
                'threshold_changed': changed,
                'adapt_reason': reason,
                'last_density': get_last_density('tetrel')
            }
        )
        finalize_funnel(
            self.instrumentation,
            pair_gen_seconds=(t_pair_end - t_pair_start),
            eval_seconds=(t_eval_end - t_eval_start),
            build_seconds=0.0
        )
        return bonds
    
    def _find_candidate_pairs(self, model) -> List[Tuple[Atom.Atom, Atom.Atom]]:
        """Find carbon-acceptor pairs within distance cutoff."""
        pairs = []
        
        # Get all carbon atoms from relevant residues
        carbon_atoms = []
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname in self.tetrel_residues:
                    for atom_name in self.tetrel_residues[resname]:
                        if atom_name in residue and residue[atom_name].element == 'C':
                            carbon_atoms.append(residue[atom_name])
        
        # Get all potential acceptor atoms
        acceptor_atoms = []
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element in self.acceptor_elements:
                        acceptor_atoms.append(atom)
        
        # Find pairs within distance cutoff
        for carbon in carbon_atoms:
            for acceptor in acceptor_atoms:
                # Skip same residue
                if carbon.get_parent() == acceptor.get_parent():
                    continue
                
                distance = self._calculate_distance(carbon, acceptor)
                if self.distance_min <= distance <= self.distance_max:
                    pairs.append((carbon, acceptor))
        
        return pairs
    
    def _has_electron_withdrawing_environment(self, carbon_atom: Atom.Atom) -> bool:
        """
        Check if carbon is bound to electron-withdrawing groups.
        Look for bonds to N, O, Cl, F or being part of C=C, C=O, aromatic systems.
        """
        residue = carbon_atom.get_parent()
        carbon_coord = carbon_atom.get_coord()
        
        # Find atoms bonded to this carbon (within 1.6 Å)
        bonded_atoms = []
        for atom in residue:
            if atom != carbon_atom and self._calculate_distance(carbon_atom, atom) < 1.6:
                bonded_atoms.append(atom)
        
        # Check for electron-withdrawing elements
        bonded_elements = [atom.element for atom in bonded_atoms]
        electron_withdrawing_elements = ['N', 'O', 'Cl', 'F', 'S']
        
        # Direct bond to electron-withdrawing atom
        if any(elem in bonded_elements for elem in electron_withdrawing_elements):
            return True
        
        # Check for aromatic/conjugated systems (multiple C bonds)
        carbon_bonds = bonded_elements.count('C')
        if carbon_bonds >= 2:  # Potentially aromatic or part of C=C
            return True
        
        # Methyl groups (especially on side chains) can have σ-hole character
        hydrogen_bonds = bonded_elements.count('H')
        if hydrogen_bonds >= 2:  # Methyl or methylene group
            return True
        
        return False
    
    def _calculate_tetrel_angles(self, carbon_atom: Atom.Atom, acceptor_atom: Atom.Atom) -> Optional[Tuple[float, Optional[float]]]:
        """
        Calculate angles for tetrel bond directionality.
        Returns (theta1, theta2) where theta is the C-Z...A angle.
        """
        try:
            residue = carbon_atom.get_parent()
            carbon_coord = np.array(carbon_atom.get_coord())
            acceptor_coord = np.array(acceptor_atom.get_coord())
            
            # Find atoms bonded to carbon (Z atoms)
            bonded_atoms = []
            for atom in residue:
                if (atom != carbon_atom and atom.element != 'H' and 
                    self._calculate_distance(carbon_atom, atom) < 1.8):
                    bonded_atoms.append(atom)
            
            if not bonded_atoms:
                return None
            
            angles = []
            for bonded_atom in bonded_atoms:
                bonded_coord = np.array(bonded_atom.get_coord())
                
                # Vector from carbon to bonded atom (C-Z)
                vec_cz = bonded_coord - carbon_coord
                # Vector from carbon to acceptor (C...A)
                vec_ca = acceptor_coord - carbon_coord
                
                # Calculate angle between C-Z bond and C...A vector
                # The σ-hole is along the extension of C-Z bond
                cos_angle = np.dot(vec_cz, vec_ca) / (np.linalg.norm(vec_cz) * np.linalg.norm(vec_ca))
                cos_angle = np.clip(cos_angle, -1.0, 1.0)
                angle = np.degrees(np.arccos(abs(cos_angle)))  # We want the acute angle
                
                # Convert to the angle from the σ-hole direction (180° - angle)
                sigma_hole_angle = 180.0 - angle
                angles.append(sigma_hole_angle)
            
            # Return the best angle(s)
            angles.sort(reverse=True)  # Sort in descending order
            theta1 = angles[0] if angles else 0.0
            theta2 = angles[1] if len(angles) > 1 else None
            
            return (theta1, theta2)
            
        except Exception as e:
            logger.warning(f"Error calculating tetrel angles: {e}")
            return None
    
    def _has_competing_hbond(self, carbon_atom: Atom.Atom, acceptor_atom: Atom.Atom) -> bool:
        """
        Check if carbon has hydrogens within vdW contact of acceptor.
        This would indicate competing conventional H-bonding.
        """
        try:
            residue = carbon_atom.get_parent()
            acceptor_coord = np.array(acceptor_atom.get_coord())
            
            # Find hydrogens bonded to carbon
            for atom in residue:
                if (atom.element == 'H' and 
                    self._calculate_distance(carbon_atom, atom) < 1.2):
                    
                    h_distance = self._calculate_distance(atom, acceptor_atom)
                    # Check if H is within vdW contact of acceptor
                    vdw_sum = self.vdw_radii.get('H', 1.2) + self.vdw_radii.get(acceptor_atom.element, 1.5)
                    
                    if h_distance <= vdw_sum:
                        return True  # Competing H-bond detected
            
            return False
            
        except Exception as e:
            logger.warning(f"Error checking H-bonding competition: {e}")
            return False
    
    def _calculate_nc_ratio(self, carbon_atom: Atom.Atom, acceptor_atom: Atom.Atom, distance: float) -> float:
        """Calculate the Nc ratio (measured distance / sum of vdW radii)."""
        vdw_carbon = self.vdw_radii.get('C', 1.70)
        vdw_acceptor = self.vdw_radii.get(acceptor_atom.element, 1.52)
        vdw_sum = vdw_carbon + vdw_acceptor
        
        return distance / vdw_sum
    
    def _calculate_distance(self, atom1: Atom.Atom, atom2: Atom.Atom) -> float:
        """Calculate distance between two atoms."""
        coord1 = atom1.get_coord()
        coord2 = atom2.get_coord()
        return math.sqrt((coord1[0] - coord2[0])**2 +
                        (coord1[1] - coord2[1])**2 +
                        (coord1[2] - coord2[2])**2)
    
    def _classify_bond_strength(self, distance: float, theta_angle: float, nc_ratio: float) -> str:
        """
        Classify tetrel bond strength based on distance, angle, and Nc ratio.
        Uses optional energy estimation for validation.
        """
        # Calculate approximate interaction energy for C...O (if applicable)
        # EDA = 3.5 × 10^7 × e^(-d_C-O/0.184)
        # This is optional but can help filter very weak interactions
        
        # Primary classification based on geometric criteria
        if distance <= 3.0 and theta_angle >= 160.0 and nc_ratio <= 1.0:
            return 'strong'
        elif distance <= 3.3 and theta_angle >= 150.0 and nc_ratio <= 1.05:
            return 'moderate'
        elif distance <= 3.6 and theta_angle >= 140.0 and nc_ratio <= 1.1:
            return 'weak'
        else:
            return 'very_weak'
    
    def to_dict_list(self, bonds: List[TetrelBond]) -> List[Dict[str, Any]]:
        """Convert bonds to list of dictionaries."""
        return [
            {
                'type': 'tetrel_bond',
                'tetrel_atom': f"{bond.tetrel_atom.get_parent().get_resname()}{bond.tetrel_atom.get_parent().id[1]}:{bond.tetrel_atom.name}",
                'acceptor_atom': f"{bond.acceptor_atom.get_parent().get_resname()}{bond.acceptor_atom.get_parent().id[1]}:{bond.acceptor_atom.name}",
                'distance': round(bond.distance, 2),
                'theta1_angle': round(bond.theta1_angle, 1),
                'theta2_angle': round(bond.theta2_angle, 1) if bond.theta2_angle else None,
                'nc_ratio': round(bond.nc_ratio, 3),
                'strength': bond.strength,
                'tetrel_residue': bond.tetrel_residue,
                'acceptor_residue': bond.acceptor_residue,
                'tetrel_chain': bond.tetrel_chain,
                'acceptor_chain': bond.acceptor_chain,
                # Standard keys for UI compatibility
                'residue1': bond.tetrel_residue,
                'residue2': bond.acceptor_residue,
                'chain1': bond.tetrel_chain,
                'chain2': bond.acceptor_chain
            }
            for bond in bonds
        ]
