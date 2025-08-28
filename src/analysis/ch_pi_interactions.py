"""
C-H···π interaction detection for protein structures.
Detects interactions between C-H bonds and π-systems.
"""

from typing import List, Dict, Any, Optional, Set, Tuple
import numpy as np
from Bio.PDB import Structure, Residue, Atom
from scipy.spatial import cKDTree

from utils.config import AppConfig

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
        interactions = []
        
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
                        interactions.extend(ch_interactions)
        
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
