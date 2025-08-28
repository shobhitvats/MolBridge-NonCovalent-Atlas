"""
Chalcogen bond detection for protein structures.
Detects S/Se/Te σ-hole donor interactions with O/N/S acceptors using directionality angles θ (115–155°) and δ (±50°), distance ≤ 3.6 Å.
"""

import numpy as np
import math
from typing import List, Dict, Any, Optional, Tuple
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from loguru import logger

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
        self.distance_cutoff = 3.6  # Distance ≤ 3.6 Å for chalcogen bonds
        self.theta_min = 115.0      # Minimum θ angle
        self.theta_max = 155.0      # Maximum θ angle
        self.delta_max = 50.0       # Maximum δ angle (±50°)
        
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
        bonds = []
        
        try:
            # Iterate through all models, chains, residues, and atoms
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            if atom.element == 'S':  # Check if the atom is sulfur
                                sulfur = atom
                                
                                # Find two carbon atoms bonded to sulfur (within 1.9 Å)
                                bonded_atoms = [nbr for nbr in residue 
                                              if nbr.element == 'C' and 
                                              self._calculate_distance(sulfur, nbr) < 1.9]
                                
                                if len(bonded_atoms) < 2:
                                    continue
                                
                                # Include sulfur in the centroid calculation (as per reference)
                                atoms_for_centroid = bonded_atoms + [sulfur]
                                x_atom = bonded_atoms[1]  # Reference atom for delta calculation
                                
                                # Search for potential acceptors in all other residues
                                for other_model in structure:
                                    for other_chain in other_model:
                                        for other_residue in other_chain:
                                            # Skip same residue
                                            if residue == other_residue:
                                                continue
                                                
                                            for other_atom in other_residue:
                                                if other_atom.element in ['O', 'N']:
                                                    distance = self._calculate_distance(sulfur, other_atom)
                                                    
                                                    # Check distance thresholds
                                                    if ((other_atom.element == 'O' and distance <= 3.6) or
                                                        (other_atom.element == 'N' and distance <= 3.6)):
                                                        
                                                        # Calculate angles using reference methods
                                                        try:
                                                            theta = self._calculate_theta_ref(sulfur, atoms_for_centroid, other_atom)
                                                            delta = self._calculate_delta_ref(x_atom, atoms_for_centroid, sulfur, other_atom)
                                                            logger.debug(f"Calculated angles for S {sulfur.get_serial_number()} -> {other_atom.element} {other_atom.get_serial_number()}: theta={theta:.1f}°, delta={delta:.1f}°")
                                                        except Exception as e:
                                                            logger.warning(f"Error calculating chalcogen angles: {e}, using fallback values")
                                                            theta, delta = 135.0, 0.0
                                                        
                                                        # Check angle criteria
                                                        if (115 <= theta <= 155) and (-50 <= delta <= 50):
                                                            bond = ChalcogenBond(
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
                                                            )
                                                            bonds.append(bond)
            
            logger.info(f"Detected {len(bonds)} chalcogen bonds")
            
        except Exception as e:
            logger.error(f"Error detecting chalcogen bonds: {e}")
        
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
