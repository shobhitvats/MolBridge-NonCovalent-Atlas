"""
Pnictogen bond detection for protein structures.
Detects interactions involving pnictogen atoms (N, P, As) as electron acceptors.
"""

import numpy as np
from typing import List, Dict, Any, Optional
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from loguru import logger

@dataclass
class PnictogenBond:
    """Represents a detected pnictogen bond."""
    pnictogen_atom: Atom.Atom
    acceptor_atom: Atom.Atom
    distance: float
    angle: float
    strength: str
    pnictogen_residue: str
    acceptor_residue: str
    pnictogen_chain: str
    acceptor_chain: str

class PnictogenBondDetector:
    """Detects pnictogen bonds in protein structures."""
    
    def __init__(self, config):
        """
        Initialize pnictogen bond detector.
        
        Args:
            config: AppConfig object containing interaction parameters
        """
        self.config = config
        self.interaction_config = config.interactions
        self.distance_cutoff = getattr(self.interaction_config, 'pnictogen_distance_cutoff', 4.0)
        self.angle_cutoff = getattr(self.interaction_config, 'pnictogen_angle_cutoff', 140.0)
        
        # Pnictogen atoms (group 15 elements)
        self.pnictogen_atoms = {'N', 'P', 'AS'}
        
        # Potential acceptor atoms (electron-rich)
        self.acceptor_atoms = {
            'O', 'N', 'S', 'SE', 'F', 'CL', 'BR', 'I'
        }
    
    def detect_pnictogen_bonds(self, structure: Structure.Structure) -> List[PnictogenBond]:
        """Detect pnictogen bonds in structure."""
        bonds = []
        
        # Get all atoms
        atoms = list(structure.get_atoms())
        
        # Find pnictogen and acceptor atoms
        pnictogen_atoms = [atom for atom in atoms if atom.element in self.pnictogen_atoms]
        acceptor_atoms = [atom for atom in atoms if atom.element in self.acceptor_atoms]
        
        logger.info(f"Found {len(pnictogen_atoms)} pnictogen atoms and {len(acceptor_atoms)} potential acceptors")
        
        for pnictogen in pnictogen_atoms:
            for acceptor in acceptor_atoms:
                # Skip if same atom or same residue
                if pnictogen == acceptor or pnictogen.get_parent() == acceptor.get_parent():
                    continue
                
                distance = self._calculate_distance(pnictogen, acceptor)
                
                if distance <= self.distance_cutoff:
                    # Calculate angle (pnictogen-bond-acceptor)
                    angle = self._calculate_pnictogen_angle(pnictogen, acceptor)
                    
                    if angle >= self.angle_cutoff:
                        strength = self._classify_bond_strength(distance)
                        
                        bond = PnictogenBond(
                            pnictogen_atom=pnictogen,
                            acceptor_atom=acceptor,
                            distance=distance,
                            angle=angle,
                            strength=strength,
                            pnictogen_residue=f"{pnictogen.get_parent().get_resname()}{pnictogen.get_parent().id[1]}",
                            acceptor_residue=f"{acceptor.get_parent().get_resname()}{acceptor.get_parent().id[1]}",
                            pnictogen_chain=pnictogen.get_parent().get_parent().id,
                            acceptor_chain=acceptor.get_parent().get_parent().id
                        )
                        bonds.append(bond)
        
        logger.info(f"Detected {len(bonds)} pnictogen bonds")
        return bonds
    
    def _calculate_distance(self, atom1: Atom.Atom, atom2: Atom.Atom) -> float:
        """Calculate distance between two atoms."""
        return np.linalg.norm(atom1.coord - atom2.coord)
    
    def _calculate_pnictogen_angle(self, pnictogen: Atom.Atom, acceptor: Atom.Atom) -> float:
        """Calculate the pnictogen bond angle."""
        # Get neighboring atoms to pnictogen
        pnictogen_residue = pnictogen.get_parent()
        neighbors = [atom for atom in pnictogen_residue.get_atoms() 
                    if atom != pnictogen and self._calculate_distance(atom, pnictogen) < 2.0]
        
        if not neighbors:
            return 180.0  # Assume linear if no neighbors found
        
        # Use the closest neighbor to define the angle
        neighbor = min(neighbors, key=lambda x: self._calculate_distance(x, pnictogen))
        
        # Calculate angle: neighbor-pnictogen-acceptor
        vec1 = neighbor.coord - pnictogen.coord
        vec2 = acceptor.coord - pnictogen.coord
        
        cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle = np.degrees(np.arccos(cos_angle))
        
        return angle
    
    def _classify_bond_strength(self, distance: float) -> str:
        """Classify bond strength based on distance."""
        if distance <= 3.0:
            return 'strong'
        elif distance <= 3.5:
            return 'moderate'
        else:
            return 'weak'
    
    def to_dict_list(self, bonds: List[PnictogenBond]) -> List[Dict[str, Any]]:
        """Convert bonds to list of dictionaries."""
        return [
            {
                'type': 'pnictogen_bond',
                'pnictogen_atom': f"{bond.pnictogen_atom.get_parent().get_resname()}{bond.pnictogen_atom.get_parent().id[1]}:{bond.pnictogen_atom.name}",
                'acceptor_atom': f"{bond.acceptor_atom.get_parent().get_resname()}{bond.acceptor_atom.get_parent().id[1]}:{bond.acceptor_atom.name}",
                'distance': round(bond.distance, 2),
                'angle': round(bond.angle, 1),
                'strength': bond.strength,
                'pnictogen_residue': bond.pnictogen_residue,
                'acceptor_residue': bond.acceptor_residue,
                'pnictogen_chain': bond.pnictogen_chain,
                'acceptor_chain': bond.acceptor_chain,
                # Standard keys for UI compatibility
                'residue1': bond.pnictogen_residue,
                'residue2': bond.acceptor_residue,
                'chain1': bond.pnictogen_chain,
                'chain2': bond.acceptor_chain
            }
            for bond in bonds
        ]
