"""
London dispersion interaction detection for protein structures.
Detects short aliphatic/aliphatic or aromatic/aromatic contacts with C...C/H...H distances ~4.1–5.4 Å.
"""

import numpy as np
from typing import List, Dict, Any, Optional
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from loguru import logger

@dataclass
class DispersionInteraction:
    """Represents a detected London dispersion interaction."""
    atom1: Atom.Atom
    atom2: Atom.Atom
    distance: float
    interaction_type: str  # 'aliphatic_aliphatic', 'aromatic_aromatic', 'aliphatic_aromatic'
    strength: str
    residue1: str
    residue2: str
    chain1: str
    chain2: str

class DispersionDetector:
    """Detects London dispersion interactions in protein structures."""
    
    def __init__(self, config):
        """
        Initialize dispersion detector.
        
        Args:
            config: AppConfig object containing interaction parameters
        """
        self.config = config
        self.interaction_config = config.interactions
        self.distance_min = 4.1  # Minimum distance for London dispersion
        self.distance_max = 5.4  # Maximum distance for London dispersion
        
        # Aliphatic carbon atoms (sp3 carbons)
        self.aliphatic_atoms = {
            'ALA': ['CB'],
            'VAL': ['CB', 'CG1', 'CG2'],
            'LEU': ['CB', 'CG', 'CD1', 'CD2'],
            'ILE': ['CB', 'CG1', 'CG2', 'CD1'],
            'MET': ['CB', 'CG', 'CE'],
            'PRO': ['CB', 'CG', 'CD'],
            'LYS': ['CB', 'CG', 'CD', 'CE'],
            'ARG': ['CB', 'CG', 'CD'],
            'GLU': ['CB', 'CG'],
            'ASP': ['CB'],
            'GLN': ['CB', 'CG'],
            'ASN': ['CB'],
            'SER': ['CB'],
            'THR': ['CB', 'CG2'],
            'CYS': ['CB']
        }
        
        # Aromatic carbon atoms (sp2 carbons in rings)
        self.aromatic_atoms = {
            'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TRP': ['CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
            'HIS': ['CG', 'CD2', 'CE1']
        }
        
        # Van der Waals radii for atoms (in Ångströms)
        self.vdw_radii = {
            'C': 1.70,
            'N': 1.55,
            'O': 1.52,
            'S': 1.80,
            'H': 1.09
        }
        
        # Atoms that typically participate in dispersion interactions
        self.dispersion_atoms = {'C', 'S', 'P', 'F', 'CL', 'BR', 'I'}
    
    def detect_dispersion_interactions(self, structure: Structure.Structure) -> List[DispersionInteraction]:
        """Detect London dispersion interactions in structure."""
        interactions = []
        
        # Get all atoms that can participate in dispersion interactions
        atoms = [atom for atom in structure.get_atoms() if atom.element in self.dispersion_atoms]
        
        logger.info(f"Found {len(atoms)} atoms for dispersion analysis")
        
        # Compare all pairs of atoms
        for i, atom1 in enumerate(atoms):
            for atom2 in atoms[i+1:]:
                # Skip if same residue or adjacent residues (backbone interactions)
                if self._should_skip_pair(atom1, atom2):
                    continue
                
                distance = self._calculate_distance(atom1, atom2)
                
                # Check if within dispersion range
                if self.min_distance <= distance <= self.distance_cutoff:
                    # Additional check: must be close to sum of van der Waals radii
                    vdw_sum = self._get_vdw_sum(atom1, atom2)
                    
                    if distance <= vdw_sum + 0.5:  # Allow some tolerance
                        strength = self._classify_interaction_strength(distance, vdw_sum)
                        
                        interaction = DispersionInteraction(
                            atom1=atom1,
                            atom2=atom2,
                            distance=distance,
                            strength=strength,
                            residue1=f"{atom1.get_parent().get_resname()}{atom1.get_parent().id[1]}",
                            residue2=f"{atom2.get_parent().get_resname()}{atom2.get_parent().id[1]}",
                            chain1=atom1.get_parent().get_parent().id,
                            chain2=atom2.get_parent().get_parent().id
                        )
                        interactions.append(interaction)
        
        logger.info(f"Detected {len(interactions)} London dispersion interactions")
        return interactions
    
    def _should_skip_pair(self, atom1: Atom.Atom, atom2: Atom.Atom) -> bool:
        """Check if atom pair should be skipped."""
        res1 = atom1.get_parent()
        res2 = atom2.get_parent()
        
        # Skip same residue
        if res1 == res2:
            return True
        
        # Skip adjacent residues in same chain (backbone interactions)
        if (res1.get_parent() == res2.get_parent() and 
            abs(res1.id[1] - res2.id[1]) <= 1):
            return True
        
        return False
    
    def _calculate_distance(self, atom1: Atom.Atom, atom2: Atom.Atom) -> float:
        """Calculate distance between two atoms."""
        return np.linalg.norm(atom1.coord - atom2.coord)
    
    def _get_vdw_sum(self, atom1: Atom.Atom, atom2: Atom.Atom) -> float:
        """Get sum of van der Waals radii for two atoms."""
        radius1 = self.vdw_radii.get(atom1.element, 1.7)  # Default to carbon
        radius2 = self.vdw_radii.get(atom2.element, 1.7)
        return radius1 + radius2
    
    def _classify_interaction_strength(self, distance: float, vdw_sum: float) -> str:
        """Classify interaction strength based on distance relative to vdW radii."""
        ratio = distance / vdw_sum
        
        if ratio <= 1.05:
            return 'strong'
        elif ratio <= 1.15:
            return 'moderate'
        else:
            return 'weak'
    
    def to_dict_list(self, interactions: List[DispersionInteraction]) -> List[Dict[str, Any]]:
        """Convert interactions to list of dictionaries."""
        return [
            {
                'type': 'dispersion',
                'atom1': f"{interaction.atom1.get_parent().get_resname()}{interaction.atom1.get_parent().id[1]}:{interaction.atom1.name}",
                'atom2': f"{interaction.atom2.get_parent().get_resname()}{interaction.atom2.get_parent().id[1]}:{interaction.atom2.name}",
                'distance': round(interaction.distance, 2),
                'strength': interaction.strength,
                'residue1': interaction.residue1,
                'residue2': interaction.residue2,
                'chain1': interaction.chain1,
                'chain2': interaction.chain2
            }
            for interaction in interactions
        ]
