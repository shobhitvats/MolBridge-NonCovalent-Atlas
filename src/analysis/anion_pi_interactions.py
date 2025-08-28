"""
Anion-π interaction detection for protein structures.
Detects electron-rich ions close to π-acidic rings (positive quadrupole); 3.5–4.5 Å; typically Asp/Glu-aromatic interactions.
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from loguru import logger

@dataclass
class AnionPiInteraction:
    """Represents a detected anion-π interaction."""
    anion_atom: Atom.Atom
    pi_center: np.ndarray
    pi_residue_obj: Residue.Residue
    distance: float
    strength: str
    anion_residue: str
    pi_residue: str
    anion_chain: str
    pi_chain: str

class AnionPiDetector:
    """Detects anion-π interactions in protein structures."""
    
    def __init__(self, config):
        """
        Initialize anion-π detector.
        
        Args:
            config: AppConfig object containing interaction parameters
        """
        self.config = config
        self.interaction_config = config.interactions
        self.distance_min = 3.5  # Minimum distance for anion-π interactions
        self.distance_max = 4.5  # Maximum distance for anion-π interactions
        
        # Anionic (electron-rich) groups in proteins
        self.anion_groups = {
            'ASP': ['OD1', 'OD2'],  # Aspartate carboxylate (typically negatively charged)
            'GLU': ['OE1', 'OE2'],  # Glutamate carboxylate (typically negatively charged)
        }
        
        # π-acidic aromatic residues (positive quadrupole moment)
        # These are electron-deficient rings that attract anions
        self.pi_acidic_residues = {
            'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],  # Benzene ring (π-acidic)
            'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],  # Phenol ring (π-acidic)
            'TRP': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],  # Indole (π-acidic)
            # Note: HIS is generally π-basic due to nitrogen lone pairs, so excluded from π-acidic interactions
        }
    
    def detect_anion_pi_interactions(self, structure: Structure.Structure) -> List[AnionPiInteraction]:
        """Detect anion-π interactions in structure."""
        interactions = []
        
        # Find anionic atoms
        anion_atoms = self._find_anion_atoms(structure)
        
        # Find aromatic rings
        aromatic_rings = self._find_aromatic_rings(structure)
        
        logger.info(f"Found {len(anion_atoms)} anion atoms and {len(aromatic_rings)} aromatic rings")
        
        for anion in anion_atoms:
            for ring_center, ring_residue in aromatic_rings:
                # Skip if same residue
                if anion.get_parent() == ring_residue:
                    continue
                
                distance = self._calculate_distance(anion.coord, ring_center)
                
                if distance <= self.distance_cutoff:
                    strength = self._classify_interaction_strength(distance)
                    
                    interaction = AnionPiInteraction(
                        anion_atom=anion,
                        pi_center=ring_center,
                        pi_residue_obj=ring_residue,
                        distance=distance,
                        strength=strength,
                        anion_residue=f"{anion.get_parent().get_resname()}{anion.get_parent().id[1]}",
                        pi_residue=f"{ring_residue.get_resname()}{ring_residue.id[1]}",
                        anion_chain=anion.get_parent().get_parent().id,
                        pi_chain=ring_residue.get_parent().id
                    )
                    interactions.append(interaction)
        
        logger.info(f"Detected {len(interactions)} anion-π interactions")
        return interactions
    
    def _find_anion_atoms(self, structure: Structure.Structure) -> List[Atom.Atom]:
        """Find all anionic atoms in the structure."""
        anion_atoms = []
        
        for residue in structure.get_residues():
            resname = residue.get_resname()
            if resname in self.anion_groups:
                for atom_name in self.anion_groups[resname]:
                    if atom_name in residue:
                        anion_atoms.append(residue[atom_name])
        
        return anion_atoms
    
    def _find_aromatic_rings(self, structure: Structure.Structure) -> List[Tuple[np.ndarray, Residue.Residue]]:
        """Find aromatic ring centers in the structure."""
        rings = []
        
        for residue in structure.get_residues():
            resname = residue.get_resname()
            if resname in self.aromatic_residues:
                ring_atoms = []
                for atom_name in self.aromatic_residues[resname]:
                    if atom_name in residue:
                        ring_atoms.append(residue[atom_name])
                
                if len(ring_atoms) >= 5:  # Minimum for aromatic ring
                    # Calculate ring center
                    coords = np.array([atom.coord for atom in ring_atoms])
                    center = np.mean(coords, axis=0)
                    rings.append((center, residue))
        
        return rings
    
    def _calculate_distance(self, coord1: np.ndarray, coord2: np.ndarray) -> float:
        """Calculate distance between two coordinates."""
        return np.linalg.norm(coord1 - coord2)
    
    def _classify_interaction_strength(self, distance: float) -> str:
        """Classify interaction strength based on distance."""
        if distance <= 3.5:
            return 'strong'
        elif distance <= 4.5:
            return 'moderate'
        else:
            return 'weak'
    
    def to_dict_list(self, interactions: List[AnionPiInteraction]) -> List[Dict[str, Any]]:
        """Convert interactions to list of dictionaries."""
        return [
            {
                'type': 'anion_pi',
                'anion_atom': f"{interaction.anion_atom.get_parent().get_resname()}{interaction.anion_atom.get_parent().id[1]}:{interaction.anion_atom.name}",
                'pi_residue': f"{interaction.pi_residue_obj.get_resname()}{interaction.pi_residue_obj.id[1]}",
                'distance': round(interaction.distance, 2),
                'strength': interaction.strength,
                'anion_residue': interaction.anion_residue,
                'pi_residue': interaction.pi_residue,
                'anion_chain': interaction.anion_chain,
                'pi_chain': interaction.pi_chain,
                # Standard keys for UI compatibility
                'residue1': interaction.anion_residue,
                'residue2': interaction.pi_residue,
                'chain1': interaction.anion_chain,
                'chain2': interaction.pi_chain
            }
            for interaction in interactions
        ]
