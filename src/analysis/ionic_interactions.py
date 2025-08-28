"""
Ionic interaction detection for protein structures.
Detects electrostatic interactions between charged residues.
"""

import numpy as np
from typing import List, Dict, Any, Optional
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from loguru import logger

@dataclass
class IonicInteraction:
    """Represents a detected ionic interaction."""
    positive_atom: Atom.Atom
    negative_atom: Atom.Atom
    distance: float
    strength: str
    positive_residue: str
    negative_residue: str
    positive_chain: str
    negative_chain: str

class IonicInteractionDetector:
    """Detects ionic interactions in protein structures."""
    
    def __init__(self, config):
        self.config = config
        self.interaction_config = config.interactions
        self.distance_cutoff = self.interaction_config.ionic_distance_cutoff
        
        # Charged residues
        self.positive_residues = {
            'ARG': ['NH1', 'NH2'],
            'LYS': ['NZ'],
            'HIS': ['ND1', 'NE2']  # Protonated form
        }
        
        self.negative_residues = {
            'ASP': ['OD1', 'OD2'],
            'GLU': ['OE1', 'OE2']
        }
    
    def detect_ionic_interactions(self, structure: Structure.Structure) -> List[IonicInteraction]:
        """Detect ionic interactions in structure."""
        interactions = []
        
        try:
            model = structure[0]
            
            positive_atoms = self._get_positive_atoms(model)
            negative_atoms = self._get_negative_atoms(model)
            
            for pos_info in positive_atoms:
                for neg_info in negative_atoms:
                    interaction = self._check_ionic_interaction(pos_info, neg_info)
                    if interaction:
                        interactions.append(interaction)
            
            logger.info(f"Detected {len(interactions)} ionic interactions")
            
        except Exception as e:
            logger.error(f"Error detecting ionic interactions: {e}")
        
        return interactions
    
    def _get_positive_atoms(self, model) -> List[Dict[str, Any]]:
        """Get positively charged atoms."""
        atoms = []
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname in self.positive_residues:
                    for atom_name in self.positive_residues[resname]:
                        if atom_name in residue:
                            atoms.append({
                                'atom': residue[atom_name],
                                'residue': residue,
                                'resname': resname,
                                'chain_id': chain.get_id(),
                                'res_id': residue.get_id()
                            })
        return atoms
    
    def _get_negative_atoms(self, model) -> List[Dict[str, Any]]:
        """Get negatively charged atoms."""
        atoms = []
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname in self.negative_residues:
                    for atom_name in self.negative_residues[resname]:
                        if atom_name in residue:
                            atoms.append({
                                'atom': residue[atom_name],
                                'residue': residue,
                                'resname': resname,
                                'chain_id': chain.get_id(),
                                'res_id': residue.get_id()
                            })
        return atoms
    
    def _check_ionic_interaction(self, pos_info: Dict, neg_info: Dict) -> Optional[IonicInteraction]:
        """Check if two charged atoms form an ionic interaction."""
        if pos_info['residue'] == neg_info['residue']:
            return None
        
        distance = np.linalg.norm(
            pos_info['atom'].get_coord() - neg_info['atom'].get_coord()
        )
        
        if distance > self.distance_cutoff:
            return None
        
        strength = 'strong' if distance < 4.0 else 'moderate' if distance < 5.0 else 'weak'
        
        return IonicInteraction(
            positive_atom=pos_info['atom'],
            negative_atom=neg_info['atom'],
            distance=distance,
            strength=strength,
            positive_residue=f"{pos_info['resname']}{pos_info['res_id'][1]}",
            negative_residue=f"{neg_info['resname']}{neg_info['res_id'][1]}",
            positive_chain=pos_info['chain_id'],
            negative_chain=neg_info['chain_id']
        )
    
    def to_dict_list(self, interactions: List[IonicInteraction]) -> List[Dict[str, Any]]:
        """Convert interactions to list of dictionaries."""
        return [
            {
                'positive_atom': i.positive_atom.get_name(),
                'negative_atom': i.negative_atom.get_name(),
                'distance': round(i.distance, 3),
                'strength': i.strength,
                'positive_residue': i.positive_residue,
                'negative_residue': i.negative_residue,
                'positive_chain': i.positive_chain,
                'negative_chain': i.negative_chain,
                'residue1': i.positive_residue,
                'residue2': i.negative_residue,
                'chain1': i.positive_chain,
                'chain2': i.negative_chain
            }
            for i in interactions
        ]
