"""
Hydrophobic contact detection for protein structures.
Detects non-polar interactions between hydrophobic residues.
"""

import numpy as np
from typing import List, Dict, Any, Optional
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from loguru import logger

@dataclass
class HydrophobicContact:
    """Represents a detected hydrophobic contact."""
    atom1: Atom.Atom
    atom2: Atom.Atom
    distance: float
    strength: str
    residue1: str
    residue2: str
    chain1: str
    chain2: str

class HydrophobicContactDetector:
    """Detects hydrophobic contacts in protein structures."""
    
    def __init__(self, config):
        self.config = config
        self.interaction_config = config.interactions
        self.distance_cutoff = self.interaction_config.hydrophobic_distance_cutoff
        
        # Hydrophobic residues and their key atoms
        self.hydrophobic_atoms = {
            'ALA': ['CB'],
            'VAL': ['CB', 'CG1', 'CG2'],
            'LEU': ['CB', 'CG', 'CD1', 'CD2'],
            'ILE': ['CB', 'CG1', 'CG2', 'CD1'],
            'PHE': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TRP': ['CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
            'TYR': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2'],
            'MET': ['CB', 'CG', 'SD', 'CE'],
            'PRO': ['CB', 'CG', 'CD']
        }
    
    def detect_hydrophobic_contacts(self, structure: Structure.Structure) -> List[HydrophobicContact]:
        """Detect hydrophobic contacts in structure."""
        contacts = []
        
        try:
            model = structure[0]
            hydrophobic_atoms = self._get_hydrophobic_atoms(model)
            
            for i, atom1_info in enumerate(hydrophobic_atoms):
                for atom2_info in hydrophobic_atoms[i+1:]:
                    contact = self._check_hydrophobic_contact(atom1_info, atom2_info)
                    if contact:
                        contacts.append(contact)
            
            logger.info(f"Detected {len(contacts)} hydrophobic contacts")
            
        except Exception as e:
            logger.error(f"Error detecting hydrophobic contacts: {e}")
        
        return contacts
    
    def _get_hydrophobic_atoms(self, model) -> List[Dict[str, Any]]:
        """Get hydrophobic atoms."""
        atoms = []
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname in self.hydrophobic_atoms:
                    for atom_name in self.hydrophobic_atoms[resname]:
                        if atom_name in residue:
                            atoms.append({
                                'atom': residue[atom_name],
                                'residue': residue,
                                'resname': resname,
                                'chain_id': chain.get_id(),
                                'res_id': residue.get_id()
                            })
        return atoms
    
    def _check_hydrophobic_contact(self, info1: Dict, info2: Dict) -> Optional[HydrophobicContact]:
        """Check if two atoms form a hydrophobic contact."""
        if info1['residue'] == info2['residue']:
            return None
        
        # Skip if residues are sequential (local contacts)
        if (info1['chain_id'] == info2['chain_id'] and 
            abs(info1['res_id'][1] - info2['res_id'][1]) < 3):
            return None
        
        distance = np.linalg.norm(
            info1['atom'].get_coord() - info2['atom'].get_coord()
        )
        
        if distance > self.distance_cutoff:
            return None
        
        strength = 'strong' if distance < 4.0 else 'moderate' if distance < 4.5 else 'weak'
        
        return HydrophobicContact(
            atom1=info1['atom'],
            atom2=info2['atom'],
            distance=distance,
            strength=strength,
            residue1=f"{info1['resname']}{info1['res_id'][1]}",
            residue2=f"{info2['resname']}{info2['res_id'][1]}",
            chain1=info1['chain_id'],
            chain2=info2['chain_id']
        )
    
    def to_dict_list(self, contacts: List[HydrophobicContact]) -> List[Dict[str, Any]]:
        """Convert contacts to list of dictionaries."""
        return [
            {
                'atom1': c.atom1.get_name(),
                'atom2': c.atom2.get_name(),
                'distance': round(c.distance, 3),
                'strength': c.strength,
                'residue1': c.residue1,
                'residue2': c.residue2,
                'chain1': c.chain1,
                'chain2': c.chain2
            }
            for c in contacts
        ]
