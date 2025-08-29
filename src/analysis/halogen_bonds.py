"""
Halogen bond detection for protein structures.
Detects X (Cl, Br, I) σ-hole donors with nucleophilic acceptors, X...acceptor < sum of vdW radii, angle ~170°.
"""

import numpy as np
from typing import List, Dict, Any, Optional
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from loguru import logger

@dataclass
class HalogenBond:
    """Represents a detected halogen bond."""
    halogen_atom: Atom.Atom
    acceptor_atom: Atom.Atom
    distance: float
    angle: float  # C-X...A angle (should be ~170°)
    strength: str
    halogen_residue: str
    acceptor_residue: str
    halogen_chain: str
    acceptor_chain: str
    halogen_type: str  # Cl, Br, I

class HalogenBondDetector:
    """Detects halogen bonds in protein structures."""
    
    def __init__(self, config):
        """
        Initialize halogen bond detector.
        
        Args:
            config: AppConfig object containing interaction parameters
        """
        self.config = config
        self.interaction_config = config.interactions
        self.angle_cutoff = 160.0  # Based on IUPAC recommendations (R-X...Y angle)
        
        # Van der Waals radii for halogen bond calculation
        self.vdw_radii = {
            'CL': 1.75,  # Chlorine
            'BR': 1.85,  # Bromine  
            'I': 1.98,   # Iodine
            'O': 1.52,   # Oxygen
            'N': 1.55,   # Nitrogen
            'S': 1.80    # Sulfur
        }
        
        # Halogen atoms that can form σ-hole bonds (not F due to high electronegativity)
        self.halogen_atoms = {
            'CL': ['CL', 'CLA', 'CLB', 'CLC'],  # Chlorine variants
            'BR': ['BR', 'BRA', 'BRB', 'BRC'],  # Bromine variants
            'I': ['I', 'IA', 'IB', 'IC']        # Iodine variants
        }
        
        # Nucleophilic acceptors for halogen bonds
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
    
    def detect_halogen_bonds(self, structure: Structure.Structure) -> List[HalogenBond]:
        """
        Detect all halogen bonds in the structure.
        
        Args:
            structure: Biopython Structure object
            
        Returns:
            List of HalogenBond objects
        """
        halogen_bonds = []
        
        try:
            model = structure[0]  # Use first model
            
            # Get all halogen atoms and acceptors
            halogens = self._get_halogen_atoms(model)
            acceptors = self._get_acceptors(model)
            
            logger.info(f"Found {len(halogens)} halogen atoms and {len(acceptors)} potential acceptors")
            
            # Check all halogen-acceptor pairs
            for halogen_info in halogens:
                for acceptor_info in acceptors:
                    halogen_bond = self._check_halogen_bond(halogen_info, acceptor_info)
                    if halogen_bond:
                        halogen_bonds.append(halogen_bond)
            
            logger.info(f"Detected {len(halogen_bonds)} halogen bonds")
            
        except Exception as e:
            logger.error(f"Error detecting halogen bonds: {e}")
        
        return halogen_bonds
    
    def _get_halogen_atoms(self, model) -> List[Dict[str, Any]]:
        """Get all halogen atoms in the structure."""
        halogens = []
        
        for chain in model:
            chain_id = chain.get_id()
            
            for residue in chain:
                resname = residue.get_resname()
                res_id = residue.get_id()
                
                for atom in residue:
                    atom_name = atom.get_name()
                    element = atom.get_element() if hasattr(atom, 'get_element') else atom_name[0]
                    
                    # Check if atom is a halogen - be more flexible with naming
                    halogen_type = None
                    for hal_type, hal_names in self.halogen_atoms.items():
                        if element.upper() == hal_type or atom_name.upper() in [name.upper() for name in hal_names] or atom_name.upper().startswith(hal_type):
                            halogen_type = hal_type
                            break
                    
                    # Also check for common halogen atom naming patterns
                    if not halogen_type:
                        if atom_name.upper() in ['CL', 'BR', 'I', 'F'] or element.upper() in ['CL', 'BR', 'I', 'F']:
                            if element.upper() == 'CL' or atom_name.upper().startswith('CL'):
                                halogen_type = 'CL'
                            elif element.upper() == 'BR' or atom_name.upper().startswith('BR'):
                                halogen_type = 'BR'
                            elif element.upper() == 'I' or atom_name.upper().startswith('I'):
                                halogen_type = 'I'
                            elif element.upper() == 'F' or atom_name.upper().startswith('F'):
                                halogen_type = 'F'  # Include fluorine for completeness
                    
                    if halogen_type:
                        # Find the carbon atom bonded to halogen for angle calculation
                        carbon_atom = self._find_bonded_carbon(residue, atom)
                        
                        halogens.append({
                            'atom': atom,
                            'carbon_atom': carbon_atom,
                            'halogen_type': halogen_type,
                            'residue': residue,
                            'resname': resname,
                            'chain_id': chain_id,
                            'res_id': res_id
                        })
        
        return halogens
    
    def _find_bonded_carbon(self, residue: Residue.Residue, halogen_atom: Atom.Atom) -> Optional[Atom.Atom]:
        """Find carbon atom bonded to halogen."""
        halogen_coord = halogen_atom.get_coord()
        min_distance = float('inf')
        closest_carbon = None
        
        # Look for carbon atoms in the same residue
        for atom in residue:
            if atom.get_element() == 'C':
                distance = np.linalg.norm(atom.get_coord() - halogen_coord)
                if distance < 2.0 and distance < min_distance:  # Typical C-X bond length
                    min_distance = distance
                    closest_carbon = atom
        
        return closest_carbon
    
    def _get_acceptors(self, model) -> List[Dict[str, Any]]:
        """Get all potential halogen bond acceptors."""
        acceptors = []
        
        for chain in model:
            chain_id = chain.get_id()
            
            for residue in chain:
                resname = residue.get_resname()
                res_id = residue.get_id()
                is_standard_aa = resname in self.acceptor_atoms
                
                for atom in residue:
                    atom_name = atom.get_name()
                    element = atom.get_element() if hasattr(atom, 'get_element') else atom_name[0]
                    
                    # Check if atom is a potential acceptor
                    is_acceptor = False
                    
                    if is_standard_aa and resname in self.acceptor_atoms:
                        if atom_name in self.acceptor_atoms[resname]:
                            is_acceptor = True
                    
                    # Backbone oxygen
                    if atom_name == 'O':
                        is_acceptor = True
                    
                    # Ligand atoms (non-standard residues)
                    if not is_standard_aa and element in ['O', 'N', 'S', 'F', 'CL', 'BR', 'I']:
                        is_acceptor = True
                    
                    if is_acceptor:
                        acceptors.append({
                            'atom': atom,
                            'residue': residue,
                            'resname': resname,
                            'chain_id': chain_id,
                            'res_id': res_id
                        })
        
        return acceptors
    
    def _check_halogen_bond(self, halogen_info: Dict[str, Any], acceptor_info: Dict[str, Any]) -> Optional[HalogenBond]:
        """Check if a halogen-acceptor pair forms a halogen bond."""
        halogen_atom = halogen_info['atom']
        acceptor_atom = acceptor_info['atom']
        carbon_atom = halogen_info['carbon_atom']
        
        # Skip if same atom or same residue
        if halogen_atom == acceptor_atom:
            return None
        
        if halogen_info['residue'] == acceptor_info['residue']:
            return None
        
        # Calculate halogen-acceptor distance
        distance = self._calculate_distance(halogen_atom, acceptor_atom)
        
        # Apply distance cutoff with van der Waals consideration
        vdw_sum = self._get_vdw_sum(halogen_info['halogen_type'], acceptor_atom)
        
        if distance > vdw_sum:
            return None
        
        # Calculate C-X...A angle if carbon is found
        angle = 180.0  # Default to linear
        if carbon_atom:
            angle = self._calculate_angle(carbon_atom, halogen_atom, acceptor_atom)
        
        # Check angle cutoff
        if angle < self.angle_cutoff:
            return None
        
        # Calculate bond strength
        strength = self._calculate_bond_strength(distance, angle, halogen_info['halogen_type'])
        
        return HalogenBond(
            halogen_atom=halogen_atom,
            acceptor_atom=acceptor_atom,
            distance=distance,
            angle=angle,
            strength=strength,
            halogen_residue=f"{halogen_info['resname']}{halogen_info['res_id'][1]}",
            acceptor_residue=f"{acceptor_info['resname']}{acceptor_info['res_id'][1]}",
            halogen_chain=halogen_info['chain_id'],
            acceptor_chain=acceptor_info['chain_id'],
            halogen_type=halogen_info['halogen_type']
        )
    
    def _get_vdw_sum(self, halogen_type: str, acceptor_atom: Atom.Atom) -> float:
        """Get sum of van der Waals radii."""
        halogen_vdw = self.vdw_radii.get(halogen_type, 1.8)
        acceptor_element = acceptor_atom.get_element() if hasattr(acceptor_atom, 'get_element') else acceptor_atom.get_name()[0]
        acceptor_vdw = self.vdw_radii.get(acceptor_element, 1.7)
        
        return halogen_vdw + acceptor_vdw
    
    def _calculate_distance(self, atom1: Atom.Atom, atom2: Atom.Atom) -> float:
        """Calculate distance between two atoms."""
        return np.linalg.norm(atom1.get_coord() - atom2.get_coord())
    
    def _calculate_angle(self, atom1: Atom.Atom, atom2: Atom.Atom, atom3: Atom.Atom) -> float:
        """Calculate angle between three atoms (atom2 is vertex)."""
        vec1 = atom1.get_coord() - atom2.get_coord()
        vec2 = atom3.get_coord() - atom2.get_coord()
        
        cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        
        return np.degrees(np.arccos(cos_angle))
    
    def _calculate_bond_strength(self, distance: float, angle: float, halogen_type: str) -> str:
        """Calculate halogen bond strength based on geometry. 
        Note: This is an approximation as the IUPAC definition does not provide a method for this.
        """
        # Halogen bond strength increases: F < Cl < Br < I
        strength_factors = {'F': 0.5, 'Cl': 1.0, 'Br': 1.5, 'I': 2.0}
        halogen_factor = strength_factors.get(halogen_type, 1.0)
        
        # Geometric score
        distance_score = max(0, (4.0 - distance) / 4.0)  # Shorter is better
        angle_score = max(0, (angle - 120) / 60)  # More linear is better
        
        overall_score = (distance_score + angle_score) * halogen_factor
        
        if overall_score > 1.5:
            return 'strong'
        elif overall_score > 0.8:
            return 'moderate'
        else:
            return 'weak'
    
    def get_statistics(self, halogen_bonds: List[HalogenBond]) -> Dict[str, Any]:
        """Get statistics about detected halogen bonds."""
        if not halogen_bonds:
            return {}
        
        stats = {
            'total_bonds': len(halogen_bonds),
            'halogen_types': {},
            'bond_strengths': {},
            'distance_stats': {
                'mean': np.mean([hb.distance for hb in halogen_bonds]),
                'std': np.std([hb.distance for hb in halogen_bonds]),
                'min': min([hb.distance for hb in halogen_bonds]),
                'max': max([hb.distance for hb in halogen_bonds])
            },
            'angle_stats': {
                'mean': np.mean([hb.angle for hb in halogen_bonds]),
                'std': np.std([hb.angle for hb in halogen_bonds]),
                'min': min([hb.angle for hb in halogen_bonds]),
                'max': max([hb.angle for hb in halogen_bonds])
            }
        }
        
        # Count by halogen type and strength
        for hbond in halogen_bonds:
            stats['halogen_types'][hbond.halogen_type] = stats['halogen_types'].get(hbond.halogen_type, 0) + 1
            stats['bond_strengths'][hbond.strength] = stats['bond_strengths'].get(hbond.strength, 0) + 1
        
        return stats
    
    def to_dict_list(self, halogen_bonds: List[HalogenBond]) -> List[Dict[str, Any]]:
        """Convert halogen bonds to list of dictionaries."""
        return [
            {
                'halogen_atom': hbond.halogen_atom.get_name(),
                'acceptor_atom': hbond.acceptor_atom.get_name(),
                'distance': round(hbond.distance, 3),
                'angle': round(hbond.angle, 1),
                'strength': hbond.strength,
                'halogen_type': hbond.halogen_type,
                'halogen_residue': hbond.halogen_residue,
                'acceptor_residue': hbond.acceptor_residue,
                'halogen_chain': hbond.halogen_chain,
                'acceptor_chain': hbond.acceptor_chain,
                'residue1': hbond.halogen_residue,
                'residue2': hbond.acceptor_residue,
                'chain1': hbond.halogen_chain,
                'chain2': hbond.acceptor_chain
            }
            for hbond in halogen_bonds
        ]
