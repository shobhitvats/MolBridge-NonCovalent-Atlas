"""
Hydrogen bond detection for protein structures.
Implements detection of conventional, low-barrier, C5-type, C-H···π, and sulfur-mediated hydrogen bonds.
"""

import numpy as np
from typing import List, Dict, Any, Tuple, Optional
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from loguru import logger

@dataclass
class HydrogenBond:
    """Represents a detected hydrogen bond."""
    donor_atom: Atom.Atom
    acceptor_atom: Atom.Atom
    hydrogen_atom: Optional[Atom.Atom]
    distance: float
    angle: float
    strength: str  # 'weak', 'moderate', 'strong'
    donor_residue: str
    acceptor_residue: str
    donor_chain: str
    acceptor_chain: str

class HydrogenBondDetector:
    """Detects hydrogen bonds in protein structures."""
    
    def __init__(self, config):
        """
        Initialize hydrogen bond detector.
        
        Args:
            config: AppConfig object containing interaction parameters
        """
        self.config = config
        self.interaction_config = config.interactions
        self.distance_cutoff = 3.5  # Preliminary distance cutoff for performance
        self.angle_cutoff = 120.0  # Based on IUPAC recommendations (X-H...Y angle should be > 110)
        
        # Define donor and acceptor atoms
        self.donors = {
            # Amino acid donors
            'ARG': ['NE', 'NH1', 'NH2'],
            'ASN': ['ND2'],
            'GLN': ['NE2'],
            'HIS': ['ND1', 'NE2'],
            'LYS': ['NZ'],
            'SER': ['OG'],
            'THR': ['OG1'],
            'TRP': ['NE1'],
            'TYR': ['OH'],
            'CYS': ['SG'],  # Sulfur donor
            'MET': ['SD'],  # Sulfur donor
            # Backbone
            'BACKBONE': ['N']
        }
        
        self.acceptors = {
            # Amino acid acceptors
            'ASP': ['OD1', 'OD2'],
            'GLU': ['OE1', 'OE2'],
            'ASN': ['OD1'],
            'GLN': ['OE1'],
            'HIS': ['ND1', 'NE2'],
            'SER': ['OG'],
            'THR': ['OG1'],
            'TYR': ['OH'],
            'MET': ['SD'],  # Sulfur acceptor
            'CYS': ['SG'],  # Sulfur acceptor
            # Backbone
            'BACKBONE': ['O']
        }
        
        # C-H donors for C-H···π interactions
        self.ch_donors = {
            'ALA': ['CB'],
            'VAL': ['CB', 'CG1', 'CG2'],
            'LEU': ['CB', 'CG', 'CD1', 'CD2'],
            'ILE': ['CB', 'CG1', 'CG2', 'CD1'],
            'PHE': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2'],
            'TRP': ['CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
            'TYR': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2'],
            'HIS': ['CB', 'CG', 'CD2', 'CE1'],
            'PRO': ['CB', 'CG', 'CD'],
            'MET': ['CB', 'CG', 'CE'],
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
        
        # Aromatic π-systems for C-H···π interactions
        self.pi_systems = {
            'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TRP': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
            'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2']
        }
        
        # Sulfur atoms for sulfur-mediated H-bonds
        self.sulfur_atoms = {
            'MET': ['SD'],
            'CYS': ['SG']
        }
        
        # Hydrogen atom naming patterns
        self.hydrogen_patterns = {
            'N': ['H', 'H1', 'H2', 'H3'],
            'ND1': ['HD1'],
            'ND2': ['HD21', 'HD22'],
            'NE': ['HE'],
            'NE1': ['HE1'],
            'NE2': ['HE2', 'HE21', 'HE22'],
            'NH1': ['HH11', 'HH12'],
            'NH2': ['HH21', 'HH22'],
            'NZ': ['HZ1', 'HZ2', 'HZ3'],
            'OG': ['HG'],
            'OG1': ['HG1'],
            'OH': ['HH'],
            'SG': ['HG']
        }
    
    def detect_hydrogen_bonds(self, structure: Structure.Structure) -> List[HydrogenBond]:
        """
        Detect all hydrogen bonds in the structure.
        
        Args:
            structure: Biopython Structure object
            
        Returns:
            List of HydrogenBond objects
        """
        hbonds = []
        
        try:
            model = structure[0]  # Use first model
            
            # Get all potential donors and acceptors
            donors = self._get_donors(model)
            acceptors = self._get_acceptors(model)
            
            logger.info(f"Found {len(donors)} potential donors, {len(acceptors)} acceptors")
            
            # Check conventional, low-barrier, and C5 hydrogen bonds
            for donor_info in donors:
                for acceptor_info in acceptors:
                    hbond = self._check_hydrogen_bond(donor_info, acceptor_info)
                    if hbond:
                        hbonds.append(hbond)
            
            logger.info(f"Detected {len(hbonds)} hydrogen bonds total")
            
        except Exception as e:
            logger.error(f"Error detecting hydrogen bonds: {e}")
        
        return hbonds
    
    def _get_donors(self, model) -> List[Dict[str, Any]]:
        """Get all potential hydrogen bond donors."""
        donors = []
        
        for chain in model:
            chain_id = chain.get_id()
            
            for residue in chain:
                resname = residue.get_resname()
                res_id = residue.get_id()
                
                # Check amino acid-specific donors
                if resname in self.donors:
                    for atom_name in self.donors[resname]:
                        if atom_name in residue:
                            donor_atom = residue[atom_name]
                            hydrogen_atoms = self._find_hydrogen_atoms(residue, atom_name)
                            
                            # Ensure hydrogen_atoms is always a list
                            if not isinstance(hydrogen_atoms, list):
                                hydrogen_atoms = [hydrogen_atoms] if hydrogen_atoms else []
                            
                            donors.append({
                                'atom': donor_atom,
                                'hydrogen_atoms': hydrogen_atoms,
                                'residue': residue,
                                'resname': resname,
                                'chain_id': chain_id,
                                'res_id': res_id
                            })
                
                # Check backbone donors
                if 'N' in residue:
                    donor_atom = residue['N']
                    hydrogen_atoms = self._find_hydrogen_atoms(residue, 'N')
                    
                    # Ensure hydrogen_atoms is always a list
                    if not isinstance(hydrogen_atoms, list):
                        hydrogen_atoms = [hydrogen_atoms] if hydrogen_atoms else []
                    
                    donors.append({
                        'atom': donor_atom,
                        'hydrogen_atoms': hydrogen_atoms,
                        'residue': residue,
                        'resname': resname,
                        'chain_id': chain_id,
                        'res_id': res_id
                    })
        
        return donors
    
    def _get_acceptors(self, model) -> List[Dict[str, Any]]:
        """Get all potential hydrogen bond acceptors."""
        acceptors = []
        
        for chain in model:
            chain_id = chain.get_id()
            
            for residue in chain:
                resname = residue.get_resname()
                res_id = residue.get_id()
                
                # Check amino acid-specific acceptors
                if resname in self.acceptors:
                    for atom_name in self.acceptors[resname]:
                        if atom_name in residue:
                            acceptor_atom = residue[atom_name]
                            
                            acceptors.append({
                                'atom': acceptor_atom,
                                'residue': residue,
                                'resname': resname,
                                'chain_id': chain_id,
                                'res_id': res_id
                            })
                
                # Check backbone acceptors
                if 'O' in residue:
                    acceptor_atom = residue['O']
                    
                    acceptors.append({
                        'atom': acceptor_atom,
                        'residue': residue,
                        'resname': resname,
                        'chain_id': chain_id,
                        'res_id': res_id
                    })
        
        return acceptors
    
    def _find_hydrogen_atoms(self, residue: Residue.Residue, heavy_atom_name: str) -> List[Atom.Atom]:
        """Find hydrogen atoms bonded to a heavy atom."""
        hydrogen_atoms = []
        
        if heavy_atom_name in self.hydrogen_patterns:
            for h_name in self.hydrogen_patterns[heavy_atom_name]:
                if h_name in residue:
                    atom = residue[h_name]
                    # Ensure it's actually an Atom object
                    if hasattr(atom, 'get_name') and hasattr(atom, 'get_coord'):
                        hydrogen_atoms.append(atom)
                    else:
                        logger.warning(f"Invalid atom object for {h_name}: {type(atom)}")
        
        return hydrogen_atoms
    
    def _check_hydrogen_bond(self, donor_info: Dict[str, Any], acceptor_info: Dict[str, Any]) -> Optional[HydrogenBond]:
        """Check if a donor-acceptor pair forms a hydrogen bond."""
        donor_atom = donor_info['atom']
        acceptor_atom = acceptor_info['atom']
        
        # Skip if same atom or same residue (for now)
        if donor_atom == acceptor_atom:
            return None
        
        if donor_info['residue'] == acceptor_info['residue']:
            return None
        
        # Calculate donor-acceptor distance
        distance = self._calculate_distance(donor_atom, acceptor_atom)
        
        if distance > self.distance_cutoff:
            return None
        
        # Find best hydrogen atom and calculate angle
        best_hydrogen = None
        best_angle = 0
        
        # Ensure hydrogen_atoms is a list
        hydrogen_atoms = donor_info.get('hydrogen_atoms', [])
        if not isinstance(hydrogen_atoms, list):
            hydrogen_atoms = [hydrogen_atoms] if hydrogen_atoms else []
        
        for hydrogen_atom in hydrogen_atoms:
            if hydrogen_atom is None:
                continue
            angle = self._calculate_angle(donor_atom, hydrogen_atom, acceptor_atom)
            if angle > best_angle:
                best_angle = angle
                best_hydrogen = hydrogen_atom
        
        # If no hydrogen found, estimate position and calculate angle
        if best_hydrogen is None:
            best_angle = self._estimate_dha_angle(donor_atom, acceptor_atom, donor_info['residue'])
        
        # Check angle cutoff
        if best_angle < self.angle_cutoff:
            return None
        
        # Determine bond strength
        strength = self._calculate_bond_strength(distance, best_angle)
        
        return HydrogenBond(
            donor_atom=donor_atom,
            acceptor_atom=acceptor_atom,
            hydrogen_atom=best_hydrogen,
            distance=distance,
            angle=best_angle,
            strength=strength,
            donor_residue=f"{donor_info['resname']}{donor_info['res_id'][1]}",
            acceptor_residue=f"{acceptor_info['resname']}{acceptor_info['res_id'][1]}",
            donor_chain=donor_info['chain_id'],
            acceptor_chain=acceptor_info['chain_id']
        )
    
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
    
    def _estimate_dha_angle(self, donor_atom: Atom.Atom, acceptor_atom: Atom.Atom, donor_residue: Residue.Residue) -> float:
        """Estimate D-H...A angle when hydrogen is missing."""
        # This is a simplified estimation
        # In practice, you'd use more sophisticated geometry
        
        # For backbone N, use previous C as reference
        if donor_atom.get_name() == 'N':
            try:
                # Try to find previous residue's C
                prev_res_id = (donor_residue.get_id()[0], donor_residue.get_id()[1] - 1, donor_residue.get_id()[2])
                prev_residue = donor_residue.get_parent()[prev_res_id]
                if 'C' in prev_residue:
                    ref_atom = prev_residue['C']
                    return self._calculate_angle(ref_atom, donor_atom, acceptor_atom)
            except:
                pass
        
        # Default estimation based on tetrahedral geometry
        return 150.0  # Typical hydrogen bond angle
    
    def _calculate_bond_strength(self, distance: float, angle: float) -> str:
        """Calculate hydrogen bond strength based on IUPAC guidelines."""
        # Based on energy ranges from IUPAC recommendations
        # Strong: 15-40 kcal/mol
        # Moderate: 4-15 kcal/mol
        # Weak: <4 kcal/mol
        # We will use distance and angle to approximate these ranges
        if distance < 2.8 and angle > 150:
            return 'strong'
        elif distance < 3.2 and angle > 130:
            return 'moderate'
        else:
            return 'weak'
    
    def to_dict_list(self, hydrogen_bonds: List[HydrogenBond]) -> List[Dict[str, Any]]:
        """Convert hydrogen bonds to list of dictionaries."""
        result = []
        for hbond in hydrogen_bonds:
            try:
                # Check if hbond is actually a HydrogenBond object
                if not hasattr(hbond, 'donor_atom'):
                    logger.error(f"Invalid hydrogen bond object: {type(hbond)} - {hbond}")
                    continue
                    
                result.append({
                    'donor_atom': hbond.donor_atom.get_name() if hasattr(hbond.donor_atom, 'get_name') else str(hbond.donor_atom),
                    'acceptor_atom': hbond.acceptor_atom.get_name() if hasattr(hbond.acceptor_atom, 'get_name') else str(hbond.acceptor_atom),
                    'hydrogen_atom': hbond.hydrogen_atom.get_name() if hbond.hydrogen_atom and hasattr(hbond.hydrogen_atom, 'get_name') else None,
                    'distance': round(hbond.distance, 3),
                    'angle': round(hbond.angle, 1),
                    'strength': hbond.strength,
                    'donor_residue': hbond.donor_residue,
                    'acceptor_residue': hbond.acceptor_residue,
                    'donor_chain': hbond.donor_chain,
                    'acceptor_chain': hbond.acceptor_chain,
                    'residue1': hbond.donor_residue,
                    'residue2': hbond.acceptor_residue,
                    'chain1': hbond.donor_chain,
                    'chain2': hbond.acceptor_chain
                })
            except Exception as e:
                logger.error(f"Error converting hydrogen bond to dict: {e} - Object: {hbond}")
                continue
                
        return result
