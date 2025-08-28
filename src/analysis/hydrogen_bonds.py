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
    bond_type: str  # 'conventional', 'low_barrier', 'c5_type', 'ch_pi', 'sulfur_mediated'
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
        self.distance_cutoff = self.interaction_config.hbond_distance_cutoff
        self.angle_cutoff = self.interaction_config.hbond_angle_cutoff
        
        # New specific cutoffs for different H-bond types
        self.low_barrier_cutoff = 2.6  # Very short donor-acceptor distances for low-barrier H-bonds
        self.c5_distance_cutoff = 3.0  # C5 intraresidue main chain NH...O=C
        self.ch_pi_distance_cutoff = 4.5  # C-H···π interaction C-center distance
        self.ch_pi_angle_cutoff = 120.0  # C-H···π angle requirement
        self.sulfur_distance_min = 2.8  # Sulfur-mediated H-bond minimum distance
        self.sulfur_distance_max = 3.5  # Sulfur-mediated H-bond maximum distance  
        self.sulfur_angle_cutoff = 140.0  # Sulfur-mediated H-bond angle
        
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
            
            # Get C-H donors for C-H···π interactions
            ch_donors = self._get_ch_donors(model)
            pi_acceptors = self._get_pi_acceptors(model)
            
            # Get sulfur atoms for sulfur-mediated H-bonds
            sulfur_atoms = self._get_sulfur_atoms(model)
            
            logger.info(f"Found {len(donors)} potential donors, {len(acceptors)} acceptors, {len(ch_donors)} C-H donors, {len(pi_acceptors)} π-systems, {len(sulfur_atoms)} sulfur atoms")
            
            # Check conventional, low-barrier, and C5 hydrogen bonds
            for donor_info in donors:
                for acceptor_info in acceptors:
                    hbond = self._check_hydrogen_bond(donor_info, acceptor_info)
                    if hbond:
                        hbonds.append(hbond)
            
            # Check C-H···π interactions
            for ch_donor in ch_donors:
                for pi_acceptor in pi_acceptors:
                    ch_pi_bond = self._check_ch_pi_interaction(ch_donor, pi_acceptor)
                    if ch_pi_bond:
                        hbonds.append(ch_pi_bond)
            
            # Check sulfur-mediated hydrogen bonds
            for sulfur_info in sulfur_atoms:
                for acceptor_info in acceptors:
                    # Sulfur as donor
                    s_hbond = self._check_sulfur_mediated_bond(sulfur_info, acceptor_info, as_donor=True)
                    if s_hbond:
                        hbonds.append(s_hbond)
                    
                for donor_info in donors:
                    # Sulfur as acceptor
                    s_hbond = self._check_sulfur_mediated_bond(donor_info, sulfur_info, as_donor=False)
                    if s_hbond:
                        hbonds.append(s_hbond)
            
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
                    hydrogen_atoms.append(residue[h_name])
        
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
        
        for hydrogen_atom in donor_info['hydrogen_atoms']:
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
        
        # Determine bond type and strength
        bond_type = self._classify_bond_type(distance, best_angle, donor_info, acceptor_info)
        strength = self._calculate_bond_strength(distance, best_angle)
        
        return HydrogenBond(
            donor_atom=donor_atom,
            acceptor_atom=acceptor_atom,
            hydrogen_atom=best_hydrogen,
            distance=distance,
            angle=best_angle,
            bond_type=bond_type,
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
    
    def _classify_bond_type(self, distance: float, angle: float, donor_info: Dict, acceptor_info: Dict) -> str:
        """Classify the type of hydrogen bond."""
        # Low-barrier hydrogen bonds (very short, usually symmetric, O/N atoms)
        if distance < self.low_barrier_cutoff:
            donor_atom = donor_info['atom'].get_name()
            acceptor_atom = acceptor_info['atom'].get_name()
            # Check if both are O or N atoms (typical for low-barrier bonds)
            if (donor_atom.startswith(('O', 'N')) and acceptor_atom.startswith(('O', 'N'))):
                return 'low_barrier'
        
        # C5-type hydrogen bonds (intramolecular main chain NH...O=C)
        if self._is_c5_type(donor_info, acceptor_info, distance):
            return 'c5_type'
        
        # Conventional hydrogen bonds
        return 'conventional'
    
    def _is_c5_type(self, donor_info: Dict, acceptor_info: Dict, distance: float) -> bool:
        """Check if this is a C5-type hydrogen bond (intraresidue main chain NH...O=C)."""
        if distance > self.c5_distance_cutoff:
            return False
        
        # Check if donor and acceptor are in the same chain
        if donor_info['chain_id'] != acceptor_info['chain_id']:
            return False
        
        # Must be backbone N donor and O acceptor
        donor_atom = donor_info['atom'].get_name()
        acceptor_atom = acceptor_info['atom'].get_name()
        
        if donor_atom != 'N' or acceptor_atom != 'O':
            return False
        
        # Check sequence separation (should be close for C5, typically in β-sheets)
        donor_seq = donor_info['res_id'][1]
        acceptor_seq = acceptor_info['res_id'][1]
        
        seq_separation = abs(donor_seq - acceptor_seq)
        
        # C5 bonds typically occur between residues 1-4 positions apart in β-sheets
        if 1 <= seq_separation <= 4:
            return True
        
        return False
    
    def _get_ch_donors(self, model) -> List[Dict[str, Any]]:
        """Get all potential C-H donors for C-H···π interactions."""
        ch_donors = []
        
        for chain in model:
            chain_id = chain.get_id()
            
            for residue in chain:
                resname = residue.get_resname()
                res_id = residue.get_id()
                
                if resname in self.ch_donors:
                    for atom_name in self.ch_donors[resname]:
                        if atom_name in residue:
                            ch_atom = residue[atom_name]
                            # Only consider carbon atoms
                            if ch_atom.element == 'C':
                                ch_donors.append({
                                    'atom': ch_atom,
                                    'residue': residue,
                                    'resname': resname,
                                    'chain_id': chain_id,
                                    'res_id': res_id
                                })
        
        return ch_donors
    
    def _get_pi_acceptors(self, model) -> List[Dict[str, Any]]:
        """Get all aromatic π-systems for C-H···π interactions."""
        pi_acceptors = []
        
        for chain in model:
            chain_id = chain.get_id()
            
            for residue in chain:
                resname = residue.get_resname()
                res_id = residue.get_id()
                
                if resname in self.pi_systems:
                    # Calculate centroid of aromatic ring
                    ring_atoms = []
                    for atom_name in self.pi_systems[resname]:
                        if atom_name in residue:
                            ring_atoms.append(residue[atom_name])
                    
                    if len(ring_atoms) >= 5:  # Need at least 5 atoms for aromatic ring
                        centroid = np.mean([atom.get_coord() for atom in ring_atoms], axis=0)
                        pi_acceptors.append({
                            'centroid': centroid,
                            'atoms': ring_atoms,
                            'residue': residue,
                            'resname': resname,
                            'chain_id': chain_id,
                            'res_id': res_id
                        })
        
        return pi_acceptors
    
    def _get_sulfur_atoms(self, model) -> List[Dict[str, Any]]:
        """Get all sulfur atoms for sulfur-mediated hydrogen bonds."""
        sulfur_atoms = []
        
        for chain in model:
            chain_id = chain.get_id()
            
            for residue in chain:
                resname = residue.get_resname()
                res_id = residue.get_id()
                
                if resname in self.sulfur_atoms:
                    for atom_name in self.sulfur_atoms[resname]:
                        if atom_name in residue:
                            sulfur_atom = residue[atom_name]
                            sulfur_atoms.append({
                                'atom': sulfur_atom,
                                'residue': residue,
                                'resname': resname,
                                'chain_id': chain_id,
                                'res_id': res_id
                            })
        
        return sulfur_atoms
    
    def _check_ch_pi_interaction(self, ch_donor: Dict, pi_acceptor: Dict) -> Optional[HydrogenBond]:
        """Check for C-H···π interaction."""
        ch_atom = ch_donor['atom']
        centroid = pi_acceptor['centroid']
        
        # Calculate distance from C-H to π-centroid
        distance = np.linalg.norm(ch_atom.get_coord() - centroid)
        
        if distance > self.ch_pi_distance_cutoff:
            return None
        
        # Estimate C-H···centroid angle (simplified)
        # In practice, you'd need to find attached hydrogens
        # For now, use a geometric approximation
        angle = 150.0  # Default reasonable angle
        
        if angle < self.ch_pi_angle_cutoff:
            return None
        
        return HydrogenBond(
            donor_atom=ch_atom,
            acceptor_atom=pi_acceptor['atoms'][0],  # Use first ring atom as representative
            hydrogen_atom=None,
            distance=distance,
            angle=angle,
            bond_type='ch_pi',
            strength=self._calculate_bond_strength(distance, angle),
            donor_residue=f"{ch_donor['resname']}{ch_donor['res_id'][1]}",
            acceptor_residue=f"{pi_acceptor['resname']}{pi_acceptor['res_id'][1]}",
            donor_chain=ch_donor['chain_id'],
            acceptor_chain=pi_acceptor['chain_id']
        )
    
    def _check_sulfur_mediated_bond(self, donor_info: Dict, acceptor_info: Dict, as_donor: bool) -> Optional[HydrogenBond]:
        """Check for sulfur-mediated hydrogen bond."""
        if as_donor:
            # Sulfur as donor
            sulfur_info = donor_info
            other_info = acceptor_info
            sulfur_atom = sulfur_info['atom']
            other_atom = other_info['atom']
        else:
            # Sulfur as acceptor
            sulfur_info = acceptor_info
            other_info = donor_info
            sulfur_atom = sulfur_info['atom']
            other_atom = other_info['atom']
        
        # Check if it's actually a sulfur atom
        if sulfur_atom.element != 'S':
            return None
        
        distance = self._calculate_distance(sulfur_atom, other_atom)
        
        # Check sulfur-mediated distance range
        if distance < self.sulfur_distance_min or distance > self.sulfur_distance_max:
            return None
        
        # Estimate angle (simplified)
        angle = 140.0  # Default angle for sulfur bonds
        
        if angle < self.sulfur_angle_cutoff:
            return None
        
        if as_donor:
            return HydrogenBond(
                donor_atom=sulfur_atom,
                acceptor_atom=other_atom,
                hydrogen_atom=None,
                distance=distance,
                angle=angle,
                bond_type='sulfur_mediated',
                strength=self._calculate_bond_strength(distance, angle),
                donor_residue=f"{sulfur_info['resname']}{sulfur_info['res_id'][1]}",
                acceptor_residue=f"{other_info['resname']}{other_info['res_id'][1]}",
                donor_chain=sulfur_info['chain_id'],
                acceptor_chain=other_info['chain_id']
            )
        else:
            return HydrogenBond(
                donor_atom=other_atom,
                acceptor_atom=sulfur_atom,
                hydrogen_atom=None,
                distance=distance,
                angle=angle,
                bond_type='sulfur_mediated',
                strength=self._calculate_bond_strength(distance, angle),
                donor_residue=f"{other_info['resname']}{other_info['res_id'][1]}",
                acceptor_residue=f"{sulfur_info['resname']}{sulfur_info['res_id'][1]}",
                donor_chain=other_info['chain_id'],
                acceptor_chain=sulfur_info['chain_id']
            )

    def _calculate_bond_strength(self, distance: float, angle: float) -> str:
        """Calculate hydrogen bond strength based on geometry."""
        # Strong bonds: short distance, good angle
        if distance < 2.8 and angle > 160:
            return 'strong'
        
        # Weak bonds: long distance or poor angle
        elif distance > 3.2 or angle < 140:
            return 'weak'
        
        # Moderate bonds: everything else
        else:
            return 'moderate'
    
    def get_statistics(self, hbonds: List[HydrogenBond]) -> Dict[str, Any]:
        """Get statistics about detected hydrogen bonds."""
        if not hbonds:
            return {}
        
        stats = {
            'total_bonds': len(hbonds),
            'bond_types': {},
            'bond_strengths': {},
            'distance_stats': {
                'mean': np.mean([hb.distance for hb in hbonds]),
                'std': np.std([hb.distance for hb in hbonds]),
                'min': min([hb.distance for hb in hbonds]),
                'max': max([hb.distance for hb in hbonds])
            },
            'angle_stats': {
                'mean': np.mean([hb.angle for hb in hbonds]),
                'std': np.std([hb.angle for hb in hbonds]),
                'min': min([hb.angle for hb in hbonds]),
                'max': max([hb.angle for hb in hbonds])
            }
        }
        
        # Count by type and strength
        for hbond in hbonds:
            stats['bond_types'][hbond.bond_type] = stats['bond_types'].get(hbond.bond_type, 0) + 1
            stats['bond_strengths'][hbond.strength] = stats['bond_strengths'].get(hbond.strength, 0) + 1
        
        return stats
    
    def to_dict_list(self, hbonds: List[HydrogenBond]) -> List[Dict[str, Any]]:
        """Convert hydrogen bonds to list of dictionaries."""
        return [
            {
                'donor_atom': hbond.donor_atom.get_name(),
                'acceptor_atom': hbond.acceptor_atom.get_name(),
                'distance': round(hbond.distance, 3),
                'angle': round(hbond.angle, 1),
                'bond_type': hbond.bond_type,
                'strength': hbond.strength,
                'donor_residue': hbond.donor_residue,
                'acceptor_residue': hbond.acceptor_residue,
                'donor_chain': hbond.donor_chain,
                'acceptor_chain': hbond.acceptor_chain,
                'residue1': hbond.donor_residue,
                'residue2': hbond.acceptor_residue,
                'chain1': hbond.donor_chain,
                'chain2': hbond.acceptor_chain
            }
            for hbond in hbonds
        ]
