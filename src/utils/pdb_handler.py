"""
PDB file handling and structure processing for Protein Interaction Explorer.
"""

import io
import warnings
import requests
from typing import Optional, Dict, Any, List, Tuple
from pathlib import Path
import numpy as np

# Suppress Bio warnings
warnings.filterwarnings("ignore")

from Bio import PDB
from Bio.PDB import Structure, Model, Chain, Residue, Atom
from Bio.PDB.PDBParser import PDBParser
try:
    from Bio.PDB.MMCIFIO import MMCIFIO
except ImportError:
    # Fallback for older BioPython versions
    try:
        from Bio.PDB import MMCIFIO
    except ImportError:
        MMCIFIO = None
from Bio.PDB.Selection import unfold_entities
try:
    import biotite.structure as struc
    import biotite.structure.io.pdb as pdb
    BIOTITE_AVAILABLE = True
except ImportError:
    BIOTITE_AVAILABLE = False
from loguru import logger

from utils.config import AppConfig

class PDBHandler:
    """Handles PDB file loading, processing, and structure manipulation."""
    
    def __init__(self, config: AppConfig):
        self.config = config
        self.parser = PDBParser(QUIET=True)
        self.cache_manager = None  # Will be set by main app
        
        # Common residue sets
        self.standard_amino_acids = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        self.standard_nucleotides = {'A', 'T', 'G', 'C', 'U', 'DA', 'DT', 'DG', 'DC'}
        
        self.water_molecules = {'HOH', 'WAT', 'H2O', 'TIP3', 'TIP4', 'SPC'}
    
    def load_structure(self, 
                      pdb_id: str, 
                      assembly: str = "biological") -> Optional[Structure.Structure]:
        """
        Load a PDB structure from file or download from RCSB.
        
        Args:
            pdb_id: 4-character PDB identifier
            assembly: "biological" or "asymmetric"
            
        Returns:
            Biopython Structure object or None if failed
        """
        try:
            # Try to get from cache first
            if hasattr(self, 'cache_manager') and self.cache_manager:
                pdb_content = self.cache_manager.get_pdb_file(pdb_id, assembly)
            else:
                pdb_content = self._download_pdb_file(pdb_id, assembly)
            
            if not pdb_content:
                logger.error(f"Failed to obtain PDB content for {pdb_id}")
                return None
            
            # Parse the structure
            pdb_io = io.StringIO(pdb_content)
            structure = self.parser.get_structure(pdb_id, pdb_io)
            
            # Process the structure based on configuration
            processed_structure = self._process_structure(structure)
            
            # Cache structure info
            if hasattr(self, 'cache_manager') and self.cache_manager:
                structure_info = self._extract_structure_info(processed_structure)
                self.cache_manager.cache_structure_info(pdb_id, structure_info)
            
            logger.info(f"Successfully loaded structure {pdb_id}")
            return processed_structure
            
        except Exception as e:
            logger.error(f"Failed to load structure {pdb_id}: {e}")
            return None
    
    def get_pdb_content(self, pdb_id: str, assembly: str = "biological") -> Optional[str]:
        """Get raw PDB content as string for visualization.
        
        Args:
            pdb_id: 4-character PDB identifier
            assembly: "biological" or "asymmetric"
            
        Returns:
            Raw PDB content as string or None if failed
        """
        try:
            # Try to get from cache first
            if hasattr(self, 'cache_manager') and self.cache_manager:
                pdb_content = self.cache_manager.get_pdb_file(pdb_id, assembly)
            else:
                pdb_content = self._download_pdb_file(pdb_id, assembly)
            
            if not pdb_content:
                logger.error(f"Failed to obtain PDB content for {pdb_id}")
                return None
                
            return pdb_content
            
        except Exception as e:
            logger.error(f"Failed to get PDB content for {pdb_id}: {e}")
            return None
    
    def _download_pdb_file(self, pdb_id: str, assembly: str) -> Optional[str]:
        """Download PDB file from RCSB."""
        try:
            if assembly == "biological":
                # Try biological assembly first
                url = f"https://files.rcsb.org/download/{pdb_id}.pdb1"
                response = requests.get(url, timeout=30)
                
                if response.status_code == 404:
                    # Fallback to asymmetric unit
                    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                    response = requests.get(url, timeout=30)
            else:
                url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                response = requests.get(url, timeout=30)
            
            response.raise_for_status()
            return response.text
            
        except requests.RequestException as e:
            logger.error(f"Failed to download {pdb_id}: {e}")
            return None
    
    def _process_structure(self, structure: Structure.Structure) -> Structure.Structure:
        """Process structure according to configuration settings."""
        # Remove waters if configured
        if self.config.exclude_waters:
            self._remove_waters(structure)
        
        # Handle ALTLOC conformers
        self._handle_altloc(structure)
        
        # Add hydrogens if missing (simplified approach)
        if self.config.processing.max_workers > 0:  # Use as flag for hydrogen addition
            self._add_missing_hydrogens(structure)
        
        return structure
    
    def _remove_waters(self, structure: Structure.Structure):
        """Remove water molecules from structure."""
        for model in structure:
            for chain in list(model):
                for residue in list(chain):
                    if residue.get_resname().strip() in self.water_molecules:
                        chain.detach_child(residue.get_id())
    
    def _handle_altloc(self, structure: Structure.Structure):
        """Handle alternative location conformers."""
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in list(residue):
                        # Keep highest occupancy conformer or first if equal
                        if atom.get_altloc() != ' ':
                            # Find all altloc variants
                            atom_name = atom.get_name()
                            altloc_atoms = [a for a in residue if a.get_name() == atom_name]
                            
                            if len(altloc_atoms) > 1:
                                # Keep highest occupancy
                                best_atom = max(altloc_atoms, key=lambda a: a.get_occupancy())
                                
                                # Remove others
                                for alt_atom in altloc_atoms:
                                    if alt_atom != best_atom:
                                        residue.detach_child(alt_atom.get_id())
                                
                                # Clear altloc identifier
                                best_atom.set_altloc(' ')
    
    def _add_missing_hydrogens(self, structure: Structure.Structure):
        """Add missing hydrogen atoms (simplified heuristic approach)."""
        # This is a simplified implementation
        # In production, you might want to use more sophisticated tools like OpenEye or RDKit
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in self.standard_amino_acids:
                        self._add_amino_acid_hydrogens(residue)
    
    def _add_amino_acid_hydrogens(self, residue: Residue.Residue):
        """Add basic hydrogen atoms to amino acid residues."""
        # This is a very simplified implementation
        # Real hydrogen addition would require proper geometry and force field calculations
        
        resname = residue.get_resname()
        
        # Add backbone amide hydrogen
        if 'N' in residue and 'H' not in residue:
            n_atom = residue['N']
            # Simple geometric placement (would need proper calculation)
            h_coord = n_atom.get_coord() + np.array([0.0, 0.0, 1.0])
            
            h_atom = Atom.Atom('H', h_coord, 0.0, 1.0, ' ', 'H', 1)
            residue.add(h_atom)
    
    def _extract_structure_info(self, structure: Structure.Structure) -> Dict[str, Any]:
        """Extract metadata and summary information from structure."""
        info = {
            "pdb_id": structure.get_id(),
            "models": len(structure),
            "chains": [],
            "total_residues": 0,
            "total_atoms": 0,
            "resolution": None,
            "organism": None,
            "method": None
        }
        
        # Count chains and residues
        for model in structure:
            for chain in model:
                chain_info = {
                    "chain_id": chain.get_id(),
                    "residues": len(chain),
                    "atoms": sum(len(residue) for residue in chain),
                    "sequence_length": len([r for r in chain if r.get_resname() in self.standard_amino_acids])
                }
                info["chains"].append(chain_info)
                info["total_residues"] += chain_info["residues"]
                info["total_atoms"] += chain_info["atoms"]
        
        return info
    
    def get_residue_atoms(self, 
                         structure: Structure.Structure, 
                         chain_id: str, 
                         residue_id: int) -> List[Atom.Atom]:
        """Get all atoms for a specific residue."""
        try:
            model = structure[0]  # Use first model
            chain = model[chain_id]
            residue = chain[residue_id]
            return list(residue.get_atoms())
        except KeyError:
            return []
    
    def get_chain_sequence(self, 
                          structure: Structure.Structure, 
                          chain_id: str) -> str:
        """Get amino acid sequence for a chain."""
        try:
            model = structure[0]
            chain = model[chain_id]
            
            # Three-letter to one-letter mapping
            aa_map = {
                'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
            }
            
            sequence = ""
            for residue in chain:
                if residue.get_resname() in aa_map:
                    sequence += aa_map[residue.get_resname()]
                else:
                    sequence += 'X'  # Unknown residue
            
            return sequence
            
        except KeyError:
            return ""
    
    def get_ligands(self, structure: Structure.Structure) -> List[Dict[str, Any]]:
        """Get all ligand molecules from structure."""
        ligands = []
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip()
                    
                    # Skip standard amino acids, nucleotides, and waters
                    if (resname not in self.standard_amino_acids and 
                        resname not in self.standard_nucleotides and
                        resname not in self.water_molecules):
                        
                        ligand_info = {
                            "name": resname,
                            "chain": chain.get_id(),
                            "residue_id": residue.get_id(),
                            "atoms": len(residue),
                            "center": self._calculate_center_of_mass(residue)
                        }
                        ligands.append(ligand_info)
        
        return ligands
    
    def _calculate_center_of_mass(self, residue: Residue.Residue) -> np.ndarray:
        """Calculate center of mass for a residue."""
        atoms = list(residue.get_atoms())
        if not atoms:
            return np.array([0.0, 0.0, 0.0])
        
        coords = np.array([atom.get_coord() for atom in atoms])
        return np.mean(coords, axis=0)
    
    def get_atom_coordinates(self, structure: Structure.Structure) -> Dict[str, np.ndarray]:
        """Get coordinates of all atoms organized by chain and residue."""
        coordinates = {}
        
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                coordinates[chain_id] = {}
                
                for residue in chain:
                    residue_id = residue.get_id()
                    coordinates[chain_id][residue_id] = {}
                    
                    for atom in residue:
                        atom_name = atom.get_name()
                        coordinates[chain_id][residue_id][atom_name] = atom.get_coord()
        
        return coordinates
    
    def calculate_distance(self, atom1: Atom.Atom, atom2: Atom.Atom) -> float:
        """Calculate distance between two atoms."""
        return np.linalg.norm(atom1.get_coord() - atom2.get_coord())
    
    def calculate_angle(self, 
                       atom1: Atom.Atom, 
                       atom2: Atom.Atom, 
                       atom3: Atom.Atom) -> float:
        """Calculate angle between three atoms (atom2 is vertex)."""
        vec1 = atom1.get_coord() - atom2.get_coord()
        vec2 = atom3.get_coord() - atom2.get_coord()
        
        cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Handle numerical errors
        
        return np.degrees(np.arccos(cos_angle))
    
    def get_nearby_atoms(self, 
                        structure: Structure.Structure,
                        center_atom: Atom.Atom,
                        radius: float = 5.0) -> List[Tuple[Atom.Atom, float]]:
        """Get all atoms within radius of center atom."""
        nearby_atoms = []
        center_coord = center_atom.get_coord()
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom != center_atom:
                            distance = np.linalg.norm(atom.get_coord() - center_coord)
                            if distance <= radius:
                                nearby_atoms.append((atom, distance))
        
        return sorted(nearby_atoms, key=lambda x: x[1])
    
    def save_structure(self, 
                      structure: Structure.Structure, 
                      filename: str,
                      format: str = "pdb") -> bool:
        """Save structure to file."""
        try:
            if format.lower() == "pdb":
                io_handler = PDB.PDBIO()
                io_handler.set_structure(structure)
                io_handler.save(filename)
            elif format.lower() == "cif":
                io_handler = MMCIFIO()
                io_handler.set_structure(structure)
                io_handler.save(filename)
            else:
                raise ValueError(f"Unsupported format: {format}")
            
            return True
            
        except Exception as e:
            logger.error(f"Failed to save structure: {e}")
            return False
    
    def validate_structure(self, structure: Structure.Structure) -> Dict[str, Any]:
        """Validate structure and return quality metrics."""
        validation = {
            "valid": True,
            "warnings": [],
            "errors": [],
            "statistics": {}
        }
        
        try:
            # Basic structure checks
            if len(structure) == 0:
                validation["errors"].append("No models found in structure")
                validation["valid"] = False
                return validation
            
            model = structure[0]
            if len(model) == 0:
                validation["errors"].append("No chains found in structure")
                validation["valid"] = False
                return validation
            
            # Check for missing atoms
            missing_atoms = 0
            total_residues = 0
            
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in self.standard_amino_acids:
                        total_residues += 1
                        required_atoms = ['N', 'CA', 'C', 'O']
                        
                        for atom_name in required_atoms:
                            if atom_name not in residue:
                                missing_atoms += 1
            
            if missing_atoms > 0:
                validation["warnings"].append(f"Missing {missing_atoms} backbone atoms")
            
            # Calculate statistics
            validation["statistics"] = {
                "total_residues": total_residues,
                "missing_atoms": missing_atoms,
                "completeness": 1 - (missing_atoms / (total_residues * 4)) if total_residues > 0 else 0
            }
            
        except Exception as e:
            validation["errors"].append(f"Validation error: {e}")
            validation["valid"] = False
        
        return validation
