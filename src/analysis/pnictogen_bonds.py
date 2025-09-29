"""
Pnictogen bond detection for protein structures.
Detects interactions involving pnictogen atoms (N, P, As) as electron acceptors.
"""

import numpy as np
from utils.angle_utils import angles_between
from typing import List, Dict, Any, Optional
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from utils.settings import get_settings
from loguru import logger
from utils.instrumentation import init_funnel, finalize_funnel

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
        self.angle_cutoff = 150.0  # Based on IUPAC recommendations (R-Pn...Y angle, allows for more deviation)
        
        # Van der Waals radii for pnictogen bond calculation
        self.vdw_radii = {
            'N': 1.55,
            'P': 1.80,
            'AS': 1.85,
            'O': 1.52,
            'S': 1.80,
            'F': 1.47,
            'CL': 1.75,
            'BR': 1.85,
            'I': 1.98
        }
        
        # Pnictogen atoms (group 15 elements)
        self.pnictogen_atoms = {'N', 'P', 'AS'}
        
        # Potential acceptor atoms (electron-rich)
        self.acceptor_atoms = {
            'O', 'N', 'S', 'SE', 'F', 'CL', 'BR', 'I'
        }
    
    def detect_pnictogen_bonds(self, structure: Structure.Structure) -> List[PnictogenBond]:
        settings = get_settings()
        if getattr(settings, 'enable_vector_geom', False):
            try:
                return self._vector_detect(structure)
            except Exception:
                return self._legacy_detect(structure)
        return self._legacy_detect(structure)

    def _legacy_detect(self, structure: Structure.Structure) -> List[PnictogenBond]:
        atoms = list(structure.get_atoms())
        pnictogen_atoms = [atom for atom in atoms if atom.element in self.pnictogen_atoms]
        acceptor_atoms = [atom for atom in atoms if atom.element in self.acceptor_atoms]
        bonds: List[PnictogenBond] = []
        for pn in pnictogen_atoms:
            for ac in acceptor_atoms:
                if pn == ac or pn.get_parent() == ac.get_parent():
                    continue
                distance = float(np.linalg.norm(pn.coord - ac.coord))
                max_vdw = self.vdw_radii.get(pn.element, 0) + self.vdw_radii.get(ac.element, 0)
                if distance > max_vdw:
                    continue
                angle = self._calculate_pnictogen_angle(pn, ac)
                if angle < self.angle_cutoff:
                    continue
                strength = self._classify_bond_strength(distance)
                bonds.append(PnictogenBond(
                    pnictogen_atom=pn,
                    acceptor_atom=ac,
                    distance=distance,
                    angle=angle,
                    strength=strength,
                    pnictogen_residue=f"{pn.get_parent().get_resname()}{pn.get_parent().id[1]}",
                    acceptor_residue=f"{ac.get_parent().get_resname()}{ac.get_parent().id[1]}",
                    pnictogen_chain=pn.get_parent().get_parent().id,
                    acceptor_chain=ac.get_parent().get_parent().id
                ))
        raw_pairs = len(pnictogen_atoms)*len(acceptor_atoms)
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=raw_pairs,
            accepted_pairs=len(bonds),
            core_pair_generation=False,
            extra={
                'pnictogens': len(pnictogen_atoms),
                'acceptors': len(acceptor_atoms)
            }
        )
        # Legacy path treats pair generation + eval as a single combined evaluation phase
        finalize_funnel(self.instrumentation, pair_gen_seconds=0.0, eval_seconds=0.0, build_seconds=0.0)
        logger.info(f"[legacy] Pnictogen bonds: {len(bonds)}")
        return bonds

    def _vector_detect(self, structure: Structure.Structure) -> List[PnictogenBond]:
        bonds: List[PnictogenBond] = []
        atoms = list(structure.get_atoms())
        pnictogen_atoms = [atom for atom in atoms if atom.element in self.pnictogen_atoms]
        acceptor_atoms = [atom for atom in atoms if atom.element in self.acceptor_atoms]
        if not pnictogen_atoms or not acceptor_atoms:
            return bonds
        p_coords = np.vstack([a.get_coord() for a in pnictogen_atoms]).astype('float32')
        a_coords = np.vstack([a.get_coord() for a in acceptor_atoms]).astype('float32')
        # Use geometry core for distance candidate pruning (upper bound on varying vdW by using max combined)
        max_combined = max(self.vdw_radii.get(pn.element, 0) for pn in pnictogen_atoms) + max(self.vdw_radii.get(ac.element, 0) for ac in acceptor_atoms)
        try:
            from geometry.core import pairwise_within_cutoff
            pi, ai = pairwise_within_cutoff(p_coords, a_coords, float(max_combined), use_kdtree=True)
            core = True
        except Exception:
            diff = p_coords[:, None, :] - a_coords[None, :, :]
            dist_mat = np.linalg.norm(diff, axis=-1)
            ii, jj = np.where(dist_mat <= max_combined)
            pi, ai = ii.astype('int32'), jj.astype('int32')
            core = False
        raw_pairs = int(p_coords.shape[0] * a_coords.shape[0])
        candidate_pairs = int(len(pi))
        for idx, (p_i, a_i) in enumerate(zip(pi.tolist(), ai.tolist())):
            pn = pnictogen_atoms[p_i]; ac = acceptor_atoms[a_i]
            if pn == ac or pn.get_parent() == ac.get_parent():
                continue
            # Exact distance and vdW validation
            d = float(np.linalg.norm(pn.coord - ac.coord))
            max_allowed = self.vdw_radii.get(pn.element, 0) + self.vdw_radii.get(ac.element, 0)
            if d > max_allowed:
                continue
            angle = self._calculate_pnictogen_angle(pn, ac)
            if angle < self.angle_cutoff:
                continue
            bonds.append(PnictogenBond(
                pnictogen_atom=pn,
                acceptor_atom=ac,
                distance=d,
                angle=angle,
                strength=self._classify_bond_strength(d),
                pnictogen_residue=f"{pn.get_parent().get_resname()}{pn.get_parent().id[1]}",
                acceptor_residue=f"{ac.get_parent().get_resname()}{ac.get_parent().id[1]}",
                pnictogen_chain=pn.get_parent().get_parent().id,
                acceptor_chain=ac.get_parent().get_parent().id
            ))
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=len(bonds),
            core_pair_generation=core,
            extra={
                'pnictogens': len(pnictogen_atoms),
                'acceptors': len(acceptor_atoms),
                'candidate_density': (candidate_pairs/raw_pairs) if raw_pairs else 0.0
            }
        )
        finalize_funnel(self.instrumentation, pair_gen_seconds=0.0, eval_seconds=0.0, build_seconds=0.0)
        logger.info(f"[vector]{'/core' if core else ''} Pnictogen bonds: {len(bonds)} (raw={raw_pairs} pruned={candidate_pairs} acc_ratio={self.instrumentation['acceptance_ratio']:.3f})")
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
        """Classify bond strength based on distance.
        Note: This is an approximation as the IUPAC definition does not provide a method for this.
        """
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
