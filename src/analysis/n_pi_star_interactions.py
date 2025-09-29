"""
n→π* interaction detection for protein structures.
Detects interactions between lone pairs and carbonyl π* orbitals with O...C ~3.0 Å, angle ~102°.
"""

import numpy as np
from typing import List, Dict, Any, Optional
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from loguru import logger
from utils.settings import get_settings
from utils.instrumentation import init_funnel, update_counts, finalize_funnel
from geometry.core import pairwise_within_cutoff, norms

@dataclass
class NPiStarInteraction:
    """Represents a detected n→π* interaction."""
    donor_atom: Atom.Atom
    acceptor_atom: Atom.Atom
    carbonyl_carbon: Atom.Atom
    distance: float
    angle: float
    strength: str
    donor_residue: str
    acceptor_residue: str
    donor_chain: str
    acceptor_chain: str

class NPiStarDetector:
    """Detects n→π* interactions in protein structures."""
    
    def __init__(self, config):
        """
        Initialize n→π* detector.
        
        Args:
            config: AppConfig object containing interaction parameters
        """
        self.config = config
        self.interaction_config = config.interactions
        self.distance_cutoff = 3.0  # O...C distance ~3.0 Å for n→π* interactions
        self.optimal_angle = 102.0  # Optimal Bürgi-Dunitz angle ~102°
        self.angle_tolerance = 15.0  # Allow ±15° around optimal angle
        
        # Lone pair donors (oxygen atoms with lone pairs)
        self.lone_pair_donors = {
            'ASP': ['OD1', 'OD2'],    # Aspartate carboxylate oxygens
            'GLU': ['OE1', 'OE2'],    # Glutamate carboxylate oxygens
            'SER': ['OG'],            # Serine hydroxyl oxygen
            'THR': ['OG1'],           # Threonine hydroxyl oxygen
            'TYR': ['OH'],            # Tyrosine hydroxyl oxygen
            'backbone': ['O']          # Backbone carbonyl oxygen
        }
        
        # Carbonyl acceptors (C=O groups with π* orbitals)
        self.carbonyl_acceptors = {
            'backbone': ['C', 'O'],    # Backbone carbonyls
            'ASN': ['CG', 'OD1'],      # Asparagine amide carbonyl
            'GLN': ['CD', 'OE1'],      # Glutamine amide carbonyl
            'ASP': ['CG', 'OD1'],      # Aspartate carboxyl (when protonated)
            'GLU': ['CD', 'OE1']       # Glutamate carboxyl (when protonated)
        }
    
    def detect_n_pi_star_interactions(self, structure: Structure.Structure) -> List[NPiStarInteraction]:
        """
        Detect n→π* interactions in structure.
        Uses Bürgi-Dunitz angle (~102°) and O...C distance ~3.0 Å.
        """
        settings = get_settings()
        if getattr(settings, 'enable_vector_geom', False):
            try:
                return self._vector_detect(structure)
            except Exception:  # pragma: no cover
                pass
        interactions: List[NPiStarInteraction] = []
        raw_pairs = 0
        candidate_pairs = 0
        import time as _t
        t_pair_start = _t.time()
        
        try:
            model = structure[0]  # Use first model
            
            # Find lone pair donor atoms
            donor_atoms = self._find_lone_pair_donors(model)
            
            # Find carbonyl groups (π* acceptors)
            carbonyl_groups = self._find_carbonyl_groups(model)
            
            logger.info(f"Found {len(donor_atoms)} lone pair donors and {len(carbonyl_groups)} carbonyl groups")
            
            raw_pairs = len(donor_atoms) * len(carbonyl_groups)
            for donor_info in donor_atoms:
                for carbonyl_info in carbonyl_groups:
                    # Skip if same residue
                    if donor_info['residue'] == carbonyl_info['residue']:
                        continue
                    
                    # Calculate distance from donor oxygen to carbonyl carbon
                    donor_atom = donor_info['atom']
                    carbon_atom = carbonyl_info['carbon']
                    oxygen_atom = carbonyl_info['oxygen']
                    
                    distance = self._calculate_distance(donor_atom, carbon_atom)
                    
                    if distance <= self.distance_cutoff:
                        # Calculate Bürgi-Dunitz angle (lone pair approaching antiparallel to C=O)
                        # This is the angle between donor...carbon vector and carbon-oxygen bond
                        angle = self._calculate_burgi_dunitz_angle(donor_atom, carbon_atom, oxygen_atom)
                        
                        # Check if angle is close to optimal 102° (within tolerance)
                        if abs(angle - self.optimal_angle) <= self.angle_tolerance:
                            candidate_pairs += 1
                            strength = self._classify_interaction_strength(distance, angle)
                            
                            interaction = NPiStarInteraction(
                                donor_atom=donor_atom,
                                acceptor_atom=oxygen_atom,  # π* acceptor oxygen
                                carbonyl_carbon=carbon_atom,
                                distance=distance,
                                angle=angle,
                                strength=strength,
                                donor_residue=f"{donor_info['resname']}{donor_info['res_id'][1]}",
                                acceptor_residue=f"{carbonyl_info['resname']}{carbonyl_info['res_id'][1]}",
                                donor_chain=donor_info['chain_id'],
                                acceptor_chain=carbonyl_info['chain_id']
                            )
                            interactions.append(interaction)
            
            logger.info(f"Detected {len(interactions)} n→π* interactions")
            
        except Exception as e:
            logger.error(f"Error detecting n→π* interactions: {e}")
        
        t_eval_end = _t.time()
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=len(interactions),
            core_pair_generation=False,
            extra={
                'donors': len(donor_atoms) if 'donor_atoms' in locals() else 0,
                'carbonyl_groups': len(carbonyl_groups) if 'carbonyl_groups' in locals() else 0,
                'candidate_density': (candidate_pairs/raw_pairs) if raw_pairs else 0.0
            }
        )
        finalize_funnel(
            self.instrumentation,
            pair_gen_seconds=0.0,
            eval_seconds=(t_eval_end - t_pair_start),
            build_seconds=0.0
        )
        return interactions

    def _vector_detect(self, structure: Structure.Structure) -> List[NPiStarInteraction]:
        interactions: List[NPiStarInteraction] = []
        import time as _t
        t_pair_start = _t.time()
        try:
            model = structure[0]
        except Exception:
            self.instrumentation = {
                'donors': 0,
                'carbonyl_groups': 0,
                'raw_pairs': 0,
                'candidate_pairs': 0,
                'accepted_pairs': 0,
                'acceptance_ratio': 0.0,
                'candidate_density': 0.0,
                'core_pair_generation': False
            }
            return interactions
        donor_atoms = self._find_lone_pair_donors(model)
        carbonyl_groups = self._find_carbonyl_groups(model)
        if not donor_atoms or not carbonyl_groups:
            self.instrumentation = {
                'donors': len(donor_atoms),
                'carbonyl_groups': len(carbonyl_groups),
                'raw_pairs': 0,
                'candidate_pairs': 0,
                'accepted_pairs': 0,
                'acceptance_ratio': 0.0,
                'candidate_density': 0.0,
                'core_pair_generation': False
            }
            return interactions
        donor_coords = np.vstack([d['atom'].get_coord() for d in donor_atoms]).astype('float32')
        carbon_coords = np.vstack([c['carbon'].get_coord() for c in carbonyl_groups]).astype('float32')
        try:
            di, ci = pairwise_within_cutoff(donor_coords, carbon_coords, self.distance_cutoff, use_kdtree=True)
            core = True
        except Exception:  # pragma: no cover
            diff = donor_coords[:, None, :] - carbon_coords[None, :, :]
            dist2 = np.sum(diff*diff, axis=-1)
            mask = dist2 <= (self.distance_cutoff ** 2)
            di, ci = np.where(mask)
            di = di.astype('int32'); ci = ci.astype('int32')
            core = False
        raw_pairs = int(donor_coords.shape[0] * carbon_coords.shape[0])
        t_pair_end = _t.time()
        t_eval_start = t_pair_end
        candidate_pairs = int(len(di))
        accepted = 0
        for idx in range(candidate_pairs):
            d_idx = int(di[idx]); c_idx = int(ci[idx])
            donor_info = donor_atoms[d_idx]
            carbonyl_info = carbonyl_groups[c_idx]
            if donor_info['residue'] == carbonyl_info['residue']:
                continue
            donor_atom = donor_info['atom']
            carbon_atom = carbonyl_info['carbon']
            oxygen_atom = carbonyl_info['oxygen']
            # Exact distance (we already know within cutoff)
            distance = float(np.linalg.norm(donor_atom.get_coord() - carbon_atom.get_coord()))
            # Bürgi-Dunitz angle
            angle = self._calculate_burgi_dunitz_angle(donor_atom, oxygen_atom, carbon_atom)
            optimal_angle = 107.0
            if abs(angle - optimal_angle) > 25.0:
                continue
            strength = self._classify_interaction_strength(distance, angle)
            interactions.append(NPiStarInteraction(
                donor_atom=donor_atom,
                acceptor_atom=oxygen_atom,
                carbonyl_carbon=carbon_atom,
                distance=distance,
                angle=angle,
                strength=strength,
                donor_residue=f"{donor_info['resname']}{donor_info['res_id'][1]}",
                acceptor_residue=f"{carbonyl_info['resname']}{carbonyl_info['res_id'][1]}",
                donor_chain=donor_info['chain_id'],
                acceptor_chain=carbonyl_info['chain_id']
            ))
            accepted += 1
        t_eval_end = _t.time()
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=accepted,
            core_pair_generation=core,
            extra={
                'donors': len(donor_atoms),
                'carbonyl_groups': len(carbonyl_groups),
                'candidate_density': (candidate_pairs / raw_pairs) if raw_pairs else 0.0
            }
        )
        finalize_funnel(
            self.instrumentation,
            pair_gen_seconds=(t_pair_end - t_pair_start),
            eval_seconds=(t_eval_end - t_eval_start),
            build_seconds=0.0
        )
        logger.info(f"n→π*[vector]{'/core' if core else ''}: {accepted} raw={raw_pairs} pruned={candidate_pairs} acc={self.instrumentation['acceptance_ratio']:.3f}")
        return interactions
    
    def _find_lone_pair_donors(self, model) -> List[Dict[str, Any]]:
        """Find all lone pair donor atoms (oxygen atoms with lone pairs)."""
        donors = []
        
        for chain in model:
            chain_id = chain.get_id()
            
            for residue in chain:
                resname = residue.get_resname()
                res_id = residue.get_id()
                
                # Check backbone oxygen
                if 'O' in residue:
                    donors.append({
                        'atom': residue['O'],
                        'residue': residue,
                        'resname': resname,
                        'chain_id': chain_id,
                        'res_id': res_id
                    })
                
                # Check side chain lone pair donors
                if resname in self.lone_pair_donors:
                    for atom_name in self.lone_pair_donors[resname]:
                        if atom_name != 'backbone' and atom_name in residue:
                            donors.append({
                                'atom': residue[atom_name],
                                'residue': residue,
                                'resname': resname,
                                'chain_id': chain_id,
                                'res_id': res_id
                            })
        
        return donors
    
    def _find_carbonyl_groups(self, model) -> List[Dict[str, Any]]:
        """Find all carbonyl groups (C=O with π* orbitals)."""
        carbonyls = []
        
        for chain in model:
            chain_id = chain.get_id()
            
            for residue in chain:
                resname = residue.get_resname()
                res_id = residue.get_id()
                
                # Check backbone carbonyl
                if 'C' in residue and 'O' in residue:
                    carbonyls.append({
                        'carbon': residue['C'],
                        'oxygen': residue['O'],
                        'residue': residue,
                        'resname': resname,
                        'chain_id': chain_id,
                        'res_id': res_id
                    })
                
                # Check side chain carbonyls
                if resname in self.carbonyl_acceptors:
                    atom_pair = self.carbonyl_acceptors[resname]
                    if atom_pair[0] != 'backbone' and len(atom_pair) == 2:
                        carbon_name, oxygen_name = atom_pair
                        if carbon_name in residue and oxygen_name in residue:
                            carbonyls.append({
                                'carbon': residue[carbon_name],
                                'oxygen': residue[oxygen_name],
                                'residue': residue,
                                'resname': resname,
                                'chain_id': chain_id,
                                'res_id': res_id
                            })
        
        return carbonyls
    
    def _calculate_burgi_dunitz_angle(self, donor_atom: Atom.Atom, carbon_atom: Atom.Atom, oxygen_atom: Atom.Atom) -> float:
        """
        Calculate the Bürgi-Dunitz angle for n→π* interaction.
        This is the angle between the donor...carbon vector and the carbon-oxygen bond.
        Optimal angle is ~102° for antiparallel approach to π* orbital.
        """
        # Vector from donor to carbon
        donor_to_carbon = carbon_atom.get_coord() - donor_atom.get_coord()
        
        # Vector from carbon to oxygen (C=O bond direction)
        carbon_to_oxygen = oxygen_atom.get_coord() - carbon_atom.get_coord()
        
        # Calculate angle between these vectors
        cos_angle = np.dot(donor_to_carbon, carbon_to_oxygen) / (
            np.linalg.norm(donor_to_carbon) * np.linalg.norm(carbon_to_oxygen)
        )
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        
        angle = np.degrees(np.arccos(cos_angle))
        
        # Convert to angle relative to antiparallel (180° - angle gives approach angle)
        return 180.0 - angle
    
    def _calculate_distance(self, atom1: Atom.Atom, atom2: Atom.Atom) -> float:
        """Calculate distance between two atoms."""
        return np.linalg.norm(atom1.get_coord() - atom2.get_coord())
    
    def _classify_interaction_strength(self, distance: float, angle: float) -> str:
        """Classify interaction strength based on distance and angle."""
        # Strong: optimal distance and angle
        if distance < 2.9 and abs(angle - self.optimal_angle) < 10:
            return 'strong'
        
        # Weak: far from optimal
        elif distance > 3.1 or abs(angle - self.optimal_angle) > 20:
            return 'weak'
        
        # Moderate: in between
        else:
            return 'moderate'
    
    def _calculate_burgi_dunitz_angle(self, donor: Atom.Atom, oxygen: Atom.Atom, carbon: Atom.Atom) -> float:
        """Calculate the Bürgi-Dunitz angle (donor-oxygen-carbon)."""
        # Vector from oxygen to donor
        vec_od = donor.coord - oxygen.coord
        # Vector from oxygen to carbon
        vec_oc = carbon.coord - oxygen.coord
        
        # Calculate angle
        cos_angle = np.dot(vec_od, vec_oc) / (np.linalg.norm(vec_od) * np.linalg.norm(vec_oc))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle = np.degrees(np.arccos(cos_angle))
        
        return angle
    
    def _classify_interaction_strength(self, distance: float, angle: float) -> str:
        """Classify interaction strength based on distance and angle."""
        # Optimal Bürgi-Dunitz angle is around 107°
        optimal_angle = 107.0
        angle_deviation = abs(angle - optimal_angle)
        
        if distance <= 3.0 and angle_deviation <= 15.0:
            return 'strong'
        elif distance <= 3.3 and angle_deviation <= 25.0:
            return 'moderate'
        else:
            return 'weak'
    
    def to_dict_list(self, interactions: List[NPiStarInteraction]) -> List[Dict[str, Any]]:
        """Convert interactions to list of dictionaries."""
        return [
            {
                'type': 'n_pi_star',
                'donor_atom': f"{interaction.donor_atom.get_parent().get_resname()}{interaction.donor_atom.get_parent().id[1]}:{interaction.donor_atom.name}",
                'acceptor_atom': f"{interaction.acceptor_atom.get_parent().get_resname()}{interaction.acceptor_atom.get_parent().id[1]}:{interaction.acceptor_atom.name}",
                'carbonyl_carbon': f"{interaction.carbonyl_carbon.get_parent().get_resname()}{interaction.carbonyl_carbon.get_parent().id[1]}:{interaction.carbonyl_carbon.name}",
                'distance': round(interaction.distance, 2),
                'angle': round(interaction.angle, 1),
                'strength': interaction.strength,
                'donor_residue': interaction.donor_residue,
                'acceptor_residue': interaction.acceptor_residue,
                'donor_chain': interaction.donor_chain,
                'acceptor_chain': interaction.acceptor_chain,
                # Standard keys for UI compatibility
                'residue1': interaction.donor_residue,
                'residue2': interaction.acceptor_residue,
                'chain1': interaction.donor_chain,
                'chain2': interaction.acceptor_chain
            }
            for interaction in interactions
        ]
