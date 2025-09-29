"""
Anion-π interaction detection for protein structures.
Detects electron-rich ions close to π-acidic rings (positive quadrupole); 3.5–4.5 Å; typically Asp/Glu-aromatic interactions.
"""

import numpy as np
import time
from typing import List, Dict, Any, Optional, Tuple
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from loguru import logger
from utils.settings import get_settings
from analysis.feature_store import get_feature_store
from geometry.core import pairwise_within_cutoff, norms, vector_fields
from utils.instrumentation import init_funnel, update_counts, finalize_funnel

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
        # Backward compatibility: some earlier code referenced self.aromatic_residues
        self.aromatic_residues = self.pi_acidic_residues
    
    def detect_anion_pi_interactions(self, structure: Structure.Structure) -> List[AnionPiInteraction]:
        settings = get_settings()
        if getattr(settings, 'enable_vector_geom', False):
            try:
                return self._vector_detect(structure)
            except Exception:  # pragma: no cover
                pass
        return self._legacy_detect(structure)

    def _legacy_detect(self, structure: Structure.Structure) -> List[AnionPiInteraction]:
        interactions: List[AnionPiInteraction] = []
        t_pair_start = time.time()
        anion_atoms = self._find_anion_atoms(structure)
        aromatic_rings = self._find_aromatic_rings(structure)
        raw_pairs = len(anion_atoms) * len(aromatic_rings)
        candidate_pairs = 0
        for anion in anion_atoms:
            for ring_center, ring_residue in aromatic_rings:
                if anion.get_parent() == ring_residue:
                    continue
                distance = self._calculate_distance(anion.coord, ring_center)
                if distance <= self.distance_max:
                    candidate_pairs += 1
                    interactions.append(AnionPiInteraction(
                        anion_atom=anion,
                        pi_center=ring_center,
                        pi_residue_obj=ring_residue,
                        distance=distance,
                        strength=self._classify_interaction_strength(distance),
                        anion_residue=f"{anion.get_parent().get_resname()}{anion.get_parent().id[1]}",
                        pi_residue=f"{ring_residue.get_resname()}{ring_residue.id[1]}",
                        anion_chain=anion.get_parent().get_parent().id,
                        pi_chain=ring_residue.get_parent().id
                    ))
        t_eval_end = time.time()
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=len(interactions),
            core_pair_generation=False,
            extra={
                'anion_atoms': len(anion_atoms),
                'rings': len(aromatic_rings),
                'candidate_density': (candidate_pairs/raw_pairs) if raw_pairs else 0.0
            }
        )
        finalize_funnel(
            self.instrumentation,
            pair_gen_seconds=(t_eval_end - t_pair_start),  # combined phase
            eval_seconds=0.0,
            build_seconds=0.0
        )
        logger.info(f"[legacy] Anion-π: {len(interactions)} (raw={raw_pairs} cand={candidate_pairs} acc_ratio={self.instrumentation['acceptance_ratio']:.3f})")
        return interactions

    def _vector_detect(self, structure: Structure.Structure) -> List[AnionPiInteraction]:
        interactions: List[AnionPiInteraction] = []
        fs = get_feature_store()
        t_pair_start = time.time()
        rings_cache = fs.ensure_rings(structure)
        anion_atoms = self._find_anion_atoms(structure)
        if not rings_cache or not anion_atoms:
            self.instrumentation = init_funnel(
                raw_pairs=0,
                candidate_pairs=0,
                accepted_pairs=0,
                core_pair_generation=False,
                extra={
                    'anion_atoms': len(anion_atoms),
                    'rings': len(rings_cache) if rings_cache else 0,
                    'candidate_density': 0.0
                }
            )
            finalize_funnel(self.instrumentation, pair_gen_seconds=0.0, eval_seconds=0.0, build_seconds=0.0)
            return interactions
        ring_centers = np.vstack([r['centroid'] for r in rings_cache]).astype('float32')
        ring_residues = [r['residue'] for r in rings_cache]
        anion_coords = np.vstack([a.get_coord() for a in anion_atoms]).astype('float32')
        try:
            ai, ri, dists = pairwise_within_cutoff(anion_coords, ring_centers, self.distance_max, use_kdtree=True, return_distances=True)
            core = True
        except Exception:  # pragma: no cover
            diff = anion_coords[:, None, :] - ring_centers[None, :, :]
            dist2 = np.sum(diff*diff, axis=-1)
            mask = dist2 <= (self.distance_max ** 2)
            ai, ri = np.where(mask)
            ai = ai.astype('int32'); ri = ri.astype('int32')
            dists = np.sqrt(dist2[ai, ri]).astype('float32')
            core = False
        raw_pairs = int(anion_coords.shape[0] * ring_centers.shape[0])
        keep = (dists >= self.distance_min) & (dists <= self.distance_max)
        ai = ai[keep]; ri = ri[keep]; dists = dists[keep]
        candidate_pairs = int(len(ai))
        t_eval_start = time.time()
        for idx in range(candidate_pairs):
            anion = anion_atoms[int(ai[idx])]
            ring_residue = ring_residues[int(ri[idx])]
            if anion.get_parent() == ring_residue:
                continue
            distance = float(dists[idx])
            interactions.append(AnionPiInteraction(
                anion_atom=anion,
                pi_center=ring_centers[int(ri[idx])],
                pi_residue_obj=ring_residue,
                distance=distance,
                strength=self._classify_interaction_strength(distance),
                anion_residue=f"{anion.get_parent().get_resname()}{anion.get_parent().id[1]}",
                pi_residue=f"{ring_residue.get_resname()}{ring_residue.id[1]}",
                anion_chain=anion.get_parent().get_parent().id,
                pi_chain=ring_residue.get_parent().id
            ))
        t_eval_end = time.time()
        self.instrumentation = init_funnel(
            raw_pairs=raw_pairs,
            candidate_pairs=candidate_pairs,
            accepted_pairs=len(interactions),
            core_pair_generation=core,
            extra={
                'anion_atoms': len(anion_atoms),
                'rings': ring_centers.shape[0],
                'candidate_density': (candidate_pairs/raw_pairs) if raw_pairs else 0.0
            }
        )
        finalize_funnel(
            self.instrumentation,
            pair_gen_seconds=(t_eval_start - t_pair_start),
            eval_seconds=(t_eval_end - t_eval_start),
            build_seconds=0.0
        )
        logger.info(f"Anion-π[vector]{'/core' if core else ''}: {len(interactions)} raw={raw_pairs} pruned={candidate_pairs} acc_ratio={self.instrumentation['acceptance_ratio']:.3f}")
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
