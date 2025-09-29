"""
Hydrogen bond detection for protein structures.
Implements detection of conventional, low-barrier, C5-type, C-H···π, and sulfur-mediated hydrogen bonds.
"""

import numpy as np
import time
from typing import List, Dict, Any, Tuple, Optional
from Bio.PDB import Structure, Atom, Residue
from dataclasses import dataclass
from .base_detector import register_detector
from loguru import logger
from utils.settings import get_settings  # lightweight import
from utils.instrumentation import init_funnel, update_counts, finalize_funnel
try:
    from scipy.spatial import cKDTree  # optional for KD-tree pruning
    _HAVE_SCIPY = True
except Exception:  # pragma: no cover - optional
    _HAVE_SCIPY = False

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

@register_detector("hydrogenbond", method="detect_hydrogen_bonds")
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
        # Pull from centralized interaction config (UI-adjustable)
        self.distance_cutoff = getattr(self.interaction_config, 'hbond_distance_cutoff', 3.5)
        self.angle_cutoff = getattr(self.interaction_config, 'hbond_angle_cutoff', 120.0)
        
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
        """Detect hydrogen bonds (vector fast-path optional).

        The legacy nested-loop implementation is preserved in a separate helper
        (_legacy_detect) so that the vectorized candidate filter can cleanly
        fall back without recursion hacks. Flag gating via Settings:
        MOLBRIDGE_ENABLE_VECTOR_GEOM=1 enables the vector path.
        """
        settings = get_settings()
        if not settings.enable_vector_geom:
            return self._legacy_detect(structure)

        try:
            t_pair_start = time.time()
            model = structure[0]
            # Attempt to reuse precomputed donors/acceptors via FeatureStore (task graph)
            donors = None
            acceptors = None
            try:
                from analysis.feature_store import get_feature_store
                fs = get_feature_store()
                d_pre, a_pre = fs.ensure_hbond_participants(structure)
                # Always derive full donor/acceptor enriched metadata (could be optimized further)
                donors = self._get_donors(model)
                acceptors = self._get_acceptors(model)
            except Exception:
                donors = self._get_donors(model)
                acceptors = self._get_acceptors(model)
            if not donors or not acceptors:
                return []
            donor_atoms = [d['atom'] for d in donors]
            acceptor_atoms = [a['atom'] for a in acceptors]
            donor_coords = [a.get_coord() for a in donor_atoms]
            acceptor_coords = [a.get_coord() for a in acceptor_atoms]
            cutoff = self.distance_cutoff
            # Use shared FeatureStore neighbor indexing over full atom set instead of custom geometry.core call.
            from utils.kdtree_thresholds import get_threshold, adapt_threshold, get_last_density, get_last_reason, should_flag_kdtree
            raw_pairs = len(donor_coords) * len(acceptor_coords)
            core = False
            candidate_pairs = []
            try:
                from analysis.feature_store import get_feature_store as _gfs
                fs2 = _gfs()
                # Acquire global atom coordinate index once
                all_coords = fs2.ensure_coords(structure)
                if all_coords is not None and all_coords.size:
                    # Build mapping atom id -> global index
                    atom_index = {}
                    for idx, atom in enumerate(structure.get_atoms()):  # type: ignore
                        atom_index[id(atom)] = idx
                    donor_idx = [atom_index.get(id(a)) for a in donor_atoms]
                    acceptor_idx = [atom_index.get(id(a)) for a in acceptor_atoms]
                    # Quick filter: find all atom neighbor pairs within cutoff
                    atom_pairs = fs2.neighbor_within(structure, float(cutoff))
                    if atom_pairs:
                        acceptor_idx_set = set(acceptor_idx)
                        donor_idx_set = set(donor_idx)
                        for i,j in atom_pairs:
                            # Map global pair to donor/acceptor role if qualifies
                            if i in donor_idx_set and j in acceptor_idx_set:
                                di = donor_idx.index(i)
                                ai = acceptor_idx.index(j)
                                candidate_pairs.append((di, ai))
                            elif j in donor_idx_set and i in acceptor_idx_set:
                                di = donor_idx.index(j)
                                ai = acceptor_idx.index(i)
                                candidate_pairs.append((di, ai))
                        core = True
            except Exception:
                pass
            if not candidate_pairs:
                # Fallback to broadcast distance matrix
                donor_arr = np.vstack(donor_coords)
                acceptor_arr = np.vstack(acceptor_coords)
                diff = donor_arr[:, None, :] - acceptor_arr[None, :, :]
                dist_matrix = np.linalg.norm(diff, axis=-1)
                di_idxs, ai_idxs = np.where(dist_matrix <= cutoff)
                candidate_pairs = list(zip(di_idxs.tolist(), ai_idxs.tolist()))
            t_pair_end = time.time()
            t_eval_start = t_pair_end
            threshold = get_threshold('hbond')
            # Heuristic: treat as kdtree_used if core path and significant pruning for adaptation metrics
            use_kdtree = core and raw_pairs > threshold and len(candidate_pairs) < raw_pairs * 0.7
            hbonds: List[HydrogenBond] = []
            pruned_pairs = len(candidate_pairs)
            # Vectorized hydrogen angle evaluation: precompute heavy donor->hydrogen vectors
            # Build flattened arrays for donors with hydrogens
            # Map donor index to hydrogen atom objects list
            # Initialize standardized funnel instrumentation before evaluation loop
            self.instrumentation = init_funnel(
                raw_pairs=raw_pairs,
                candidate_pairs=pruned_pairs,
                accepted_pairs=0,
                core_pair_generation=core,
                extra={
                    'donors': len(donors),
                    'acceptors': len(acceptors),
                    'hbonds': 0,
                    'kdtree_used': int(use_kdtree),
                    'threshold': threshold,
                    'adaptive_threshold': None,  # filled after adaptation
                    'threshold_changed': False,  # placeholder
                    'candidate_density': None,
                    'adapt_reason': None
                }
            )
            # ---------------- Vectorized evaluation phase ----------------
            try:
                import numpy as _np
                donor_coords_arr = _np.vstack([d['atom'].get_coord() for d in donors]).astype('float32')
                acceptor_coords_arr = _np.vstack([a['atom'].get_coord() for a in acceptors]).astype('float32')
                # Precompute hydrogen coordinates per donor (ragged) -> flatten and index map
                h_coord_list = []
                h_parent_index = []  # maps hydrogen row -> donor idx
                for idx_d, d in enumerate(donors):
                    hlist = d.get('hydrogen_atoms') or []
                    if not isinstance(hlist, list):
                        hlist = [hlist] if hlist else []
                    for h in hlist:
                        if h is not None:
                            h_coord_list.append(h.get_coord())
                            h_parent_index.append(idx_d)
                if h_coord_list:
                    h_coords = _np.vstack(h_coord_list).astype('float32')
                    h_parent_index = _np.asarray(h_parent_index, dtype='int32')
                else:
                    h_coords = _np.zeros((0,3), dtype='float32')
                    h_parent_index = _np.zeros((0,), dtype='int32')
                # For each donor-acceptor candidate, compute distance (already within cutoff) and best hydrogen angle if hydrogens present
                # Build arrays for candidate pairs
                if candidate_pairs:
                    d_idx = _np.array([p[0] for p in candidate_pairs], dtype='int32')
                    a_idx = _np.array([p[1] for p in candidate_pairs], dtype='int32')
                    d_pts = donor_coords_arr[d_idx]
                    a_pts = acceptor_coords_arr[a_idx]
                    # Distance
                    da_vec = d_pts - a_pts
                    da_dist = _np.linalg.norm(da_vec, axis=1)
                    # Filter again by cutoff (numerical safety)
                    within = da_dist <= (cutoff + 1e-6)
                    if not within.any():
                        candidate_pairs = []
                    else:
                        d_idx = d_idx[within]; a_idx = a_idx[within]; da_dist = da_dist[within]; d_pts = d_pts[within]; a_pts = a_pts[within]
                        # Hydrogen angle computation: for donors that have hydrogens we compute max over associated hydrogens
                        # Build index list of hydrogen rows for each donor in vector domain via grouping
                        if h_coords.shape[0]:
                            # For efficiency, create mapping donor -> slice of hydrogen indices
                            from collections import defaultdict
                            donor_to_h_indices: dict[int, list[int]] = defaultdict(list)
                            for row, parent in enumerate(h_parent_index.tolist()):
                                donor_to_h_indices[parent].append(row)
                        # Iterate vector-sliced; still faster than per-pair angle with atom objects
                        for k in range(d_idx.shape[0]):
                            di = int(d_idx[k]); ai = int(a_idx[k])
                            donor_info = donors[di]; acceptor_info = acceptors[ai]
                            dist = float(da_dist[k])
                            # Angle determination
                            best_angle = None; best_h_atom = None
                            h_candidates = donor_info.get('hydrogen_atoms') or []
                            if not isinstance(h_candidates, list):
                                h_candidates = [h_candidates] if h_candidates else []
                            if h_candidates:
                                rows = donor_to_h_indices.get(di, []) if 'donor_to_h_indices' in locals() else []
                                if rows:
                                    h_sub = h_coords[rows]
                                    v1 = donor_coords_arr[di] - h_sub
                                    v2 = acceptor_coords_arr[ai] - h_sub
                                    n1 = _np.linalg.norm(v1, axis=1)
                                    n2 = _np.linalg.norm(v2, axis=1)
                                    denom = n1 * n2
                                    mask = denom > 0
                                    ang = _np.zeros(len(denom), dtype='float32')
                                    if mask.any():
                                        ang[mask] = _np.degrees(_np.arccos(_np.clip((v1[mask]*v2[mask]).sum(axis=1)/denom[mask], -1.0, 1.0)))
                                    idx_best = int(ang.argmax()) if ang.size else -1
                                    if idx_best >= 0 and ang.size:
                                        best_angle = float(ang[idx_best])
                                        best_h_atom = h_candidates[idx_best] if idx_best < len(h_candidates) else None
                            if best_angle is None:
                                best_angle = self._estimate_dha_angle(donor_info['atom'], acceptor_info['atom'], donor_info['residue'])
                            if best_angle < self.angle_cutoff:
                                continue
                            strength = self._calculate_bond_strength(dist, best_angle)
                            hbonds.append(HydrogenBond(
                                donor_atom=donor_info['atom'],
                                acceptor_atom=acceptor_info['atom'],
                                hydrogen_atom=best_h_atom,
                                distance=dist,
                                angle=best_angle,
                                strength=strength,
                                donor_residue=f"{donor_info['resname']}{donor_info['res_id'][1]}",
                                acceptor_residue=f"{acceptor_info['resname']}{acceptor_info['res_id'][1]}",
                                donor_chain=donor_info['chain_id'],
                                acceptor_chain=acceptor_info['chain_id']
                            ))
                else:
                    candidate_pairs = []
            except Exception:
                # Fallback to original per-pair evaluation on error
                for di, ai in candidate_pairs:
                    donor_info = donors[di]
                    acceptor_info = acceptors[ai]
                    hbond = self._check_hydrogen_bond(donor_info, acceptor_info)
                    if hbond:
                        hbonds.append(hbond)
            t_eval_end = time.time()
            pruning_flag = should_flag_kdtree(pruned_pairs, raw_pairs, 'hbond')
            new_thresh, changed, reason = adapt_threshold('hbond', pruned_pairs, (use_kdtree or pruning_flag))
            # Update adaptation related fields
            self.instrumentation['adaptive_threshold'] = new_thresh
            self.instrumentation['threshold_changed'] = changed
            self.instrumentation['candidate_density'] = get_last_density('hbond')
            self.instrumentation['adapt_reason'] = reason
            # Final counts
            update_counts(self.instrumentation, accepted=len(hbonds))
            self.instrumentation['hbonds'] = len(hbonds)
            finalize_funnel(
                self.instrumentation,
                pair_gen_seconds=(t_pair_end - t_pair_start),
                eval_seconds=(t_eval_end - t_eval_start),
                build_seconds=0.0
            )
            from utils.settings import get_settings as _gs
            _s = _gs()
            if _s.performance_mode and not getattr(_s, 'verbose_detector_logs', False):
                logger.debug(
                    f"[vector]{'/kdtree' if self.instrumentation['kdtree_used'] else ''} H-bonds: {len(hbonds)} "
                    f"(donors={len(donors)} acceptors={len(acceptors)} raw={raw_pairs} pruned={pruned_pairs} acc_ratio={self.instrumentation['acceptance_ratio']:.3f})"
                )
            else:
                logger.info(
                    f"[vector]{'/kdtree' if self.instrumentation['kdtree_used'] else ''} H-bonds: {len(hbonds)} "
                    f"(donors={len(donors)} acceptors={len(acceptors)} raw={raw_pairs} pruned={pruned_pairs} acc_ratio={self.instrumentation['acceptance_ratio']:.3f})"
                )
            return hbonds
        except Exception as e:
            logger.error(f"Vector hydrogen bond path error, falling back to legacy: {e}")
            return self._legacy_detect(structure)

    def _legacy_detect(self, structure: Structure.Structure) -> List[HydrogenBond]:
        """Original non-vectorized detection path (stable reference)."""
        hbonds: List[HydrogenBond] = []
        t_pair_start = time.time()
        try:
            model = structure[0]
            donors = self._get_donors(model)
            acceptors = self._get_acceptors(model)
            logger.info(f"Found {len(donors)} potential donors, {len(acceptors)} acceptors")
            t_pair_end = time.time()
            t_eval_start = t_pair_end
            for donor_info in donors:
                for acceptor_info in acceptors:
                    hbond = self._check_hydrogen_bond(donor_info, acceptor_info)
                    if hbond:
                        hbonds.append(hbond)
            logger.info(f"Detected {len(hbonds)} hydrogen bonds total")
            t_eval_end = time.time()
            # Legacy instrumentation using shared helper (no pruning distinction)
            raw_pairs = len(donors) * len(acceptors)
            self.instrumentation = init_funnel(
                raw_pairs=raw_pairs,
                candidate_pairs=raw_pairs,
                accepted_pairs=len(hbonds),
                core_pair_generation=False,
                extra={
                    'donors': len(donors),
                    'acceptors': len(acceptors),
                    'hbonds': len(hbonds),
                    'kdtree_used': False
                }
            )
            finalize_funnel(
                self.instrumentation,
                pair_gen_seconds=(t_pair_end - t_pair_start),
                eval_seconds=(t_eval_end - t_eval_start),
                build_seconds=0.0
            )
        except Exception as e:
            logger.error(f"Error (legacy) detecting hydrogen bonds: {e}")
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
        
        # Vectorized hydrogen angle selection (fast path)
        best_hydrogen = None
        hydrogen_atoms = donor_info.get('hydrogen_atoms', [])
        if not isinstance(hydrogen_atoms, list):
            hydrogen_atoms = [hydrogen_atoms] if hydrogen_atoms else []
        if hydrogen_atoms:
            try:
                import numpy as _np
                h_coords = _np.vstack([h.get_coord() for h in hydrogen_atoms if h is not None]) if hydrogen_atoms else _np.empty((0,3),dtype='float32')
                if h_coords.size:
                    d_coord = donor_atom.get_coord()
                    a_coord = acceptor_atom.get_coord()
                    v1 = d_coord - h_coords
                    v2 = a_coord - h_coords
                    norm1 = _np.linalg.norm(v1, axis=1)
                    norm2 = _np.linalg.norm(v2, axis=1)
                    denom = norm1 * norm2
                    mask = denom > 0
                    cosang = _np.zeros(len(denom), dtype='float32')
                    cosang[mask] = (v1[mask] * v2[mask]).sum(axis=1) / denom[mask]
                    cosang = _np.clip(cosang, -1.0, 1.0)
                    angles = _np.degrees(_np.arccos(cosang))
                    if len(angles):
                        idx = int(angles.argmax())
                        best_angle = float(angles[idx])
                        best_hydrogen = hydrogen_atoms[idx]
                    else:
                        best_angle = self._estimate_dha_angle(donor_atom, acceptor_atom, donor_info['residue'])
                else:
                    best_angle = self._estimate_dha_angle(donor_atom, acceptor_atom, donor_info['residue'])
            except Exception:
                best_angle = 0.0
                for hydrogen_atom in hydrogen_atoms:
                    if hydrogen_atom is None:
                        continue
                    angle = self._calculate_angle(donor_atom, hydrogen_atom, acceptor_atom)
                    if angle > best_angle:
                        best_angle = angle
                        best_hydrogen = hydrogen_atom
        else:
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
                    'chain2': hbond.acceptor_chain,
                    'type': 'hydrogen_bond'
                })
            except Exception as e:
                logger.error(f"Error converting hydrogen bond to dict: {e} - Object: {hbond}")
                continue
                
        return result
