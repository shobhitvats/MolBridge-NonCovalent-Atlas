"""
PDB file handling and structure processing for Protein Interaction Explorer.
"""

import io
import warnings
import requests
from typing import Optional, Dict, Any, List, Tuple
from pathlib import Path
import numpy as np
import asyncio
try:
    import httpx  # for async multi-fetch
    _HTTPX_AVAILABLE = True
except Exception:  # pragma: no cover
    _HTTPX_AVAILABLE = False

# Suppress Bio warnings
warnings.filterwarnings("ignore")

from Bio import PDB
from Bio.PDB import Structure, Model, Chain, Residue, Atom
from Bio.PDB.PDBParser import PDBParser
try:
    try:  # some BioPython builds may not expose MMCIFIO submodule path identically
        from Bio.PDB.MMCIFIO import MMCIFIO  # type: ignore
    except Exception:  # pragma: no cover
        try:
            from Bio.PDB import MMCIFIO  # type: ignore
        except Exception:
            MMCIFIO = None  # type: ignore
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
from utils.structure_hash import compute_structure_coord_hash
import numpy as np

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
            
            # Lightweight surrogate cache path (avoids reparsing heavy hierarchy)
            surrogate_key = f"surrogate_{pdb_id}_{assembly}"
            structure = None
            if hasattr(self, 'cache_manager') and self.cache_manager:
                surrogate = self.cache_manager.cache.get(surrogate_key)
                if surrogate is not None:
                    try:
                        # Rehydrate minimal structure surrogate into Bio.PDB Structure
                        structure = self._rehydrate_surrogate(pdb_id, surrogate)
                    except Exception:
                        structure = None
            if structure is None:
                # Parse the structure (original path)
                # Attempt fast-path reconstruction from cached binary coordinate snapshot (.npz)
                fast_rehydrated = False
                if hasattr(self, 'cache_manager') and self.cache_manager:
                    try:
                        # Direct index key (new) to avoid scanning full key set
                        idx_key = f"npz_index_{pdb_id}_{assembly}"
                        npz_key = self.cache_manager.cache.get(idx_key)
                        if npz_key is None:
                            # Backwards compatibility: fall back to heuristic scan once, then store index key
                            hits = []
                            for key in list(self.cache_manager.cache.iterkeys())[:500]:  # heuristic cap
                                if isinstance(key, str) and key.startswith(f"npz_coords_{pdb_id}_"):
                                    hits.append(key)
                            if hits:
                                npz_key = hits[-1]
                                # Cache the mapping for next time
                                try:
                                    self.cache_manager.cache.set(idx_key, npz_key, expire=7*24*3600)
                                except Exception:
                                    pass
                        if npz_key:
                            blob = self.cache_manager.cache.get(npz_key)
                            if blob:
                                import io as _io, numpy as _np
                                buf = _io.BytesIO(blob)
                                npz = _np.load(buf)
                                coords = npz['coords'] if 'coords' in npz.files else None
                                surrogate = None
                                try:
                                    surrogate = self.cache_manager.cache.get(f"surrogate_{pdb_id}_{assembly}")
                                except Exception:
                                    surrogate = None
                                if surrogate and isinstance(coords, _np.ndarray) and coords.shape[0] == surrogate.get('atom_count'):
                                    surrogate = dict(surrogate)
                                    surrogate['coords'] = coords.astype('float32')
                                    structure = self._rehydrate_surrogate(pdb_id, surrogate)
                                    fast_rehydrated = structure is not None
                    except Exception:
                        fast_rehydrated = False
                if not fast_rehydrated:
                    pdb_io = io.StringIO(pdb_content)
                    structure = self.parser.get_structure(pdb_id, pdb_io)
                # Build & cache surrogate for future fast loads
                if hasattr(self, 'cache_manager') and self.cache_manager:
                    try:
                        surrogate = self._build_surrogate(structure)
                        self.cache_manager.cache.set(surrogate_key, surrogate, expire=7*24*3600)
                    except Exception:  # pragma: no cover - non-critical
                        pass
                # Attempt binary coordinate cache (.npz) persistence keyed by hash
                try:
                    struct_hash = compute_structure_coord_hash(structure)
                    if hasattr(self, 'cache_manager') and self.cache_manager:
                        npz_key = f"npz_coords_{pdb_id}_{struct_hash}"
                        existing = self.cache_manager.cache.get(npz_key)
                        if not existing:
                            arr = np.asarray([a.get_coord() for a in structure.get_atoms()], dtype='float32')  # type: ignore
                            meta = np.asarray([0, arr.shape[0]], dtype='int32')  # simple header placeholder
                            import io as _io
                            buf = _io.BytesIO()
                            # Store coordinates in compressed npz
                            np.savez_compressed(buf, coords=arr, meta=meta)
                            self.cache_manager.cache.set(npz_key, buf.getvalue(), expire=7*24*3600)
                            # Store index key for direct lookup next time
                            try:
                                idx_key = f"npz_index_{pdb_id}_{assembly}"
                                self.cache_manager.cache.set(idx_key, npz_key, expire=7*24*3600)
                            except Exception:
                                pass
                except Exception:  # pragma: no cover
                    pass
            
            # Process the structure based on configuration
            processed_structure = self._process_structure(structure)
            
            # Cache structure info
            if hasattr(self, 'cache_manager') and self.cache_manager:
                structure_info = self._extract_structure_info(processed_structure)
                self.cache_manager.cache_structure_info(pdb_id, structure_info)
            
            logger.info(f"Successfully loaded structure {pdb_id}")
            # Fire-and-forget background prewarm (KD-tree, rings) to hide first-detector cost
            try:
                import threading
                def _prewarm(struct_ref):  # pragma: no cover - best-effort thread
                    try:
                        from analysis.feature_store import get_feature_store
                        fs = get_feature_store()
                        fs.ensure_coords(struct_ref)
                        fs.ensure_kdtree(struct_ref)
                        fs.ensure_rings(struct_ref)
                        fs.ensure_classification(struct_ref)
                    except Exception:
                        pass
                threading.Thread(target=_prewarm, args=(processed_structure,), daemon=True).start()
            except Exception:
                pass
            return processed_structure
            
        except Exception as e:
            logger.error(f"Failed to load structure {pdb_id}: {e}")
            return None

    # ---------------- Asynchronous multi-PDB fetch with ETag -----------------
    async def fetch_multiple_pdbs(self, pdb_ids: List[str], assembly: str = "biological", concurrency: int = 5, use_etag: bool = True) -> Dict[str, Optional[str]]:
        """Fetch multiple PDB files concurrently using httpx with optional ETag conditional requests.

        Returns mapping pdb_id -> content (None if failed or 304 not modified and cache manager provides stored version).
        """
        results: Dict[str, Optional[str]] = {}
        if not _HTTPX_AVAILABLE:
            # Fallback sequential using existing synchronous path
            for pid in pdb_ids:
                try:
                    results[pid] = self.cache_manager.get_pdb_file(pid, assembly) if self.cache_manager else self._download_pdb_file(pid, assembly)
                except Exception:
                    results[pid] = None
            return results
        sem = asyncio.Semaphore(concurrency)
        base_url_primary = "https://files.rcsb.org/download/{pid}.pdb1" if assembly == "biological" else "https://files.rcsb.org/download/{pid}.pdb"
        base_url_fallback = "https://files.rcsb.org/download/{pid}.pdb"

        async def _fetch(pid: str):
            async with sem:
                url = base_url_primary.format(pid=pid)
                headers = {}
                if use_etag and self.cache_manager:
                    et = self.cache_manager.get_pdb_etag(pid, assembly)
                    if et:
                        headers['If-None-Match'] = et
                try:
                    async with httpx.AsyncClient(timeout=30) as client:
                        r = await client.get(url, headers=headers)
                        if r.status_code == 404 and assembly == 'biological':
                            # retry fallback
                            url2 = base_url_fallback.format(pid=pid)
                            r = await client.get(url2, headers=headers)
                        if r.status_code == 304 and self.cache_manager:
                            # Not modified; retrieve cached
                            cached = self.cache_manager.get_pdb_file(pid, assembly)
                            results[pid] = cached
                            return
                        if r.status_code != 200:
                            results[pid] = None
                            return
                        text = r.text
                        results[pid] = text
                        # Cache content & update ETag
                        if self.cache_manager:
                            try:
                                self.cache_manager.cache.set(f"pdb_{pid}_{assembly}", text, expire=7*24*3600)
                                etag = r.headers.get('ETag')
                                if etag:
                                    self.cache_manager.set_pdb_etag(pid, assembly, etag)
                            except Exception:
                                pass
                except Exception:
                    results[pid] = None
        await asyncio.gather(*[_fetch(pid) for pid in pdb_ids])
        return results

    # ---------------- Surrogate Build / Rehydrate -----------------
    def _build_surrogate(self, structure: Structure.Structure) -> Dict[str, Any]:
        """Extract minimal arrays sufficient for downstream coordinate analyses.

        Returns dict with keys: chains -> list; residues -> list of tuples;
        atoms -> ndarray float32; atom_meta -> list of (resname, resid, chain, atom_name, element)
        """
        coords = []
        meta = []
        try:
            for chain in structure.get_chains():
                cid = chain.get_id()
                for residue in chain:
                    resname = residue.get_resname()
                    resid = residue.get_id()[1]
                    for atom in residue:
                        coords.append(atom.get_coord())
                        meta.append((resname, resid, cid, atom.get_name(), getattr(atom, 'element', '')))
            arr = np.asarray(coords, dtype=np.float32) if coords else np.zeros((0,3), dtype=np.float32)
            return {
                'coords': arr,
                'meta': meta,
                'atom_count': arr.shape[0]
            }
        except Exception as e:  # pragma: no cover
            logger.debug(f"Surrogate build failed: {e}")
            return {'coords': np.zeros((0,3), dtype=np.float32), 'meta': [], 'atom_count': 0}

    def _rehydrate_surrogate(self, pdb_id: str, surrogate: Dict[str, Any]) -> Optional[Structure.Structure]:
        """Reconstruct a minimal Bio.PDB Structure with atom coordinates and metadata.

        Hierarchy: Structure -> Model(0) -> Chain -> Residue -> Atom.
        Only fields required by detectors (resname, id, chain id, element, coord) are populated.
        """
        try:
            builder = PDB.StructureBuilder.StructureBuilder()  # type: ignore
            builder.init_structure(pdb_id)
            builder.init_model(0)
            current_chain = None
            current_res_key = None
            for (resname, resid, cid, atom_name, element), coord in zip(surrogate.get('meta', []), surrogate.get('coords', [])):
                if current_chain != cid:
                    builder.init_chain(cid)
                    current_chain = cid
                res_key = (cid, resid, resname)
                if current_res_key != res_key:
                    builder.init_seg(' ')
                    builder.init_residue(resname, ' ', resid, ' ')
                    current_res_key = res_key
                builder.init_atom(atom_name, coord, 1.0, 1.0, ' ', atom_name, 0, element or atom_name[0])
            # Attach cached coordinate hash if available to avoid recomputation
            try:
                setattr(builder.get_structure(), '_molbridge_coord_hash', (16, compute_structure_coord_hash(builder.get_structure()), surrogate.get('atom_count', 0)))
            except Exception:
                pass
            return builder.get_structure()
        except Exception as e:  # pragma: no cover
            logger.debug(f"Surrogate rehydrate failed: {e}")
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

    # ------------------------------------------------------------------
    # New Variant Fetching & Filtering API
    # ------------------------------------------------------------------
    def fetch_structure_variant(self,
                                pdb_id: str,
                                assembly: str = "biological",
                                include_ligands: bool = True,
                                exclude_waters: bool = True,
                                ligand_whitelist: Optional[List[str]] = None) -> Optional[Structure.Structure]:
        """Fetch a structure variant with assembly selection and optional filtering.

        This method wraps load_structure but applies ligand/water filtering AFTER
        initial processing while caching variants to avoid repeated parsing.

        Cache key format: f"{pdb_id}|{assembly}|L{int(include_ligands)}|W{int(exclude_waters)}"
        """
        try:
            if ligand_whitelist is None:
                ligand_whitelist = []
            if not hasattr(self, '_variant_cache'):
                self._variant_cache = {}
            cache_key = f"{pdb_id}|{assembly}|L{int(include_ligands)}|W{int(exclude_waters)}"
            if cache_key in self._variant_cache:
                return self._variant_cache[cache_key]

            base_structure = self.load_structure(pdb_id, assembly=assembly)
            if base_structure is None:
                return None

            # Deep copy the structure (BioPython structures are mutable); a light clone
            import copy
            structure = copy.deepcopy(base_structure)

            # Apply dynamic water removal override (overriding global config if necessary)
            if exclude_waters:
                self._remove_waters(structure)

            if not include_ligands:
                self._remove_ligands(structure, whitelist=ligand_whitelist)

            self._variant_cache[cache_key] = structure
            return structure
        except Exception as e:
            logger.error(f"Failed to fetch structure variant for {pdb_id}: {e}")
            return None

    def _remove_ligands(self, structure: Structure.Structure, whitelist: Optional[List[str]] = None):
        """Remove non-standard residues considered ligands unless whitelisted."""
        if whitelist is None:
            whitelist = []
        whitelist = {w.upper() for w in whitelist}
        for model in structure:
            for chain in list(model):
                for residue in list(chain):
                    resname = residue.get_resname().strip().upper()
                    if (resname not in self.standard_amino_acids and
                        resname not in self.standard_nucleotides and
                        resname not in self.water_molecules and
                        resname not in whitelist):
                        chain.detach_child(residue.get_id())
    
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

    # ------------------------------------------------------------------
    # Minimal internal representation extraction for downstream heuristics
    # ------------------------------------------------------------------
    def to_internal_representation(self, structure: Structure.Structure) -> List[Dict[str, Any]]:
        """Convert a Bio.PDB structure to the lightweight format expected by
        structural extensions (list with one dict containing 'residues').

        Each residue entry contains:
            {
              'id': f"{chain_id}:{resseq}",
              'name': resname,
              'seq': resseq,
              'chain': chain_id,
              'atoms': [ {'name': atom_name, 'coord': (x,y,z)} ... ]
            }
        """
        residues = []
        try:
            model = structure[0]
        except Exception:
            return []
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                resname = residue.get_resname().strip()
                # Skip waters / non-standard if desired? Keep all – filters happen later.
                hetflag, resseq, icode = residue.get_id()
                # Build atom list with Cartesian coords
                atoms = []
                for atom in residue:
                    try:
                        c = atom.get_coord()
                        atoms.append({'name': atom.get_name(), 'coord': (float(c[0]), float(c[1]), float(c[2]))})
                    except Exception:
                        continue
                residues.append({
                    'id': f"{chain_id}:{resseq}",
                    'name': resname,
                    'seq': resseq,
                    'chain': chain_id,
                    'atoms': atoms
                })
        return [{ 'residues': residues }]
    
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
