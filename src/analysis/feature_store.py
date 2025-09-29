"""Ephemeral per-structure feature cache (geometry primitives scaffold).

Extended: now includes aromatic ring centroid & normal extraction; future
additions will cover charged group centers & spatial indexes.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional
import weakref
import math
try:
    from scipy.spatial import cKDTree  # optional
    _HAVE_SCIPY = True
except Exception:  # pragma: no cover
    cKDTree = None  # type: ignore
    _HAVE_SCIPY = False


@dataclass
class FeatureBundle:
    coords: Optional[Any] = None  # numpy ndarray (N,3)
    # Ring features
    rings: Optional[Any] = None   # legacy list-of-dicts (lazy constructed if accessed)
    ring_centroids: Optional[Any] = None  # np.ndarray (R,3)
    ring_normals: Optional[Any] = None    # np.ndarray (R,3)
    ring_resnames: Optional[Any] = None   # list[str]
    ring_chain_ids: Optional[Any] = None  # list[str]
    ring_residues: Optional[Any] = None   # list[Residue]
    charged_group_centroids: Optional[Any] = None
    acidic_group_centroids: Optional[Any] = None
    hb_acceptors: Optional[Any] = None
    hb_donors: Optional[Any] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    # Spatial index cache
    kdtree: Optional[Any] = None
    kdtree_size: int = 0
    kdtree_build_time: float | None = None
    # Simple neighbor query stats
    neighbor_queries: int = 0
    neighbor_hits: int = 0
    # Cached neighbor pair lists for frequently requested radii (global atom mode) with rudimentary LRU metadata
    neighbor_cache: Dict[float, Any] | None = None  # radius -> {'pairs': list[(i,j)], 'age': int}
    neighbor_cache_clock: int = 0
    neighbor_cache_max_entries: int = 8  # avoid unbounded memory growth
    # Residue / atom classification caches
    residue_flags: Optional[Any] = None  # numpy uint8 flags per residue
    atom_flags: Optional[Any] = None     # numpy uint8 flags per atom
    residue_index_map: Optional[Dict[int,int]] = None  # residue id(structure) -> index
    atom_index_map: Optional[Dict[int,int]] = None      # atom id(structure) -> index
    residue_meta: Optional[Any] = None  # list of (chain,resname,resid)


class FeatureStore:
    def __init__(self):
        self._store: "weakref.WeakValueDictionary[int, FeatureBundle]" = weakref.WeakValueDictionary()

    def get_bundle(self, structure) -> FeatureBundle:
        sid = id(structure)
        bundle = self._store.get(sid)
        if bundle is None:
            bundle = FeatureBundle(metadata={"structure_id": sid})
            self._store[sid] = bundle
        return bundle

    def ensure_coords(self, structure):  # placeholder; vectorize later
        b = self.get_bundle(structure)
        if b.coords is None:
            try:
                import numpy as np
                coords = [atom.get_coord() for atom in structure.get_atoms()]  # type: ignore[attr-defined]
                b.coords = np.asarray(coords, dtype='float32') if coords else None
            except Exception:
                b.coords = None
        return b.coords

    # ---------------- Classification Caching -----------------
    def ensure_classification(self, structure):
        """Assign lightweight bitmask flags to residues and atoms once per structure.

        Bit layout (uint8):
          bit0: standard_amino_acid
          bit1: aromatic_ring_residue
          bit2: charged_positive (LYS/ARG/HIS)
          bit3: charged_negative (ASP/GLU)
          bit4: sulfur_containing (CYS/MET)
        """
        b = self.get_bundle(structure)
        if b.residue_flags is not None and b.atom_flags is not None:
            return b.residue_flags, b.atom_flags
        import numpy as np
        residue_flags: list[int] = []
        residue_meta: list[tuple[str,str,int]] = []
        residue_index_map: dict[int,int] = {}
        atom_flags: list[int] = []
        atom_index_map: dict[int,int] = {}
        aromatic = {'PHE','TYR','TRP','HIS'}
        pos = {'LYS','ARG','HIS'}
        neg = {'ASP','GLU'}
        sulfur = {'CYS','MET'}
        std = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
        try:
            model = structure[0]
        except Exception:
            b.residue_flags = np.zeros((0,), dtype='uint8')
            b.atom_flags = np.zeros((0,), dtype='uint8')
            b.residue_index_map = {}
            b.atom_index_map = {}
            b.residue_meta = []
            return b.residue_flags, b.atom_flags
        r_idx = 0
        a_idx = 0
        for chain in model:
            cid = chain.get_id()
            for residue in chain:
                resname = residue.get_resname().strip()
                res_id = residue.get_id()[1]
                flag = 0
                if resname in std: flag |= 1
                if resname in aromatic: flag |= (1<<1)
                if resname in pos: flag |= (1<<2)
                if resname in neg: flag |= (1<<3)
                if resname in sulfur: flag |= (1<<4)
                residue_flags.append(flag)
                residue_meta.append((cid,resname,res_id))
                residue_index_map[id(residue)] = r_idx
                r_idx += 1
                for atom in residue:
                    # Basic atom flags derivation from residue flags for now (could refine with atom roles)
                    atom_flags.append(flag)
                    atom_index_map[id(atom)] = a_idx
                    a_idx += 1
        b.residue_flags = np.asarray(residue_flags, dtype='uint8')
        b.atom_flags = np.asarray(atom_flags, dtype='uint8')
        b.residue_index_map = residue_index_map
        b.atom_index_map = atom_index_map
        b.residue_meta = residue_meta
        return b.residue_flags, b.atom_flags

    # ---------------- Spatial Index / Neighbor API -----------------
    def ensure_kdtree(self, structure, rebuild: bool = False):
        """Build (or reuse) a shared KD-tree over atomic coordinates.

        Rebuild automatically if atom count changed or explicit rebuild flag.
        """
        if not _HAVE_SCIPY:
            return None
        import time as _t
        b = self.get_bundle(structure)
        coords = self.ensure_coords(structure)
        if coords is None or coords.size == 0:
            return None
        n = coords.shape[0]
        if (not rebuild) and b.kdtree is not None and b.kdtree_size == n:
            return b.kdtree
        t0 = _t.time()
        try:
            b.kdtree = cKDTree(coords, balanced_tree=False, compact_nodes=True)  # type: ignore[arg-type]
            b.kdtree_size = n
            b.kdtree_build_time = _t.time() - t0
        except Exception:  # pragma: no cover
            b.kdtree = None
        return b.kdtree

    def neighbor_within(self, structure, radius: float, points=None):
        """Return pairs (i,j) within ``radius``.

        Two modes:
          1. Global atom mode (points is None): pairs of atom indices over the full
             structure coordinate array (original behavior).
          2. Subset mode (points is np.ndarray (M,3)): pairs of *subset-local* indices
             using only the provided coordinate set. This allows detectors to prune
             candidates for a filtered atom/centroid set (e.g., hydrophobic atoms,
             ring centroids) while reusing the same adaptive KD-tree heuristics.

        Strategy:
          * If SciPy available and size >= 512 -> use cKDTree ball query
          * Else fall back to vectorized brute force O(M^2)

        Returns list[tuple[int,int]] with i<j.
        """
        import numpy as np
        b = self.get_bundle(structure)
        b.neighbor_queries += 1
        # Subset / custom points mode ---------------------------------
        if points is not None:
            arr = np.asarray(points, dtype='float32')
            m = arr.shape[0] if arr.ndim == 2 else 0
            if m == 0:
                return []
            if not _HAVE_SCIPY or m < 512:
                diff = arr[:, None, :] - arr[None, :, :]
                dist2 = (diff * diff).sum(axis=-1)
                r2 = float(radius) * float(radius)
                ii, jj = np.where((dist2 <= r2) & (dist2 > 0))
                mask = ii < jj
                pairs = list(zip(ii[mask].tolist(), jj[mask].tolist()))
                b.neighbor_hits += len(pairs)
                return pairs
            # KD-tree subset (ephemeral – not cached on bundle to avoid polluting global tree)
            try:
                tree = cKDTree(arr, balanced_tree=False, compact_nodes=True)  # type: ignore[arg-type]
                lists = tree.query_ball_tree(tree, r=float(radius))
                pairs: list[tuple[int,int]] = []
                for i, nbrs in enumerate(lists):
                    for j in nbrs:
                        if j > i:
                            pairs.append((i, j))
                b.neighbor_hits += len(pairs)
                return pairs
            except Exception:  # pragma: no cover
                return []
        # Global atom mode ---------------------------------------------
        coords = self.ensure_coords(structure)
        if coords is None or coords.size == 0:
            return []
        n = coords.shape[0]
        # Global mode radius cache (only when not using SciPy path to avoid large memory duplication of sparse sets)
        if points is None:
            b = self.get_bundle(structure)
            if b.neighbor_cache is None:
                b.neighbor_cache = {}
            entry = b.neighbor_cache.get(float(radius)) if b.neighbor_cache else None
            if entry and isinstance(entry, dict) and 'pairs' in entry:
                entry['age'] = b.neighbor_cache_clock = b.neighbor_cache_clock + 1
                pairs_cached = entry['pairs']
                b.neighbor_hits += len(pairs_cached)
                return pairs_cached
        if not _HAVE_SCIPY or n < 512:
            diff = coords[:, None, :] - coords[None, :, :]
            dist2 = (diff * diff).sum(axis=-1)
            r2 = float(radius) * float(radius)
            ii, jj = np.where((dist2 <= r2) & (dist2 > 0))
            mask = ii < jj
            pairs = list(zip(ii[mask].tolist(), jj[mask].tolist()))
            b.neighbor_hits += len(pairs)
            if points is None and b.neighbor_cache is not None and len(pairs) < 2_000_000:  # cap to avoid huge memory
                try:
                    # Insert with age counter & LRU eviction if capacity exceeded
                    b.neighbor_cache_clock += 1
                    b.neighbor_cache[float(radius)] = {'pairs': pairs, 'age': b.neighbor_cache_clock}
                    if len(b.neighbor_cache) > b.neighbor_cache_max_entries:
                        # Evict oldest age
                        oldest_r = None
                        oldest_age = 10**18
                        for r, meta in b.neighbor_cache.items():
                            age = meta.get('age', 0) if isinstance(meta, dict) else 0
                            if age < oldest_age:
                                oldest_age = age
                                oldest_r = r
                        if oldest_r is not None:
                            b.neighbor_cache.pop(oldest_r, None)
                except Exception:
                    pass
            return pairs
        tree = self.ensure_kdtree(structure)
        if tree is None:
            return []
        try:
            lists = tree.query_ball_tree(tree, r=radius)
            pairs = []
            for i, nbrs in enumerate(lists):
                for j in nbrs:
                    if j > i:
                        pairs.append((i, j))
            b.neighbor_hits += len(pairs)
            return pairs
        except Exception:  # pragma: no cover
            return []

    def neighbor_between(self, structure, radius: float, points_a, points_b):
        """Return cross-set pairs (i,j) with i over points_a, j over points_b within radius.

        Uses cKDTree if available & sufficiently large, else vectorized brute force.
        Indices are local to provided arrays (0..len(points_a)-1, 0..len(points_b)-1).
        """
        import numpy as np
        a = np.asarray(points_a, dtype='float32')
        b = np.asarray(points_b, dtype='float32')
        if a.ndim != 2 or b.ndim != 2 or a.shape[0] == 0 or b.shape[0] == 0:
            return []
        r = float(radius)
        if not _HAVE_SCIPY or (a.shape[0] + b.shape[0]) < 512:
            diff = a[:, None, :] - b[None, :, :]
            dist2 = (diff * diff).sum(axis=-1)
            r2 = r * r
            ii, jj = np.where(dist2 <= r2)
            return list(zip(ii.tolist(), jj.tolist()))
        try:
            tree_a = cKDTree(a, balanced_tree=False, compact_nodes=True)  # type: ignore[arg-type]
            tree_b = cKDTree(b, balanced_tree=False, compact_nodes=True)  # type: ignore[arg-type]
            lists = tree_a.query_ball_tree(tree_b, r=r)
            pairs: list[tuple[int,int]] = []
            for i, nbrs in enumerate(lists):
                for j in nbrs:
                    pairs.append((i, j))
            return pairs
        except Exception:  # pragma: no cover
            return []

    # ---------------- Shared snapshot (for potential multiprocess sharing) -----------------
    def snapshot_bundle(self, structure):
        """Produce a lightweight snapshot of frequently reused arrays.

        This is a shallow copy (arrays reused) intended for read-only sharing in
        multiprocess scenarios. Current consumer stubs may attach this to a task
        payload; no in-place mutation should be performed on returned data.
        """
        b = self.get_bundle(structure)
        snap = {
            'coords': b.coords,
            'ring_centroids': b.ring_centroids,
            'ring_normals': b.ring_normals,
            'charged_centroids': getattr(b, 'charged_group_centroids', None),
            'acidic_centroids': getattr(b, 'acidic_group_centroids', None),
            'residue_flags': b.residue_flags,
            'atom_flags': b.atom_flags,
        }
        return snap

    # Basic hydrogen bond donor/acceptor extraction (precompute skeleton) -- can be expanded
    def ensure_hbond_participants(self, structure):
        b = self.get_bundle(structure)
        if b.hb_acceptors is not None and b.hb_donors is not None:
            return b.hb_donors, b.hb_acceptors
        from collections import defaultdict
        hb_donors = []
        hb_acceptors = []
        try:
            model = structure[0]
        except Exception:
            b.hb_donors, b.hb_acceptors = [], []
            return b.hb_donors, b.hb_acceptors
        donor_elements = {'N','O','S'}  # simplistic; refine later
        acceptor_elements = {'O','N','S'}
        for chain in model:
            for residue in chain:
                for atom in residue:
                    el = getattr(atom, 'element', '').upper()
                    if el in donor_elements:
                        hb_donors.append(atom)
                    if el in acceptor_elements:
                        hb_acceptors.append(atom)
        b.hb_donors = hb_donors
        b.hb_acceptors = hb_acceptors
        return b.hb_donors, b.hb_acceptors

    def ensure_rings(self, structure):
        b = self.get_bundle(structure)
        # If SoA already built, optionally synthesize legacy list if first legacy access
        if b.ring_centroids is not None and b.rings is not None:
            return b.rings
        import numpy as np
        ring_defs = {
            'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TRP': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
            'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
        }
        centroids = []
        normals = []
        resnames = []
        chain_ids = []
        residues = []
        ring_dicts = []
        try:
            model = structure[0]
        except Exception:
            b.ring_centroids = np.zeros((0,3), dtype='float32')
            b.ring_normals = np.zeros((0,3), dtype='float32')
            b.ring_resnames = []
            b.ring_chain_ids = []
            b.ring_residues = []
            b.rings = []
            return b.rings
        for chain in model:
            cid = chain.get_id()
            for residue in chain:
                rname = residue.get_resname()
                if rname not in ring_defs:
                    continue
                atoms = [residue[an] for an in ring_defs[rname] if an in residue]
                if len(atoms) < 3:
                    continue
                coords = np.vstack([a.get_coord() for a in atoms])
                centroid = coords.mean(axis=0)
                try:
                    cov = np.cov(coords.T)
                    eigvals, eigvecs = np.linalg.eigh(cov)
                    normal = eigvecs[:, eigvals.argmin()]
                except Exception:
                    normal = np.array([0.0, 0.0, 1.0])
                norm_unit = normal / (np.linalg.norm(normal) or 1.0)
                centroids.append(centroid.astype('float32'))
                normals.append(norm_unit.astype('float32'))
                resnames.append(rname)
                chain_ids.append(cid)
                residues.append(residue)
                ring_dicts.append({
                    'residue': residue,
                    'resname': rname,
                    'chain_id': cid,
                    'centroid': centroid,
                    'normal': norm_unit,
                    'atom_ids': [id(a) for a in atoms]
                })
        b.ring_centroids = np.vstack(centroids) if centroids else np.zeros((0,3), dtype='float32')
        b.ring_normals = np.vstack(normals) if normals else np.zeros((0,3), dtype='float32')
        b.ring_resnames = resnames
        b.ring_chain_ids = chain_ids
        b.ring_residues = residues
        # Lazy: only build list-of-dicts if existing code expects it
        b.rings = ring_dicts
        return b.rings

    # ---- Charged group centroids (for cation-π, salt bridges etc.) ----
    def ensure_charged_centers(self, structure):
        """Cache centroids for positively charged side-chain groups.

        Returns list of dicts: {residue, resname, chain_id, centroid}
        Currently supports LYS (NZ/CE/CD), ARG (CZ/NH1/NH2), HIS (ring centroid heuristic).
        """
        b = self.get_bundle(structure)
        if getattr(b, 'charged_group_centroids', None) is not None:
            return b.charged_group_centroids
        import numpy as np
        try:
            model = structure[0]
        except Exception:
            b.charged_group_centroids = []
            return b.charged_group_centroids
        results = []
        for chain in model:
            cid = chain.get_id()
            for residue in chain:
                rname = residue.get_resname()
                centroid = None
                if rname == 'LYS':
                    atom_names = [n for n in ['NZ','CE','CD'] if n in residue]
                elif rname == 'ARG':
                    atom_names = [n for n in ['CZ','NH1','NH2'] if n in residue]
                elif rname == 'HIS':
                    atom_names = [n for n in ['CG','ND1','CD2','CE1','NE2'] if n in residue]
                else:
                    continue
                if not atom_names:
                    continue
                coords = [residue[n].get_coord() for n in atom_names]
                if coords:
                    centroid = np.mean(np.vstack(coords), axis=0)
                if centroid is not None:
                    results.append({
                        'residue': residue,
                        'resname': rname,
                        'chain_id': cid,
                        'centroid': centroid
                    })
        b.charged_group_centroids = results
        return b.charged_group_centroids

    # ---- Acidic (negatively charged) group centroids (ASP/GLU) ----
    def ensure_acidic_centers(self, structure):
        """Cache centroids for negatively charged side-chain groups (carboxylates).

        Returns list of dicts: {residue, resname, chain_id, centroid}
        ASP: OD1/OD2/CG, GLU: OE1/OE2/CD; fall back to OD1/OE1 only if others missing.
        """
        b = self.get_bundle(structure)
        if getattr(b, 'acidic_group_centroids', None) is not None:
            return b.acidic_group_centroids
        import numpy as np
        try:
            model = structure[0]
        except Exception:
            b.acidic_group_centroids = []
            return b.acidic_group_centroids
        results = []
        for chain in model:
            cid = chain.get_id()
            for residue in chain:
                rname = residue.get_resname()
                if rname == 'ASP':
                    atom_names = [n for n in ['OD1','OD2','CG'] if n in residue]
                elif rname == 'GLU':
                    atom_names = [n for n in ['OE1','OE2','CD'] if n in residue]
                else:
                    continue
                if not atom_names:
                    continue
                coords = [residue[n].get_coord() for n in atom_names]
                if coords:
                    centroid = np.mean(np.vstack(coords), axis=0)
                    results.append({
                        'residue': residue,
                        'resname': rname,
                        'chain_id': cid,
                        'centroid': centroid
                    })
        b.acidic_group_centroids = results
        return b.acidic_group_centroids


_GLOBAL_FEATURE_STORE = FeatureStore()


def get_feature_store() -> FeatureStore:
    return _GLOBAL_FEATURE_STORE
