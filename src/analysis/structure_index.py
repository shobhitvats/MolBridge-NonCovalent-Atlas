"""Lightweight structure preprocessing indexes.

Currently provides an AromaticRingIndex used by vectorized detectors (π-π, cation–π, CH–π, etc.).

Design goals:
  * Zero external deps beyond numpy/Bio.PDB.
  * Pure-functional build helper returning small dataclass (cheap to discard / rebuild).
  * Non-breaking: detectors can still construct rings ad‑hoc if flag disabled.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Sequence
import numpy as np
from Bio.PDB import Residue, Atom


@dataclass
class AromaticRingRecord:
    residue: Residue.Residue
    ring_type: str
    atoms: List[Atom.Atom]
    center: np.ndarray  # (3,)
    normal: np.ndarray  # (3,)


@dataclass
class AromaticRingIndex:
    rings: List[AromaticRingRecord]
    centers: np.ndarray  # (N,3)
    normals: np.ndarray  # (N,3)

    def __len__(self):  # pragma: no cover - trivial
        return len(self.rings)


AROMATIC_ATOMS = {
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TRP': {
        'indole': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    },
    'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2']
}


def _ring_geometry(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (center, normal) for a ring set of coordinates.

    Uses SVD plane fit; ensures normal has positive z for consistency (avoids
    sign flips causing large angle deviations when comparing nearly parallel planes).
    """
    center = coords.mean(axis=0)
    centered = coords - center
    try:
        # For tiny rings (<3 atoms) fallback (shouldn't happen if caller guards)
        if centered.shape[0] < 3:
            normal = np.array([0.0, 0.0, 1.0])
        else:
            _, _, vh = np.linalg.svd(centered)
            normal = vh[-1]
    except Exception:  # pragma: no cover - extremely unlikely
        normal = np.array([0.0, 0.0, 1.0])
    if normal[2] < 0:
        normal = -normal
    # Normalize
    n = np.linalg.norm(normal) + 1e-12
    normal = normal / n
    return center, normal


def build_aromatic_ring_index(model) -> AromaticRingIndex:
    """Extract aromatic rings from a Bio.PDB model and return index.

    Current policy: treat TRP indole as a single fused ring (indole system) for speed.
    Could be extended later to split sub-rings with additional metadata.
    """
    rings: List[AromaticRingRecord] = []
    for chain in model:
        for residue in chain:
            resname = residue.get_resname()
            if resname in ('PHE', 'TYR', 'HIS'):
                atom_names: Sequence[str] = AROMATIC_ATOMS[resname]  # type: ignore
                atoms = [residue[name] for name in atom_names if name in residue]
                if len(atoms) < 3:
                    continue
                coords = np.array([a.get_coord() for a in atoms])
                center, normal = _ring_geometry(coords)
                rings.append(AromaticRingRecord(residue, resname, atoms, center, normal))
            elif resname == 'TRP':
                atom_names = AROMATIC_ATOMS['TRP']['indole']  # type: ignore
                atoms = [residue[name] for name in atom_names if name in residue]
                if len(atoms) < 6:
                    continue
                coords = np.array([a.get_coord() for a in atoms])
                center, normal = _ring_geometry(coords)
                rings.append(AromaticRingRecord(residue, 'TRP_indole', atoms, center, normal))
    if not rings:
        return AromaticRingIndex([], np.zeros((0, 3)), np.zeros((0, 3)))
    centers = np.vstack([r.center for r in rings])
    normals = np.vstack([r.normal for r in rings])
    return AromaticRingIndex(rings, centers, normals)
