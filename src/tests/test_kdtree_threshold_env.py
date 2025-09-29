"""Tests environment-variable tunable KD-tree thresholds for ionic & hydrophobic detectors.

The goal is to ensure kdtree_used toggles when threshold is forced low for a tiny
synthetic structure (forcing KD-tree) vs high (disabling it). Uses small mock
coordinate duplication to simulate size.
"""
from __future__ import annotations

import os
import importlib
import types
from utils.config import AppConfig
from utils.settings import get_settings
from analysis.ionic_interactions import IonicInteractionDetector
from analysis.hydrophobic_contacts import HydrophobicContactDetector
from Bio.PDB import StructureBuilder
import pytest


def build_synthetic(repeats: int = 50):
    sb = StructureBuilder.StructureBuilder()
    sb.init_structure('synthetic')
    sb.init_model(0)
    sb.init_chain('A')
    # We'll add alternating residues to create positive/negative & hydrophobic pattern
    # Use ALA (hydrophobic) and ASP/LYS (ionic donors/acceptors) for simplicity.
    from Bio.PDB import Atom
    import numpy as np
    for i in range(repeats):
        # ALA residue with CB atom
        sb.init_residue('ALA', ' ', i*3+1, ' ')
        coord = np.array([float(i), 0.0, 0.0])
        atom = Atom.Atom('CB', coord, 1.0, 1.0, ' ', 'CB', i*10+1, element='C')
        sb.structure[0]['A'][( ' ', i*3+1, ' ')].add(atom)
        # LYS residue with NZ atom
        sb.init_residue('LYS', ' ', i*3+2, ' ')
        coord2 = np.array([float(i), 0.5, 0.0])
        atom2 = Atom.Atom('NZ', coord2, 1.0, 1.0, ' ', 'NZ', i*10+2, element='N')
        sb.structure[0]['A'][( ' ', i*3+2, ' ')].add(atom2)
        # ASP residue with OD1 atom
        sb.init_residue('ASP', ' ', i*3+3, ' ')
        coord3 = np.array([float(i), 1.0, 0.0])
        atom3 = Atom.Atom('OD1', coord3, 1.0, 1.0, ' ', 'OD1', i*10+3, element='O')
        sb.structure[0]['A'][( ' ', i*3+3, ' ')].add(atom3)
    return sb.get_structure()


@pytest.mark.parametrize("detector_cls, env_name", [
    (IonicInteractionDetector, 'MOLBRIDGE_KDTREE_IONIC_THRESHOLD'),
    (HydrophobicContactDetector, 'MOLBRIDGE_KDTREE_HYDRO_THRESHOLD'),
])
def test_kdtree_env_threshold_toggle(detector_cls, env_name):
    structure = build_synthetic(repeats=20)  # ensures above minimal threshold values
    from utils.settings import get_settings
    # Force vector path
    os.environ['MOLBRIDGE_ENABLE_VECTOR_GEOM'] = '1'
    get_settings.cache_clear()  # type: ignore
    cfg = AppConfig()

    # High threshold disables KD-tree
    os.environ[env_name] = '999999'
    get_settings.cache_clear()  # type: ignore
    det = detector_cls(cfg)
    det.detect_ionic_interactions(structure) if detector_cls is IonicInteractionDetector else det.detect_hydrophobic_contacts(structure)
    instr1 = getattr(det, 'instrumentation', {})

    # Low threshold enables KD-tree
    os.environ[env_name] = '1'
    get_settings.cache_clear()  # type: ignore
    det2 = detector_cls(cfg)
    det2.detect_ionic_interactions(structure) if detector_cls is IonicInteractionDetector else det2.detect_hydrophobic_contacts(structure)
    instr2 = getattr(det2, 'instrumentation', {})

    assert instr1.get('kdtree_used') != instr2.get('kdtree_used'), 'KD-tree usage flag should toggle with threshold extremes'
