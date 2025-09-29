"""Parity test for π-π detector vector vs legacy paths.

We build a minimal synthetic structure containing two stacked phenyl rings
so that a face_to_face interaction is detected. The vector fast-path is
gated by the MOLBRIDGE_ENABLE_VECTOR_GEOM env flag. We flip it on/off and
compare the resulting interaction sets (ignoring ordering and floating
precision on continuous fields).
"""
import os
import math
import numpy as np
from Bio.PDB.StructureBuilder import StructureBuilder
from analysis.pi_pi_stacking import PiPiDetector
from utils.settings import get_settings


RING_ATOMS = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']


def _hexagon(radius: float = 1.4, z: float = 0.0):
    coords = []
    for k in range(6):
        theta = k * math.pi / 3.0
        coords.append((radius * math.cos(theta), radius * math.sin(theta), z))
    return coords


def build_two_phe_structure(distance: float = 3.6):
    sb = StructureBuilder()
    sb.init_structure("TEST")
    sb.init_model(0)
    sb.init_chain("A")
    sb.init_seg("    ")
    # Residue 1
    sb.init_residue("PHE", " ", 1, " ")
    coords1 = _hexagon(z=0.0)
    for idx, name in enumerate(RING_ATOMS):
        x, y, z = coords1[idx]
        sb.init_atom(name, np.array([x, y, z], dtype=float), 1.0, 1.0, " ", name, idx + 1, element=name[0])
    # Residue 2 stacked above
    sb.init_residue("PHE", " ", 2, " ")
    coords2 = _hexagon(z=distance)
    for idx, name in enumerate(RING_ATOMS):
        x, y, z = coords2[idx]
        sb.init_atom(name, np.array([x, y, z], dtype=float), 1.0, 1.0, " ", name, idx + 100, element=name[0])
    return sb.get_structure()


def _extract_pairs(result):
    return {tuple(sorted((i.ring1_residue, i.ring2_residue))) for i in result}


def test_pipi_vector_parity():
    structure = build_two_phe_structure()
    from utils.config import load_config
    config = load_config()

    # Legacy path
    os.environ['MOLBRIDGE_ENABLE_VECTOR_GEOM'] = '0'
    get_settings.cache_clear()  # type: ignore
    legacy_detector = PiPiDetector(config)
    legacy = legacy_detector.detect_pi_pi_interactions(structure)
    assert legacy, "Legacy path failed to detect expected π-π interaction"

    # Vector path
    os.environ['MOLBRIDGE_ENABLE_VECTOR_GEOM'] = '1'
    get_settings.cache_clear()  # type: ignore
    vector_detector = PiPiDetector(config)
    vector = vector_detector.detect_pi_pi_interactions(structure)
    assert vector, "Vector path failed to detect expected π-π interaction"

    assert _extract_pairs(legacy) == _extract_pairs(vector)
    # Basic sanity on geometry similarity (distance diff < 1e-6)
    assert abs(legacy[0].distance - vector[0].distance) < 1e-6
