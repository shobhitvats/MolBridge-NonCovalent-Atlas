"""Parity test for cation–π detector vector vs legacy paths.

Creates a minimal structure containing one LYS (as cation) and one PHE ring
within distance cutoff to ensure detection in both paths.
"""
import os
import math
import numpy as np
from Bio.PDB.StructureBuilder import StructureBuilder
from analysis.cation_pi_interactions import CationPiDetector
from utils.settings import get_settings
from utils.config import load_config

RING_ATOMS = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']


def _hexagon(radius: float = 1.4, center=(0.0, 0.0, 0.0)):
    cx, cy, cz = center
    coords = []
    for k in range(6):
        theta = k * math.pi / 3.0
        coords.append((cx + radius * math.cos(theta), cy + radius * math.sin(theta), cz))
    return coords


def build_lys_phe(distance: float = 4.0):
    sb = StructureBuilder()
    sb.init_structure("TEST")
    sb.init_model(0)
    sb.init_chain("A")
    sb.init_seg("    ")
    # PHE ring base at z=0
    sb.init_residue("PHE", " ", 1, " ")
    for idx, name in enumerate(RING_ATOMS):
        x, y, z = _hexagon()[idx]
        sb.init_atom(name, np.array([x, y, z], dtype=float), 1.0, 1.0, " ", name, idx + 1, element=name[0])
    # LYS side chain simplified: place NZ at specified distance along z, with CE/CD as nearby to mimic centroid
    sb.init_residue("LYS", " ", 2, " ")
    nz_coord = np.array([0.0, 0.0, distance])
    ce_coord = np.array([0.5, 0.0, distance - 0.5])
    cd_coord = np.array([-0.5, 0.0, distance - 0.5])
    sb.init_atom("NZ", nz_coord, 1.0, 1.0, " ", "NZ", 101, element='N')
    sb.init_atom("CE", ce_coord, 1.0, 1.0, " ", "CE", 102, element='C')
    sb.init_atom("CD", cd_coord, 1.0, 1.0, " ", "CD", 103, element='C')
    return sb.get_structure()


def _pairs(result):
    return {tuple(sorted((r['residue1'], r['residue2']))) for r in result}


def test_cation_pi_vector_parity():
    structure = build_lys_phe()
    config = load_config()
    os.environ['MOLBRIDGE_ENABLE_VECTOR_GEOM'] = '0'
    get_settings.cache_clear()  # type: ignore
    legacy_detector = CationPiDetector(config)
    legacy = legacy_detector.detect_cation_pi(structure)
    assert legacy, "Legacy path failed to detect cation–π"

    os.environ['MOLBRIDGE_ENABLE_VECTOR_GEOM'] = '1'
    get_settings.cache_clear()  # type: ignore
    vector_detector = CationPiDetector(config)
    vector = vector_detector.detect_cation_pi(structure)
    assert vector, "Vector path failed to detect cation–π"

    assert _pairs(legacy) == _pairs(vector)