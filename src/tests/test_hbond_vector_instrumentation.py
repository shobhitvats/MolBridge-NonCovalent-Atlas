import os
import pytest
from utils.config import AppConfig
from utils.pdb_handler import PDBHandler
from analysis.hydrogen_bonds import HydrogenBondDetector
from utils.settings import get_settings

@pytest.mark.skipif(not os.getenv('MOLBRIDGE_TEST_PDB'), reason='Set MOLBRIDGE_TEST_PDB to run this instrumentation test')
def test_hbond_vector_instrumentation(monkeypatch):
    monkeypatch.setenv('MOLBRIDGE_ENABLE_VECTOR_GEOM', '1')
    pdb_id = os.getenv('MOLBRIDGE_TEST_PDB')
    cfg = AppConfig()
    s = PDBHandler(cfg).load_structure(pdb_id)
    det = HydrogenBondDetector(cfg)
    bonds = det.detect_hydrogen_bonds(s)
    assert isinstance(bonds, list)
    assert hasattr(det, 'instrumentation')
    instr = det.instrumentation
    assert 'donors' in instr and 'acceptors' in instr and 'candidate_pairs' in instr
