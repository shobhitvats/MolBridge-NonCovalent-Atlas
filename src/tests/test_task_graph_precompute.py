import os
import pytest
from analysis.feature_store import get_feature_store
from performance.parallel_processor import HighPerformanceProcessor
from utils.config import AppConfig
from utils.pdb_handler import PDBHandler

@pytest.mark.skipif(not os.getenv('MOLBRIDGE_TEST_PDB'), reason="Set MOLBRIDGE_TEST_PDB to a small PDB ID to run this")
def test_task_graph_precompute_smoke(monkeypatch):
    pdb_id = os.getenv('MOLBRIDGE_TEST_PDB')
    cfg = AppConfig()
    handler = PDBHandler(cfg)
    structure = handler.load_structure(pdb_id)
    fs = get_feature_store()
    coords = fs.ensure_coords(structure)
    rings = fs.ensure_rings(structure)
    charged = fs.ensure_charged_centers(structure)
    hb = fs.ensure_hbond_participants(structure)
    assert coords is None or coords.ndim == 2
    assert isinstance(rings, list)
    assert isinstance(charged, list)
    assert isinstance(hb, tuple) and len(hb) == 2

    # Ensure processor doesn't crash with task graph enabled
    monkeypatch.setenv('MOLBRIDGE_TASK_GRAPH', '1')
    proc = HighPerformanceProcessor(config=cfg)
    # No detectors pass (empty) just to exercise path
    out = proc.process_interactions_parallel(structure, [], config=cfg)
    assert isinstance(out, dict)
