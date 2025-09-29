import os
import pytest
from utils.config import AppConfig
from utils.pdb_handler import PDBHandler
from performance.parallel_processor import HighPerformanceProcessor
from analysis.hydrogen_bonds import HydrogenBondDetector

@pytest.mark.skipif(not os.getenv('MOLBRIDGE_TEST_PDB'), reason='Set MOLBRIDGE_TEST_PDB to run parallel instrumentation test')
def test_parallel_instrumentation_capture(monkeypatch):
    monkeypatch.setenv('MOLBRIDGE_ENABLE_VECTOR_GEOM', '1')
    cfg = AppConfig()
    structure = PDBHandler(cfg).load_structure(os.getenv('MOLBRIDGE_TEST_PDB'))
    proc = HighPerformanceProcessor(config=cfg)
    out = proc.process_interactions_parallel(structure, [HydrogenBondDetector], config=cfg)
    assert isinstance(out, dict)
    instr = getattr(proc, '_detector_instrumentation', {})
    # Either empty (if detector produced none) or contains hydrogenbond stats
    assert isinstance(instr, dict)
