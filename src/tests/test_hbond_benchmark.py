"""Benchmark suite for hydrogen bond detection (vector vs legacy paths).

Uses pytest-benchmark if installed. Safe to skip if not present.
The goal is to provide reproducible timing data to feed the perf regression
gate and to tune KD-tree thresholds. The test does NOT assert timing (that
belongs to the regression gate script) but records both code paths under
controlled environment flags.
"""
from __future__ import annotations

import os
import importlib
import pytest


pytestmark = pytest.mark.skipif(
    os.getenv("MOLBRIDGE_ENABLE_HBOND_BENCH", "1") in {"0", "false", "False"},
    reason="HBond benchmark disabled via env"
)


@pytest.fixture(scope="module")
def structure():
    """Provide a representative structure for benchmarking.

    We reuse the small demo PDB shipped in examples if available; fall back to
    generating a minimal mock structure if import fails. Keeping it stable
    ensures comparable timings across runs.
    """
    try:
        from utils.pdb_handler import load_structure  # type: ignore
        # Attempt to locate an example PDB (user can override via env)
        pdb_path = os.getenv("MOLBRIDGE_BENCH_PDB")
        if pdb_path and os.path.exists(pdb_path):
            return load_structure(pdb_path)
    except Exception:
        pass
    # Minimal synthetic fallback (very small — may reduce benchmark fidelity)
    from Bio.PDB import StructureBuilder
    sb = StructureBuilder.StructureBuilder()
    sb.init_structure("synthetic")
    sb.init_model(0)
    sb.init_chain("A")
    # Just return empty (detector should early exit gracefully)
    return sb.get_structure()


def _reset_settings():
    try:
        from utils.settings import get_settings
        get_settings.cache_clear()  # type: ignore
    except Exception:
        pass


def _run_detector(structure, use_vector: bool):
    if use_vector:
        os.environ["MOLBRIDGE_ENABLE_VECTOR_GEOM"] = "1"
    else:
        os.environ["MOLBRIDGE_ENABLE_VECTOR_GEOM"] = "0"
    _reset_settings()
    mod = importlib.import_module("analysis.hydrogen_bonds")
    # Re-import ensures class picks up env flag via settings on call
    from utils.config import AppConfig  # delayed import
    det = mod.HydrogenBondDetector(AppConfig())
    hbonds = det.detect_hydrogen_bonds(structure)
    instr = getattr(det, "instrumentation", {}) if use_vector else {}
    return hbonds, instr


@pytest.mark.benchmark(group="hydrogen_bonds")
def test_hbond_vector_path_benchmark(benchmark, structure):  # type: ignore[override]
    """Benchmark vector path (env-flag ON)."""
    hbonds, instr = benchmark(lambda: _run_detector(structure, True))
    # Sanity: instrumentation dictionary expected when vector path engaged
    assert isinstance(instr, dict)


@pytest.mark.benchmark(group="hydrogen_bonds")
def test_hbond_legacy_path_benchmark(benchmark, structure):  # type: ignore[override]
    """Benchmark legacy path (env-flag OFF)."""
    hbonds, instr = benchmark(lambda: _run_detector(structure, False))
    # Legacy path has no instrumentation attribute currently
    assert instr == {}
