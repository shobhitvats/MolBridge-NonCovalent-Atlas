import math
from utils.config import AppConfig
from analysis.registry import DETECTOR_REGISTRY
from performance.parallel_processor import HighPerformanceProcessor

# Basic invariant tests for funnel metrics across detectors
# Ensures schema presence and logical constraints (candidate <= raw, accepted <= candidate, ratios within [0,1])

def test_funnel_invariants_basic():
    config = AppConfig()
    # Use a small synthetic-like structure by reusing an existing builder path if available
    # For simplicity, construct from first available detector's expected usage
    # We'll pick hydrogenbond detector structure expectation: create minimal mock using Bio.PDB if present.
    from Bio.PDB import StructureBuilder
    sb = StructureBuilder.StructureBuilder()
    sb.init_structure('X')
    sb.init_model(0)
    sb.init_chain('A')
    # build 10 residues with basic backbone atoms to feed various detectors
    import numpy as np
    for i in range(10):
        sb.init_seg('    ')
        sb.init_residue('GLY', ' ', i+1, ' ')
        base = np.array([i*1.5, 0.0, 0.0])
        for name, elem, offset in [('N','N',0.0),('O','O',0.6),('C','C',1.0),('CA','C',0.3)]:
            sb.init_atom(name, base + np.array([0.0, offset, 0.0]), 1.0, 1.0, ' ', name, i*10 + len(name), elem)
    structure = sb.get_structure()

    # Select a subset of detectors for quick invariants (all registry entries if lightweight)
    detector_classes = []
    for key in ['hydrogenbond','pipi','ionicinteraction','hydrophobiccontact','halogenbond','chpi']:
        cls = DETECTOR_REGISTRY.get(key)
        if cls:
            detector_classes.append(cls)
    processor = HighPerformanceProcessor(config=config)
    processor.process_interactions_parallel(structure, detector_classes, config=config)
    instr = getattr(processor, '_detector_instrumentation', {})
    assert instr, 'No instrumentation captured'
    for name, metrics in instr.items():
        if not isinstance(metrics, dict):
            continue
        raw_pairs = metrics.get('raw_pairs', 0)
        candidate_pairs = metrics.get('candidate_pairs', 0)
        accepted_pairs = metrics.get('accepted_pairs', 0)
        acceptance_ratio = metrics.get('acceptance_ratio', 0)
        # Logical bounds
        assert raw_pairs >= 0
        assert candidate_pairs >= 0
        assert accepted_pairs >= 0
        assert candidate_pairs <= raw_pairs or raw_pairs == 0, f"candidate_pairs > raw_pairs for {name}"  # allow zero raw edge
        assert accepted_pairs <= candidate_pairs, f"accepted_pairs > candidate_pairs for {name}"
        if candidate_pairs > 0:
            assert 0.0 <= acceptance_ratio <= 1.0 + 1e-6, f"acceptance_ratio out of range for {name}"  # slight float tolerance
        # If ratio provided, recompute close
        if candidate_pairs > 0 and 'acceptance_ratio' in metrics:
            comp = accepted_pairs / candidate_pairs if candidate_pairs else 0.0
            assert math.isclose(comp, acceptance_ratio, rel_tol=1e-3, abs_tol=1e-3), f"acceptance_ratio mismatch for {name}"  
