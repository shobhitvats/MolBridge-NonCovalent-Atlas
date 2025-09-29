"""Unit tests for the shared pruning utilization heuristic `should_flag_kdtree`.

Ensures stable behavior around ratio and target boundaries so that
adaptive threshold tuning does not oscillate due to minor numerical noise.
"""
from importlib import reload

def test_should_flag_kdtree_basic_threshold(monkeypatch):
    # Force clean reload to ensure fresh state
    import utils.kdtree_thresholds as kt
    reload(kt)
    # Provide synthetic raw & candidate counts for different detectors
    # Using a detector with target density (e.g., hydro target=0.8)
    # Case 1: Heavy pruning (ratio well below 0.75 * target) should flag True
    target = kt._DENSITY_TARGETS['hydro']  # type: ignore[attr-defined]
    raw_pairs = 10000
    candidate_pairs = int(raw_pairs * target * 0.4)  # clearly below 0.75 band
    assert kt.should_flag_kdtree(candidate_pairs, raw_pairs, 'hydro') is True

    # Case 2: Slight pruning just above 0.75 * target should NOT flag
    boundary = target * 0.75
    candidate_pairs = int(raw_pairs * boundary * 1.05)  # a touch above boundary
    assert kt.should_flag_kdtree(candidate_pairs, raw_pairs, 'hydro') is False

    # Case 3: No pruning (candidate==raw) -> False
    assert kt.should_flag_kdtree(raw_pairs, raw_pairs, 'hydro') is False

    # Case 4: Unknown detector key falls back to generic 0.6 ratio rule
    raw_pairs = 5000
    pruned = int(raw_pairs * 0.55)  # below 0.6 -> True
    assert kt.should_flag_kdtree(pruned, raw_pairs, 'unknown_new_detector') is True
    unpruned = int(raw_pairs * 0.7)  # above 0.6 -> False
    assert kt.should_flag_kdtree(unpruned, raw_pairs, 'unknown_new_detector') is False

def test_should_flag_kdtree_zero_and_negative():
    import utils.kdtree_thresholds as kt
    # Zero raw pairs -> False (avoid division error / meaningless)
    assert kt.should_flag_kdtree(0, 0, 'hydro') is False
    # Negative raw (should not happen, defensive) -> False
    assert kt.should_flag_kdtree(10, -5, 'hydro') is False
