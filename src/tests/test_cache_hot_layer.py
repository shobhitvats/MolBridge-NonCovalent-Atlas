"""Smoke tests for CacheManager hot cache and batch functionality."""
from pathlib import Path
from utils.cache import CacheManager


def test_hot_cache_roundtrip(tmp_path: Path):
    cm = CacheManager(cache_dir=tmp_path, hot_cache_size=4)

    params = {"cutoff": 5.0}
    data = {"value": 123, "items": [1,2,3]}
    cm.cache_analysis_result("1ABC", "hydrogen_bonds", params, data, compress=False)

    # First fetch should populate hot cache
    r1 = cm.get_analysis_result("1ABC", "hydrogen_bonds", params)
    assert r1 == data

    # Internal hot cache should now have the key
    key_count_before = cm.get_cache_stats()["hot_cache_items"]
    assert key_count_before >= 1

    # Fetch again, should come directly from hot cache (behavioral - just ensure same result)
    r2 = cm.get_analysis_result("1ABC", "hydrogen_bonds", params)
    assert r2 == data


def test_batch_context(tmp_path: Path):
    cm = CacheManager(cache_dir=tmp_path, hot_cache_size=0)
    params = {"p": 1}
    # Use batch to queue writes
    with cm.batch():
        for i in range(5):
            cm.cache_analysis_result(f"P{i}", "pi_stack", params, {"i": i}, compress=False)
    # After context, entries should exist
    for i in range(5):
        assert cm.get_analysis_result(f"P{i}", "pi_stack", params) == {"i": i}
