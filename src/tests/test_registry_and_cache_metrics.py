import os
from pathlib import Path
from utils.cache import CacheManager
from analysis.registry import DETECTOR_REGISTRY, get_detector


def test_registry_contains_expected_core_keys():
    expected = {"hydrogenbond", "halogenbond", "chpi", "pipi", "ionicinteraction"}
    assert expected.issubset(set(DETECTOR_REGISTRY.keys()))


def test_get_detector_returns_instance_and_method():
    inst, method = get_detector("hydrogenbond")
    assert inst is not None
    assert isinstance(method, str)
    assert method


def test_cache_metrics_increment(tmp_path: Path):
    cm = CacheManager(tmp_path)
    params = {"a": 1}
    # Miss first
    assert cm.get_analysis_result("X", "hydrogenbond", params) is None
    m1 = cm.get_cache_metrics()
    assert m1["analysis_misses"] == 1
    # Set
    cm.cache_analysis_result("X", "hydrogenbond", params, [1,2,3])
    # Hit
    assert cm.get_analysis_result("X", "hydrogenbond", params) is not None
    stats = cm.get_cache_stats()
    assert stats["analysis_hits"] >= 1
    assert "analysis_hit_ratio" in stats
