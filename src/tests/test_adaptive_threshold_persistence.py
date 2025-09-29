import os, json, tempfile
from importlib import reload

def test_adaptive_threshold_persistence(monkeypatch):
    # Create temp cache file path
    fd, path = tempfile.mkstemp(prefix='adaptive', suffix='.json')
    os.close(fd)
    monkeypatch.setenv('MOLBRIDGE_ADAPTIVE_CACHE_PATH', path)
    # Fresh import to trigger load with empty file
    import utils.kdtree_thresholds as kt
    reload(kt)
    # Simulate threshold change by calling adapt_threshold with forced kdtree usage scenario
    base = kt.get_threshold('hydro')
    new_val, changed, reason = kt.adapt_threshold('hydro', int(base*0.1), True)  # should raise threshold
    assert changed
    # Ensure file persisted
    with open(path, 'r') as f:
        data = json.load(f)
    assert 'hydro' in data and isinstance(data['hydro'], int)
    # Reload module to ensure value restored from cache
    reload(kt)
    restored = kt.get_threshold('hydro')
    assert restored == data['hydro']
