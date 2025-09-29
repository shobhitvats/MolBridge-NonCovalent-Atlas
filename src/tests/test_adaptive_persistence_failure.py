import os, tempfile
from importlib import reload


def test_adaptive_persistence_failure(monkeypatch):
    # Point to unwritable directory (simulate by using a file path within a removed dir)
    bad_dir = tempfile.mkdtemp()
    # Remove directory to force failure
    os.rmdir(bad_dir)
    bad_path = os.path.join(bad_dir, 'adaptive.json')
    monkeypatch.setenv('MOLBRIDGE_ADAPTIVE_CACHE_PATH', bad_path)
    import utils.kdtree_thresholds as kt
    reload(kt)
    base = kt.get_threshold('hydro')
    # Should not raise even if persistence fails internally
    kt.adapt_threshold('hydro', int(base * 0.1), True)
