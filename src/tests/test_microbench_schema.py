import json, tempfile, os
from importlib import reload
import subprocess, sys, pathlib


def test_microbench_schema(monkeypatch):
    monkeypatch.setenv('PYTHONPATH', 'src')
    out_file = tempfile.NamedTemporaryFile(delete=False, suffix='.json')
    out_file.close()
    cmd = [sys.executable, '-m', 'performance.microbench']
    # Run with PYTHONPATH injection
    env = dict(os.environ)
    env['PYTHONPATH'] = 'src'
    subprocess.run(cmd, stdout=open(out_file.name, 'w'), stderr=subprocess.DEVNULL, env=env, check=False)
    data = json.loads(pathlib.Path(out_file.name).read_text())
    assert 'runs' in data and isinstance(data['runs'], list)
    assert 'aggregate' in data and isinstance(data['aggregate'], dict)
