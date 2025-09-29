"""Test that metrics emission writes a JSON line when metrics file env is set."""
from __future__ import annotations

import os, json, tempfile
from utils.logging_config import emit_metrics

def test_emit_metrics_file_output():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, 'metrics.log')
        os.environ['MOLBRIDGE_METRICS_FILE'] = path
        # Reconfigure metrics channel (emit_metrics will work regardless)
        emit_metrics({'custom': 'value', 'count': 3})
        # Read back
        with open(path, 'r') as f:
            lines = [l.strip() for l in f if l.strip()]
        assert lines, 'No metrics lines written'
        loaded = json.loads(lines[-1])
        assert loaded.get('custom') == 'value'
        assert loaded.get('count') == 3
