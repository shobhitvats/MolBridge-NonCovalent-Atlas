import json, tempfile, subprocess, sys, os


def test_bench_regression_gate_exit_code():
    # Create synthetic baseline and current results with higher p95 to trigger failure
    baseline = {"aggregate": {"hydrogenbond": {"p95": 1.0, "mean": 0.9, "n": 2, "mean_acceptance_ratio": 0.1, "warnings": []}}}
    current = {"aggregate": {"hydrogenbond": {"p95": 2.0, "mean": 1.8, "n": 2, "mean_acceptance_ratio": 0.1, "warnings": []}}}
    b_file = tempfile.NamedTemporaryFile(delete=False, suffix='.json'); b_file.write(json.dumps(baseline).encode()); b_file.close()
    c_file = tempfile.NamedTemporaryFile(delete=False, suffix='.json'); c_file.write(json.dumps(current).encode()); c_file.close()
    cmd = [sys.executable, '-m', 'performance.bench_regression_gate', c_file.name, '--baseline', b_file.name, '--max-p95-growth', '1.5']
    env = dict(os.environ)
    env['PYTHONPATH'] = 'src'
    proc = subprocess.run(cmd, env=env)
    assert proc.returncode == 1, "Expected non-zero exit due to p95 regression"