#!/usr/bin/env python3
"""Utility to refresh the committed microbenchmark baseline.

Runs the microbench harness (synthetic structures by default) and overwrites
`microbench_baseline.json` with merged timing & funnel aggregates.
Intended to be invoked manually after intentional performance optimizations.

Usage:
  python scripts/update_microbench_baseline.py --repeat 2 --seed 42
"""
from __future__ import annotations
import json, argparse, os, sys, time
from pathlib import Path


def run_microbench(seed: int | None, repeat: int) -> dict:
    env = dict(os.environ)
    env['PYTHONPATH'] = 'src'
    if seed is not None:
        env['MOLBRIDGE_BENCH_SEED'] = str(seed)
    env['MOLBRIDGE_BENCH_REPEAT'] = str(repeat)
    import subprocess, tempfile
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.json')
    tmp.close()
    cmd = [sys.executable, '-m', 'performance.microbench']
    subprocess.run(cmd, env=env, stdout=open(tmp.name, 'w'), check=False)
    try:
        data = json.loads(Path(tmp.name).read_text())
    except Exception:
        data = {}
    Path(tmp.name).unlink(missing_ok=True)
    return data


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=42, help='Deterministic synthetic structure seed')
    ap.add_argument('--repeat', type=int, default=2, help='Repeats per size tier')
    ap.add_argument('--output', default='microbench_baseline.json')
    args = ap.parse_args(argv)
    data = run_microbench(args.seed, args.repeat)
    meta = {
        'refreshed_at': time.time(),
        'seed': args.seed,
        'repeat': args.repeat
    }
    data['baseline_meta'] = meta
    Path(args.output).write_text(json.dumps(data, indent=2))
    print(f"Updated baseline written to {args.output}")

if __name__ == '__main__':  # pragma: no cover
    raise SystemExit(main())
