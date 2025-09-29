"""MIT License

Copyright (c) 2025 MolBridge Contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

CI regression gate for microbenchmark results.

Usage (after running microbench):
  python -m performance.microbench > bench.json
  python -m performance.bench_regression_gate bench.json --max-p95-growth 1.25 --baseline baseline.json

Baseline file optional; if provided, compares detector p95 values and fails (exit 1) if growth exceeds factor.
Acceptance ratio sanity warnings can escalate to failure if --strict-accept passed.
"""
from __future__ import annotations
import json, sys, argparse, math
from pathlib import Path

def load_json(path: str):
    try:
        return json.loads(Path(path).read_text())
    except Exception as e:  # pragma: no cover
        print(f"ERROR loading {path}: {e}", file=sys.stderr)
        return {}

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('results', help='Path to microbench JSON output')
    ap.add_argument('--baseline', help='Baseline JSON (prior run)')
    ap.add_argument('--max-p95-growth', type=float, default=1.5, help='Max allowed p95 slowdown factor vs baseline')
    ap.add_argument('--min-accept-ratio', type=float, default=0.0, help='Minimum mean acceptance ratio (warn if below)')
    ap.add_argument('--strict-accept', action='store_true', help='Fail build on acceptance ratio warnings')
    args = ap.parse_args(argv)

    data = load_json(args.results)
    baseline = load_json(args.baseline) if args.baseline else {}
    agg = data.get('aggregate', {})
    base_agg = baseline.get('aggregate', {}) if baseline else {}
    failures = []
    warnings = []
    for det, stats in agg.items():
        p95 = stats.get('p95')
        if p95 is None:
            continue
        b = base_agg.get(det, {})
        b_p95 = b.get('p95')
        if b_p95 and b_p95 > 0:
            factor = p95 / b_p95
            if factor > args.max_p95_growth:
                failures.append(f"{det}: p95 {p95:.4f}s exceeds {args.max_p95_growth:.2f}x baseline ({b_p95:.4f}s)")
        mean_accept = stats.get('mean_acceptance_ratio', 0.0)
        if mean_accept < args.min_accept_ratio:
            warnings.append(f"{det}: mean acceptance ratio {mean_accept:.6f} < {args.min_accept_ratio}")
        if stats.get('warnings'):
            warnings.append(f"{det}: instrumentation warnings -> {stats['warnings']}")
    if warnings and args.strict_accept:
        failures.extend([f"(strict) {w}" for w in warnings])
    if failures:
        print("BENCH REGRESSION FAILURES:")
        for f in failures:
            print(" -", f)
        if warnings and not args.strict_accept:
            print("WARNINGS:")
            for w in warnings:
                print(" -", w)
        sys.exit(1)
    if warnings:
        print("WARNINGS (non-fatal):")
        for w in warnings:
            print(" -", w)
    print("Benchmark regression gate passed.")
    return 0

if __name__ == '__main__':  # pragma: no cover
    main()
