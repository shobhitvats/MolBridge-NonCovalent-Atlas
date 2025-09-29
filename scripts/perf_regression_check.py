#!/usr/bin/env python3
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

Performance regression gate script.

Runs a comparison against an existing golden baseline directory.
Exit codes:
 0 - success, no regressions beyond tolerance
 1 - regressions detected
 2 - usage / configuration error

Example:
  python scripts/perf_regression_check.py --baseline golden_baseline --pdb 1CRN 4HHB --tolerance 0.05
"""
from __future__ import annotations
import argparse
import sys
from utils.config import AppConfig
from utils.golden_dataset import compare_to_baseline


def parse_args(argv=None):
    p = argparse.ArgumentParser("perf-regression-check")
    p.add_argument('--baseline', required=True, help='Golden baseline directory path')
    p.add_argument('--pdb', nargs='+', required=True, help='Subset of PDB IDs (must exist in baseline)')
    p.add_argument('--tolerance', type=float, default=0.05, help='Fractional tolerance for interaction count drift')
    p.add_argument('--time-tolerance', type=float, default=0.30, help='Fractional tolerance for slowdown in processing time')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    cfg = AppConfig()
    try:
        out = compare_to_baseline(
            args.baseline,
            args.pdb,
            cfg,
            tolerance=args.tolerance,
            time_tolerance=args.time_tolerance,
        )
    except AssertionError as e:
        print(f"Configuration error: {e}")
        return 2
    regressions = out.get('regressions', {})
    timing_regressions = out.get('timing_regressions', {})
    had_count = any(v for v in regressions.values())
    had_time = bool(timing_regressions)
    if had_count or had_time:
        print("Regressions detected:")
        if had_count:
            print(" Interaction count drift:")
            for pid, data in regressions.items():
                print(f"  {pid}: {data}")
        if had_time:
            print(" Timing slowdowns:")
            for pid, data in timing_regressions.items():
                print(f"  {pid}: {data}")
        return 1
    print("No regressions beyond tolerances (counts/timing)")
    return 0

if __name__ == '__main__':  # pragma: no cover
    raise SystemExit(main(sys.argv[1:]))
