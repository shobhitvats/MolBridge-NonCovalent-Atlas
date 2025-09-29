"""CLI utilities for golden dataset snapshot & comparison.

Usage:
  python -m src.cli_golden snapshot --pdb 1CRN 4HHB --out golden_baseline
  python -m src.cli_golden compare --baseline golden_baseline --pdb 1CRN 4HHB --tolerance 0.05
"""
from __future__ import annotations
import argparse
import sys
from utils.config import AppConfig
from utils.golden_dataset import build_golden_snapshot, compare_to_baseline

def parse_args(argv=None):
    p = argparse.ArgumentParser("golden-dataset")
    sub = p.add_subparsers(dest='cmd', required=True)

    snap = sub.add_parser('snapshot', help='Generate golden baseline snapshots')
    snap.add_argument('--pdb', nargs='*', help='PDB IDs to process')
    snap.add_argument('--pdb-file', help='Path to file containing newline-delimited PDB IDs')
    snap.add_argument('--out', required=True, help='Output directory for baseline')

    cmp_p = sub.add_parser('compare', help='Compare current run to baseline')
    cmp_p.add_argument('--baseline', required=True, help='Baseline directory (with manifest)')
    cmp_p.add_argument('--pdb', nargs='*', help='PDB IDs to process & compare')
    cmp_p.add_argument('--pdb-file', help='Path to file containing newline-delimited PDB IDs')
    cmp_p.add_argument('--tolerance', type=float, default=0.05, help='Fractional tolerance for interaction count drift')

    return p.parse_args(argv)

def main(argv=None):
    args = parse_args(argv)
    cfg = AppConfig()
    if args.cmd == 'snapshot':
        pdb_ids = args.pdb or []
        if args.pdb_file:
            from pathlib import Path
            file_ids = [l.strip() for l in Path(args.pdb_file).read_text().splitlines() if l.strip()]
            pdb_ids.extend(file_ids)
        if not pdb_ids:
            print('No PDB IDs supplied', file=sys.stderr)
            return 2
        results = build_golden_snapshot(pdb_ids, args.out, cfg)
        print(f"Snapshot complete for {len(results)} PDB IDs -> {args.out}")
    elif args.cmd == 'compare':
        pdb_ids = args.pdb or []
        if args.pdb_file:
            from pathlib import Path
            file_ids = [l.strip() for l in Path(args.pdb_file).read_text().splitlines() if l.strip()]
            pdb_ids.extend(file_ids)
        if not pdb_ids:
            print('No PDB IDs supplied', file=sys.stderr)
            return 2
        out = compare_to_baseline(args.baseline, pdb_ids, cfg, tolerance=args.tolerance)
        if out['regressions']:
            print("Regressions detected:")
            for pid, data in out['regressions'].items():
                print(f"  {pid}: {data}")
        else:
            print("No regressions detected")
    else:
        raise SystemExit(1)

if __name__ == '__main__':  # pragma: no cover
    main(sys.argv[1:])
