"""CLI for metrics aggregation utilities.

Usage:
  python -m cli_metrics summarize --file metrics.log --out summary.json
  python -m cli_metrics csv --file metrics.log --out metrics.csv

Accepts JSON lines produced by emit_metrics (and per-detector structured records).
"""
from __future__ import annotations
import argparse, json, statistics, sys
from pathlib import Path
from typing import List, Dict, Any


def load_lines(path: Path) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
                if isinstance(rec, dict):
                    records.append(rec)
            except Exception:
                continue
    return records


def _p95(values: List[float]) -> float:
    if not values:
        return 0.0
    if len(values) == 1:
        return values[0]
    vs = sorted(values)
    idx = int(0.95 * (len(vs)-1))
    return vs[idx]


def summarize(records: List[Dict[str, Any]]) -> Dict[str, Any]:
    batch = [r for r in records if r.get('event') == 'batch_metrics']
    det_complete = [r for r in records if r.get('event') == 'detector_complete']
    out: Dict[str, Any] = {
        'total_records': len(records),
        'batches': len(batch),
        'detector_completions': len(det_complete),
    }
    if batch:
        durations = [b.get('duration_sec', 0) for b in batch]
        out['batch_duration_mean'] = statistics.mean(durations)
        out['batch_duration_p95'] = statistics.quantiles(durations, n=20)[-1] if len(durations) > 1 else durations[0]
    # Aggregate detector timing if present in batches
    det_times: Dict[str, List[float]] = {}
    for b in batch:
        for k, v in (b.get('detector_timings') or {}).items():
            det_times.setdefault(k, []).append(v)
    out['detector_timing_mean'] = {k: statistics.mean(vs) for k, vs in det_times.items() if vs}
    out['detector_timing_p95'] = {k: _p95(vs) for k, vs in det_times.items() if vs and len(vs) > 1}
    # Adaptive threshold traces (if instrumentation embedded in batch metrics)
    adaptive_info: Dict[str, Dict[str, Any]] = {}
    for b in batch:
        instr = b.get('instrumentation') or {}
        for det, meta in instr.items():
            if not isinstance(meta, dict):
                continue
            # Expect keys like adapted_threshold, original_threshold
            rec = adaptive_info.setdefault(det, {'adapted': [], 'original': []})
            if 'adapted_threshold' in meta:
                rec['adapted'].append(meta.get('adapted_threshold'))
            if 'original_threshold' in meta:
                rec['original'].append(meta.get('original_threshold'))
    if adaptive_info:
        out['adaptive_thresholds'] = {
            det: {
                'original_mean': statistics.mean(v['original']) if v['original'] else None,
                'adapted_mean': statistics.mean(v['adapted']) if v['adapted'] else None,
                'adapted_p95': _p95([x for x in v['adapted'] if isinstance(x,(int,float))]) if v['adapted'] else None
            } for det, v in adaptive_info.items()
        }
    return out


def write_csv(records: List[Dict[str, Any]], out_path: Path) -> None:
    import csv
    # Flatten known fields; dynamic keys dropped
    fieldnames = ['event','strategy','detector','count','duration_sec','queued_delay_sec','timestamp']
    with out_path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in records:
            row = {k: r.get(k) for k in fieldnames}
            w.writerow(row)


def main():  # pragma: no cover - CLI integration
    ap = argparse.ArgumentParser(description='Metrics log utilities')
    sub = ap.add_subparsers(dest='cmd', required=True)
    s1 = sub.add_parser('summarize')
    s1.add_argument('--file', required=True)
    s1.add_argument('--out', required=False)
    s2 = sub.add_parser('export-csv')
    s2.add_argument('--file', required=True)
    s2.add_argument('--out', required=True)
    args = ap.parse_args()
    path = Path(args.file)
    if not path.exists():
        print('Input file not found', file=sys.stderr)
        sys.exit(2)
    records = load_lines(path)
    if args.cmd == 'summarize':
        summary = summarize(records)
        if args.out:
            Path(args.out).write_text(json.dumps(summary, indent=2))
        else:
            print(json.dumps(summary, indent=2))
    else:  # export-csv
        write_csv(records, Path(args.out))

if __name__ == '__main__':  # pragma: no cover
    main()
