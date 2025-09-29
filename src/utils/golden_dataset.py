"""Golden dataset management utilities.

Provides helpers to:
  * Collect a curated set of PDB IDs (or local files) and run all detectors
  * Persist canonical interaction snapshots + summary metrics to JSON
  * Compare a fresh run against stored baseline (regression detection)

Usage (script style):
    from utils.golden_dataset import build_golden_snapshot, compare_to_baseline
    build_golden_snapshot(['1CRN','4HHB'], out_dir='golden_baseline')
    regressions = compare_to_baseline('golden_baseline', ['1CRN','4HHB'])
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Any
import json
import time

from utils.config import AppConfig
from analysis.high_performance_batch import HighPerformanceBatchProcessor

@dataclass
class GoldenResult:
    pdb_id: str
    interactions: Dict[str, Any]
    summary: Dict[str, Any]
    processing_time: float


def build_golden_snapshot(pdb_ids: List[str], out_dir: str, config: AppConfig | None = None) -> List[GoldenResult]:
    cfg = config or AppConfig()
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    proc = HighPerformanceBatchProcessor(cfg, use_parallel=True)
    results: List[GoldenResult] = []
    for pid in pdb_ids:
        start = time.time()
        res = proc.process_single_protein(pid)
        processing_time = time.time() - start
        interactions = res.get('interactions', {})
        summary = res.get('summary', {})
        golden = GoldenResult(pid, interactions, summary, processing_time)
        # Write per-PDB snapshot
        with (out_path / f"{pid}.json").open('w') as f:
            json.dump({
                'pdb_id': pid,
                'interactions': interactions,
                'summary': summary,
                'processing_time': processing_time
            }, f, indent=2)
        results.append(golden)
    # Write manifest
    manifest = {
        'pdb_ids': pdb_ids,
        'count': len(pdb_ids),
        'generated_at': time.time()
    }
    with (out_path / 'manifest.json').open('w') as f:
        json.dump(manifest, f, indent=2)
    return results


def compare_to_baseline(baseline_dir: str, pdb_ids: List[str], config: AppConfig | None = None, tolerance: float = 0.05, time_tolerance: float = 0.30) -> Dict[str, Any]:
    """Compare current run to baseline.

    tolerance: fractional allowed deviation in interaction counts per detector.
    time_tolerance: fractional allowed slowdown in total processing time per PDB.
    """
    cfg = config or AppConfig()
    base = Path(baseline_dir)
    assert base.exists(), f"Baseline dir {baseline_dir} does not exist"
    processor = HighPerformanceBatchProcessor(cfg, use_parallel=True)
    regressions = {}
    timing_regressions = {}
    for pid in pdb_ids:
        baseline_file = base / f"{pid}.json"
        if not baseline_file.exists():
            regressions[pid] = {'error': 'baseline missing'}
            continue
        with baseline_file.open() as f:
            baseline = json.load(f)
        current = processor.process_single_protein(pid)
        base_counts = baseline.get('summary', {}).get('interaction_counts', {})
        cur_counts = current.get('summary', {}).get('interaction_counts', {})
        drift = {}
        for k, v in base_counts.items():
            cur = cur_counts.get(k, 0)
            if v == 0 and cur == 0:
                continue
            frac = abs(cur - v) / max(1, v)
            if frac > tolerance:
                drift[k] = {'baseline': v, 'current': cur, 'fractional_change': round(frac, 3)}
        if drift:
            regressions[pid] = {'drift': drift}
        # Timing comparison
        b_time = baseline.get('processing_time') or baseline.get('summary', {}).get('processing_time')
        c_time = current.get('processing_time') or current.get('summary', {}).get('processing_time')
        if b_time and c_time:
            if c_time > b_time * (1 + time_tolerance):
                timing_regressions[pid] = {'baseline_sec': round(b_time,3), 'current_sec': round(c_time,3), 'slowdown_factor': round(c_time / b_time, 3)}
    return {'regressions': regressions, 'timing_regressions': timing_regressions, 'tolerance': tolerance, 'time_tolerance': time_tolerance}
