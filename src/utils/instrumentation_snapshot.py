"""Aggregation helpers for detector instrumentation.

Provides utilities to collect instrumentation dicts from instantiated
(detector, interactions) results and produce a summarized snapshot
suitable for UI panels, logging, or JSON export.

The snapshot schema:
{
  'detectors': {
       <key>: { <instrumentation fields> }
  },
  'summary': {
       'total_raw_pairs': int,
       'total_candidate_pairs': int,
       'total_accepted_pairs': int,
       'mean_acceptance_ratio': float,
       'detector_count': int,
       'timestamp': iso8601
  }
}
"""
from __future__ import annotations
from typing import Dict, Any, Iterable, Mapping
from datetime import datetime
from utils.detector_meta import list_detector_keys, get_meta

# Keys we attempt to aggregate if present
AGG_KEYS = ('raw_pairs', 'candidate_pairs', 'accepted_pairs', 'acceptance_ratio')


def collect_snapshot(detectors: Iterable[Any]) -> Dict[str, Any]:
    """Collect instrumentation from detector objects.

    detectors: iterable of detector instances (each may have `instrumentation`).
    Missing instrumentation is skipped.
    """
    out: Dict[str, Any] = {'detectors': {}, 'summary': {}}
    total_raw = total_cand = total_acc = 0
    ratios = []
    for det in detectors:
        key = getattr(det, 'registry_key', None) or getattr(det, 'key', None)
        # Fallback: try class name lower
        if key is None:
            key = det.__class__.__name__.replace('Detector', '').lower()
        instr = getattr(det, 'instrumentation', None)
        if not isinstance(instr, dict):
            continue
        out['detectors'][key] = instr
        total_raw += int(instr.get('raw_pairs') or 0)
        total_cand += int(instr.get('candidate_pairs') or 0)
        total_acc += int(instr.get('accepted_pairs') or 0)
        ratio = instr.get('acceptance_ratio')
        if isinstance(ratio, (int, float)):
            ratios.append(float(ratio))
    det_count = len(out['detectors'])
    out['summary'] = {
        'total_raw_pairs': total_raw,
        'total_candidate_pairs': total_cand,
        'total_accepted_pairs': total_acc,
        'mean_acceptance_ratio': (sum(ratios)/len(ratios)) if ratios else 0.0,
        'detector_count': det_count,
        'timestamp': datetime.utcnow().isoformat() + 'Z'
    }
    return out

__all__ = ['collect_snapshot']
