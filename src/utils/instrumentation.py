"""Shared instrumentation helpers to reduce duplication across detectors.

Provides standardized creation and updating of funnel metrics so each detector
can focus on domain–specific counting logic without rewriting boilerplate.

All helpers operate on plain dicts to keep them serializable and JSON‑friendly.
"""
from __future__ import annotations
from typing import Dict, Optional

BASE_KEYS = (
    'raw_pairs', 'candidate_pairs', 'accepted_pairs', 'acceptance_ratio',
    'core_pair_generation', 'phase_pair_gen_ms', 'phase_eval_ms', 'phase_build_ms'
)

def init_funnel(raw_pairs: int = 0,
                candidate_pairs: int = 0,
                accepted_pairs: int = 0,
                core_pair_generation: bool = False,
                extra: Optional[Dict[str, object]] = None) -> Dict[str, object]:
    """Initialize a funnel instrumentation dictionary with standard keys.

    extra: optional dict merged in after base keys (does not override unless intentional).
    """
    instr: Dict[str, object] = {
        'raw_pairs': int(raw_pairs),
        'candidate_pairs': int(candidate_pairs),
        'accepted_pairs': int(accepted_pairs),
        'acceptance_ratio': (float(accepted_pairs) / candidate_pairs) if candidate_pairs else 0.0,
        'core_pair_generation': bool(core_pair_generation),
        'phase_pair_gen_ms': None,  # filled later
        'phase_eval_ms': None,
        'phase_build_ms': None,
    }
    if extra:
        instr.update(extra)
    return instr

def update_counts(instr: Dict[str, object], *, raw: Optional[int] = None,
                  candidate: Optional[int] = None, accepted: Optional[int] = None) -> None:
    """Update count fields and recompute acceptance ratio in place."""
    if raw is not None:
        instr['raw_pairs'] = int(raw)
    if candidate is not None:
        instr['candidate_pairs'] = int(candidate)
    if accepted is not None:
        instr['accepted_pairs'] = int(accepted)
    cand = int(instr.get('candidate_pairs') or 0)
    acc = int(instr.get('accepted_pairs') or 0)
    instr['acceptance_ratio'] = (acc / cand) if cand else 0.0


def finalize_funnel(instr: Dict[str, object], *,
                    pair_gen_seconds: float,
                    eval_seconds: float,
                    build_seconds: float = 0.0) -> None:
    """Finalize timing phases (store in ms, rounded) leaving existing counts untouched."""
    instr['phase_pair_gen_ms'] = round(pair_gen_seconds * 1000.0, 3)
    instr['phase_eval_ms'] = round(eval_seconds * 1000.0, 3)
    instr['phase_build_ms'] = round(build_seconds * 1000.0, 3)

__all__ = ['init_funnel', 'update_counts', 'finalize_funnel', 'BASE_KEYS']
