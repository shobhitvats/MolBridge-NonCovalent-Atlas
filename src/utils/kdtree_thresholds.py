"""Central KD-tree threshold configuration helper.

Provides a single place to parse environment variables controlling when
KD-tree pruning is enabled for different detector types.

Variables:
  MOLBRIDGE_KDTREE_IONIC_THRESHOLD
  MOLBRIDGE_KDTREE_HYDRO_THRESHOLD
  MOLBRIDGE_KDTREE_HBOND_THRESHOLD

Defaults chosen based on preliminary benchmarks; can be tuned without code changes.
"""
from __future__ import annotations
import os
from functools import lru_cache
from typing import Dict, Tuple
import json
import tempfile
from pathlib import Path

_DEFAULTS = {
    'ionic': 1200,
    'hydro': 600,
    'hbond': 5000,
    'pipi': 800,
    'chalcogen': 1200,
    'tetrel': 1500,
}

_ENV_MAP = {
    'ionic': 'MOLBRIDGE_KDTREE_IONIC_THRESHOLD',
    'hydro': 'MOLBRIDGE_KDTREE_HYDRO_THRESHOLD',
    'hbond': 'MOLBRIDGE_KDTREE_HBOND_THRESHOLD',
    'pipi': 'MOLBRIDGE_KDTREE_PIPI_THRESHOLD',
    'chalcogen': 'MOLBRIDGE_KDTREE_CHALCOGEN_THRESHOLD',
    'tetrel': 'MOLBRIDGE_KDTREE_TETREL_THRESHOLD'
}

# Adaptive runtime overrides (in-memory; not persisted across processes)
_ADAPTIVE: Dict[str, int] = {}
# EMA trackers for smoothing (value & utilization density) and change cooldown counters
_EMA_THRESHOLD: Dict[str, float] = {}
_COOLDOWN: Dict[str, int] = {}  # runs until another change allowed
_EMA_ALPHA = float(os.getenv('MOLBRIDGE_KDTREE_EMA_ALPHA', '0.3'))
_COOLDOWN_RUNS = int(os.getenv('MOLBRIDGE_KDTREE_COOLDOWN', '2'))
_DENSITY_TARGETS: Dict[str, float] = {
    # target candidate pairs per atom (donor*acceptor scale) for stable runtime
    'hbond': 0.9,
    'ionic': 1.2,
    'hydro': 0.8,
    'pipi': 0.7,
    'chalcogen': 0.9,
    'tetrel': 1.0
}
_LAST_DENSITY: Dict[str, float] = {}
_REASONS: Dict[str, str] = {}
_ADAPTIVE_BOUNDS = {
    'ionic': (200, 20000),
    'hydro': (50, 5000),
    'hbond': (500, 20000),
    'pipi': (100, 10000),
    'chalcogen': (200, 15000),
    'tetrel': (300, 18000)
}
_WARM_START_DONE = False  # guard to avoid repeated warm-start scaling

# Optional persistence (opt-in via ENV path); loads on import if present
_CACHE_ENV = 'MOLBRIDGE_ADAPTIVE_CACHE_PATH'
_CACHE_PATH = os.getenv(_CACHE_ENV)
if _CACHE_PATH:
    try:  # pragma: no cover - load best effort
        p = Path(_CACHE_PATH)
        if p.is_file():
            data = json.loads(p.read_text())
            if isinstance(data, dict):
                for k, v in data.items():
                    if isinstance(v, int):
                        _ADAPTIVE[k] = v
    except Exception:
        pass

def _persist_adaptive_state():  # pragma: no cover (I/O side effect)
    if not _CACHE_PATH:
        return
    try:
        p = Path(_CACHE_PATH)
        p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_suffix('.tmp')
        tmp.write_text(json.dumps(_ADAPTIVE, indent=2))
        tmp.replace(p)
    except Exception:
        pass


@lru_cache(maxsize=None)
def get_threshold(detector_key: str) -> int:
    key = detector_key.lower()
    env_var = _ENV_MAP.get(key)
    # Environment variable (when present) has highest precedence so tests that
    # toggle thresholds at runtime see immediate effect regardless of prior adaptation.
    if env_var and os.getenv(env_var) is not None:
        try:
            return int(os.getenv(env_var))
        except Exception:
            pass
    if key in _ADAPTIVE:
        return _ADAPTIVE[key]
    if env_var:
        try:
            return int(os.getenv(env_var, str(_DEFAULTS.get(key, 10**9))))
        except Exception:
            return _DEFAULTS.get(key, 10**9)
    return _DEFAULTS.get(key, 10**9)


def adapt_threshold(detector_key: str, candidate_complexity: int, kdtree_used: bool) -> Tuple[int, bool, str]:
    """Adaptive tuner adjusting threshold based on observed usage.

    Returns (new_threshold, changed, reason_code)
    reason_code examples:
      'raise_low_util', 'lower_high_miss', 'density_adjust_up', 'density_adjust_down', 'unchanged'
    """
    key = detector_key.lower()
    current = get_threshold(key)  # will read adaptive if already set
    low, high = _ADAPTIVE_BOUNDS.get(key, (100, 10**7))
    changed = False
    reason = 'unchanged'
    # Cooldown enforcement
    cd = _COOLDOWN.get(key, 0)
    if cd > 0:
        _COOLDOWN[key] = cd - 1
        # Still update density tracking but skip structural changes
        pass_change = False
    else:
        pass_change = True
    # Heuristics (conservative adjustments)
    if pass_change and kdtree_used and candidate_complexity < current * 0.4:
        # using tree when maybe not needed -> raise threshold moderately
        new_val = min(high, int(current * 1.25) + 1)
        if new_val != current:
            _ADAPTIVE[key] = new_val
            changed = True
            reason = 'raise_low_util'
    elif pass_change and (not kdtree_used) and candidate_complexity > current * 1.6:
        # missed opportunity -> lower threshold
        new_val = max(low, int(candidate_complexity * 0.8))
        if new_val != current:
            _ADAPTIVE[key] = new_val
            changed = True
            reason = 'lower_high_miss'

    # Candidate density adjustment (only if not just changed above)
    try:
        if key in _DENSITY_TARGETS and candidate_complexity > 0:
            target = _DENSITY_TARGETS[key]
            # Approximate atom count scale: sqrt(candidate_complexity) for donor*acceptor pairing domain
            import math
            approx_atoms = max(1, int(math.sqrt(candidate_complexity)))
            density = candidate_complexity / approx_atoms
            _LAST_DENSITY[key] = density
            # Adjust within +/- 15% band; skip if previous change to avoid oscillation
            if not changed and pass_change:
                if density > target * 1.4 and current > low:
                    new_val = max(low, int(current * 0.9))
                    if new_val != current:
                        _ADAPTIVE[key] = new_val
                        changed = True
                        reason = 'density_adjust_down'
                elif density < target * 0.6 and current < high:
                    new_val = min(high, int(current * 1.1))
                    if new_val != current:
                        _ADAPTIVE[key] = new_val
                        changed = True
                        reason = 'density_adjust_up'
    except Exception:
        pass
    # EMA smoothing: update smoothed threshold view (does not override real threshold immediately)
    try:
        ema_prev = _EMA_THRESHOLD.get(key, float(current))
        ema_new = (1 - _EMA_ALPHA) * ema_prev + _EMA_ALPHA * float(_ADAPTIVE.get(key, current))
        _EMA_THRESHOLD[key] = ema_new
        # If difference between ema_new and stored adaptive > 15% and cooldown expired, snap adaptive to ema
        if pass_change and not changed:
            adaptive_val = _ADAPTIVE.get(key, current)
            if abs(ema_new - adaptive_val)/max(1.0, adaptive_val) > 0.15:
                snap = int(min(high, max(low, ema_new)))
                if snap != adaptive_val:
                    _ADAPTIVE[key] = snap
                    changed = True
                    reason = 'ema_snap'
    except Exception:
        pass
    if changed:
        _COOLDOWN[key] = _COOLDOWN_RUNS
    _REASONS[key] = reason
    if changed:
        try:
            _persist_adaptive_state()
        except Exception:
            pass
    return _ADAPTIVE.get(key, current), changed, _REASONS.get(key, reason)

def warm_start_thresholds(atom_count: int):
    """Heuristically warm-start adaptive thresholds once per process.

    Rationale: The initial defaults can be far from ideal for very small or
    very large structures leading to unnecessary early KD-tree usage (or
    missed opportunities) before enough adaptation cycles accumulate.

    Strategy: scale each default by a factor derived from atom_count:
        scale = clamp( 0.6 + 0.4 * (atom_count / 3000), 0.6, 1.6 )
    This gently raises thresholds for large structures (delaying tree usage
    where brute force vectorized paths are still competitive) and lowers for
    tiny structures so KD-tree isn't considered at all.

    Only applies to detector keys that do not already have an adaptive value
    and only runs once per process lifetime.
    """
    global _WARM_START_DONE
    if _WARM_START_DONE:
        return
    try:
        scale = 0.6 + 0.4 * min(1.0, max(0.0, atom_count / 3000.0))
        scale = max(0.6, min(1.6, scale))
        for k, base in _DEFAULTS.items():
            # Respect explicit environment overrides: skip warm-start for those
            env_var = _ENV_MAP.get(k)
            if env_var and os.getenv(env_var) is not None:
                continue
            if k not in _ADAPTIVE:  # do not override previously persisted/adapted values
                low, high = _ADAPTIVE_BOUNDS.get(k, (int(base*0.25), int(base*4)))
                val = int(base * scale)
                val = min(high, max(low, val))
                _ADAPTIVE[k] = val
        _WARM_START_DONE = True
    except Exception:
        _WARM_START_DONE = True  # prevent retry storms


def get_last_density(detector_key: str) -> float:
    return _LAST_DENSITY.get(detector_key.lower(), 0.0)

def get_last_reason(detector_key: str) -> str:
    return _REASONS.get(detector_key.lower(), 'unknown')


__all__ = ['get_threshold', 'adapt_threshold', 'get_last_density', 'get_last_reason']

def should_flag_kdtree(candidate_pairs: int, raw_pairs: int, detector_key: str) -> bool:
    """Heuristic to decide if KD-tree (or analogous pruning) was effectively used.

    Uses ratio of candidate/raw vs density target to infer pruning benefit.
    Returns True if pruning appears beneficial (i.e., candidate density below target band) so that
    adapt_threshold can consider raising threshold (less need) or track utilization.
    """
    if raw_pairs <= 0:
        return False
    ratio = candidate_pairs / raw_pairs
    target = _DENSITY_TARGETS.get(detector_key.lower())
    if target is None:
        # Fallback generic: treat below 0.6 as pruned
        return ratio < 0.6
    # If ratio is significantly below target, consider tree utilized
    return ratio < (target * 0.75)

__all__.append('should_flag_kdtree')