"""Automatic performance flag heuristic application.

Applies recommended defaults unless the user explicitly overrides a flag or selects a manual profile.

Profiles (via MOLBRIDGE_PERF_PROFILE):
  auto    -> Heuristic (recommended)
  full    -> Enable all acceleration features
  minimal -> Disable expensive parallel/accel features, keep safe fast paths
  manual  -> Respect explicit checkboxes / env values only

This module is intentionally lightweight and side-effectful only when it decides to adjust flags.
"""
from __future__ import annotations
import os
from typing import Sequence, Any

def _rust_available() -> bool:
    try:  # pragma: no cover
        from utils import rust_geom  # noqa: F401
        return True
    except Exception:
        return False

def _explicit_overridden(var: str) -> bool:
    return var in os.environ

_MANAGED_VARS = [
    'MOLBRIDGE_ENABLE_VECTOR_GEOM',
    'MOLBRIDGE_USE_PROCESS_POOL',
    'MOLBRIDGE_USE_SHM',
    'MOLBRIDGE_USE_NUMBA',
    'MOLBRIDGE_TASK_GRAPH'
]

def apply_auto_flags(structure, detector_classes: Sequence[Any]):
    profile = os.getenv('MOLBRIDGE_PERF_PROFILE', 'auto').lower()
    # Always ensure vector path ON for all profiles except explicit manual override off
    if profile != 'manual' and not _explicit_overridden('MOLBRIDGE_ENABLE_VECTOR_GEOM'):
        os.environ['MOLBRIDGE_ENABLE_VECTOR_GEOM'] = '1'

    if profile == 'manual':
        return {'profile': profile, 'applied': False}

    # FULL profile: turn on everything aggressively
    if profile == 'full':
        os.environ.setdefault('MOLBRIDGE_USE_PROCESS_POOL', '1')
        os.environ.setdefault('MOLBRIDGE_USE_SHM', '1')
        os.environ.setdefault('MOLBRIDGE_USE_NUMBA', '1')
        os.environ.setdefault('MOLBRIDGE_TASK_GRAPH', '1')
        return {'profile': profile, 'applied': True, 'mode': 'full'}

    # MINIMAL profile: disable heavier layers (except vector geometry which we keep ON)
    if profile == 'minimal':
        os.environ.setdefault('MOLBRIDGE_USE_PROCESS_POOL', '0')
        os.environ.setdefault('MOLBRIDGE_USE_SHM', '0')
        os.environ.setdefault('MOLBRIDGE_USE_NUMBA', '0')
        os.environ.setdefault('MOLBRIDGE_TASK_GRAPH', '1')  # still beneficial
        return {'profile': profile, 'applied': True, 'mode': 'minimal'}

    # AUTO heuristic (default)
    # Skip if user explicitly set any managed var (treat as manual micro overrides for those vars)
    explicit_any = any(_explicit_overridden(v) for v in _MANAGED_VARS if v != 'MOLBRIDGE_ENABLE_VECTOR_GEOM')
    if explicit_any:
        return {'profile': profile, 'applied': False, 'reason': 'user explicit vars present'}

    # Compute workload characteristics
    try:
        n_atoms = sum(1 for _ in structure.get_atoms())
    except Exception:
        n_atoms = 0
    n_det = len(detector_classes) if detector_classes else 0
    est_work = n_atoms * max(1, n_det)

    # Heuristics thresholds (empirically chosen):
    heavy = n_atoms > 12000 and n_det >= 6
    very_heavy = est_work > 90000
    moderate = est_work > 60000

    # Decide process pool + SHM
    if heavy or very_heavy:
        os.environ.setdefault('MOLBRIDGE_USE_PROCESS_POOL', '1')
        os.environ.setdefault('MOLBRIDGE_USE_SHM', '1')
    else:
        os.environ.setdefault('MOLBRIDGE_USE_PROCESS_POOL', '0')
        os.environ.setdefault('MOLBRIDGE_USE_SHM', '0')

    # Numba if moderate & no rust
    if moderate and not _rust_available():
        os.environ.setdefault('MOLBRIDGE_USE_NUMBA', '1')
    else:
        os.environ.setdefault('MOLBRIDGE_USE_NUMBA', '0')

    os.environ.setdefault('MOLBRIDGE_TASK_GRAPH', '1')

    return {
        'profile': profile,
        'applied': True,
        'mode': 'auto',
        'n_atoms': n_atoms,
        'n_detectors': n_det,
        'est_work': est_work,
        'process_pool': os.environ.get('MOLBRIDGE_USE_PROCESS_POOL'),
        'shm': os.environ.get('MOLBRIDGE_USE_SHM'),
        'numba': os.environ.get('MOLBRIDGE_USE_NUMBA')
    }
