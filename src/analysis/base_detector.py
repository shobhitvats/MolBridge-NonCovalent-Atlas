"""Common detector interfaces and registration decorator.

This module introduces a lightweight protocol / base class pattern plus a
decorator-based registration pathway so that new interaction detectors can be
added with a single annotation instead of editing multiple mapping dicts.

The existing eager registry in ``analysis.registry`` remains the source of
truth; the decorator simply appends / overrides entries there to avoid
drift. Backwards compatibility: all previously imported detector classes
continue to work as they are already present in ``DETECTOR_REGISTRY``.
"""
from __future__ import annotations

from typing import Protocol, Iterable, Any, Callable, runtime_checkable, List, Dict

try:  # local import – safe circular avoidance (function level if needed)
    from .registry import DETECTOR_REGISTRY
except Exception:  # pragma: no cover - during early import bootstrap
    DETECTOR_REGISTRY = {}


@runtime_checkable
class DetectorProtocol(Protocol):  # pragma: no cover - structural typing
    """Minimum surface all detectors should expose going forward.

    A detector implementation MAY still provide a specialised method name
    (e.g. ``detect_hydrogen_bonds``). The batch / parallel processors resolve
    that via the registry metadata. For new detectors, prefer implementing a
    single ``detect(structure)`` plus ``to_dict_list`` to normalize output.
    """

    def detect(self, structure: Any) -> Any:  # broad typing – structure is a Bio.PDB structure
        ...

    def to_dict_list(self, interactions: Any) -> List[Dict[str, Any]]:  # noqa: D401
        ...


def register_detector(key: str, method: str | None = None) -> Callable[[type], type]:
    """Decorator to register a detector class.

    Args:
        key: canonical interaction key (lowercase, snake / compact style)
        method: optional explicit method name; if omitted 'detect' assumed.
    """

    def _decorator(cls: type) -> type:
        mname = method or 'detect'
        # Late import to avoid circular import at module load
        try:
            from . import registry as _reg  # type: ignore
            _reg.DETECTOR_REGISTRY[key] = (cls, mname)
        except Exception:
            # Fallback to local dict if registry import not yet resolved
            DETECTOR_REGISTRY[key] = (cls, mname)
        return cls

    return _decorator


def resolve_detector_key(detector_cls: type) -> tuple[str | None, str | None]:
    """Resolve registry key + method for an already known detector class.

    Returns (key, method) or (None, None) if not registered.
    """
    for k, (cls, method) in DETECTOR_REGISTRY.items():
        if cls is detector_cls:
            return k, method
    return None, None
