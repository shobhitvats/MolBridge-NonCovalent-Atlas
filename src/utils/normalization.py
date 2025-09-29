"""Normalization helpers for detector outputs (non-breaking layer).

These functions are **not** yet wired into processors; they can be adopted incrementally.
"""
from __future__ import annotations
from typing import Any, Dict, List, Mapping

try:  # optional import
    from analysis.base import to_interaction_dict  # type: ignore
except Exception:  # pragma: no cover
    def to_interaction_dict(obj: Any, fallback_kind: str):  # type: ignore
        return obj if isinstance(obj, dict) else {"type": fallback_kind}


def normalize_detector_output(kind: str, raw: Any, detector: Any | None = None) -> List[Dict[str, Any]]:
    """Return list[dict] with at least a 'type' key.

    Accepted raw forms:
      * list[dict]
      * list[dataclass/objects]
      * single dict/object
    """
    if raw is None:
        return []
    if isinstance(raw, list):
        out: List[Dict[str, Any]] = []
        for item in raw:
            if isinstance(item, dict):
                if 'type' not in item:
                    item['type'] = kind
                out.append(item)
            else:
                out.append(to_interaction_dict(item, kind))
        return out
    if isinstance(raw, dict):
        return [raw if 'type' in raw else {**raw, 'type': kind}]
    return [to_interaction_dict(raw, kind)]


def normalize_interaction_map(raw_map: Mapping[str, Any]) -> Dict[str, List[Dict[str, Any]]]:
    return {k: normalize_detector_output(k, v) for k, v in raw_map.items()}
