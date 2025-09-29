"""Unified serialization helpers for interaction results.

Non-breaking: Accepts legacy dict lists, dataclass objects, or new Interaction objects.
Provides export utilities to JSON/CSV/Markdown.
"""
from __future__ import annotations
from typing import Any, Iterable, List, Dict, Tuple
import json
import os
import csv
import io
from utils.settings import get_settings
import threading

# Optional orjson fast path
try:  # pragma: no cover (speed path)
    import orjson  # type: ignore
    _HAVE_ORJSON = True
except Exception:  # pragma: no cover
    orjson = None  # type: ignore
    _HAVE_ORJSON = False

try:
    from analysis.base import Interaction, to_interaction_dict  # type: ignore
except Exception:  # pragma: no cover
    Interaction = None  # type: ignore
    def to_interaction_dict(obj: Any, fallback_kind: str):  # type: ignore
        return obj if isinstance(obj, dict) else {"type": fallback_kind}

NORMALIZED_KEYS = [
    "type","distance","score","residue1","residue2","chain1","chain2","angle","theta_angle","delta_angle"
]

def _normalize_item(item: Any, fallback_kind: str) -> Dict[str, Any]:
    if isinstance(item, dict):
        if 'type' not in item:
            item['type'] = fallback_kind
        return item
    return to_interaction_dict(item, fallback_kind)

def normalize_interaction_collection(raw: Any, fallback_kind: str) -> List[Dict[str, Any]]:
    if raw is None:
        return []
    if isinstance(raw, list):
        return [_normalize_item(r, fallback_kind) for r in raw]
    # Single object
    return [_normalize_item(raw, fallback_kind)]

_SER_CACHE: Dict[Tuple[str,int], str] = {}
_SER_LOCK = threading.Lock()

def _float_compact(v: Any):  # best-effort trimming for floats
    if isinstance(v, float):
        return float(f"{v:.3f}")
    if isinstance(v, list):
        return [_float_compact(x) for x in v]
    if isinstance(v, dict):
        return {k: _float_compact(val) for k, val in v.items()}
    return v

def export_json(interactions_map: Dict[str, Any], cache_key: str | None = None, *, already_columnar: bool = False) -> str:
    settings = get_settings()
    use_msgpack = os.getenv('MOLBRIDGE_USE_MSGPACK', '0') in {'1','true','True'}
    # Fast path: if interactions_map is already in columnar payload form (lists of aligned fields)
    # we wrap it under a marker key to remain distinguishable if needed downstream.
    if already_columnar and settings.performance_mode:
        # Expect structure: {type: {residue1: [...], ...}}
        payload = interactions_map  # trust caller
        if use_msgpack:
            try:
                import msgpack, base64  # type: ignore
                packed = msgpack.packb({'_columnar': True, 'format_version': 1, 'schema': 'columnar_v1', 'data': payload}, use_bin_type=True)
                data = 'mp:' + base64.b64encode(packed).decode('ascii')
            except Exception:
                if _HAVE_ORJSON:
                    data = orjson.dumps({'_columnar': True, 'format_version': 1, 'schema': 'columnar_v1', 'data': payload}, option=orjson.OPT_SERIALIZE_NUMPY).decode('utf-8')
                else:
                    data = json.dumps({'_columnar': True, 'format_version': 1, 'schema': 'columnar_v1', 'data': payload}, separators=(',', ':'))
        else:
            if _HAVE_ORJSON:
                data = orjson.dumps({'_columnar': True, 'format_version': 1, 'schema': 'columnar_v1', 'data': payload}, option=orjson.OPT_SERIALIZE_NUMPY).decode('utf-8')
            else:
                data = json.dumps({'_columnar': True, 'format_version': 1, 'schema': 'columnar_v1', 'data': payload}, separators=(',', ':'))
        if cache_key:
            with _SER_LOCK:
                _SER_CACHE[(cache_key, len(payload))] = data
        return data
    payload = {k: normalize_interaction_collection(v, k) for k, v in interactions_map.items()}
    if settings.performance_mode:
        # Compact numeric representation & no indentation for bandwidth
        compact = {k: [_float_compact(item) for item in coll] for k, coll in payload.items()}
        if use_msgpack:
            try:
                import msgpack, base64  # type: ignore
                packed = msgpack.packb(compact, use_bin_type=True)
                data = 'mp:' + base64.b64encode(packed).decode('ascii')
            except Exception:
                if _HAVE_ORJSON:
                    data = orjson.dumps(compact, option=orjson.OPT_SERIALIZE_NUMPY).decode('utf-8')
                else:
                    data = json.dumps(compact, separators=(',', ':'))
        else:
            if _HAVE_ORJSON:
                data = orjson.dumps(compact, option=orjson.OPT_SERIALIZE_NUMPY).decode('utf-8')
            else:
                data = json.dumps(compact, separators=(',', ':'))
        if cache_key:
            with _SER_LOCK:
                _SER_CACHE[(cache_key, len(compact))] = data
        return data
    # Non-performance path retains readability
    if _HAVE_ORJSON:
        return orjson.dumps(payload, option=orjson.OPT_INDENT_2 | orjson.OPT_SERIALIZE_NUMPY).decode('utf-8')
    return json.dumps(payload, indent=2)

def export_csv(interactions_map: Dict[str, Any]) -> str:
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(["group"] + NORMALIZED_KEYS)
    for k, collection in interactions_map.items():
        norm = normalize_interaction_collection(collection, k)
        for item in norm:
            row = [k] + [item.get(col, "") for col in NORMALIZED_KEYS]
            writer.writerow(row)
    return output.getvalue()

def export_markdown_summary(interactions_map: Dict[str, Any]) -> str:
    lines = ["# Interaction Summary","", "| Type | Count |", "|------|-------|"]
    for k, collection in interactions_map.items():
        count = len(normalize_interaction_collection(collection, k))
        lines.append(f"| {k} | {count} |")
    return "\n".join(lines) + "\n"

def export_funnel_metrics_json(instrumentation: Dict[str, Dict[str, Any]]) -> str:
    """Export detector funnel metrics (raw/candidate/accepted) to JSON.

    Accepts the HighPerformanceProcessor _detector_instrumentation mapping.
    Only includes known funnel keys to keep payload concise.
    """
    if not instrumentation:
        return '{}'
    filtered: Dict[str, Any] = {}
    keys = {'raw_pairs','candidate_pairs','accepted_pairs','acceptance_ratio','candidate_density','adapt_reason','threshold_changed'}
    for name, instr in instrumentation.items():
        if not isinstance(instr, dict):
            continue
        filtered[name] = {k: instr.get(k) for k in keys if k in instr}
    settings = get_settings()
    if _HAVE_ORJSON:
        if settings.performance_mode:
            return orjson.dumps(filtered, option=orjson.OPT_SERIALIZE_NUMPY).decode('utf-8')
        return orjson.dumps(filtered, option=orjson.OPT_INDENT_2 | orjson.OPT_SERIALIZE_NUMPY).decode('utf-8')
    return json.dumps(filtered, indent=0 if settings.performance_mode else 2)

__all__ = ['export_json','export_csv','export_markdown_summary','normalize_interaction_collection','export_funnel_metrics_json']
