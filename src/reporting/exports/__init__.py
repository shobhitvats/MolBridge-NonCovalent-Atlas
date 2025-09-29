"""Central export facade building on serializer (non-breaking)."""
from __future__ import annotations
from typing import Dict, Any
from reporting.serializer import export_json, export_csv, export_markdown_summary

_FORMATTERS = {
    "json": export_json,
    "csv": export_csv,
    "md": export_markdown_summary,
    "markdown": export_markdown_summary,
}

def export(format: str, interactions_map: Dict[str, Any]) -> str:
    fmt = format.lower()
    fn = _FORMATTERS.get(fmt)
    if not fn:
        raise ValueError(f"Unsupported export format: {format}")
    return fn(interactions_map)
