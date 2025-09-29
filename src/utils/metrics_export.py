"""Structured metrics export (JSON) for performance & detector timings."""
from __future__ import annotations
import json, os, time
from pathlib import Path
from typing import Any, Dict

def export_metrics(data: Dict[str, Any], out_dir: str | Path | None = None, prefix: str = "metrics") -> str | None:
    if out_dir is None:
        out_dir = os.getenv("MOLBRIDGE_METRICS_DIR", "./metrics")
    p = Path(out_dir)
    try:
        p.mkdir(parents=True, exist_ok=True)
        ts = int(time.time())
        fp = p / f"{prefix}_{ts}.json"
        with fp.open('w') as f:
            json.dump(data, f, indent=2)
        return str(fp)
    except Exception:
        return None

__all__ = ["export_metrics"]