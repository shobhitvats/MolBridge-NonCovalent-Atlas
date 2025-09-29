"""Performance baseline script (non-breaking).

Run to produce a JSON snapshot of detector counts & timing for a PDB ID set.
Usage example:
    python -m performance.baseline --pdb 1ABC 2XYZ --out docs/performance/baseline.json
"""
from __future__ import annotations
import argparse
import json
import time
from pathlib import Path
from typing import List

from analysis.registry import DETECTOR_REGISTRY, get_detector
from utils.config import AppConfig
from utils.pdb_handler import PDBHandler


def _parse() -> argparse.Namespace:
    p = argparse.ArgumentParser("baseline")
    p.add_argument("--pdb", nargs="+", required=True, help="PDB IDs")
    p.add_argument("--out", required=True, help="Output JSON path")
    return p.parse_args()


def main():  # pragma: no cover
    args = _parse()
    config = AppConfig()
    ph = PDBHandler(config)
    baseline = {"pdb_ids": args.pdb, "detectors": list(DETECTOR_REGISTRY.keys()), "results": []}
    for pdb_id in args.pdb:
        structure = ph.load_structure(pdb_id, assembly=config.default_assembly)
        if not structure:
            baseline["results"].append({"pdb_id": pdb_id, "error": "load_failed"})
            continue
        entry = {"pdb_id": pdb_id, "detectors": {}}
        for key in DETECTOR_REGISTRY:
            det, method_name = get_detector(key, config)
            method = getattr(det, method_name, None)
            if not method:
                continue
            t0 = time.time()
            try:
                raw = method(structure)
                count = len(raw) if isinstance(raw, list) else 0
            except Exception:
                count = -1
            entry["detectors"][key] = {"count": count, "time": round(time.time() - t0, 4)}
        baseline["results"].append(entry)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(baseline, indent=2))

if __name__ == "__main__":  # pragma: no cover
    main()
