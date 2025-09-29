"""Command-line interface for MolBridge (initial minimal version).

Example:
    python -m molbridge --pdb 1ABC --interactions hydrogenbond,pi_pi

Future: add subcommands (analyze, cache, manifest, benchmark).
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from typing import List

import sys, os
# Ensure project root (parent of this file's directory) is on sys.path when executed directly
_THIS_DIR = Path(__file__).resolve().parent
_ROOT = _THIS_DIR.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from utils.config import load_config, get_interaction_types
from utils.cache import CacheManager
from analysis.registry import DETECTOR_REGISTRY, get_detector
from utils.settings import get_settings


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser("molbridge")
    parser.add_argument("--pdb", required=True, help="PDB ID to analyze")
    parser.add_argument("--interactions", help="Comma-separated interaction keys (default: all)")
    parser.add_argument("--preset", default=None, help="Preset name (default: env/ settings)")
    parser.add_argument("--json", action="store_true", help="Emit JSON result only")
    parser.add_argument("--instrumentation", action="store_true", help="Emit instrumentation snapshot (detector funnel metrics)")
    return parser.parse_args()


def main():  # pragma: no cover - thin wrapper
    args = _parse_args()
    settings = get_settings()
    config = load_config()
    cache = CacheManager(config.cache_dir)

    interaction_keys = get_interaction_types()
    if args.interactions:
        requested = [k.strip() for k in args.interactions.split(",") if k.strip()]
        interaction_keys = [k for k in requested if k in DETECTOR_REGISTRY]
        if not interaction_keys:
            print("No valid interaction keys provided.")
            return 1

    preset = args.preset or settings.cli_default_preset
    if preset not in config.presets:
        preset = "literature_default"

    # Apply preset simply by mutating config.interactions
    for pname, val in config.presets[preset].items():
        if hasattr(config.interactions, pname):
            setattr(config.interactions, pname, val)

    # Load structure
    from utils.pdb_handler import PDBHandler
    ph = PDBHandler(config)
    structure = ph.load_structure(args.pdb)
    if not structure:
        print(f"Failed to load structure {args.pdb}")
        return 2

    results = {}
    detectors_used = []
    for key in interaction_keys:
        detector, method_name = get_detector(key, config)
        if not detector:
            continue
        method = getattr(detector, method_name, None)
        if not method:
            continue
        try:
            detected = method(structure)
            if hasattr(detector, 'to_dict_list'):
                if isinstance(detected, list):
                    results[key] = detector.to_dict_list(detected)
                else:
                    results[key] = detector.to_dict_list()
            else:
                results[key] = detected if isinstance(detected, list) else []
            detectors_used.append(detector)
        except Exception as exc:  # pragma: no cover
            results[key] = []
            print(f"Error in detector {key}: {exc}")

    summary = {k: len(v) for k, v in results.items()}
    payload = {"pdb_id": args.pdb, "preset": preset, "summary": summary, "interactions": results}
    # Coercion utility for entire payload (instrumentation may contain numpy scalars)
    def _coerce_all(obj):
        try:
            import numpy as _np
        except Exception:  # pragma: no cover
            _np = None
        if isinstance(obj, dict):
            return {k: _coerce_all(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_coerce_all(v) for v in obj]
        if _np is not None and isinstance(obj, (_np.generic,)):
            return obj.item()
        return obj

    if args.instrumentation:
        try:
            from utils.instrumentation_snapshot import collect_snapshot
            snap = collect_snapshot(detectors_used)
            payload['instrumentation'] = _coerce_all(snap)
        except Exception as _exc:  # pragma: no cover
            payload['instrumentation_error'] = str(_exc)
    payload = _coerce_all(payload)

    if args.json:
        print(json.dumps(payload))
    else:
        print(json.dumps(payload, indent=2))
    return 0

if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
