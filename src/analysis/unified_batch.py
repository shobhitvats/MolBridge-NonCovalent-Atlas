"""Unified batch processor facade (non-breaking).

Does NOT modify existing batch_processor classes. Provides an optional cleaner API
relying solely on the detector registry.
"""
from __future__ import annotations
import time
from typing import Any, Dict, List, Optional

from analysis.registry import DETECTOR_REGISTRY, get_detector
from utils.config import AppConfig
from utils.pdb_handler import PDBHandler
from utils.normalization import normalize_detector_output
from utils.settings import get_settings
from utils.provenance import build_provenance  # lightweight (hash only)

class UnifiedBatchProcessor:
    def __init__(self, config: Optional[AppConfig] = None):
        self.config = config or AppConfig()
        self.pdb_handler = PDBHandler(self.config)

    def analyze(self, pdb_id: str, interaction_keys: Optional[List[str]] = None) -> Dict[str, Any]:
        t0 = time.time()
        interaction_keys = interaction_keys or list(DETECTOR_REGISTRY.keys())
        structure = self.pdb_handler.load_structure(pdb_id, assembly=self.config.default_assembly)
        if not structure:
            return {"pdb_id": pdb_id, "success": False, "error": "structure_load_failed"}
        results: Dict[str, List[dict]] = {}
        settings = get_settings()
        normalized_results: Dict[str, Any] = {}
        for key in interaction_keys:
            detector, method_name = get_detector(key, self.config)
            if not detector:
                continue
            method = getattr(detector, method_name, None)
            if not method:
                continue
            try:
                raw = method(structure)
                if hasattr(detector, 'to_dict_list') and isinstance(raw, list):
                    # many existing detectors provide to_dict_list(list)
                    results[key] = detector.to_dict_list(raw)  # type: ignore[attr-defined]
                elif hasattr(detector, 'to_dict_list') and not isinstance(raw, list):
                    # sometimes they store internal list
                    results[key] = detector.to_dict_list()  # type: ignore
                else:
                    results[key] = normalize_detector_output(key, raw, detector)
                if settings.enable_normalization:
                    # produce canonical normalized variant (does not replace legacy list)
                    try:
                        normalized_results[key] = normalize_detector_output(key, results[key], detector)
                    except Exception:
                        pass
            except Exception:
                results[key] = []
        elapsed = time.time() - t0
        summary = {k: len(v) for k, v in results.items()}
        response: Dict[str, Any] = {"pdb_id": pdb_id, "success": True, "processing_time": elapsed, "summary": summary, "interactions": results}
        if settings.enable_normalization and normalized_results:
            response["interactions_normalized"] = normalized_results
        if settings.enable_provenance:
            # collect minimal param dict (distance/angle cutoffs etc.) from config.interactions
            try:
                params_dict = {k: getattr(self.config.interactions, k) for k in dir(self.config.interactions) if not k.startswith('_') and not callable(getattr(self.config.interactions, k))}
            except Exception:
                params_dict = {}
            response["provenance"] = build_provenance(pdb_id, params_dict, assembly=self.config.default_assembly)
        return response


def get_batch_processor(unified: bool = False, **kwargs):  # factory helper
    if not unified:
        # allow env flag override
        try:
            from utils.settings import get_settings
            if get_settings().enable_unified_processor:
                unified = True
        except Exception:
            pass
    if unified:
        return UnifiedBatchProcessor(**kwargs)
    from analysis.batch_processor import HighPerformanceBatchProcessor  # type: ignore
    return HighPerformanceBatchProcessor(**kwargs)
