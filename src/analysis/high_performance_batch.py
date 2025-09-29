"""High performance batch processing (refactored out of previous monolithic file).

Primary responsibilities:
  * Orchestrate per-structure multi-detector analysis using the unified
    ``DETECTOR_REGISTRY`` mapping.
  * Optional parallel execution through ``HighPerformanceProcessor``.
  * Produce normalized key mapping, detector timing & counts without re-running
    detectors (timing captured inline).
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional
import time
import os
import logging
import traceback

from .registry import DETECTOR_REGISTRY
from .keys import resolve_key
from utils.normalization import normalize_detector_output
from utils.settings import get_settings
from utils.provenance import build_provenance
from utils.metrics_export import export_metrics
from utils.cache import CacheManager
from utils.pdb_handler import PDBHandler
from utils.config import AppConfig
from utils.param_hash import interaction_params_signature
from performance.parallel_processor import HighPerformanceProcessor, get_global_processor
from reporting.serializer import export_json
from utils.structure_hash import compute_structure_coord_hash
from analysis.columnar import get_columnar_store
from utils.hash_utils import stable_hash
import hashlib
import numpy as np

logger = logging.getLogger(__name__)


class HighPerformanceBatchProcessor:
    """Refactored high performance batch processor relying on DETECTOR_REGISTRY.

    Backwards compatible public methods:
      * process_single_protein
      * process_multiple_proteins
    """

    def __init__(self, config: AppConfig | None = None, use_parallel: bool = True, max_workers: Optional[int] = None):
        self.config = config or AppConfig()
        self.use_parallel = use_parallel
        self.pdb_handler = PDBHandler(self.config)
        if use_parallel:
            self.processor = HighPerformanceProcessor(max_workers=max_workers, config=self.config)
        else:
            self.processor = get_global_processor(config=self.config)
            if max_workers:
                self.processor.max_workers = max_workers
        # Lazy cache manager (point at global cache dir from config)
        try:
            self.cache_manager = CacheManager(self.config.cache_dir)
        except Exception:
            self.cache_manager = None
        logger.info("HighPerformanceBatchProcessor init: detectors=%d parallel=%s workers=%s", len(DETECTOR_REGISTRY), use_parallel, self.processor.max_workers)
        self._adaptive_parallel_disabled = False
        self._last_detector_mean_ms: float | None = None
        # Expose last detector objects for external instrumentation snapshot endpoint
        self.last_detectors: list[Any] = []
        # Rolling acceptance ratio history per detector key
        self._rolling_acceptance: dict[str, list[float]] = {}
        from utils.settings import get_settings
        try:
            self._rolling_window = max(1, int(get_settings().rolling_acceptance_window))
        except Exception:
            self._rolling_window = 10
        # Adaptive parallel threshold learning state
        # Samples: (atom_count, ring_count, charged_count, mean_detector_ms)
        self._adaptive_samples: list[tuple[int,int,int,float]] = []
        # Linear model coefficients for features (a_atom, a_ring, a_charged, intercept)
        self._adaptive_model: tuple[float,float,float,float] | None = None
        self._adaptive_persist_path = os.path.join(os.getenv('MOLBRIDGE_STATE_DIR', '/tmp'), 'adaptive_parallel_threshold.json')
        self._load_adaptive_model()

    # ---------------- Internal helpers -----------------
    def _filter_registry(self, interaction_filters: Optional[List[str]]) -> Dict[str, tuple[type, str]]:
        if not interaction_filters:
            return DETECTOR_REGISTRY
        wanted = set(interaction_filters)
        return {k: v for k, v in DETECTOR_REGISTRY.items() if k in wanted}

    # --------------- Adaptive Parallel Threshold Learning ---------------
    def _load_adaptive_model(self):  # pragma: no cover - simple IO
        import json
        try:
            if os.path.isfile(self._adaptive_persist_path):
                with open(self._adaptive_persist_path, 'r') as fh:
                    data = json.load(fh)
                # Backward compatibility: either old 2-param or new 4-param model
                if isinstance(data, dict):
                    if all(k in data for k in ('a_atom','a_ring','a_charged','intercept')):
                        self._adaptive_model = (float(data['a_atom']), float(data['a_ring']), float(data['a_charged']), float(data['intercept']))
                    elif 'slope' in data and 'intercept' in data:  # legacy
                        self._adaptive_model = (float(data['slope']), 0.0, 0.0, float(data['intercept']))
        except Exception:
            self._adaptive_model = None

    def _save_adaptive_model(self):  # pragma: no cover - simple IO
        import json, math
        if not self._adaptive_model:
            return
        a_atom, a_ring, a_charged, intercept = self._adaptive_model
        payload = {'a_atom': a_atom, 'a_ring': a_ring, 'a_charged': a_charged, 'intercept': intercept}
        try:
            os.makedirs(os.path.dirname(self._adaptive_persist_path), exist_ok=True)
            with open(self._adaptive_persist_path, 'w') as fh:
                json.dump(payload, fh)
        except Exception:
            pass

    def _record_adaptive_sample(self, atom_count: int, ring_count: int, charged_count: int, mean_ms: float):
        # Maintain bounded sample list
        self._adaptive_samples.append((atom_count, ring_count, charged_count, mean_ms))
        if len(self._adaptive_samples) > 128:
            self._adaptive_samples = self._adaptive_samples[-128:]
        # Train simple linear model y = a1*atoms + a2*rings + a3*charged + b (all normalized)
        if len(self._adaptive_samples) >= 10:
            import numpy as np
            atoms = np.array([s[0] for s in self._adaptive_samples], dtype='float64')
            rings = np.array([s[1] for s in self._adaptive_samples], dtype='float64')
            charged = np.array([s[2] for s in self._adaptive_samples], dtype='float64')
            y = np.array([s[3] for s in self._adaptive_samples], dtype='float64')
            # Normalize each feature
            def _norm(arr):
                m = arr.max() or 1.0
                return arr / m
            A = np.vstack([_norm(atoms), _norm(rings), _norm(charged), np.ones_like(y)]).T
            try:
                sol, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
                a_atom, a_ring, a_charged, intercept = [float(x) for x in sol]
                self._adaptive_model = (a_atom, a_ring, a_charged, intercept)
                self._save_adaptive_model()
            except Exception:
                pass

    def _predict_threshold(self, atom_count: int, ring_count: int, charged_count: int, default_threshold: float) -> float:
        if not self._adaptive_model or atom_count <= 0:
            return float(default_threshold)
        a_atom, a_ring, a_charged, intercept = self._adaptive_model
        atoms_max = max([s[0] for s in self._adaptive_samples] + [atom_count]) or 1.0
        rings_max = max([s[1] for s in self._adaptive_samples] + [ring_count]) or 1.0
        charged_max = max([s[2] for s in self._adaptive_samples] + [charged_count]) or 1.0
        predicted = (a_atom * (atom_count/atoms_max) + a_ring * (ring_count/rings_max) + a_charged * (charged_count/charged_max) + intercept)
        return max(1.0, min(default_threshold, predicted * 0.75))

    # ---------------- Public API -----------------
    def process_single_protein(self, pdb_id: str, interaction_filters: Optional[List[str]] = None, structure=None) -> Dict[str, Any]:
        start = time.time()
        settings = get_settings()
        try:
            reg = self._filter_registry(interaction_filters)
            if not reg:
                return {"pdb_id": pdb_id, "success": False, "error": "No valid interaction types selected", "processing_time": 0.0}
            if structure is None:
                structure = self.pdb_handler.load_structure(pdb_id)
            if not structure:
                return {"pdb_id": pdb_id, "success": False, "error": f"Failed to load structure {pdb_id}", "processing_time": 0.0}

            # Stable interaction parameter signature (version + hash)
            param_hash = interaction_params_signature(self.config.interactions)

            # Pre-compute coordinate hash for potential compact blob reuse (only if performance_mode)
            structure_hash = None
            if settings.performance_mode:
                try:
                    structure_hash = compute_structure_coord_hash(structure)
                except Exception:  # pragma: no cover
                    structure_hash = None

            # Try bundle + (optionally) compact blob cache first
            if self.cache_manager:
                # First attempt composite package (if structure_hash known)
                if structure_hash:
                    pkg = self.cache_manager.get_interaction_package(pdb_id, param_hash, structure_hash)
                    if pkg and 'canonical' in pkg:
                        cached_bundle = pkg.get('canonical', {})
                        processing_time = time.time() - start
                        summary_counts = {k: len(v) for k, v in cached_bundle.items() if isinstance(v, list)}
                        result_hit = {
                            'pdb_id': pdb_id,
                            'success': True,
                            'interactions': cached_bundle,
                            'summary': {
                                'total_interactions': sum(summary_counts.values()),
                                'interaction_counts': summary_counts,
                                'processing_time': processing_time,
                                'performance_metrics': self.processor.get_performance_summary(),
                                'detector_timings': {k: 0.0 for k in summary_counts},
                                'detector_counts': summary_counts,
                                'analyzed_interactions': interaction_filters or list(reg.keys()),
                                'cache_hit': True,
                                'param_signature': param_hash,
                                'batch_phases': {'detect_ms': 0.0, 'normalize_ms': 0.0, 'serialize_ms': 0.0}
                            },
                            'processing_time': processing_time,
                            'cache_hit': True,
                            'param_signature': param_hash,
                            'structure_hash': pkg.get('structure_hash', structure_hash)
                        }
                        compact_blob = pkg.get('compact')
                        if compact_blob:
                            result_hit['serialized_compact'] = compact_blob
                        return result_hit
                cached_bundle = self.cache_manager.get_interaction_bundle(pdb_id, param_hash)
                compact_blob: bytes | None = None
                if settings.performance_mode and structure_hash:
                    compact_blob = self.cache_manager.get_compact_blob(pdb_id, param_hash, structure_hash)
                if cached_bundle is not None:
                    processing_time = time.time() - start
                    summary_counts = {k: len(v) for k, v in cached_bundle.items() if isinstance(v, list)}
                    result_hit = {
                        'pdb_id': pdb_id,
                        'success': True,
                        'interactions': cached_bundle,
                        'summary': {
                            'total_interactions': sum(summary_counts.values()),
                            'interaction_counts': summary_counts,
                            'processing_time': processing_time,
                            'performance_metrics': self.processor.get_performance_summary(),
                            'detector_timings': {k: 0.0 for k in summary_counts},  # cache hit timings not measured
                            'detector_counts': summary_counts,
                            'analyzed_interactions': interaction_filters or list(reg.keys()),
                            'cache_hit': True,
                            'param_signature': param_hash
                        },
                        'processing_time': processing_time,
                        'cache_hit': True,
                        'param_signature': param_hash
                    }
                    if compact_blob is not None:
                        # Provide passthrough of compact blob (already serialized)
                        try:
                            result_hit['serialized_compact'] = compact_blob.decode('utf-8')
                            result_hit['structure_hash'] = structure_hash
                        except Exception:  # pragma: no cover
                            pass
                    return result_hit

            # Execute detection with timing capture
            detector_timings: Dict[str, float] = {}
            detector_instrumentation: Dict[str, Dict[str, int]] = {}
            # Warm-start KD-tree thresholds once we know atom count (best effort)
            try:
                from utils.kdtree_thresholds import warm_start_thresholds
                atom_est = sum(1 for _ in structure.get_atoms())  # type: ignore
                warm_start_thresholds(atom_est)
            except Exception:
                pass
            detect_t0 = time.time()
            self.last_detectors = []  # reset per invocation
            if self.use_parallel:
                # Collect all detector classes and process in a single parallel invocation to allow
                # shared feature precompute utilization inside processor.
                classes = [cls for cls, _m in reg.values()]
                t0_all = time.time()
                try:
                    res_map = self.processor.process_interactions_parallel(structure, classes, config=self.config)
                except Exception as e:
                    logger.error("Parallel batch processing failed, falling back sequential inside parallel branch: %s", e)
                    res_map = {}
                batch_duration = time.time() - t0_all
                # Map back by registry key using class name heuristic
                # processor returns dict keyed by detector_class base name lowercased w/o 'Detector'
                reverse_index = {}
                for key, (cls, method_name) in reg.items():
                    canonical_name = cls.__name__.replace('Detector','').lower()
                    reverse_index[canonical_name] = key
                raw_results = {}
                for k, interactions in (res_map or {}).items():
                    reg_key = reverse_index.get(k, k)
                    raw_results[reg_key] = interactions
                # Capture timings per detector if processor recorded them; else approximate even split
                proc_timings = getattr(self.processor, '_detector_timings', {})
                total_timing = sum(proc_timings.values()) or batch_duration
                for key, (cls, _m) in reg.items():
                    cname = cls.__name__.replace('Detector','').lower()
                    detector_timings[key] = proc_timings.get(cname, total_timing / max(1,len(reg)))
                # Collect instrumentation from processor (process-mode currently only)
                proc_instr = getattr(self.processor, '_detector_instrumentation', {}) or {}
                proc_detectors = getattr(self.processor, '_detector_instances', None)
                if isinstance(proc_detectors, dict):
                    for cname, det_obj in proc_detectors.items():  # class-name key
                        reg_key = reverse_index.get(cname, cname)
                        if det_obj not in self.last_detectors:
                            self.last_detectors.append(det_obj)
                for k_proc, stats in proc_instr.items():
                    # Map class-name key back to registry key
                    reg_key = reverse_index.get(k_proc, k_proc)
                    if isinstance(stats, dict):
                        detector_instrumentation[reg_key] = {sk: int(sv) if isinstance(sv, (int,float)) else sv for sk, sv in stats.items()}
            else:
                raw_results: Dict[str, List] = {}
                for key, (cls, method_name) in reg.items():
                    t0 = time.time()
                    try:
                        detector = cls(self.config)
                        method = getattr(detector, method_name, None)
                        if method is None:
                            logger.error("Detector %s missing method %s", key, method_name)
                            raw_results[key] = []
                        else:
                            interactions = method(structure)
                            # If detector exposes instrumentation stats, capture them
                            stats = getattr(detector, 'instrumentation', None)
                            if isinstance(stats, dict):
                                detector_instrumentation[key] = {k: int(v) for k, v in stats.items() if isinstance(v, (int, float))}
                            if hasattr(detector, 'to_dict_list') and callable(getattr(detector, 'to_dict_list')):
                                raw_results[key] = detector.to_dict_list(interactions)
                            elif isinstance(interactions, list):
                                raw_results[key] = interactions
                            else:
                                raw_results[key] = []
                            # Track detector object for external snapshot
                            self.last_detectors.append(detector)
                    except Exception as e:  # noqa: BLE001
                        logger.error("Error detecting %s: %s", key, e)
                        raw_results[key] = []
                    finally:
                        self.processor.performance_metrics.append(  # reuse metrics holder for sequential path
                            # lightweight synthetic metric record
                            type(self.processor.performance_metrics[0]) if self.processor.performance_metrics else None  # type: ignore
                        ) if False else None  # placeholder; real granular timing done below
                    raw_results[key] = raw_results.get(key, [])

            detect_ms = (time.time() - detect_t0) * 1000.0
            # Canonical key map; optionally also build normalized in same pass (single-pass fusion)
            canonical: Dict[str, List] = {}
            normalized: Dict[str, List[dict]] = {}
            # Auto-enable columnar if performance mode and force_columnar_in_performance flag set
            use_columnar = getattr(settings, 'enable_columnar', False) or (settings.performance_mode and getattr(settings, 'force_columnar_in_performance', True))
            columnar_store = get_columnar_store() if use_columnar else None
            if columnar_store:
                # Reset global store per invocation (simple strategy)
                columnar_store.blocks.clear(); columnar_store.finalized = False
            # Decide normalization strategy early
            runtime_include = getattr(settings, '_runtime_include_normalized', False)
            do_normalize = (settings.enable_normalization and (runtime_include or not settings.performance_mode)) or runtime_include
            partial_norm: Dict[str, List[dict]] | None = None
            if not do_normalize and settings.performance_mode and getattr(settings, 'partial_normalization', False):
                partial_norm = {}
            for k, v in raw_results.items():
                ckey = resolve_key(k)
                items = v if isinstance(v, list) else [v]
                if columnar_store:
                    for item in items:
                        if isinstance(item, dict):
                            columnar_store.add_interaction(ckey, item)
                        else:
                            columnar_store.add_interaction(ckey, {'type': ckey})
                canonical.setdefault(ckey, []).extend(items)
                # Inline normalization build (single-pass) if enabled
                if do_normalize:
                    try:
                        normalized[ckey] = normalize_detector_output(ckey, items)
                    except Exception:
                        normalized[ckey] = []
                elif partial_norm is not None:
                    out_min: List[dict] = []
                    for itm in items:
                        if isinstance(itm, dict):
                            out_min.append({
                                'type': itm.get('type', ckey),
                                'residue1': itm.get('residue1'),
                                'residue2': itm.get('residue2'),
                                'chain1': itm.get('chain1'),
                                'chain2': itm.get('chain2'),
                                'distance': itm.get('distance')
                            })
                    partial_norm[ckey] = out_min
            columnar_payload = None
            if columnar_store:
                columnar_store.finalize()
                if settings.performance_mode and getattr(settings, 'direct_columnar_serialization', False):
                    # Keep canonical map unchanged for in-memory usage, but prepare columnar payload for serialization fast path
                    columnar_payload = columnar_store.to_columnar_payload()
                else:
                    # Replace canonical map with reconstructed dict lists (legacy outward contract)
                    canonical = columnar_store.to_interactions_map()

            norm_t0 = time.time()
            # In fused path, normalization work already performed during canonical loop; just measure zero-cost placeholder.
            normalize_ms = (time.time() - norm_t0) * 1000.0

            # Build timing / counts directly from result sizes (avoid recomputation)
            detector_counts = {k: len(v) for k, v in canonical.items()}
            # Map timings onto canonical keys (registry keys already canonical form mostly)
            detector_timings = {resolve_key(k): detector_timings.get(k, 0.0) for k in reg.keys() if resolve_key(k) in canonical}
            # Merge instrumentation (counts) keyed by canonical name
            detector_stats = {}
            for k, stats in detector_instrumentation.items():
                ckey = resolve_key(k)
                detector_stats[ckey] = stats
                # Rolling acceptance ratio update
                try:
                    acc_ratio = float(stats.get('acceptance_ratio', 0.0))
                    hist = self._rolling_acceptance.setdefault(ckey, [])
                    hist.append(acc_ratio)
                    if len(hist) > self._rolling_window:
                        del hist[0:len(hist)-self._rolling_window]
                    detector_stats[ckey]['rolling_acceptance_mean'] = round(sum(hist)/len(hist), 6)
                except Exception:
                    pass
            total_interactions = sum(detector_counts.values())
            processing_time = time.time() - start
            perf_summary = self.processor.get_performance_summary()
            result: Dict[str, Any] = {
                'pdb_id': pdb_id,
                'success': True,
                'interactions': canonical,
                'summary': {
                    'total_interactions': total_interactions,
                    'interaction_counts': detector_counts,
                    'processing_time': processing_time,
                    'performance_metrics': perf_summary,
                    'detector_timings': detector_timings,
                    'detector_counts': detector_counts,
                    'detector_stats': detector_stats,
                    'analyzed_interactions': interaction_filters or list(reg.keys()),
                    'param_signature': param_hash
                },
                'processing_time': processing_time,
                'param_signature': param_hash
            }
            if do_normalize:
                result['interactions_normalized'] = normalized
            # Lightweight metrics history ring buffer (in-memory, process scoped)
            try:
                hist_size = getattr(settings, 'metrics_history_size', 500)
                if not hasattr(self, '_metrics_history'):
                    self._metrics_history: list[dict] = []  # type: ignore
                entry = {
                    'pdb_id': pdb_id,
                    'ts': time.time(),
                    'total_interactions': total_interactions,
                    'detector_counts': detector_counts,
                    'detector_stats': {k: {'accepted_pairs': v.get('accepted_pairs'), 'acceptance_ratio': v.get('acceptance_ratio')} for k, v in detector_stats.items()},
                    'mean_detector_ms': self._last_detector_mean_ms,
                }
                self._metrics_history.append(entry)
                if len(self._metrics_history) > hist_size:
                    # Drop oldest half to amortize list slicing cost
                    del self._metrics_history[: len(self._metrics_history) - hist_size]
            except Exception:
                pass
            # Attach compressed JSON blob in performance_mode for downstream quick reuse
            serialize_t0 = time.time()
            if settings.performance_mode:
                try:
                    if structure_hash is None:
                        structure_hash = compute_structure_coord_hash(structure)
                    cache_key = f"{pdb_id}:{param_hash}:{structure_hash}"
                    if columnar_payload is not None:
                        compact_json = export_json(columnar_payload, cache_key=cache_key, already_columnar=True)
                    else:
                        compact_json = export_json(canonical, cache_key=cache_key)
                    result['serialized_compact'] = compact_json
                    result['structure_hash'] = structure_hash
                    # Persist compact blob for reuse
                    if self.cache_manager and compact_json:
                        self.cache_manager.cache_compact_blob(pdb_id, param_hash, structure_hash, compact_json)
                except Exception:  # pragma: no cover
                    pass
            serialize_ms = (time.time() - serialize_t0) * 1000.0
            if settings.enable_provenance:
                try:
                    params_dict = {attr: getattr(self.config.interactions, attr) for attr in dir(self.config.interactions) if not attr.startswith('_') and not callable(getattr(self.config.interactions, attr))}
                except Exception:
                    params_dict = {}
                result['provenance'] = build_provenance(pdb_id, params_dict)
            # Store composite package if possible, else fallback to legacy bundle + compact
            if self.cache_manager:
                try:
                    self.cache_manager.cache_interaction_package(pdb_id, param_hash, structure_hash, canonical, result.get('serialized_compact'))
                except Exception:
                    try:
                        self.cache_manager.cache_interaction_bundle(pdb_id, param_hash, canonical)
                    except Exception:
                        pass
            # Optional metrics export (env flag)
            if settings.enable_normalization and os.getenv('MOLBRIDGE_EXPORT_METRICS', '0') in {'1','true','True'}:
                try:
                    export_metrics({
                        'pdb_id': pdb_id,
                        'processing_time': processing_time,
                        'detector_counts': detector_counts,
                        'detector_timings': detector_timings,
                        'detector_stats': detector_stats,
                        'total_interactions': total_interactions,
                        'cache_hit': result.get('cache_hit', False)
                    })
                except Exception:
                    pass
            if settings.performance_mode and not getattr(settings, 'verbose_detector_logs', False):
                logger.debug("[perf] Completed %s interactions=%d time=%.2fs cache=%s", pdb_id, total_interactions, processing_time, bool(self.cache_manager))
            else:
                logger.info("Completed %s with %d interactions in %.2fs (cache_store=%s)", pdb_id, total_interactions, processing_time, bool(self.cache_manager))
            # Batch phase instrumentation
            phases = {
                'detect_ms': round(detect_ms, 2),
                'normalize_ms': round(normalize_ms, 2),
                'serialize_ms': round(serialize_ms, 2)
            }
            result['summary']['batch_phases'] = phases
            # Aggregate detector phase timings if present in instrumentation
            if detector_instrumentation:
                agg = {'phase_pair_gen_ms': 0.0, 'phase_eval_ms': 0.0, 'phase_build_ms': 0.0}
                counted = 0
                for stats in detector_instrumentation.values():
                    if not isinstance(stats, dict):
                        continue
                    updated = False
                    for k in list(agg.keys()):
                        if k in stats and isinstance(stats[k], (int, float)):
                            agg[k] += float(stats[k])
                            updated = True
                    if updated:
                        counted += 1
                if counted:
                    for k in agg:
                        agg[k] = round(agg[k], 3)
                    result['summary']['detector_phase_totals'] = agg
            if partial_norm is not None:
                result['interactions_partial_norm'] = partial_norm

            # Adaptive parallelism: if parallel was used and mean per-detector time low, disable
            if self.use_parallel and not self._adaptive_parallel_disabled:
                n_det = max(1, len(detector_counts))
                mean_ms = detect_ms / n_det
                self._last_detector_mean_ms = mean_ms
                # Learn structure size (atom count) if available
                atom_count = 0
                try:
                    atom_count = sum(1 for _ in structure.get_atoms())  # type: ignore
                except Exception:
                    pass
                ring_count = 0; charged_count = 0
                try:
                    from analysis.feature_store import get_feature_store
                    fs = get_feature_store(); b = fs.get_bundle(structure)
                    if b.ring_centroids is not None:
                        ring_count = getattr(b.ring_centroids, 'shape', [0])[0] if hasattr(b.ring_centroids, 'shape') else len(b.ring_centroids or [])
                    if b.charged_group_centroids is not None:
                        charged_count = len(b.charged_group_centroids)
                except Exception:
                    pass
                if atom_count:
                    self._record_adaptive_sample(atom_count, ring_count, charged_count, mean_ms)
                # Predict threshold if model trained else fallback to settings
                base_threshold = getattr(get_settings(), 'adaptive_parallel_threshold_ms', 5)
                dynamic_threshold = self._predict_threshold(atom_count, ring_count, charged_count, base_threshold)
                if mean_ms < dynamic_threshold:
                    self._adaptive_parallel_disabled = True
                    self.use_parallel = False
                    self.processor = get_global_processor(config=self.config)
                    logger.info("AdaptiveParallel: Disabled thread pool (mean %.2f ms < %.2f ms threshold, atoms=%d)", mean_ms, dynamic_threshold, atom_count)
            return result
        except Exception as e:  # noqa: BLE001
            logger.error("Failure processing %s: %s", pdb_id, e)
            logger.debug(traceback.format_exc())
            return {'pdb_id': pdb_id, 'success': False, 'error': str(e), 'processing_time': time.time() - start}

    def process_multiple_proteins(self, pdb_ids: List[str], interaction_filters: Optional[List[str]] = None) -> List[Dict[str, Any]]:
        out: List[Dict[str, Any]] = []
        for pid in pdb_ids:
            out.append(self.process_single_protein(pid, interaction_filters))
        return out
