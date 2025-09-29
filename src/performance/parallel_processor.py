"""
High-performance parallel processing module for protein interaction analysis.
Optimizes computation using multiprocessing, threading, and various performance techniques.
"""

import os
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from typing import Dict, List, Any, Callable, Tuple, Optional
import time
import logging
from functools import partial
import numpy as np
from dataclasses import dataclass
import asyncio
import threading
from queue import Queue, Empty
try:
    from multiprocessing import shared_memory
    _HAVE_SHM = True
except Exception:  # pragma: no cover
    _HAVE_SHM = False
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False
    print("psutil not available - performance monitoring will be limited")

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def _run_detector_process_task(detector_class, config, structure, detector_name: str, start_ts: float, shared_meta: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Isolated process task runner.
    Avoids capturing outer self; builds detector, runs appropriate method, returns serializable dict.
    """
    t0 = time.time()
    # Minimal mapping (duplicated to keep function self-contained for pickling)
    detector_methods = {
        'hydrogenbond': 'detect_hydrogen_bonds',
        'halogenbond': 'detect_halogen_bonds',
        'chpi': 'detect_ch_pi_interactions',
        'pipi': 'detect_pi_pi_interactions',
        'ionicinteraction': 'detect_ionic_interactions',
        'hydrophobiccontact': 'detect_hydrophobic_contacts',
        'chalcogenbond': 'detect_chalcogen_bonds',
        'pnictogenbond': 'detect_pnictogen_bonds',
        'tetrelbond': 'detect_tetrel_bonds',
        'anionpi': 'detect_anion_pi_interactions',
        'npistar': 'detect_n_pi_star_interactions',
        'dispersion': 'detect_dispersion_interactions'
    }
    try:
        # Inject shared precomputed features into FeatureStore (best-effort)
        if shared_meta:
            try:
                from analysis.feature_store import get_feature_store  # local import
                fs = get_feature_store()
                bundle = fs.get_bundle(structure)
                # Direct Python object copies
                if 'rings' in shared_meta and shared_meta['rings'] is not None:
                    bundle.rings = shared_meta['rings']
                if 'charged_centers' in shared_meta and shared_meta['charged_centers'] is not None:
                    bundle.charged_group_centroids = shared_meta['charged_centers']
                # Decode shared memory packed ring arrays if present
                if 'rings_shm' in shared_meta and shared_meta['rings_shm']:
                    meta = shared_meta['rings_shm']
                    try:  # pragma: no cover
                        from multiprocessing import shared_memory as _shm
                        import numpy as _np
                        shm_block = _shm.SharedMemory(name=meta['name'])
                        arr = _np.ndarray(tuple(meta['shape']), dtype=meta['dtype'], buffer=shm_block.buf)
                        # Each row: centroid(0:3) + normal(3:6)
                        ring_dicts = []
                        for row in arr:
                            centroid = row[0:3].astype('float32')
                            normal = row[3:6].astype('float32')
                            ring_dicts.append({'centroid': centroid, 'normal': normal})
                        # Only populate if not already from full Python objects
                        if not bundle.rings:
                            bundle.rings = ring_dicts
                        # We don't unlink here (producer cleans up); just close
                        shm_block.close()
                    except Exception:
                        pass
                # Decode charged centers shared memory
                if 'charged_shm' in shared_meta and shared_meta['charged_shm']:
                    meta = shared_meta['charged_shm']
                    try:  # pragma: no cover
                        from multiprocessing import shared_memory as _shm
                        import numpy as _np
                        ch_block = _shm.SharedMemory(name=meta['name'])
                        ch_arr = _np.ndarray(tuple(meta['shape']), dtype=meta['dtype'], buffer=ch_block.buf)
                        if bundle.charged_group_centroids is None:
                            bundle.charged_group_centroids = [
                                {'centroid': row.copy()} for row in ch_arr
                            ]
                        ch_block.close()
                    except Exception:
                        pass
                # Decode acidic centers shared memory
                if 'acidic_shm' in shared_meta and shared_meta['acidic_shm']:
                    meta = shared_meta['acidic_shm']
                    try:  # pragma: no cover
                        from multiprocessing import shared_memory as _shm
                        import numpy as _np
                        ac_block = _shm.SharedMemory(name=meta['name'])
                        ac_arr = _np.ndarray(tuple(meta['shape']), dtype=meta['dtype'], buffer=ac_block.buf)
                        if bundle.acidic_group_centroids is None:
                            bundle.acidic_group_centroids = [
                                {'centroid': row.copy()} for row in ac_arr
                            ]
                        ac_block.close()
                    except Exception:
                        pass
                # Decode HBond donor/acceptor participant coordinate arrays
                if 'hbond_donor_shm' in shared_meta and shared_meta['hbond_donor_shm']:
                    meta = shared_meta['hbond_donor_shm']
                    try:  # pragma: no cover
                        from multiprocessing import shared_memory as _shm
                        import numpy as _np
                        d_block = _shm.SharedMemory(name=meta['name'])
                        d_arr = _np.ndarray(tuple(meta['shape']), dtype=meta['dtype'], buffer=d_block.buf)
                        # Store as lightweight pseudo atoms (only coord retained) if bundle lacks donors
                        if not getattr(bundle, 'hb_donors', None):
                            bundle.hb_donors = [ {'coord': row.copy()} for row in d_arr ]
                        d_block.close()
                    except Exception:
                        pass
                if 'hbond_acceptor_shm' in shared_meta and shared_meta['hbond_acceptor_shm']:
                    meta = shared_meta['hbond_acceptor_shm']
                    try:  # pragma: no cover
                        from multiprocessing import shared_memory as _shm
                        import numpy as _np
                        a_block = _shm.SharedMemory(name=meta['name'])
                        a_arr = _np.ndarray(tuple(meta['shape']), dtype=meta['dtype'], buffer=a_block.buf)
                        if not getattr(bundle, 'hb_acceptors', None):
                            bundle.hb_acceptors = [ {'coord': row.copy()} for row in a_arr ]
                        a_block.close()
                    except Exception:
                        pass
            except Exception:  # pragma: no cover
                pass
        if config:
            detector = detector_class(config)
        else:
            try:
                from utils.config import AppConfig  # local import for process
                detector = detector_class(AppConfig())
            except Exception:
                detector = detector_class()
        if hasattr(detector, 'set_performance_mode'):
            try:
                detector.set_performance_mode(high_performance=True)
            except Exception:
                pass
        method_name = detector_methods.get(detector_name, 'detect')
        if not hasattr(detector, method_name):
            return {"detector": detector_name, "interactions": [], "duration": 0.0, "error": f"missing {method_name}"}
        method = getattr(detector, method_name)
        interactions = method(structure)
        instrumentation = getattr(detector, 'instrumentation', None)
        if hasattr(detector, 'to_dict_list') and callable(getattr(detector, 'to_dict_list')):
            interactions_list = detector.to_dict_list(interactions)
        elif isinstance(interactions, list):
            interactions_list = interactions
        else:
            interactions_list = []
        result_payload = {
            "detector": detector_name,
            "interactions": interactions_list,
            "duration": time.time() - t0,
            "queued_delay": t0 - start_ts,
            "instrumentation": instrumentation if isinstance(instrumentation, dict) else None
        }
        # Structured streaming log (process path)
        try:  # pragma: no cover
            import json
            result_payload['_log'] = json.dumps({
                'event': 'detector_complete',
                'detector': detector_name,
                'count': len(interactions_list),
                'duration_sec': round(result_payload['duration'], 6),
                'queued_delay_sec': round(result_payload['queued_delay'], 6),
                'strategy': 'process',
                'timestamp': time.time()
            })
        except Exception:
            pass
        return result_payload
    except Exception as e:  # noqa: BLE001
        return {"detector": detector_name, "interactions": [], "duration": time.time() - t0, "error": str(e)}

@dataclass
class PerformanceMetrics:
    """Track performance metrics for optimization."""
    start_time: float
    end_time: float
    cpu_usage: float
    memory_usage: float
    threads_used: int
    tasks_completed: int
    
    @property
    def duration(self) -> float:
        return self.end_time - self.start_time
    
    @property
    def tasks_per_second(self) -> float:
        return self.tasks_completed / self.duration if self.duration > 0 else 0

class OptimalThreadCalculator:
    """Calculate optimal number of threads based on system resources and task type."""
    
    @staticmethod
    def get_optimal_threads(task_type: str = "io_bound") -> int:
        """Calculate optimal thread count based on task type and system resources."""
        cpu_count = os.cpu_count() or 2
        
        if PSUTIL_AVAILABLE:
            try:
                available_memory_gb = psutil.virtual_memory().available / (1024**3)
            except:
                available_memory_gb = 4  # Default assumption
        else:
            available_memory_gb = 4  # Default assumption
        
        if task_type == "cpu_bound":
            # For CPU-bound tasks, use number of cores
            return min(cpu_count, 8)  # Cap at 8 to prevent overwhelming
        elif task_type == "io_bound":
            # For I/O-bound tasks, can use more threads
            return min(cpu_count * 4, 32)  # More threads for I/O
        elif task_type == "mixed":
            # For mixed workloads
            return min(cpu_count * 2, 16)
        else:
            return cpu_count

class HighPerformanceProcessor:
    """High-performance processor for protein interaction analysis."""
    
    def __init__(self, max_workers: Optional[int] = None, config=None):
        """Initialize with adaptive worker selection.
        Selection order (first non-None wins):
          1. Explicit argument `max_workers`
          2. ENV override MOLBRIDGE_MAX_WORKERS
          3. config.processing.max_workers
          4. Heuristic (cores-1 for >=3 cores else 1)
        A soft cap (12) is applied to avoid resource thrash on dev machines.
        If psutil is available we also clamp further if available memory < 2GB.
        """
        env_override = None
        try:
            if os.getenv('MOLBRIDGE_MAX_WORKERS'):
                env_override = int(os.getenv('MOLBRIDGE_MAX_WORKERS'))
        except Exception:
            env_override = None
        config_limit = None
        if config is not None and hasattr(config, 'processing') and hasattr(config.processing, 'max_workers'):
            try:
                config_limit = int(config.processing.max_workers)
            except Exception:
                config_limit = None
        if max_workers is not None:
            chosen = max_workers
            source = 'arg'
        elif env_override is not None:
            chosen = env_override
            source = 'env'
        elif config_limit is not None:
            chosen = config_limit
            source = 'config'
        else:
            cores = os.cpu_count() or 2
            heuristic = max(1, cores - 1)  # leave one core free
            chosen = min(heuristic, OptimalThreadCalculator.get_optimal_threads("mixed"))
            source = 'heuristic'
        soft_capped = min(chosen, 12)
        if PSUTIL_AVAILABLE:
            try:
                avail_gb = psutil.virtual_memory().available / (1024**3)
                if avail_gb < 2 and soft_capped > 4:
                    soft_capped = 4  # memory protection
            except Exception:
                pass
        self.max_workers = max(1, soft_capped)
        self.config = config
        # Per-run metrics + simple result cache
        self.performance_metrics: List[PerformanceMetrics] = []
        self.cache: Dict[str, Any] = {}
        self._detector_timings: Dict[str, float] = {}
        # Strategy selection (env-gated experimental process pool)
        self.use_process_pool = os.getenv('MOLBRIDGE_USE_PROCESS_POOL', '0') in {'1', 'true', 'True'}
        # We keep a feature flag so we can permanently disable after first failure
        self._process_pool_supported = True
        # Optional detector instance pool (stateless detectors reused across runs)
        self._detector_pool: Dict[str, List[Any]] = {}
        try:
            from utils.settings import get_settings  # local import to avoid early import cycles
            s = get_settings()
            self._pool_enabled = getattr(s, 'enable_detector_pool', False)
            self._pool_cap = int(getattr(s, 'detector_pool_max_size', 64))
        except Exception:
            self._pool_enabled = False
            self._pool_cap = 0
        if self._pool_enabled:
            logger.info("Detector pooling enabled (cap=%d)", self._pool_cap)

        # Set optimal process/thread limits (defensive caps for native libs)
        os.environ['OMP_NUM_THREADS'] = str(self.max_workers)
        os.environ['NUMEXPR_MAX_THREADS'] = str(self.max_workers)

        logger.info(
            f"Initialized HighPerformanceProcessor workers={self.max_workers} (source={source}, raw_choice={chosen})"
        )
        # Ring buffer cap for performance metrics (avoid unbounded growth)
        self._perf_cap = 1000
    
    def process_interactions_parallel(self, 
                                    structure, 
                                    detector_classes: List[Any],
                                    config=None,
                                    chunk_size: Optional[int] = None,
                                    dry_run: bool = False,
                                    progress_callback: Optional[Callable[[Dict[str, Any]], None]] = None) -> Dict[str, List]:
        """Process all interaction types in parallel with optimal performance.

        Parameters
        ----------
        dry_run : bool
            If True, detectors are instantiated but their heavy acceptance / final geometry
            filtering path is skipped when detector provides a `preview`/`preview_candidates`
            method or instrumentation exposing raw/candidate counts. Fallback: run full but
            discard heavy detail list to reduce memory transfer.
        progress_callback : callable(dict) -> None
            Invoked after each detector completion with a payload containing:
              {detector, count, duration, strategy, dry_run, timestamp, funnel?}
        """
        start_time = time.time()
        phase_times = {}
        
        if PSUTIL_AVAILABLE:
            try:
                start_cpu = psutil.cpu_percent()
                start_memory = psutil.virtual_memory().percent
            except:
                start_cpu = 0
                start_memory = 0
        else:
            start_cpu = 0
            start_memory = 0
        
        # Auto performance profile application (heuristic or user-selected profile)
        try:
            from performance.auto_tune import apply_auto_flags  # local import to avoid hard dependency earlier
            auto_meta = apply_auto_flags(structure, detector_classes)
        except Exception as _e:  # pragma: no cover
            auto_meta = {'error': str(_e)}

        # Calculate optimal chunk size based on structure size
        if chunk_size is None:
            residue_count = len(list(structure.get_residues()))
            chunk_size = max(10, residue_count // (self.max_workers * 2))
        
        # Optional task graph precompute (shared feature extraction)
        task_graph_enabled = os.getenv('MOLBRIDGE_TASK_GRAPH', '1') in {'1', 'true', 'True'}
        shared_features: Dict[str, Any] = {}
        if task_graph_enabled:
            try:  # best-effort; failures just skip shared opt
                _t0 = time.time()
                from analysis.feature_store import get_feature_store
                fs = get_feature_store()
                shared_features['coords'] = fs.ensure_coords(structure)
                shared_features['rings'] = fs.ensure_rings(structure)
                shared_features['charged_centers'] = fs.ensure_charged_centers(structure)
                try:
                    shared_features['acidic_centers'] = fs.ensure_acidic_centers(structure)
                except Exception:
                    shared_features['acidic_centers'] = []
                shared_features['hb_participants'] = fs.ensure_hbond_participants(structure)
                phase_times['feature_extract'] = time.time() - _t0
            except Exception as e:  # pragma: no cover
                logger.debug(f"Task graph precompute skipped: {e}")

        # Create task queue (capture start times for timing)
        tasks = []
        for entry in detector_classes:
            # Support either a bare detector class or a (class, method_name) tuple as some tests supply
            if isinstance(entry, tuple):
                detector_class, explicit_method = entry[0], (entry[1] if len(entry) > 1 else None)
            else:
                detector_class, explicit_method = entry, None
            try:
                detector_name = detector_class.__name__.replace('Detector', '').lower()
            except AttributeError:
                # Skip invalid entries defensively
                continue
            tasks.append((detector_name, detector_class, structure, config, time.time(), shared_features))

        # Shared memory distribution (coords + ring centroids/normals + charged/acidic + HBond donors/acceptors) if enabled and using process pool
        shm_meta = None
        ring_shm_meta = None
        charged_shm_meta = None
        acidic_shm_meta = None
        hbond_donor_shm_meta = None
        hbond_acceptor_shm_meta = None
        use_shm = os.getenv('MOLBRIDGE_USE_SHM', '0') in {'1','true','True'} and self.use_process_pool and _HAVE_SHM
        if use_shm:
            try:
                coords = shared_features.get('coords')
                if coords is not None:
                    buf = coords.astype('float32', copy=False)
                    shm_block = shared_memory.SharedMemory(create=True, size=buf.nbytes)
                    shm_buf = np.ndarray(buf.shape, dtype='float32', buffer=shm_block.buf)
                    shm_buf[:] = buf[:]
                    shm_meta = {'name': shm_block.name,'shape': buf.shape,'dtype':'float32'}
                    shared_features['coords_shm'] = shm_meta
                rings = shared_features.get('rings') or []
                if rings:
                    centroids = np.vstack([r['centroid'] for r in rings]).astype('float32', copy=False)
                    normals = np.vstack([r['normal'] for r in rings]).astype('float32', copy=False)
                    packed = np.concatenate([centroids, normals], axis=1)
                    ring_block = shared_memory.SharedMemory(create=True, size=packed.nbytes)
                    ring_buf = np.ndarray(packed.shape, dtype='float32', buffer=ring_block.buf)
                    ring_buf[:] = packed[:]
                    ring_shm_meta = {'name': ring_block.name,'shape': packed.shape,'dtype':'float32','layout':'centroid[0:3]|normal[3:6]'}
                    shared_features['rings_shm'] = ring_shm_meta
                charged = shared_features.get('charged_centers') or []
                if charged:
                    c_coords = np.vstack([c['centroid'] for c in charged]).astype('float32', copy=False)
                    ch_block = shared_memory.SharedMemory(create=True, size=c_coords.nbytes)
                    ch_buf = np.ndarray(c_coords.shape, dtype='float32', buffer=ch_block.buf)
                    ch_buf[:] = c_coords[:]
                    charged_shm_meta = {'name': ch_block.name,'shape': c_coords.shape,'dtype':'float32'}
                    shared_features['charged_shm'] = charged_shm_meta
                acidic = shared_features.get('acidic_centers') or []
                if acidic:
                    a_coords = np.vstack([a['centroid'] for a in acidic]).astype('float32', copy=False)
                    ac_block = shared_memory.SharedMemory(create=True, size=a_coords.nbytes)
                    ac_buf = np.ndarray(a_coords.shape, dtype='float32', buffer=ac_block.buf)
                    ac_buf[:] = a_coords[:]
                    acidic_shm_meta = {'name': ac_block.name,'shape': a_coords.shape,'dtype':'float32'}
                    shared_features['acidic_shm'] = acidic_shm_meta
                hb_participants = shared_features.get('hb_participants')
                if hb_participants:
                    donors, acceptors = hb_participants
                    if donors:
                        d_coords = np.vstack([d.get_coord() for d in donors]).astype('float32', copy=False)
                        d_block = shared_memory.SharedMemory(create=True, size=d_coords.nbytes)
                        d_buf = np.ndarray(d_coords.shape, dtype='float32', buffer=d_block.buf)
                        d_buf[:] = d_coords[:]
                        hbond_donor_shm_meta = {'name': d_block.name,'shape': d_coords.shape,'dtype':'float32'}
                        shared_features['hbond_donor_shm'] = hbond_donor_shm_meta
                    if acceptors:
                        a_coords2 = np.vstack([a.get_coord() for a in acceptors]).astype('float32', copy=False)
                        a_block2 = shared_memory.SharedMemory(create=True, size=a_coords2.nbytes)
                        a_buf2 = np.ndarray(a_coords2.shape, dtype='float32', buffer=a_block2.buf)
                        a_buf2[:] = a_coords2[:]
                        hbond_acceptor_shm_meta = {'name': a_block2.name,'shape': a_coords2.shape,'dtype':'float32'}
                        shared_features['hbond_acceptor_shm'] = hbond_acceptor_shm_meta
            except Exception as e:  # pragma: no cover
                logger.debug(f"Shared memory setup failed: {e}")
                shm_meta = ring_shm_meta = charged_shm_meta = acidic_shm_meta = None
                hbond_donor_shm_meta = hbond_acceptor_shm_meta = None

        results = {}
        strategy = 'thread'
        # Attempt process pool if enabled
        collected_instrumentation: Dict[str, Dict[str, Any]] = {}
        detect_phase_start = time.time()
        if self.use_process_pool and self._process_pool_supported:
            try:
                strategy = 'process'
                with ProcessPoolExecutor(max_workers=self.max_workers) as proc_exec:
                    future_to_detector = {}
                    for detector_name, detector_class, structure, config, start_ts, shared in tasks:
                        # Prepare lightweight serializable shared meta (rings & charged centers + ring_shm meta)
                        shared_meta = {
                            'rings': shared.get('rings'),
                            'charged_centers': shared.get('charged_centers'),
                            'acidic_centers': shared.get('acidic_centers'),
                            'rings_shm': shared.get('rings_shm'),
                            'charged_shm': shared.get('charged_shm'),
                            'acidic_shm': shared.get('acidic_shm'),
                            'hbond_donor_shm': shared.get('hbond_donor_shm'),
                            'hbond_acceptor_shm': shared.get('hbond_acceptor_shm')
                        } if isinstance(shared, dict) else None
                        fut = proc_exec.submit(_run_detector_process_task, detector_class, config, structure, detector_name, start_ts, shared_meta)
                        future_to_detector[fut] = detector_name
                    for fut in as_completed(future_to_detector):
                        dname = future_to_detector[fut]
                        try:
                            info = fut.result(timeout=60)
                            interactions_list = info.get('interactions', [])
                            results[dname] = interactions_list
                            if not hasattr(self, '_detector_timings'):
                                self._detector_timings = {}
                            self._detector_timings[dname] = info.get('duration', 0.0)
                            if info.get('instrumentation'):
                                collected_instrumentation[dname] = info['instrumentation']
                            # Emit structured log if provided
                            if '_log' in info:
                                try:
                                    logger.info(info['_log'])
                                except Exception:
                                    pass
                            logger.info(f"[process] Completed {dname}: {len(interactions_list)} in {info.get('duration',0.0):.3f}s")
                        except Exception as e:  # noqa: BLE001
                            logger.error(f"Process task failed for {dname}: {e}")
                            results[dname] = []
            except Exception as e:
                logger.warning(f"Process pool fallback due to error: {e}; reverting to threads")
                self._process_pool_supported = False
                results = {}
        if not results:  # thread strategy or fallback
            strategy = 'thread'
            try:
                with ThreadPoolExecutor(max_workers=self.max_workers) as thread_executor:
                    future_to_detector = {}
                    for detector_name, detector_class, structure, config, start_ts, shared in tasks:
                        if dry_run and hasattr(detector_class, 'preview'):  # detector implements lightweight candidate preview
                            fut = thread_executor.submit(self._preview_wrapper, detector_name, detector_class, structure, config, start_ts)
                        else:
                            fut = thread_executor.submit(self._timed_detect_wrapper, detector_name, detector_class, structure, config, start_ts)
                        future_to_detector[fut] = detector_name
                    for future in as_completed(future_to_detector):
                        detector_name = future_to_detector[future]
                        try:
                            info = future.result(timeout=30)
                            interactions_list = info.get('interactions', [])
                            if dry_run:
                                # shrink payload - we only need counts for preview; discard list if large
                                interactions_list = []
                            results[detector_name] = interactions_list
                            if not hasattr(self, '_detector_timings'):
                                self._detector_timings = {}
                            self._detector_timings[detector_name] = info.get('duration', 0.0)
                            instr = info.get('instrumentation')
                            if instr:
                                collected_instrumentation[detector_name] = instr
                            if progress_callback:
                                try:
                                    progress_callback({
                                        'detector': detector_name,
                                        'count': len(info.get('interactions', [])),
                                        'duration': info.get('duration', 0.0),
                                        'strategy': 'thread',
                                        'dry_run': dry_run,
                                        'funnel': instr,
                                        'timestamp': time.time()
                                    })
                                except Exception:
                                    pass
                            # Timed wrapper currently doesn’t return instrumentation; could be extended
                            # Human-readable log
                            logger.info(f"[thread] Completed {detector_name}: {len(interactions_list)} interactions in {info.get('duration', 0.0):.3f}s (queue_delay={info.get('queued_delay', 0.0):.3f}s)")
                            # Structured streaming log (JSON) for downstream aggregation
                            try:  # pragma: no cover
                                import json, time as _t
                                logger.info(json.dumps({
                                    'event': 'detector_complete',
                                    'detector': detector_name,
                                    'count': len(interactions_list),
                                    'duration_sec': round(info.get('duration', 0.0), 6),
                                    'queued_delay_sec': round(info.get('queued_delay', 0.0), 6),
                                    'strategy': 'thread',
                                    'timestamp': _t.time()
                                }))
                            except Exception:
                                pass
                        except Exception as e:
                            logger.error(f"Error in {detector_name}: {str(e)}")
                            results[detector_name] = []
            except Exception as e:
                logger.error(f"Error in parallel processing: {e}")
                results = self._fallback_sequential_processing(tasks)
        
        phase_times['detection'] = time.time() - detect_phase_start
        # Record performance metrics
        end_time = time.time()
        phase_times['total'] = end_time - start_time
        if 'feature_extract' not in phase_times:
            phase_times['feature_extract'] = 0.0
        phase_times['overhead'] = max(0.0, phase_times['total'] - (phase_times['feature_extract'] + phase_times['detection']))
        
        if PSUTIL_AVAILABLE:
            try:
                end_cpu = psutil.cpu_percent()
                end_memory = psutil.virtual_memory().percent
            except:
                end_cpu = 0
                end_memory = 0
        else:
            end_cpu = 0
            end_memory = 0
        
        metrics = PerformanceMetrics(
            start_time=start_time,
            end_time=end_time,
            cpu_usage=(start_cpu + end_cpu) / 2,
            memory_usage=(start_memory + end_memory) / 2,
            threads_used=self.max_workers,
            tasks_completed=len(tasks)
        )
        self.performance_metrics.append(metrics)

        logger.info(
            f"Parallel processing ({strategy}) completed in {metrics.duration:.2f}s ({metrics.tasks_per_second:.1f} tasks/sec)"
        )
        if auto_meta and isinstance(auto_meta, dict):
            try:  # pragma: no cover
                import json as _json, time as _t
                logger.info(_json.dumps({
                    'event': 'auto_perf_profile',
                    'meta': auto_meta,
                    'timestamp': _t.time()
                }))
            except Exception:
                pass
        # Structured log emission (simple JSON line)
        try:
            import json
            structured_record = {
                'event': 'batch_complete',
                'strategy': strategy,
                'duration_sec': round(metrics.duration, 4),
                'tasks': len(tasks),
                'workers': self.max_workers,
                'phase_times': {k: round(v,4) for k,v in phase_times.items()},
                'detector_timings': getattr(self, '_detector_timings', {}),
                'funnel_summary': self._build_funnel_summary(collected_instrumentation),
                'timestamp': time.time()
            }
            logger.info(json.dumps(structured_record))
        except Exception:  # pragma: no cover
            pass

        # Attach instrumentation for callers that want to integrate into summaries
        self._detector_instrumentation = collected_instrumentation
        # Metrics channel emission (best-effort)
        try:  # pragma: no cover
            from utils.logging_config import emit_metrics
            emit_metrics({
                'event': 'batch_metrics',
                'strategy': strategy,
                'duration_sec': round(metrics.duration, 6),
                'workers': self.max_workers,
                'detector_timings': getattr(self, '_detector_timings', {}),
                'instrumentation': collected_instrumentation
            })
        except Exception:
            pass
        # Cleanup shared memory if allocated
        if 'shm_block' in locals():
            try:
                shm_block.close()
                shm_block.unlink()
            except Exception:
                pass
        if 'ring_block' in locals():
            try:
                ring_block.close()
                ring_block.unlink()
            except Exception:
                pass
        if 'ch_block' in locals():
            try:
                ch_block.close()
                ch_block.unlink()
            except Exception:
                pass
        if 'ac_block' in locals():
            try:
                ac_block.close()
                ac_block.unlink()
            except Exception:
                pass
        return results

    def _build_funnel_summary(self, instrumentation: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """Aggregate funnel metrics (raw -> candidate -> accepted) across detectors.

        Returns per-detector acceptance ratios and a global summary. Safe on partial data.
        """
        summary = {'detectors': {}, 'global': {}}
        total_raw = total_candidate = total_accepted = 0
        for name, instr in (instrumentation or {}).items():
            try:
                raw_pairs = int(instr.get('raw_pairs') or instr.get('pairs_considered') or 0)
                candidate_pairs = int(instr.get('candidate_pairs') or instr.get('pairs_within_cutoff') or 0)
                accepted_pairs = int(instr.get('accepted_pairs') or instr.get('hbonds') or 0)
                ratio = (accepted_pairs / candidate_pairs) if candidate_pairs else 0.0
                summary['detectors'][name] = {
                    'raw': raw_pairs,
                    'candidate': candidate_pairs,
                    'accepted': accepted_pairs,
                    'acceptance_ratio': round(ratio, 6)
                }
                total_raw += raw_pairs
                total_candidate += candidate_pairs
                total_accepted += accepted_pairs
            except Exception:
                continue
        summary['global'] = {
            'raw': total_raw,
            'candidate': total_candidate,
            'accepted': total_accepted,
            'acceptance_ratio': round((total_accepted/ total_candidate) if total_candidate else 0.0, 6)
        }
        return summary
    
    def _detect_interactions_optimized(self, detector_name: str, detector_class: Any, structure, config=None) -> List[Dict]:
        """Optimized interaction detection with caching and performance improvements."""
        
        # Check cache first
        cache_key = f"{detector_name}_{id(structure)}"
        if cache_key in self.cache:
            return self.cache[cache_key]
        
        # Mapping of detector names to their method names
        detector_methods = {
            'hydrogenbond': 'detect_hydrogen_bonds',
            'halogenbond': 'detect_halogen_bonds', 
            'chpi': 'detect_ch_pi_interactions',
            'pipi': 'detect_pi_pi_interactions',
            'ionicinteraction': 'detect_ionic_interactions',
            'hydrophobiccontact': 'detect_hydrophobic_contacts',
            'chalcogenbond': 'detect_chalcogen_bonds',
            'pnictogenbond': 'detect_pnictogen_bonds',
            'tetrelbond': 'detect_tetrel_bonds',
            'anionpi': 'detect_anion_pi_interactions',
            'npistar': 'detect_n_pi_star_interactions',
            'dispersion': 'detect_dispersion_interactions'
        }
        
        try:
            # Initialize detector with config if available
            if config:
                detector = detector_class(config)
            else:
                # Try to create a minimal config for backwards compatibility
                try:
                    from utils.config import AppConfig
                    detector = detector_class(AppConfig())
                except:
                    # Some detectors might not need config - try without it
                    try:
                        detector = detector_class()
                    except:
                        logger.error(f"Could not initialize {detector_name} detector")
                        return []
            
            # Apply performance optimizations
            if hasattr(detector, 'set_performance_mode'):
                detector.set_performance_mode(high_performance=True)
            
            # Get the correct method name for this detector
            method_name = detector_methods.get(detector_name, 'detect')
            if not hasattr(detector, method_name):
                logger.error(f"Detector {detector_name} does not have method {method_name}")
                return []
            
            # Detect interactions using the correct method
            detect_method = getattr(detector, method_name)
            interactions = detect_method(structure)
            
            # Convert to optimized format
            if hasattr(detector, 'to_dict_list') and callable(getattr(detector, 'to_dict_list')):
                result = detector.to_dict_list(interactions)
            elif isinstance(interactions, list):
                result = interactions
            else:
                result = []
            
            # Cache results for reuse
            self.cache[cache_key] = result
            
            return result
            
        except Exception as e:
            logger.error(f"Error in {detector_name} detection: {str(e)}")
            return []

    def _timed_detect_wrapper(self, detector_name: str, detector_class: Any, structure, config, start_ts: float) -> Dict[str, Any]:
        t0 = time.time()
        # Instantiate detector directly here so we can capture instrumentation
        instrumentation = None
        interactions_list: List[Any] = []
        try:
            detector = None
            pooled_list = None
            if getattr(self, '_pool_enabled', False):
                pooled_list = self._detector_pool.get(detector_name)
                if pooled_list:
                    try:
                        detector = pooled_list.pop()
                    except Exception:
                        detector = None
            if detector is None:
                if config:
                    detector = detector_class(config)
                else:
                    try:
                        from utils.config import AppConfig
                        detector = detector_class(AppConfig())
                    except Exception:
                        detector = detector_class()
            # Note: detectors should be largely stateless between runs. If stateful fields accumulate
            # (e.g., instrumentation dicts), detectors should clear them internally; pooling assumes
            # safe reuse. If a detector is not safe for pooling, it can expose an attribute
            # `__no_pool__ = True` to opt-out in future enhancement.
            if hasattr(detector, 'set_performance_mode'):
                try:
                    detector.set_performance_mode(high_performance=True)
                except Exception:
                    pass
            # Resolve method mapping similarly to process path
            method_map = {
                'hydrogenbond': 'detect_hydrogen_bonds',
                'halogenbond': 'detect_halogen_bonds',
                'chpi': 'detect_ch_pi_interactions',
                'pipi': 'detect_pi_pi_interactions',
                'ionicinteraction': 'detect_ionic_interactions',
                'hydrophobiccontact': 'detect_hydrophobic_contacts',
                'chalcogenbond': 'detect_chalcogen_bonds',
                'pnictogenbond': 'detect_pnictogen_bonds',
                'tetrelbond': 'detect_tetrel_bonds',
                'anionpi': 'detect_anion_pi_interactions',
                'npistar': 'detect_n_pi_star_interactions',
                'dispersion': 'detect_dispersion_interactions'
            }
            mname = method_map.get(detector_name, 'detect')
            method = getattr(detector, mname, None)
            if method is None:
                interactions = []
            else:
                interactions = method(structure)
            if hasattr(detector, 'to_dict_list') and callable(getattr(detector, 'to_dict_list')):
                interactions_list = detector.to_dict_list(interactions)
            elif isinstance(interactions, list):
                interactions_list = interactions
            instrumentation = getattr(detector, 'instrumentation', None)
            # Return detector to pool if enabled and pool not over capacity
            if self._pool_enabled and detector is not None:
                bucket = self._detector_pool.setdefault(detector_name, [])
                if len(bucket) < self._pool_cap:
                    bucket.append(detector)
        except Exception as e:
            logging.getLogger(__name__).debug(f"Detector wrapper error for {detector_name}: {e}")
            interactions_list = []
        return {"detector": detector_name, "interactions": interactions_list, "duration": time.time() - t0, "queued_delay": t0 - start_ts, "instrumentation": instrumentation if isinstance(instrumentation, dict) else None}

    def _preview_wrapper(self, detector_name: str, detector_class, structure, config, start_ts: float):
        """Attempt a lightweight detector preview for dry_run mode.
        Looks for `preview` or `preview_candidates` method returning a dict with counts.
        Fallback: run full detection but keep only length meta.
        """
        t0 = time.time()
        try:
            detector = detector_class(config) if config else detector_class()
            preview_fn = None
            for cand in ('preview', 'preview_candidates'):
                if hasattr(detector, cand) and callable(getattr(detector, cand)):
                    preview_fn = getattr(detector, cand)
                    break
            if preview_fn:
                preview_data = preview_fn(structure)
                instrumentation = getattr(detector, 'instrumentation', None)
                count = 0
                if isinstance(preview_data, dict):
                    count = preview_data.get('accepted_pairs') or preview_data.get('candidates') or preview_data.get('raw_pairs') or 0
                return {
                    'detector': detector_name,
                    'interactions': [None] * int(count),
                    'duration': time.time() - t0,
                    'queued_delay': t0 - start_ts,
                    'instrumentation': instrumentation if isinstance(instrumentation, dict) else None
                }
            # fallback to full detection
            return self._timed_detect_wrapper(detector_name, detector_class, structure, config, start_ts)
        except Exception as e:  # noqa: BLE001
            return {
                'detector': detector_name,
                'interactions': [],
                'duration': time.time() - t0,
                'queued_delay': t0 - start_ts,
                'error': str(e)
            }
    
    def _fallback_sequential_processing(self, tasks: List[Tuple]) -> Dict[str, List]:
        """Fallback sequential processing if parallel processing fails."""
        logger.warning("Falling back to sequential processing")
        
        results = {}
        for task in tasks:
            if len(task) == 4:  # With config
                detector_name, detector_class, structure, config = task
            else:  # Without config (backwards compatibility)
                detector_name, detector_class, structure = task
                config = None
                
            try:
                results[detector_name] = self._detect_interactions_optimized(detector_name, detector_class, structure, config)
            except Exception as e:
                logger.error(f"Error in sequential {detector_name}: {str(e)}")
                results[detector_name] = []
        
        return results
    
    def process_batch_structures(self, structures: List[Any], detector_classes: List[Any]) -> List[Dict]:
        """Process multiple structures in parallel batches."""
        
        if not structures:
            return []
        
        # Calculate optimal batch size
        batch_size = max(1, len(structures) // self.max_workers)
        if batch_size == 0:
            batch_size = 1
        
        # Create batches
        batches = [structures[i:i + batch_size] for i in range(0, len(structures), batch_size)]
        
        results = []
        
        # Process batches in parallel
        with ProcessPoolExecutor(max_workers=min(len(batches), self.max_workers)) as executor:
            batch_futures = [
                executor.submit(self._process_structure_batch, batch, detector_classes)
                for batch in batches
            ]
            
            for future in as_completed(batch_futures):
                try:
                    batch_results = future.result(timeout=120)  # 2 minute timeout per batch
                    results.extend(batch_results)
                except Exception as e:
                    logger.error(f"Error processing batch: {str(e)}")
        
        return results
    
    def _process_structure_batch(self, structure_batch: List[Any], detector_classes: List[Any]) -> List[Dict]:
        """Process a batch of structures."""
        batch_results = []
        
        for structure in structure_batch:
            try:
                result = self.process_interactions_parallel(structure, detector_classes, self.config)
                batch_results.append(result)
            except Exception as e:
                logger.error(f"Error processing structure in batch: {str(e)}")
                batch_results.append({})
        
        return batch_results
    
    def get_performance_summary(self) -> Dict[str, Any]:
        """Get performance summary of recent operations."""
        # Enforce ring buffer cap
        if len(self.performance_metrics) > self._perf_cap:
            # Drop oldest half to amortize cost
            del self.performance_metrics[:len(self.performance_metrics) - self._perf_cap]
        if not self.performance_metrics:
            return {"message": "No performance data available"}
        
        recent_metrics = self.performance_metrics[-10:]  # Last 10 operations
        
        avg_duration = np.mean([m.duration for m in recent_metrics])
        avg_cpu = np.mean([m.cpu_usage for m in recent_metrics])
        avg_memory = np.mean([m.memory_usage for m in recent_metrics])
        avg_throughput = np.mean([m.tasks_per_second for m in recent_metrics])
        
        return {
            "average_duration": f"{avg_duration:.2f}s",
            "average_cpu_usage": f"{avg_cpu:.1f}%",
            "average_memory_usage": f"{avg_memory:.1f}%",
            "average_throughput": f"{avg_throughput:.1f} tasks/sec",
            "total_operations": len(self.performance_metrics),
            "max_workers": self.max_workers,
            "cache_size": len(self.cache)
        }
    
    def clear_cache(self):
        """Clear the results cache to free memory."""
        self.cache.clear()
        logger.info("Performance cache cleared")

class AsyncInteractionProcessor:
    """Asynchronous processor for real-time interaction updates."""
    
    def __init__(self, max_concurrent: int = 10):
        self.max_concurrent = max_concurrent
        self.active_tasks = set()
        self.result_queue = Queue()
    
    async def process_interaction_async(self, structure, detector_class) -> Dict[str, Any]:
        """Process single interaction type asynchronously."""
        
        loop = asyncio.get_event_loop()
        
        # Run CPU-bound task in thread pool
        with ThreadPoolExecutor(max_workers=2) as executor:
            future = loop.run_in_executor(
                executor, 
                self._sync_detect_interactions, 
                structure, 
                detector_class
            )
            
            result = await future
            return result
    
    def _sync_detect_interactions(self, structure, detector_class) -> Dict[str, Any]:
        """Synchronous wrapper for async processing."""
        try:
            detector = detector_class()
            interactions = detector.detect(structure)
            
            if hasattr(interactions, 'to_dict_list'):
                return {
                    'detector': detector_class.__name__,
                    'interactions': interactions.to_dict_list(),
                    'count': len(interactions.to_dict_list())
                }
            else:
                return {
                    'detector': detector_class.__name__,
                    'interactions': [],
                    'count': 0
                }
        except Exception as e:
            logger.error(f"Async detection error: {str(e)}")
            return {
                'detector': detector_class.__name__,
                'interactions': [],
                'count': 0,
                'error': str(e)
            }

# Global processor instance
_global_processor = None

def get_global_processor(config=None) -> HighPerformanceProcessor:
    """Get or create global processor instance."""
    global _global_processor
    if _global_processor is None:
        _global_processor = HighPerformanceProcessor(config=config)
    return _global_processor

def optimize_system_for_performance():
    """Apply system-level optimizations for better performance."""
    
    # Set optimal environment variables
    os.environ['PYTHONHASHSEED'] = '0'  # Reproducible hashing
    os.environ['OMP_DYNAMIC'] = 'false'  # Disable dynamic thread adjustment
    
    # Set process priority (if possible)
    if PSUTIL_AVAILABLE:
        try:
            p = psutil.Process(os.getpid())
            if hasattr(p, 'nice'):
                p.nice(-5)  # Higher priority (requires permissions)
        except:
            pass
    
    # Optimize garbage collection
    import gc
    gc.set_threshold(700, 10, 10)  # More aggressive GC
    
    logger.info("System optimizations applied")

# Apply optimizations on import (env gated)
if os.getenv('MOLBRIDGE_AUTO_TUNE', '1') not in {'0', 'false', 'False'}:
    optimize_system_for_performance()
