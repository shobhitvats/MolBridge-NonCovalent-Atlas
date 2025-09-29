"""Experimental distributed execution scaffold (Ray / Dask optional).

Usage:
  Set MOLBRIDGE_DISTRIBUTED=ray (or dask) to enable. Falls back silently
  if libraries are missing. This is a thin abstraction over the existing
  HighPerformanceProcessor to prepare for future scale-out validation.
"""
from __future__ import annotations

import os
from typing import List, Any, Dict
from .parallel_processor import HighPerformanceProcessor

try:  # Optional imports
    import ray  # type: ignore
    _HAVE_RAY = True
except Exception:  # pragma: no cover
    _HAVE_RAY = False

try:
    from dask.distributed import Client, LocalCluster  # type: ignore
    _HAVE_DASK = True
except Exception:  # pragma: no cover
    _HAVE_DASK = False


class DistributedProcessor:
    """Facade providing a future-compatible distributed API.

    Current implementation keeps logic minimal: distributes whole-structure
    processing tasks across workers. Fine-grained detector partitioning will
    be explored once baseline benefit is demonstrated.
    """
    def __init__(self, backend: str | None = None, config=None):
        self.backend = backend or os.getenv("MOLBRIDGE_DISTRIBUTED") or ""
        self.config = config
        self._ray_inited = False
        self._dask_client = None

        if self.backend == "ray" and _HAVE_RAY:
            if not ray.is_initialized():  # pragma: no cover - ray not in default deps
                ray.init(ignore_reinit_error=True, include_dashboard=False, log_to_driver=False)
            self._ray_inited = True
        elif self.backend == "dask" and _HAVE_DASK:
            self._dask_cluster = LocalCluster(processes=True, n_workers=min(4, os.cpu_count() or 2))
            self._dask_client = Client(self._dask_cluster)

    def shutdown(self):  # pragma: no cover - not critical path for unit tests
        if self._ray_inited:
            import ray
            try:
                ray.shutdown()
            except Exception:
                pass
        if self._dask_client:
            try:
                self._dask_client.close()
                self._dask_cluster.close()
            except Exception:
                pass

    def process_structures(self, structures: List[Any], detector_classes: List[Any]) -> List[Dict]:
        """Distribute structure processing based on selected backend.

        Fallback chain: ray -> dask -> local high performance processor.
        """
        if self.backend == "ray" and self._ray_inited:
            import ray

            @ray.remote
            def _process_one(serialized, det_classes):  # simple remote wrapper
                # Reconstruct pickled structure (ray handles serialization)
                from .parallel_processor import HighPerformanceProcessor
                proc = HighPerformanceProcessor()
                return proc.process_interactions_parallel(serialized, det_classes)

            futures = [_process_one.remote(s, detector_classes) for s in structures]
            results = ray.get(futures)
            return results
        elif self.backend == "dask" and self._dask_client:
            def _process_one(s):
                proc = HighPerformanceProcessor()
                return proc.process_interactions_parallel(s, detector_classes)
            futures = [self._dask_client.submit(_process_one, s) for s in structures]
            return self._dask_client.gather(futures)
        else:
            # Local fallback
            local = HighPerformanceProcessor(config=self.config)
            return [local.process_interactions_parallel(s, detector_classes) for s in structures]


__all__ = ["DistributedProcessor"]
