"""
Cache management for Protein Interaction Explorer.
Handles caching of PDB files, analysis results, and computed data.
"""

import os
import json
import pickle
import gzip
from collections import OrderedDict
from contextlib import contextmanager
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Optional, Dict, Union, Iterable, Tuple
import requests
import diskcache as dc
from loguru import logger
from utils.hash_utils import stable_hash
import time

class CacheManager:
    """Manages caching of PDB files and analysis results.

    Optimization additions (2025-09):
      * Unified key normalization (lowercase, no spaces).
      * Transparent lazy decompression for analysis & batch results.
      * Optional in-memory hot cache (small) for most recent lookups.
      * Streaming manifest sizing (avoid full re-pickling large objects).
      * Batch write context manager to group cache.set calls.
    """

    def __init__(self, cache_dir: Path, size_limit: int = 2**30, hot_cache_size: int = 32, schema_version: str | None = None):  # 1GB default
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Initialize disk cache
        self.cache = dc.Cache(str(self.cache_dir / "diskcache"), size_limit=size_limit)

        # Separate directories for different types of cached data
        self.pdb_dir = self.cache_dir / "pdb_files"
        self.results_dir = self.cache_dir / "analysis_results"
        self.temp_dir = self.cache_dir / "temp"

        for directory in [self.pdb_dir, self.results_dir, self.temp_dir]:
            directory.mkdir(exist_ok=True)

        # LRU hot cache to bypass disk for very recent lookups
        env_hot = os.getenv("MOLBRIDGE_CACHE_HOT_SIZE")
        if env_hot:
            try:
                hot_cache_size = int(env_hot)
            except ValueError:
                logger.warning(f"Invalid MOLBRIDGE_CACHE_HOT_SIZE={env_hot}, falling back to provided default {hot_cache_size}")
        self._hot_cache: "OrderedDict[str, Any]" = OrderedDict()
        self._hot_cache_size = max(0, hot_cache_size)

        # Internal flag to coalesce writes inside batch() context
        self._batch_active = False
        self._batch_ops: list[Tuple[str, Any, Optional[int]]] = []

        # In-memory runtime metrics (reset each process start)
        self._metrics = {
            "analysis_hits": 0,
            "analysis_misses": 0,
            "analysis_decompressions": 0,
            "analysis_sets": 0,
            "hot_hits": 0,
            "hot_misses": 0,
            "approx_bytes": 0,
            "evictions_soft": 0,
            "aged_evictions": 0,
        }
        # Last-access tracking (in-memory only) and optional aging
        self._last_access: dict[str, float] = {}
        try:
            self._max_age_seconds = float(os.getenv("MOLBRIDGE_CACHE_MAX_AGE_SEC", "0"))
        except Exception:
            self._max_age_seconds = 0.0
        # Schema / key versioning (v1 default). Environment override: MOLBRIDGE_CACHE_SCHEMA
        env_schema = os.getenv("MOLBRIDGE_CACHE_SCHEMA")
        self.schema_version = (schema_version or env_schema or "v1").strip()
        # Soft memory limit (MB) for composite/compact layers (optional)
        try:
            self._soft_limit_mb = float(os.getenv("MOLBRIDGE_CACHE_SOFT_LIMIT_MB", "0"))
        except Exception:
            self._soft_limit_mb = 0.0

    # --------------- Memory Accounting Helpers ---------------
    def _approx_size(self, value: Any) -> int:
        try:
            if isinstance(value, (bytes, bytearray)):
                return len(value)
            if isinstance(value, str):
                return len(value.encode('utf-8'))
            # Lightweight pickle size approximation
            import pickle as _p
            return len(_p.dumps(value, protocol=4))
        except Exception:
            return 0

    def _bump_size(self, delta: int):
        if delta <= 0:
            return
        self._metrics["approx_bytes"] += delta
        if self._soft_limit_mb > 0:
            if (self._metrics["approx_bytes"] / (1024*1024)) > self._soft_limit_mb:
                self._enforce_soft_limit()

    def _enforce_soft_limit(self):
        """Best-effort soft eviction: drop oldest non-PDB composite / compact entries until under limit.

        We iterate keys and remove those with prefixes in target_prefixes. diskcache already enforces hard size_limit, this is an extra layer.
        """
        target_prefixes = (f"{self.schema_version}:pkg_", f"{self.schema_version}:compact_", f"{self.schema_version}:bundle_") if self.schema_version else ("pkg_","compact_","bundle_")
        bytes_before = self._metrics.get("approx_bytes", 0)
        removed = 0
        for key in list(self.cache.iterkeys()):
            if key.startswith(target_prefixes):
                try:
                    val = self.cache.pop(key, default=None)
                    removed += 1
                    if val is not None:
                        self._metrics["approx_bytes"] -= self._approx_size(val)
                    if (self._metrics["approx_bytes"] / (1024*1024)) <= self._soft_limit_mb * 0.90:  # leave headroom
                        break
                except Exception:
                    continue
        if removed:
            self._metrics["evictions_soft"] += removed
            # Clamp minimum
            if self._metrics["approx_bytes"] < 0:
                self._metrics["approx_bytes"] = 0
            logger.info(f"SoftCacheEvict: removed={removed} bytes~{bytes_before-self._metrics['approx_bytes']} new_total={self._metrics['approx_bytes']/1024/1024:.2f}MB")

    # ---------------- Internal Utilities -----------------
    def _normalize_key(self, key: str) -> str:
        base = key.strip().lower().replace(" ", "_")
        # Only version selected namespaces; skip pdb_ for compatibility.
        versioned_prefixes = ("analysis_", "structure_info_", "batch_", "bundle_", "compact_")
        if base.startswith(versioned_prefixes) and self.schema_version:
            return f"{self.schema_version}:{base}"
        return base

    def _legacy_key(self, key: str) -> str:
        """Return legacy (pre-versioned) normalized key for backward read attempts."""
        return key.strip().lower().replace(" ", "_")

    def _hot_cache_put(self, key: str, value: Any):
        if self._hot_cache_size <= 0:
            return
        # Promote / insert
        if key in self._hot_cache:
            try:
                self._hot_cache.pop(key)
            except KeyError:
                pass
        self._hot_cache[key] = value
        # Evict oldest if size exceeded
        while len(self._hot_cache) > self._hot_cache_size:
            self._hot_cache.popitem(last=False)

    def _hot_cache_get(self, key: str) -> Optional[Any]:
        if key not in self._hot_cache:
            self._metrics["hot_misses"] += 1
            return None
        # Move to MRU position
        try:
            value = self._hot_cache.pop(key)
        except KeyError:
            return None
        self._hot_cache[key] = value
        self._metrics["hot_hits"] += 1
        self._last_access[key] = time.time()
        return value

    @contextmanager
    def batch(self):
        """Group multiple cache writes to reduce disk metadata overhead.

        Usage:
            with cache_manager.batch():
                cache_manager.cache_analysis_result(...)
                cache_manager.cache_analysis_result(...)
        """
        self._batch_active = True
        try:
            yield
            # Flush queued ops
            for key, value, expire in self._batch_ops:
                try:
                    self.cache.set(key, value, expire=expire)
                except Exception as e:
                    logger.error(f"Batch write failed for {key}: {e}")
        finally:
            self._batch_active = False
            self._batch_ops.clear()
    
    def _generate_key(self, *args) -> str:
        """Generate a cache key from arguments using unified stable_hash (v2)."""
        key_string = "|".join(str(arg) for arg in args)
        return f"v2_{stable_hash(key_string, length=16)}"
    
    def get_pdb_file(self, pdb_id: str, assembly: str = "biological") -> Optional[str]:
        """Get PDB file content from cache or download if not cached."""
        cache_key = f"pdb_{pdb_id}_{assembly}"
        
        # Check cache first
        cached_content = self.cache.get(cache_key)
        if cached_content:
            logger.info(f"Retrieved {pdb_id} from cache")
            self._last_access[cache_key] = time.time()
            return cached_content
        
        # Download if not in cache
        content = self._download_pdb(pdb_id, assembly)
        if content:
            # Cache the content with 7-day expiration
            self.cache.set(cache_key, content, expire=7*24*3600)
            logger.info(f"Downloaded and cached {pdb_id}")
            self._last_access[cache_key] = time.time()
            return content
        
        return None
    
    def _download_pdb(self, pdb_id: str, assembly: str = "biological") -> Optional[str]:
        """Download PDB file from RCSB."""
        try:
            if assembly == "biological":
                url = f"https://files.rcsb.org/download/{pdb_id}.pdb1"
            else:
                url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            return response.text
            
        except requests.RequestException as e:
            logger.error(f"Failed to download {pdb_id}: {e}")
            return None
    
    def cache_analysis_result(self,
                              pdb_id: str,
                              interaction_type: str,
                              parameters: Dict[str, Any],
                              result: Any,
                              compress: bool = True) -> None:
        """Cache analysis result."""
        # Unified parameter hash (v2)
        param_hash = stable_hash(parameters, length=16)
        cache_key = self._normalize_key(f"analysis_{pdb_id}_{interaction_type}_{param_hash}")
        
        try:
            if compress:
                # Compress large results
                data = pickle.dumps(result)
                if len(data) > 10000:  # Compress if larger than 10KB
                    data = gzip.compress(data)
                    cache_key += "_compressed"
            else:
                data = result
            
            # Cache with 7-day expiration
            if self._batch_active:
                self._batch_ops.append((cache_key, data, 7*24*3600))
            else:
                self.cache.set(cache_key, data, expire=7*24*3600)
            self._last_access[cache_key] = time.time()
            # Put small (uncompressed) items into hot cache for quicker re-fetch
            if not cache_key.endswith("_compressed"):
                if isinstance(data, (bytes, bytearray)):
                    if len(data) < 8000:
                        try:
                            obj = pickle.loads(data)
                            self._hot_cache_put(cache_key, obj)
                        except Exception:
                            pass
                else:
                    # Directly cache small python objects too
                    try:
                        # Heuristic: only cache if pickle size < 8KB to avoid large RAM hits
                        approx = pickle.dumps(data)
                        if len(approx) < 8000:
                            self._hot_cache_put(cache_key, data)
                    except Exception:
                        pass
            # Memory accounting (approx)
            try:
                self._bump_size(self._approx_size(data))
            except Exception:
                pass
            self._metrics["analysis_sets"] += 1
            logger.debug(f"Cached analysis result for {pdb_id}:{interaction_type}")
            
        except Exception as e:
            logger.error(f"Failed to cache analysis result: {e}")
    
    def get_analysis_result(self,
                            pdb_id: str,
                            interaction_type: str,
                            parameters: Dict[str, Any]) -> Optional[Any]:
        """Get cached analysis result."""
        param_hash = stable_hash(parameters, length=16)
        versioned_key = self._normalize_key(f"analysis_{pdb_id}_{interaction_type}_{param_hash}")
        cache_key = versioned_key
        
        # Hot cache first
        hot = self._hot_cache_get(cache_key)
        if hot is not None:
            self._metrics["analysis_hits"] += 1
            return hot

        # Try regular key (could be pickled bytes or already object if small)
        result = self.cache.get(cache_key)
        if result is not None:
            self._metrics["analysis_hits"] += 1
            # If bytes, attempt unpickle
            if isinstance(result, (bytes, bytearray)):
                try:
                    obj = pickle.loads(result)
                    self._hot_cache_put(cache_key, obj)
                    self._last_access[cache_key] = time.time()
                    return obj
                except Exception:
                    return result
            # For plain python objects, populate hot cache as well
            try:
                self._hot_cache_put(cache_key, result)
            except Exception:
                pass
            self._last_access[cache_key] = time.time()
            return result
        
        # Try compressed key
        compressed_key = cache_key + "_compressed"
        compressed_result = self.cache.get(compressed_key)
        if compressed_result is not None:
            try:
                decompressed = gzip.decompress(compressed_result)
                obj = pickle.loads(decompressed)
                self._hot_cache_put(compressed_key, obj)
                self._metrics["analysis_hits"] += 1
                self._metrics["analysis_decompressions"] += 1
                return obj
            except Exception as e:
                logger.error(f"Failed to decompress cached result: {e}")
                # Remove corrupted cache entry
                try:
                    del self.cache[compressed_key]
                except Exception:
                    pass
        # Backward compatibility read (legacy unversioned key) if miss so far and schema versioning enabled
        if self.schema_version:
            legacy_base = self._legacy_key(f"analysis_{pdb_id}_{interaction_type}_{param_hash}")
            if legacy_base != cache_key:
                legacy_result = self.cache.get(legacy_base)
                if legacy_result is not None:
                    try:
                        if isinstance(legacy_result, (bytes, bytearray)):
                            obj = pickle.loads(legacy_result)
                        else:
                            obj = legacy_result
                        # Write-through upgrade to versioned key
                        self.cache.set(cache_key, legacy_result, expire=7*24*3600)
                        self._hot_cache_put(cache_key, obj)
                        self._last_access[cache_key] = time.time()
                        self._metrics["analysis_hits"] += 1
                        return obj
                    except Exception:
                        pass
        self._metrics["analysis_misses"] += 1
        return None

    # ---------------- PDB ETag helpers (for async multi-fetch) -----------------
    def get_pdb_etag(self, pdb_id: str, assembly: str = "biological") -> Optional[str]:
        try:
            key = self._normalize_key(f"pdb_etag_{pdb_id}_{assembly}")
            val = self.cache.get(key)
            if isinstance(val, str):
                return val
            return None
        except Exception:
            return None

    def set_pdb_etag(self, pdb_id: str, assembly: str, etag: str):
        try:
            key = self._normalize_key(f"pdb_etag_{pdb_id}_{assembly}")
            self.cache.set(key, etag, expire=7*24*3600)
        except Exception:
            pass
    
    def cache_structure_info(self, pdb_id: str, info: Dict[str, Any]) -> None:
        """Cache structure metadata."""
        cache_key = self._normalize_key(f"structure_info_{pdb_id}")
        if self._batch_active:
            self._batch_ops.append((cache_key, info, 30*24*3600))
        else:
            self.cache.set(cache_key, info, expire=30*24*3600)  # 30 days
    
    def get_structure_info(self, pdb_id: str) -> Optional[Dict[str, Any]]:
        """Get cached structure metadata."""
        cache_key = self._normalize_key(f"structure_info_{pdb_id}")
        hot = self._hot_cache_get(cache_key)
        if hot is not None:
            return hot
        info = self.cache.get(cache_key)
        if isinstance(info, dict):
            self._hot_cache_put(cache_key, info)
        return info
    
    def cache_batch_results(self, batch_id: str, results: Dict[str, Any]) -> None:
        """Cache batch analysis results."""
        cache_key = self._normalize_key(f"batch_{batch_id}")
        
        try:
            # Compress batch results as they can be large
            data = pickle.dumps(results)
            compressed_data = gzip.compress(data)
            
            if self._batch_active:
                self._batch_ops.append((cache_key, compressed_data, 24*3600))
            else:
                self.cache.set(cache_key, compressed_data, expire=24*3600)  # 1 day
            logger.info(f"Cached batch results for {batch_id}")
            
        except Exception as e:
            logger.error(f"Failed to cache batch results: {e}")
    
    def get_batch_results(self, batch_id: str) -> Optional[Dict[str, Any]]:
        """Get cached batch results."""
        cache_key = self._normalize_key(f"batch_{batch_id}")
        compressed_data = self.cache.get(cache_key)
        
        if compressed_data:
            try:
                data = gzip.decompress(compressed_data)
                return pickle.loads(data)
            except Exception as e:
                logger.error(f"Failed to decompress batch results: {e}")
                del self.cache[cache_key]
        
        return None
    
    def save_temp_file(self, data: Union[str, bytes], extension: str = ".tmp") -> str:
        """Save temporary file and return path."""
        import tempfile
        
        with tempfile.NamedTemporaryFile(
            dir=self.temp_dir, 
            delete=False, 
            suffix=extension
        ) as tmp:
            if isinstance(data, str):
                tmp.write(data.encode())
            else:
                tmp.write(data)
            
            return tmp.name
    
    def cleanup_temp_files(self, max_age_hours: int = 24) -> None:
        """Clean up temporary files older than specified hours."""
        cutoff_time = datetime.now() - timedelta(hours=max_age_hours)
        
        for temp_file in self.temp_dir.glob("*"):
            try:
                if datetime.fromtimestamp(temp_file.stat().st_mtime) < cutoff_time:
                    temp_file.unlink()
            except Exception:
                continue
    
    def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache statistics including aging sweep if enabled."""
        try:
            aged = 0
            if self._max_age_seconds > 0 and self._last_access:
                now = time.time()
                for k, ts in list(self._last_access.items()):
                    if (now - ts) > self._max_age_seconds:
                        try:
                            self.cache.pop(k, default=None)
                            self._last_access.pop(k, None)
                            aged += 1
                        except Exception:
                            pass
                if aged:
                    self._metrics['aged_evictions'] += aged
            stats = {
                "total_items": len(self.cache),
                "cache_size_mb": sum(
                    os.path.getsize(f) for f in self.cache_dir.rglob("*") if f.is_file()
                ) / (1024 * 1024),
                "pdb_files_cached": len(list(self.cache.iterkeys("pdb_*"))),
                "analysis_results_cached": len(list(self.cache.iterkeys("analysis_*"))),
                "batch_results_cached": len(list(self.cache.iterkeys("batch_*"))),
                "structure_info_cached": len(list(self.cache.iterkeys("structure_info_*"))),
                "hot_cache_items": len(self._hot_cache),
                "tracked_access_keys": len(self._last_access),
                "max_age_seconds": self._max_age_seconds,
                "aged_evicted_in_sweep": aged,
            }
            stats.update(self._metrics)
            hits = stats.get("analysis_hits", 0)
            misses = stats.get("analysis_misses", 0)
            total = hits + misses or 1
            stats["analysis_hit_ratio"] = round(hits / total, 4)
            hh = stats.get("hot_hits", 0)
            hm = stats.get("hot_misses", 0)
            htotal = hh + hm or 1
            stats["hot_hit_ratio"] = round(hh / htotal, 4)
            return stats
        except Exception as e:
            logger.error(f"Failed to get cache stats: {e}")
            return {}

    def get_cache_metrics(self) -> Dict[str, Any]:
        """Return raw metrics counters (copy)."""
        return dict(self._metrics)
    
    def clear_cache(self, pattern: Optional[str] = None) -> int:
        """Clear cache entries matching pattern (or all if no pattern)."""
        if pattern:
            count = 0
            for key in list(self.cache.iterkeys(pattern)):
                try:
                    del self.cache[key]
                    count += 1
                except Exception:
                    pass
            return count
        count = len(self.cache)
        self.cache.clear()
        self._hot_cache.clear()
        return count
    
    def clear_expired_entries(self) -> int:
        """Clear expired cache entries."""
        initial_count = len(self.cache)
        self.cache.expire()
        final_count = len(self.cache)
        return initial_count - final_count
    
    def preload_common_pdbs(self, pdb_list: list) -> None:
        """Preload commonly used PDB files."""
        logger.info(f"Preloading {len(pdb_list)} PDB files...")
        
        for pdb_id in pdb_list:
            if not self.cache.get(f"pdb_{pdb_id}_biological"):
                self.get_pdb_file(pdb_id, "biological")
            if not self.cache.get(f"pdb_{pdb_id}_asymmetric"):
                self.get_pdb_file(pdb_id, "asymmetric")
    
    def export_cache_manifest(self, limit: int = 200, pattern: str | None = None) -> Dict[str, Any]:
        """Export cache manifest for debugging/monitoring.

        Args:
            limit: Max number of entries to list (avoid huge JSON payloads). Can be overridden
                   by env MOLBRIDGE_MANIFEST_LIMIT.
            pattern: Optional glob-like prefix filter (simple startswith match) e.g. 'analysis_'.
        """
        # ENV override for limit
        env_lim = os.getenv("MOLBRIDGE_MANIFEST_LIMIT")
        if env_lim:
            try:
                limit = max(1, int(env_lim))
            except ValueError:
                logger.warning(f"Invalid MOLBRIDGE_MANIFEST_LIMIT={env_lim}, ignoring")
        manifest = {
            "generated_at": datetime.now().isoformat(),
            "cache_stats": self.get_cache_stats(),
            "entries": [],
            "truncated": False,
            "filter_pattern": pattern,
            "limit": limit,
        }
        SIMPLE_SIZER_TYPES = (str, int, float, bool)
        try:
            matched = 0
            for key in self.cache:
                if pattern and not key.startswith(pattern):
                    continue
                if matched >= limit:
                    manifest["truncated"] = True
                    break
                try:
                    val = self.cache[key]
                    # Fast size heuristics to avoid re-pickling large containers repeatedly.
                    if isinstance(val, (bytes, bytearray)):
                        size_bytes = len(val)
                    elif isinstance(val, SIMPLE_SIZER_TYPES):
                        size_bytes = len(str(val))
                    elif isinstance(val, (list, tuple)) and len(val) <= 50:
                        # approximate via pickling small sequences only
                        size_bytes = len(pickle.dumps(val))
                    elif isinstance(val, dict) and len(val) <= 50:
                        size_bytes = len(pickle.dumps({k: val[k] for k in list(val)[:25]}))
                    else:
                        # For large / unknown types defer to len of pickle of a lightweight header
                        try:
                            size_bytes = len(pickle.dumps(val))
                        except Exception:
                            size_bytes = -1
                    manifest["entries"].append({
                        "key": key,
                        "size_bytes": size_bytes,
                        "type": key.split("_")[0] if "_" in key else "unknown"
                    })
                    matched += 1
                except Exception as inner_e:
                    logger.debug(f"Manifest entry failed for {key}: {inner_e}")
        except Exception as e:
            logger.error(f"Failed to generate cache manifest: {e}")
        return manifest

    # ------------- Bundle Interaction Caching (aggregate) -------------
    def cache_interaction_bundle(self, pdb_id: str, param_hash: str, bundle: Dict[str, Any]) -> None:
        """Store an aggregate bundle of interactions for all types for a structure.

        Keyed as bundle_<pdb_id>_<param_hash>. Compressed with gzip.
        """
        key = self._normalize_key(f"bundle_{pdb_id}_{param_hash}")
        try:
            raw = pickle.dumps(bundle)
            comp = gzip.compress(raw)
            self.cache.set(key, comp, expire=7*24*3600)
            self._bump_size(self._approx_size(comp))
        except Exception as e:  # pragma: no cover
            logger.error(f"Failed to cache bundle {pdb_id}: {e}")

    def get_interaction_bundle(self, pdb_id: str, param_hash: str) -> Optional[Dict[str, Any]]:
        key = self._normalize_key(f"bundle_{pdb_id}_{param_hash}")
        data = self.cache.get(key)
        if not data:
            return None
        try:
            return pickle.loads(gzip.decompress(data))
        except Exception:
            try:
                del self.cache[key]
            except Exception:
                pass
            return None

    # -------- Composite package (canonical + compact) --------
    def cache_interaction_package(self, pdb_id: str, param_hash: str, structure_hash: str | None, canonical: Dict[str, Any], compact_blob: str | None) -> None:
        """Store both canonical interactions and optional compact blob as one entry.

        Layout: gzip-compressed JSON of {canonical: {...}, compact: str|None, structure_hash, version}
        """
        key = self._normalize_key(f"pkg_{pdb_id}_{param_hash}_{structure_hash or 'na'}")
        try:
            import gzip, json
            payload = {
                'version': 1,
                'canonical': canonical,
                'compact': compact_blob,
                'structure_hash': structure_hash
            }
            raw = json.dumps(payload, separators=(',', ':')).encode('utf-8')
            comp = gzip.compress(raw) if len(raw) > 4096 else raw
            self.cache.set(key, comp, expire=7*24*3600)
            self._bump_size(self._approx_size(comp))
        except Exception as e:  # pragma: no cover
            logger.debug(f"Failed to cache composite package: {e}")

    def get_interaction_package(self, pdb_id: str, param_hash: str, structure_hash: str | None):
        key = self._normalize_key(f"pkg_{pdb_id}_{param_hash}_{structure_hash or 'na'}")
        data = self.cache.get(key)
        if not data:
            return None
        import json, gzip
        try:
            if isinstance(data, (bytes, bytearray)):
                buf = bytes(data)
                if buf[:2] == b'\x1f\x8b':  # gzip magic
                    buf = gzip.decompress(buf)
                payload = json.loads(buf.decode('utf-8'))
                return payload
        except Exception:  # pragma: no cover
            try:
                del self.cache[key]
            except Exception:
                pass
        return None

    # ------------- Compact Serialized Blob Caching -------------
    def cache_compact_blob(self, pdb_id: str, param_hash: str, structure_hash: str, blob: bytes | str) -> None:
        """Persist a compact JSON (or similar) blob for quick retrieval.

        The blob should already be a serialized representation of canonical
        interactions (NOT pickled). We store it as raw bytes (utf-8) possibly
        gzip-compressed if large ( >32KB ) to reduce disk IO. Keyed by
        compact_<pdb_id>_<param_hash>_<structure_hash>.
        """
        try:
            key = self._normalize_key(f"compact_{pdb_id}_{param_hash}_{structure_hash}")
            if isinstance(blob, str):
                data_bytes = blob.encode("utf-8")
            else:
                data_bytes = blob
            if len(data_bytes) > 32_000:  # compress larger blobs
                data_bytes = gzip.compress(data_bytes)
                key += "_gz"
            self.cache.set(key, data_bytes, expire=7*24*3600)
            self._bump_size(self._approx_size(data_bytes))
        except Exception as e:  # pragma: no cover
            logger.debug(f"Failed to cache compact blob {pdb_id}: {e}")

    def get_compact_blob(self, pdb_id: str, param_hash: str, structure_hash: str) -> Optional[bytes]:
        """Return previously cached compact blob (raw bytes) if present."""
        base_key = self._normalize_key(f"compact_{pdb_id}_{param_hash}_{structure_hash}")
        # Try uncompressed first
        data = self.cache.get(base_key)
        if data is not None:
            if isinstance(data, (bytes, bytearray)):
                return bytes(data)
            if isinstance(data, str):
                return data.encode("utf-8")
        # Try compressed variant
        gz_key = base_key + "_gz"
        data_gz = self.cache.get(gz_key)
        if data_gz is not None:
            try:
                if isinstance(data_gz, (bytes, bytearray)):
                    return gzip.decompress(data_gz)
            except Exception:
                try:
                    del self.cache[gz_key]
                except Exception:  # pragma: no cover
                    pass
        return None

    # ------------- Parameter Hashing Utilities -------------
    @staticmethod
    def stable_param_hash(parameters: Dict[str, Any], include: Optional[Iterable[str]] = None) -> str:
        """Produce a stable hash for a parameter subset (v2 unified hashing).

        Args:
            parameters: Full parameter dict.
            include: Optional iterable restricting keys considered.
        Returns:
            16-char truncated stable hash string.
        """
        if include is not None:
            sub = {k: parameters[k] for k in sorted(include) if k in parameters}
        else:
            sub = {k: parameters[k] for k in sorted(parameters)}
        return stable_hash(sub, length=16)

    # Convenience high-level helper
    def get_or_compute_analysis(self, pdb_id: str, interaction_type: str, parameters: Dict[str, Any], compute_fn) -> Any:
        """Retrieve cached analysis result or compute & cache.

        compute_fn: Callable returning the result when cache miss.
        """
        cached = self.get_analysis_result(pdb_id, interaction_type, parameters)
        if cached is not None:
            return cached
        result = compute_fn()
        self.cache_analysis_result(pdb_id, interaction_type, parameters, result)
        return result
    
    def __del__(self):
        """Cleanup on object destruction."""
        try:
            self.cleanup_temp_files()
            self.cache.close()
        except Exception:
            pass
