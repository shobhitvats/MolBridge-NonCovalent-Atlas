"""
Cache management for Protein Interaction Explorer.
Handles caching of PDB files, analysis results, and computed data.
"""

import os
import json
import pickle
import hashlib
import gzip
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Optional, Dict, Union
import requests
import diskcache as dc
from loguru import logger

class CacheManager:
    """Manages caching of PDB files and analysis results."""
    
    def __init__(self, cache_dir: Path, size_limit: int = 2**30):  # 1GB default
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
    
    def _generate_key(self, *args) -> str:
        """Generate a cache key from arguments."""
        key_string = "|".join(str(arg) for arg in args)
        return hashlib.md5(key_string.encode()).hexdigest()
    
    def get_pdb_file(self, pdb_id: str, assembly: str = "biological") -> Optional[str]:
        """Get PDB file content from cache or download if not cached."""
        cache_key = f"pdb_{pdb_id}_{assembly}"
        
        # Check cache first
        cached_content = self.cache.get(cache_key)
        if cached_content:
            logger.info(f"Retrieved {pdb_id} from cache")
            return cached_content
        
        # Download if not in cache
        content = self._download_pdb(pdb_id, assembly)
        if content:
            # Cache the content with 7-day expiration
            self.cache.set(cache_key, content, expire=7*24*3600)
            logger.info(f"Downloaded and cached {pdb_id}")
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
        # Create cache key from PDB ID, interaction type, and parameters
        param_hash = hashlib.md5(json.dumps(parameters, sort_keys=True).encode()).hexdigest()
        cache_key = f"analysis_{pdb_id}_{interaction_type}_{param_hash}"
        
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
            self.cache.set(cache_key, data, expire=7*24*3600)
            logger.debug(f"Cached analysis result for {pdb_id}:{interaction_type}")
            
        except Exception as e:
            logger.error(f"Failed to cache analysis result: {e}")
    
    def get_analysis_result(self, 
                          pdb_id: str, 
                          interaction_type: str, 
                          parameters: Dict[str, Any]) -> Optional[Any]:
        """Get cached analysis result."""
        param_hash = hashlib.md5(json.dumps(parameters, sort_keys=True).encode()).hexdigest()
        cache_key = f"analysis_{pdb_id}_{interaction_type}_{param_hash}"
        
        # Try regular key first
        result = self.cache.get(cache_key)
        if result is not None:
            return result
        
        # Try compressed key
        compressed_key = cache_key + "_compressed"
        compressed_result = self.cache.get(compressed_key)
        if compressed_result is not None:
            try:
                decompressed = gzip.decompress(compressed_result)
                return pickle.loads(decompressed)
            except Exception as e:
                logger.error(f"Failed to decompress cached result: {e}")
                # Remove corrupted cache entry
                del self.cache[compressed_key]
        
        return None
    
    def cache_structure_info(self, pdb_id: str, info: Dict[str, Any]) -> None:
        """Cache structure metadata."""
        cache_key = f"structure_info_{pdb_id}"
        self.cache.set(cache_key, info, expire=30*24*3600)  # 30 days
    
    def get_structure_info(self, pdb_id: str) -> Optional[Dict[str, Any]]:
        """Get cached structure metadata."""
        cache_key = f"structure_info_{pdb_id}"
        return self.cache.get(cache_key)
    
    def cache_batch_results(self, batch_id: str, results: Dict[str, Any]) -> None:
        """Cache batch analysis results."""
        cache_key = f"batch_{batch_id}"
        
        try:
            # Compress batch results as they can be large
            data = pickle.dumps(results)
            compressed_data = gzip.compress(data)
            
            self.cache.set(cache_key, compressed_data, expire=24*3600)  # 1 day
            logger.info(f"Cached batch results for {batch_id}")
            
        except Exception as e:
            logger.error(f"Failed to cache batch results: {e}")
    
    def get_batch_results(self, batch_id: str) -> Optional[Dict[str, Any]]:
        """Get cached batch results."""
        cache_key = f"batch_{batch_id}"
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
        """Get cache statistics."""
        try:
            stats = {
                "total_items": len(self.cache),
                "cache_size_mb": sum(
                    os.path.getsize(f) for f in self.cache_dir.rglob("*") if f.is_file()
                ) / (1024 * 1024),
                "pdb_files_cached": len(list(self.cache.iterkeys("pdb_*"))),
                "analysis_results_cached": len(list(self.cache.iterkeys("analysis_*"))),
                "batch_results_cached": len(list(self.cache.iterkeys("batch_*"))),
                "structure_info_cached": len(list(self.cache.iterkeys("structure_info_*"))),
            }
            return stats
        except Exception as e:
            logger.error(f"Failed to get cache stats: {e}")
            return {}
    
    def clear_cache(self, pattern: Optional[str] = None) -> int:
        """Clear cache entries matching pattern (or all if no pattern)."""
        if pattern:
            count = 0
            for key in list(self.cache.iterkeys(pattern)):
                del self.cache[key]
                count += 1
            return count
        else:
            count = len(self.cache)
            self.cache.clear()
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
    
    def export_cache_manifest(self) -> Dict[str, Any]:
        """Export cache manifest for debugging/monitoring."""
        manifest = {
            "generated_at": datetime.now().isoformat(),
            "cache_stats": self.get_cache_stats(),
            "entries": []
        }
        
        try:
            for key in self.cache:
                entry_info = {
                    "key": key,
                    "size_bytes": len(pickle.dumps(self.cache[key])),
                    "type": key.split("_")[0] if "_" in key else "unknown"
                }
                manifest["entries"].append(entry_info)
        except Exception as e:
            logger.error(f"Failed to generate cache manifest: {e}")
        
        return manifest
    
    def __del__(self):
        """Cleanup on object destruction."""
        try:
            self.cleanup_temp_files()
            self.cache.close()
        except Exception:
            pass
