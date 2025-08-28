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
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False
    print("psutil not available - performance monitoring will be limited")

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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
        """Initialize with optimal worker count and configuration."""
        self.max_workers = max_workers or OptimalThreadCalculator.get_optimal_threads("mixed")
        self.config = config
        self.performance_metrics = []
        self.cache = {}
        
        # Set optimal process/thread limits
        os.environ['OMP_NUM_THREADS'] = str(self.max_workers)
        os.environ['NUMEXPR_MAX_THREADS'] = str(self.max_workers)
        
        logger.info(f"Initialized HighPerformanceProcessor with {self.max_workers} workers")
    
    def process_interactions_parallel(self, 
                                    structure, 
                                    detector_classes: List[Any],
                                    config=None,
                                    chunk_size: Optional[int] = None) -> Dict[str, List]:
        """Process all interaction types in parallel with optimal performance."""
        
        start_time = time.time()
        
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
        
        # Calculate optimal chunk size based on structure size
        if chunk_size is None:
            residue_count = len(list(structure.get_residues()))
            chunk_size = max(10, residue_count // (self.max_workers * 2))
        
        # Create task queue for different interaction types
        tasks = []
        for detector_class in detector_classes:
            detector_name = detector_class.__name__.replace('Detector', '').lower()
            tasks.append((detector_name, detector_class, structure, config))
        
        # Process using ThreadPoolExecutor for I/O bound tasks (PDB parsing, caching)
        # and ProcessPoolExecutor for CPU bound tasks (calculations)
        results = {}
        
        try:
            # Use threading for I/O operations and light computations
            with ThreadPoolExecutor(max_workers=self.max_workers) as thread_executor:
                # Submit all detection tasks
                future_to_detector = {
                    thread_executor.submit(self._detect_interactions_optimized, detector_name, detector_class, structure, config): detector_name
                    for detector_name, detector_class, structure, config in tasks
                }
                
                # Collect results as they complete
                for future in as_completed(future_to_detector):
                    detector_name = future_to_detector[future]
                    try:
                        interactions = future.result(timeout=30)  # 30 second timeout per detector
                        results[detector_name] = interactions
                        logger.info(f"Completed {detector_name}: {len(interactions)} interactions")
                    except Exception as e:
                        logger.error(f"Error in {detector_name}: {str(e)}")
                        results[detector_name] = []
        
        except Exception as e:
            logger.error(f"Error in parallel processing: {str(e)}")
            # Fallback to sequential processing
            results = self._fallback_sequential_processing(tasks)
        
        # Record performance metrics
        end_time = time.time()
        
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
        
        logger.info(f"Parallel processing completed in {metrics.duration:.2f}s "
                   f"({metrics.tasks_per_second:.1f} tasks/sec)")
        
        return results
    
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

# Apply optimizations on import
optimize_system_for_performance()
