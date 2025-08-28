"""
Batch processing module for protein interaction analysis.
Enhanced with high-performance parallel processing capabilities.
"""

import logging
from typing import Dict, List, Any, Optional, Tuple
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import sys

# Add src to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.pdb_handler import PDBHandler
from performance.parallel_processor import HighPerformanceProcessor, get_global_processor

# Import all detector classes
from analysis.hydrogen_bonds import HydrogenBondDetector
from analysis.halogen_bonds import HalogenBondDetector
from analysis.ch_pi_interactions import CHPiDetector
from analysis.pi_pi_stacking import PiPiDetector
from analysis.ionic_interactions import IonicInteractionDetector
from analysis.hydrophobic_contacts import HydrophobicContactDetector
from analysis.chalcogen_bonds import ChalcogenBondDetector
from analysis.pnictogen_bonds import PnictogenBondDetector
from analysis.tetrel_bonds import TetrelBondDetector
from analysis.anion_pi_interactions import AnionPiDetector
from analysis.n_pi_star_interactions import NPiStarDetector
from analysis.london_dispersion import DispersionDetector

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class HighPerformanceBatchProcessor:
    """
    Enhanced batch processor with parallel processing and performance optimizations.
    """
    
    def __init__(self, config=None, use_parallel: bool = True, max_workers: Optional[int] = None):
        """
        Initialize the batch processor.
        
        Args:
            config: Application configuration (optional for backwards compatibility)
            use_parallel: Whether to use parallel processing
            max_workers: Maximum number of worker threads (auto-detected if None)
        """
        self.use_parallel = use_parallel
        self.config = config
        
        # Initialize PDB handler with or without config
        if config:
            self.pdb_handler = PDBHandler(config)
        else:
            # Create a minimal config or use defaults
            try:
                from utils.config import AppConfig
                self.pdb_handler = PDBHandler(AppConfig())
            except:
                # Fallback - try without config
                self.pdb_handler = PDBHandler()
        
        if use_parallel:
            # Initialize high-performance processor with config
            self.processor = HighPerformanceProcessor(max_workers=max_workers, config=config)
        else:
            # Use the global high-performance processor with config
            self.processor = get_global_processor(config=config)
            if max_workers:
                self.processor.max_workers = max_workers        # All available detectors for comprehensive analysis
        self.detector_classes = [
            HydrogenBondDetector,
            HalogenBondDetector,
            CHPiDetector,
            PiPiDetector,
            IonicInteractionDetector,
            HydrophobicContactDetector,
            ChalcogenBondDetector,
            PnictogenBondDetector,
            TetrelBondDetector,
            AnionPiDetector,
            NPiStarDetector,
            DispersionDetector
        ]
        
        logger.info(f"Initialized HighPerformanceBatchProcessor with {len(self.detector_classes)} detectors")
        logger.info(f"Parallel processing: {'Enabled' if use_parallel else 'Disabled'}")
        logger.info(f"Max workers: {self.processor.max_workers}")
    
    def _filter_detectors(self, interaction_filters: Optional[List[str]] = None) -> List:
        """Filter detector classes based on selected interaction types."""
        if not interaction_filters:
            return self.detector_classes
        
        # Map interaction types to detector classes
        detector_mapping = {
            'hydrogenbond': HydrogenBondDetector,
            'halogenbond': HalogenBondDetector,
            'ionicinteraction': IonicInteractionDetector,
            'hydrophobiccontact': HydrophobicContactDetector,
            'pipi': PiPiDetector,
            'chpi': CHPiDetector,
            'chalcogenbond': ChalcogenBondDetector,
            'pnictogenbond': PnictogenBondDetector,
            'tetrelbond': TetrelBondDetector,
            'anionpi': AnionPiDetector,
            'npistar': NPiStarDetector,
            'dispersion': DispersionDetector
        }
        
        # Filter detector classes based on selected interactions
        filtered_detectors = []
        for interaction_type in interaction_filters:
            if interaction_type in detector_mapping:
                detector_class = detector_mapping[interaction_type]
                if detector_class in self.detector_classes:
                    filtered_detectors.append(detector_class)
        
        return filtered_detectors

    def process_single_protein(self, pdb_id: str, interaction_filters: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Process a single protein with selected interaction detectors using high-performance processing.
        
        Args:
            pdb_id: PDB identifier to process
            interaction_filters: List of interaction types to analyze (None for all)
            
        Returns:
            Dictionary containing all detected interactions and metadata
        """
        start_time = time.time()
        
        try:
            logger.info(f"Starting high-performance analysis of {pdb_id}")
            
            # Filter detector classes based on selected interactions
            detector_classes_to_use = self._filter_detectors(interaction_filters)
            
            if not detector_classes_to_use:
                return {
                    'pdb_id': pdb_id,
                    'success': False,
                    'error': 'No valid interaction types selected',
                    'processing_time': time.time() - start_time
                }
            
            logger.info(f"Using {len(detector_classes_to_use)} detector classes for {pdb_id}")
            
            # Load protein structure
            structure = self.pdb_handler.load_structure(pdb_id)
            if not structure:
                return {
                    'pdb_id': pdb_id,
                    'success': False,
                    'error': f'Failed to load structure for {pdb_id}',
                    'processing_time': time.time() - start_time
                }
            
            # Process selected interactions in parallel
            if self.use_parallel:
                interactions = self.processor.process_interactions_parallel(
                    structure, 
                    detector_classes_to_use,
                    self.config
                )
            else:
                interactions = self._process_sequential(structure, detector_classes_to_use)
            
            # Calculate summary statistics
            total_interactions = sum(len(interaction_list) for interaction_list in interactions.values())
            
            processing_time = time.time() - start_time
            
            # Get performance metrics
            perf_summary = self.processor.get_performance_summary()
            
            result = {
                'pdb_id': pdb_id,
                'success': True,
                'interactions': interactions,
                'summary': {
                    'total_interactions': total_interactions,
                    'interaction_counts': {k: len(v) for k, v in interactions.items()},
                    'processing_time': processing_time,
                    'performance_metrics': perf_summary,
                    'analyzed_interactions': interaction_filters or self._get_all_interaction_types()
                },
                'processing_time': processing_time
            }
            
            logger.info(f"Completed {pdb_id}: {total_interactions} total interactions in {processing_time:.2f}s")
            return result
            
        except Exception as e:
            error_msg = f"Error processing {pdb_id}: {str(e)}"
            logger.error(error_msg)
            logger.error(traceback.format_exc())
            
            return {
                'pdb_id': pdb_id,
                'success': False,
                'error': error_msg,
                'processing_time': time.time() - start_time
            }
    
    def process_multiple_proteins(self, pdb_ids: List[str], interaction_filters: Optional[List[str]] = None) -> List[Dict[str, Any]]:
        """
        Process multiple proteins in parallel with optimal batching and selected interactions.
        
        Args:
            pdb_ids: List of PDB identifiers to process
            interaction_filters: List of interaction types to analyze (None for all)
            
        Returns:
            List of analysis results for each protein
        """
        if not pdb_ids:
            return []
        
        logger.info(f"Starting batch processing of {len(pdb_ids)} proteins with {len(interaction_filters) if interaction_filters else 'all'} interaction types")
        start_time = time.time()
        
        results = []
        
        if self.use_parallel and len(pdb_ids) > 1:
            # Use parallel processing for multiple proteins
            max_workers = min(len(pdb_ids), self.processor.max_workers)
            
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                # Submit all tasks with interaction filters
                future_to_pdb = {
                    executor.submit(self.process_single_protein, pdb_id, interaction_filters): pdb_id 
                    for pdb_id in pdb_ids
                }
                
                # Collect results as they complete
                for future in as_completed(future_to_pdb):
                    pdb_id = future_to_pdb[future]
                    try:
                        result = future.result(timeout=300)  # 5 minute timeout per protein
                        results.append(result)
                        logger.info(f"Completed batch processing for {pdb_id}")
                    except Exception as e:
                        error_result = {
                            'pdb_id': pdb_id,
                            'success': False,
                            'error': f'Batch processing timeout or error: {str(e)}',
                            'processing_time': 0
                        }
                        results.append(error_result)
                        logger.error(f"Error in batch processing {pdb_id}: {str(e)}")
        else:
            # Sequential processing for single protein or if parallel is disabled
            for pdb_id in pdb_ids:
                result = self.process_single_protein(pdb_id, interaction_filters)
                results.append(result)
        
        total_time = time.time() - start_time
        successful_results = [r for r in results if r.get('success', False)]
        
        logger.info(f"Batch processing completed: {len(successful_results)}/{len(pdb_ids)} successful in {total_time:.2f}s")
        
        return results
    
    def _process_sequential(self, structure, detector_classes=None) -> Dict[str, List]:
        """Fallback sequential processing method."""
        if detector_classes is None:
            detector_classes = self.detector_classes
            
        interactions = {}
        
        for detector_class in detector_classes:
            detector_name = detector_class.__name__.replace('Detector', '').lower()
            try:
                # Initialize detector with config if available
                if hasattr(self, 'config') and self.config:
                    detector = detector_class(self.config)
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
                            interactions[detector_name] = []
                            continue
                
                result = detector.detect(structure)
                
                if hasattr(result, 'to_dict_list'):
                    interactions[detector_name] = result.to_dict_list()
                else:
                    interactions[detector_name] = result if isinstance(result, list) else []
                    
            except Exception as e:
                logger.error(f"Error in sequential {detector_name}: {str(e)}")
                interactions[detector_name] = []
        
        return interactions
    
    def get_detector_info(self) -> List[Dict[str, str]]:
        """Get information about available detectors."""
        detector_info = []
        
        for detector_class in self.detector_classes:
            info = {
                'name': detector_class.__name__,
                'short_name': detector_class.__name__.replace('Detector', ''),
                'description': getattr(detector_class, '__doc__', 'No description available').strip()
            }
            detector_info.append(info)
        
        return detector_info
    
    def benchmark_performance(self, pdb_id: str = "1A2B", iterations: int = 3) -> Dict[str, Any]:
        """
        Benchmark the performance of the batch processor.
        
        Args:
            pdb_id: PDB ID to use for benchmarking
            iterations: Number of iterations to run
            
        Returns:
            Performance benchmark results
        """
        logger.info(f"Starting performance benchmark with {iterations} iterations")
        
        # Clear cache for fair testing
        self.processor.clear_cache()
        
        times = []
        results_data = []
        
        for i in range(iterations):
            start_time = time.time()
            result = self.process_single_protein(pdb_id)
            end_time = time.time()
            
            iteration_time = end_time - start_time
            times.append(iteration_time)
            results_data.append(result)
            
            logger.info(f"Benchmark iteration {i+1}/{iterations}: {iteration_time:.2f}s")
        
        # Calculate statistics
        avg_time = sum(times) / len(times)
        min_time = min(times)
        max_time = max(times)
        
        # Get system info
        cpu_info = {
            'cpu_count': os.cpu_count(),
            'memory_total': 8,  # Default assumption
            'memory_available': 4  # Default assumption
        }
        
        try:
            import psutil
            cpu_info.update({
                'cpu_percent': psutil.cpu_percent(interval=1),
                'memory_total': psutil.virtual_memory().total // (1024**3),  # GB
                'memory_available': psutil.virtual_memory().available // (1024**3)  # GB
            })
        except ImportError:
            cpu_info['cpu_percent'] = 0
        
        benchmark_result = {
            'benchmark_summary': {
                'iterations': iterations,
                'average_time': f"{avg_time:.2f}s",
                'min_time': f"{min_time:.2f}s",
                'max_time': f"{max_time:.2f}s",
                'throughput': f"{1/avg_time:.2f} proteins/sec",
                'parallel_enabled': self.use_parallel,
                'max_workers': self.processor.max_workers
            },
            'system_info': cpu_info,
            'performance_metrics': self.processor.get_performance_summary(),
            'detailed_results': results_data
        }
        
        logger.info(f"Benchmark completed - Average: {avg_time:.2f}s, Throughput: {1/avg_time:.2f} proteins/sec")
        
        return benchmark_result

# Legacy compatibility class
class BatchProcessor(HighPerformanceBatchProcessor):
    """Legacy compatibility wrapper for the enhanced batch processor."""
    
    def __init__(self):
        super().__init__(use_parallel=True)
        logger.warning("Using legacy BatchProcessor class. Consider upgrading to HighPerformanceBatchProcessor.")
    
    def process_protein(self, pdb_id: str) -> Dict[str, Any]:
        """Legacy method name compatibility."""
        return self.process_single_protein(pdb_id)

import time
import uuid
import concurrent.futures
from typing import List, Dict, Any, Optional, Callable
from pathlib import Path
import multiprocessing as mp
from dataclasses import asdict
from loguru import logger

from Bio.PDB import Structure

from utils.config import AppConfig, get_interaction_types
from utils.pdb_handler import PDBHandler
from analysis.hydrogen_bonds import HydrogenBondDetector
from analysis.halogen_bonds import HalogenBondDetector
from analysis.pi_pi_stacking import PiPiDetector
from analysis.ionic_interactions import IonicInteractionDetector
from analysis.hydrophobic_contacts import HydrophobicContactDetector
from analysis.ch_pi_interactions import CHPiDetector
from analysis.chalcogen_bonds import ChalcogenBondDetector
from analysis.pnictogen_bonds import PnictogenBondDetector
from analysis.tetrel_bonds import TetrelBondDetector
from analysis.anion_pi_interactions import AnionPiDetector
from analysis.n_pi_star_interactions import NPiStarDetector
from analysis.london_dispersion import DispersionDetector

class BatchProcessor:
    """Processes multiple protein structures for interaction analysis."""
    
    def __init__(self, config: AppConfig):
        self.config = config
        self.pdb_handler = PDBHandler(config)
        
        # Initialize all interaction detectors
        self.detectors = self._initialize_detectors()
        
        # Cache for processed structures
        self.structure_cache = {}
    
    def _initialize_detectors(self) -> Dict[str, Any]:
        """Initialize all interaction detection modules."""
        detectors = {}
        
        # Get current parameters (could be from preset or custom)
        params = self.config.interactions
        
        detectors['hydrogen_bond'] = HydrogenBondDetector(self.config)
        
        detectors['halogen_bond'] = HalogenBondDetector(self.config)
        
        detectors['pi_pi'] = PiPiDetector(self.config)
        
        detectors['ionic'] = IonicInteractionDetector(self.config)
        
        detectors['hydrophobic'] = HydrophobicContactDetector(self.config)
        
        detectors['ch_pi'] = CHPiDetector(self.config)
        
        detectors['chalcogen_bond'] = ChalcogenBondDetector(self.config)
        
        detectors['pnictogen_bond'] = PnictogenBondDetector(self.config)
        
        detectors['tetrel_bond'] = TetrelBondDetector(self.config)
        
        detectors['anion_pi'] = AnionPiDetector(self.config)
        
        detectors['n_pi_star'] = NPiStarDetector(self.config)
        
        detectors['dispersion'] = DispersionDetector(self.config)
        
        return detectors
    
    def process_batch(self, 
                     pdb_ids: List[str],
                     preset: str = "literature_default",
                     interaction_types: Optional[List[str]] = None,
                     progress_callback: Optional[Callable] = None) -> Dict[str, Dict[str, Any]]:
        """
        Process a batch of PDB structures.
        
        Args:
            pdb_ids: List of PDB identifiers
            preset: Parameter preset to use
            interaction_types: List of interaction types to analyze (None = all)
            progress_callback: Function to call with progress updates
            
        Returns:
            Dictionary mapping PDB IDs to analysis results
        """
        batch_id = str(uuid.uuid4())
        start_time = time.time()
        
        logger.info(f"Starting batch analysis {batch_id[:8]} with {len(pdb_ids)} structures")
        
        # Apply preset parameters
        self._apply_preset(preset)
        
        # Set default interaction types if not specified
        if interaction_types is None:
            interaction_types = get_interaction_types()
        
        results = {}
        failed_structures = []
        
        # Process structures based on configuration
        if self.config.processing.max_workers == 1:
            # Sequential processing
            for i, pdb_id in enumerate(pdb_ids):
                try:
                    if progress_callback:
                        progress_callback(i, len(pdb_ids), pdb_id)
                    
                    result = self._process_single_structure(pdb_id, interaction_types, preset)
                    if result:
                        results[pdb_id] = result
                    else:
                        failed_structures.append(pdb_id)
                        
                except Exception as e:
                    logger.error(f"Failed to process {pdb_id}: {e}")
                    failed_structures.append(pdb_id)
        
        else:
            # Parallel processing
            results, failed_structures = self._process_parallel(
                pdb_ids, interaction_types, preset, progress_callback
            )
        
        # Final progress update
        if progress_callback:
            progress_callback(len(pdb_ids), len(pdb_ids), "Complete")
        
        total_time = time.time() - start_time
        
        logger.info(f"Batch analysis complete: {len(results)} successful, {len(failed_structures)} failed, {total_time:.2f}s")
        
        # Cache batch results
        if hasattr(self, 'cache_manager'):
            batch_summary = {
                'results': results,
                'failed_structures': failed_structures,
                'batch_id': batch_id,
                'processing_time': total_time,
                'preset': preset,
                'interaction_types': interaction_types
            }
            self.cache_manager.cache_batch_results(batch_id, batch_summary)
        
        return results
    
    def _process_parallel(self, 
                         pdb_ids: List[str],
                         interaction_types: List[str],
                         preset: str,
                         progress_callback: Optional[Callable] = None) -> tuple:
        """Process structures in parallel."""
        results = {}
        failed_structures = []
        completed_count = 0
        
        # Create chunks for processing
        chunk_size = min(self.config.processing.chunk_size, len(pdb_ids))
        chunks = [pdb_ids[i:i+chunk_size] for i in range(0, len(pdb_ids), chunk_size)]
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=self.config.processing.max_workers) as executor:
            # Submit chunks for processing
            future_to_chunk = {
                executor.submit(
                    self._process_chunk, 
                    chunk, 
                    interaction_types, 
                    preset
                ): chunk for chunk in chunks
            }
            
            # Collect results as they complete
            for future in concurrent.futures.as_completed(future_to_chunk):
                chunk = future_to_chunk[future]
                
                try:
                    chunk_results, chunk_failed = future.result(timeout=self.config.processing.timeout_seconds)
                    results.update(chunk_results)
                    failed_structures.extend(chunk_failed)
                    
                    completed_count += len(chunk)
                    
                    if progress_callback:
                        progress_callback(completed_count, len(pdb_ids), f"Processed {len(chunk)} structures")
                        
                except Exception as e:
                    logger.error(f"Chunk processing failed: {e}")
                    failed_structures.extend(chunk)
        
        return results, failed_structures
    
    def _process_chunk(self, 
                      pdb_ids: List[str],
                      interaction_types: List[str],
                      preset: str) -> tuple:
        """Process a chunk of PDB IDs (for parallel processing)."""
        # This method runs in a separate process
        chunk_results = {}
        chunk_failed = []
        
        # Reinitialize components for this process
        local_processor = BatchProcessor(self.config)
        local_processor._apply_preset(preset)
        
        for pdb_id in pdb_ids:
            try:
                result = local_processor._process_single_structure(pdb_id, interaction_types, preset)
                if result:
                    chunk_results[pdb_id] = result
                else:
                    chunk_failed.append(pdb_id)
            except Exception as e:
                logger.error(f"Failed to process {pdb_id} in chunk: {e}")
                chunk_failed.append(pdb_id)
        
        return chunk_results, chunk_failed
    
    def analyze_single_structure(self, 
                                pdb_id: str,
                                structure: Structure.Structure,
                                preset: str = "literature_default",
                                interaction_types: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Analyze a single structure that's already loaded.
        
        Args:
            pdb_id: PDB identifier
            structure: Loaded Biopython Structure
            preset: Parameter preset
            interaction_types: List of interaction types to analyze
            
        Returns:
            Analysis results dictionary
        """
        start_time = time.time()
        
        # Apply preset parameters
        self._apply_preset(preset)
        
        # Set default interaction types
        if interaction_types is None:
            interaction_types = get_interaction_types()
        
        # Run analysis
        interactions = {}
        metadata = {
            'pdb_id': pdb_id,
            'analysis_time': 0,
            'preset': preset,
            'interaction_types_analyzed': interaction_types,
            'total_interactions': 0
        }
        
        # Detect each interaction type
        for interaction_type in interaction_types:
            if interaction_type in self.detectors:
                try:
                    detector = self.detectors[interaction_type]
                    detected_interactions = self._run_detector(detector, structure, interaction_type)
                    interactions[interaction_type] = detected_interactions
                    metadata['total_interactions'] += len(detected_interactions)
                    
                except Exception as e:
                    logger.error(f"Failed to detect {interaction_type} in {pdb_id}: {e}")
                    interactions[interaction_type] = []
        
        # Calculate analysis time
        metadata['analysis_time'] = time.time() - start_time
        
        # Structure validation
        validation = self.pdb_handler.validate_structure(structure)
        
        # Get structure information
        structure_info = self.pdb_handler._extract_structure_info(structure)
        
        # Compile results
        result = {
            'interactions': interactions,
            'metadata': metadata,
            'structure_info': structure_info,
            'validation': validation,
            'hotspots': self._identify_hotspots(interactions),
            'statistics': self._calculate_statistics(interactions)
        }
        
        return result
    
    def _process_single_structure(self, 
                                 pdb_id: str,
                                 interaction_types: List[str],
                                 preset: str) -> Optional[Dict[str, Any]]:
        """Process a single PDB structure."""
        try:
            # Load structure
            structure = self.pdb_handler.load_structure(pdb_id, assembly=self.config.default_assembly)
            if not structure:
                return None
            
            return self.analyze_single_structure(pdb_id, structure, preset, interaction_types)
            
        except Exception as e:
            logger.error(f"Error processing {pdb_id}: {e}")
            return None
    
    def _run_detector(self, detector, structure: Structure.Structure, interaction_type: str) -> List[Dict[str, Any]]:
        """Run a specific interaction detector."""
        if interaction_type == 'hydrogen_bond':
            bonds = detector.detect_hydrogen_bonds(structure)
            return detector.to_dict_list(bonds)
        
        elif interaction_type == 'halogen_bond':
            bonds = detector.detect_halogen_bonds(structure)
            return detector.to_dict_list(bonds)
        
        elif interaction_type == 'pi_pi':
            interactions = detector.detect_pi_pi_interactions(structure)
            return detector.to_dict_list(interactions)
        
        elif interaction_type == 'ionic':
            interactions = detector.detect_ionic_interactions(structure)
            return detector.to_dict_list(interactions)
        
        elif interaction_type == 'hydrophobic':
            contacts = detector.detect_hydrophobic_contacts(structure)
            return detector.to_dict_list(contacts)
        
        elif interaction_type == 'ch_pi':
            interactions = detector.detect_ch_pi_interactions(structure)
            return interactions  # Already returns dict list
        
        elif interaction_type == 'chalcogen_bond':
            bonds = detector.detect_chalcogen_bonds(structure)
            return detector.to_dict_list(bonds)
        
        elif interaction_type == 'pnictogen_bond':
            bonds = detector.detect_pnictogen_bonds(structure)
            return detector.to_dict_list(bonds)
        
        elif interaction_type == 'tetrel_bond':
            bonds = detector.detect_tetrel_bonds(structure)
            return detector.to_dict_list(bonds)
        
        elif interaction_type == 'anion_pi':
            interactions = detector.detect_anion_pi_interactions(structure)
            return detector.to_dict_list(interactions)
        
        elif interaction_type == 'n_pi_star':
            interactions = detector.detect_n_pi_star_interactions(structure)
            return detector.to_dict_list(interactions)
        
        elif interaction_type == 'dispersion':
            interactions = detector.detect_dispersion_interactions(structure)
            return detector.to_dict_list(interactions)
        
        # Should not reach here as all detectors are now implemented
        else:
            logger.warning(f"Unknown interaction type: {interaction_type}")
            return []
    
    def _apply_preset(self, preset: str):
        """Apply parameter preset to detectors."""
        if preset not in self.config.presets:
            logger.warning(f"Unknown preset {preset}, using literature_default")
            preset = "literature_default"
        
        preset_params = self.config.presets[preset]
        
        # Update interaction configuration
        for param_name, value in preset_params.items():
            if hasattr(self.config.interactions, param_name):
                setattr(self.config.interactions, param_name, value)
        
        # Reinitialize detectors with new parameters
        self.detectors = self._initialize_detectors()
    
    def _identify_hotspots(self, interactions: Dict[str, List[Dict]]) -> List[Dict[str, Any]]:
        """Identify interaction hotspots (residues with many interactions)."""
        residue_counts = {}
        
        # Count interactions per residue
        for interaction_type, interaction_list in interactions.items():
            for interaction in interaction_list:
                res1 = f"{interaction.get('chain1', '')}{interaction.get('residue1', '')}"
                res2 = f"{interaction.get('chain2', '')}{interaction.get('residue2', '')}"
                
                residue_counts[res1] = residue_counts.get(res1, 0) + 1
                residue_counts[res2] = residue_counts.get(res2, 0) + 1
        
        # Sort by interaction count and identify hotspots
        sorted_residues = sorted(residue_counts.items(), key=lambda x: x[1], reverse=True)
        
        hotspots = []
        for residue, count in sorted_residues[:20]:  # Top 20 hotspots
            if count >= 3:  # Minimum threshold
                hotspots.append({
                    'residue': residue,
                    'interaction_count': count,
                    'normalized_score': count / max(residue_counts.values()) if residue_counts else 0
                })
        
        return hotspots
    
    def _calculate_statistics(self, interactions: Dict[str, List[Dict]]) -> Dict[str, Any]:
        """Calculate summary statistics for all interactions."""
        stats = {
            'total_interactions': sum(len(interaction_list) for interaction_list in interactions.values()),
            'interaction_counts': {},
            'chain_interactions': {},
            'average_distances': {}
        }
        
        # Count by interaction type
        for interaction_type, interaction_list in interactions.items():
            stats['interaction_counts'][interaction_type] = len(interaction_list)
            
            # Calculate average distances if available
            distances = [i.get('distance', 0) for i in interaction_list if 'distance' in i]
            if distances:
                stats['average_distances'][interaction_type] = sum(distances) / len(distances)
        
        return stats
    
    def get_batch_statistics(self, batch_results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """Calculate statistics across a batch of structures."""
        if not batch_results:
            return {}
        
        batch_stats = {
            'total_structures': len(batch_results),
            'successful_structures': len([r for r in batch_results.values() if r]),
            'total_interactions': 0,
            'interaction_type_counts': {},
            'average_interactions_per_structure': 0,
            'hotspot_residues': {},
            'processing_times': []
        }
        
        # Aggregate statistics
        for pdb_id, result in batch_results.items():
            if not result:
                continue
                
            batch_stats['total_interactions'] += result['statistics']['total_interactions']
            batch_stats['processing_times'].append(result['metadata']['analysis_time'])
            
            # Count interaction types
            for int_type, count in result['statistics']['interaction_counts'].items():
                batch_stats['interaction_type_counts'][int_type] = \
                    batch_stats['interaction_type_counts'].get(int_type, 0) + count
            
            # Aggregate hotspots
            for hotspot in result.get('hotspots', []):
                residue = hotspot['residue']
                if residue not in batch_stats['hotspot_residues']:
                    batch_stats['hotspot_residues'][residue] = {
                        'structures': [],
                        'total_interactions': 0
                    }
                
                batch_stats['hotspot_residues'][residue]['structures'].append(pdb_id)
                batch_stats['hotspot_residues'][residue]['total_interactions'] += hotspot['interaction_count']
        
        # Calculate averages
        if batch_stats['successful_structures'] > 0:
            batch_stats['average_interactions_per_structure'] = \
                batch_stats['total_interactions'] / batch_stats['successful_structures']
        
        if batch_stats['processing_times']:
            import numpy as np
            batch_stats['average_processing_time'] = np.mean(batch_stats['processing_times'])
            batch_stats['total_processing_time'] = sum(batch_stats['processing_times'])
        
        return batch_stats

    def _filter_detectors(self, interaction_filters: Optional[List[str]]) -> List:
        """Filter detector classes based on selected interaction types."""
        if not interaction_filters:
            return self.detector_classes
        
        # Map interaction type names to detector classes
        detector_mapping = {
            'hydrogenbond': HydrogenBondDetector,
            'halogenbond': HalogenBondDetector,
            'chpi': CHPiDetector,
            'pipi': PiPiDetector,
            'ionicinteraction': IonicInteractionDetector,
            'hydrophobiccontact': HydrophobicContactDetector,
            'chalcogenbond': ChalcogenBondDetector,
            'pnictogenbond': PnictogenBondDetector,
            'tetrelbond': TetrelBondDetector,
            'anionpi': AnionPiDetector,
            'npistar': NPiStarDetector,
            'dispersion': DispersionDetector
        }
        
        filtered_detectors = []
        for interaction_type in interaction_filters:
            if interaction_type in detector_mapping:
                filtered_detectors.append(detector_mapping[interaction_type])
            else:
                logger.warning(f"Unknown interaction type: {interaction_type}")
        
        return filtered_detectors
    
    def _get_all_interaction_types(self) -> List[str]:
        """Get list of all supported interaction types."""
        return [
            'hydrogenbond', 'halogenbond', 'chpi', 'pipi', 'ionicinteraction',
            'hydrophobiccontact', 'chalcogenbond', 'pnictogenbond', 'tetrelbond',
            'anionpi', 'npistar', 'dispersion'
        ]
