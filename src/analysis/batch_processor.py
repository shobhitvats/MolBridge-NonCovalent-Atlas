"""Backwards compatibility shim for refactored batch processing modules.

Re-exports:
    - HighPerformanceBatchProcessor
    - BatchProcessor (legacy shim)
"""

from .high_performance_batch import HighPerformanceBatchProcessor  # noqa: F401
from .legacy_batch import BatchProcessor  # noqa: F401

__all__ = ["HighPerformanceBatchProcessor", "BatchProcessor"]
