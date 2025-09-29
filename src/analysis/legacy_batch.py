"""Legacy batch processor shim.

Provides a thin wrapper maintaining the original ``BatchProcessor`` public
interface expected elsewhere. Internally delegates to
``HighPerformanceBatchProcessor`` for consistency.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional
from .high_performance_batch import HighPerformanceBatchProcessor
from utils.config import AppConfig
import logging

logger = logging.getLogger(__name__)


class BatchProcessor(HighPerformanceBatchProcessor):  # pragma: no cover - shim
    def __init__(self, config: AppConfig | None = None):
        logger.warning("BatchProcessor is deprecated; use HighPerformanceBatchProcessor instead.")
        super().__init__(config=config or AppConfig(), use_parallel=True)

    # Backwards alias
    def process_protein(self, pdb_id: str) -> Dict[str, Any]:
        return self.process_single_protein(pdb_id)
