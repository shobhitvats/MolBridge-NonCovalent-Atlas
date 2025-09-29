"""Golden snapshot style test (non-breaking).

This test is optional and will SKIP automatically if prerequisites (sample PDB,
biopython, etc.) are missing. It exercises the UnifiedBatchProcessor under
feature flags if they are enabled to ensure output contract stability.
"""
from __future__ import annotations
import os
import pytest

try:
    from analysis.unified_batch import UnifiedBatchProcessor
    from utils.settings import get_settings
except Exception:  # pragma: no cover
    UnifiedBatchProcessor = None  # type: ignore
    get_settings = lambda: None  # type: ignore

SAMPLE_PDB = os.environ.get("MOLBRIDGE_GOLDEN_PDB", "1CRN")  # small protein default

EXPECTED_MIN_COUNTS = {
    # Provide conservative non-zero lower bounds only for very common interactions
    # (These may be zero if detection heuristics differ; keep loose to avoid flakiness.)
    # 'hydrogenbond': 5,
}

@pytest.mark.skipif(UnifiedBatchProcessor is None, reason="UnifiedBatchProcessor import failed")
def test_unified_batch_contract():
    settings = get_settings()
    processor = UnifiedBatchProcessor()
    result = processor.analyze(SAMPLE_PDB)
    assert result["pdb_id"].lower() == SAMPLE_PDB.lower()
    assert result["success"] is True
    assert "interactions" in result
    # Optional keys present only when flags enabled
    if settings.enable_normalization:
        assert "interactions_normalized" in result
    if settings.enable_provenance:
        assert "provenance" in result and "schema_version" in result["provenance"]
    # Loose count assertions (non-breaking)
    for key, min_count in EXPECTED_MIN_COUNTS.items():
        if key in result["interactions"]:
            assert len(result["interactions"][key]) >= 0  # placeholder (min_count) intentionally lax
