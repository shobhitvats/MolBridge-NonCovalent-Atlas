"""Tests for compact serialized blob caching & structure coordinate hash utility."""
from __future__ import annotations

from utils.structure_hash import compute_structure_coord_hash
from utils.cache import CacheManager
from reporting.serializer import export_json
import types


class _DummyAtom:
    def __init__(self, coord):
        self._c = coord
    def get_coord(self):  # pragma: no cover - simple accessor
        return self._c

class _DummyStructure:
    def __init__(self, coords):
        self._atoms = [_DummyAtom(c) for c in coords]
    def get_atoms(self):  # yields atoms
        for a in self._atoms:
            yield a


def test_structure_hash_and_compact_blob(tmp_path):
    # Two simple coordinate sets (3 atoms)
    coords = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
    ]
    s = _DummyStructure(coords)
    h = compute_structure_coord_hash(s)
    assert len(h) == 16

    # Prepare canonical interactions dict (minimal)
    canonical = {"hydrogenbond": [{"a": 1, "b": 2}]}
    blob = export_json(canonical, cache_key=f"test:{h}")
    assert isinstance(blob, str) and "hydrogenbond" in blob

    cm = CacheManager(tmp_path)
    cm.cache_compact_blob("TEST", "phash", h, blob)
    retrieved = cm.get_compact_blob("TEST", "phash", h)
    assert retrieved is not None
    assert b"hydrogenbond" in retrieved
