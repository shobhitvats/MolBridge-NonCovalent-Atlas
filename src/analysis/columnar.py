"""Columnar interaction storage and adapters.

Provides a lightweight Structure-of-Arrays (SoA) representation for interactions
to reduce Python object overhead. Activated when settings.enable_columnar=1.

Design goals:
  * Immutable-ish after finalize() so it can be safely cached / reused.
  * Minimal per-field lists; converted to numpy float32/int32 where feasible.
  * Adapter to list-of-dicts only at serialization time.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Any, Iterable
import numpy as np

@dataclass
class ColumnarBlock:
    type_name: str
    # Core scalar/string fields kept as python lists for flexibility; numeric -> numpy
    residue1: List[str] = field(default_factory=list)
    residue2: List[str] = field(default_factory=list)
    chain1: List[str] = field(default_factory=list)
    chain2: List[str] = field(default_factory=list)
    subtype: List[str] = field(default_factory=list)
    # Numeric arrays captured as python list then compacted
    distance: List[float] = field(default_factory=list)
    angle: List[float] = field(default_factory=list)
    score: List[float] = field(default_factory=list)
    _finalized: bool = False
    # Numpy materializations
    _distance_arr: Any | None = None
    _angle_arr: Any | None = None
    _score_arr: Any | None = None

    # Internal preallocated buffers (optional growth strategy)
    _cap: int = 0
    _len: int = 0
    _distance_buf: Any | None = None
    _angle_buf: Any | None = None
    _score_buf: Any | None = None

    def _ensure_capacity(self):
        if self._finalized:
            return
        if self._len < self._cap:
            return
        import numpy as _np
        new_cap = 32 if self._cap == 0 else min(self._cap * 2, self._cap + 65536)
        # Grow python lists for categorical fields (append cost amortized; no need prealloc) only track cap for numeric arrays
        if self._distance_buf is None:
            self._distance_buf = _np.empty(new_cap, dtype='float32')
            self._angle_buf = _np.empty(new_cap, dtype='float32')
            self._score_buf = _np.empty(new_cap, dtype='float32')
        else:
            self._distance_buf = _np.resize(self._distance_buf, new_cap)
            self._angle_buf = _np.resize(self._angle_buf, new_cap)
            self._score_buf = _np.resize(self._score_buf, new_cap)
        self._cap = new_cap

    def add(self, item: Dict[str, Any]):
        if self._finalized:
            raise RuntimeError("Block already finalized")
        self._ensure_capacity()
        self.residue1.append(item.get('residue1',''))
        self.residue2.append(item.get('residue2',''))
        self.chain1.append(item.get('chain1',''))
        self.chain2.append(item.get('chain2',''))
        self.subtype.append(item.get('subtype',''))
        dist = float(item.get('distance', 0.0))
        ang = float(item.get('angle', item.get('theta_angle', 0.0)))
        sco = float(item.get('score', 0.0))
        if self._distance_buf is not None:
            self._distance_buf[self._len] = dist
            self._angle_buf[self._len] = ang
            self._score_buf[self._len] = sco
        else:  # fallback to list accumulation (first few adds before ensure_capacity sets buffers)
            self.distance.append(dist)
            self.angle.append(ang)
            self.score.append(sco)
        self._len += 1

    def finalize(self):
        if self._finalized:
            return
        if self._distance_buf is not None:
            # Slice to logical length
            self._distance_arr = self._distance_buf[:self._len].copy()
            self._angle_arr = self._angle_buf[:self._len].copy()
            self._score_arr = self._score_buf[:self._len].copy()
            # Drop temp buffers
            self._distance_buf = self._angle_buf = self._score_buf = None
            self.distance = []  # type: ignore
            self.angle = []     # type: ignore
            self.score = []     # type: ignore
        else:
            self._distance_arr = np.asarray(self.distance, dtype='float32')
            self._angle_arr = np.asarray(self.angle, dtype='float32')
            self._score_arr = np.asarray(self.score, dtype='float32')
            self.distance = []  # type: ignore
            self.angle = []     # type: ignore
            self.score = []     # type: ignore
        self._finalized = True

    def to_dict_list(self) -> List[Dict[str, Any]]:
        # Reconstruct list-of-dicts (used at serialization boundary)
        n = len(self.residue1)
        if self._finalized:
            dist = self._distance_arr.tolist() if self._distance_arr is not None else []
            ang = self._angle_arr.tolist() if self._angle_arr is not None else []
            sco = self._score_arr.tolist() if self._score_arr is not None else []
        else:
            if self._distance_buf is not None:
                dist = self._distance_buf[:self._len].tolist()
                ang = self._angle_buf[:self._len].tolist()
                sco = self._score_buf[:self._len].tolist()
            else:
                dist = self.distance
                ang = self.angle
                sco = self.score
        out = []
        for i in range(n):
            out.append({
                'type': self.type_name,
                'residue1': self.residue1[i],
                'residue2': self.residue2[i],
                'chain1': self.chain1[i],
                'chain2': self.chain2[i],
                'subtype': self.subtype[i],
                'distance': dist[i] if i < len(dist) else 0.0,
                'angle': ang[i] if i < len(ang) else 0.0,
                'score': sco[i] if i < len(sco) else 0.0,
            })
        return out


class ColumnarStore:
    def __init__(self):
        self.blocks: Dict[str, ColumnarBlock] = {}
        self.finalized = False
        # Reusable buffer pools keyed by interaction type -> (distance, angle, score arrays, cap, len)
        self._buffer_pool: Dict[str, Dict[str, Any]] = {}

    def _acquire_block(self, type_name: str) -> ColumnarBlock:
        blk = ColumnarBlock(type_name)
        pool_entry = self._buffer_pool.get(type_name)
        if pool_entry:
            # Reuse numpy arrays by assigning as buffers (will be resized if needed)
            blk._distance_buf = pool_entry.get('distance')
            blk._angle_buf = pool_entry.get('angle')
            blk._score_buf = pool_entry.get('score')
            if blk._distance_buf is not None:
                blk._cap = blk._distance_buf.shape[0]
        return blk

    def _release_block(self, blk: ColumnarBlock):
        if not blk._finalized or blk._distance_arr is None:
            return
        # Keep a modest cap for reuse (avoid unbounded growth): shrink if oversized
        import numpy as _np
        max_keep = 65536
        if blk._distance_arr.shape[0] > max_keep:
            return
        # Store fresh empty buffers sized to next power-of-two >= current length for amortization
        n = blk._distance_arr.shape[0]
        cap = 1
        while cap < max(32, n):
            cap <<= 1
        self._buffer_pool[blk.type_name] = {
            'distance': _np.empty(cap, dtype='float32'),
            'angle': _np.empty(cap, dtype='float32'),
            'score': _np.empty(cap, dtype='float32')
        }

    def add_interaction(self, type_name: str, item: Dict[str, Any]):
        if type_name not in self.blocks:
            self.blocks[type_name] = self._acquire_block(type_name)
        self.blocks[type_name].add(item)

    def finalize(self):
        if self.finalized:
            return
        for blk in self.blocks.values():
            blk.finalize()
            # Return empty reusable buffers for this type for next invocation
            self._release_block(blk)
        self.finalized = True

    def to_interactions_map(self) -> Dict[str, List[Dict[str, Any]]]:
        return {k: blk.to_dict_list() for k, blk in self.blocks.items()}

    def to_columnar_payload(self) -> Dict[str, Dict[str, Any]]:
        """Return a columnar JSON-serializable payload without reconstructing list-of-dicts.

        Structure:
          {
            interaction_type: {
               'residue1': [...], 'residue2': [...], 'chain1': [...], 'chain2': [...],
               'subtype': [...], 'distance': <list|None>, 'angle': <list|None>, 'score': <list|None>
            },
            ...
          }
        Numeric numpy arrays converted to python lists (float32 → short lists) only once.
        This path is intended for performance_mode direct serialization.
        """
        payload: Dict[str, Dict[str, Any]] = {}
        for name, blk in self.blocks.items():
            if blk._finalized:
                dist = blk._distance_arr.tolist() if blk._distance_arr is not None else []
                ang = blk._angle_arr.tolist() if blk._angle_arr is not None else []
                sco = blk._score_arr.tolist() if blk._score_arr is not None else []
            else:
                dist = blk.distance
                ang = blk.angle
                sco = blk.score
            payload[name] = {
                'residue1': blk.residue1,
                'residue2': blk.residue2,
                'chain1': blk.chain1,
                'chain2': blk.chain2,
                'subtype': blk.subtype,
                'distance': dist,
                'angle': ang,
                'score': sco,
            }
        return payload

_global_columnar_store: ColumnarStore | None = None

def get_columnar_store() -> ColumnarStore:
    global _global_columnar_store
    if _global_columnar_store is None:
        _global_columnar_store = ColumnarStore()
    return _global_columnar_store

__all__ = ['ColumnarStore','ColumnarBlock','get_columnar_store']
