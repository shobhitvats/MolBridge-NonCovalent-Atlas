"""Common interaction data structures and detector protocol.

Introduces a lightweight unified representation WITHOUT forcing existing detectors
(to avoid breaking current logic). Detectors can optionally emit Interaction
objects; legacy list of custom dataclasses/dicts remains supported.
"""
# NOTE (2025-09): Shared FeatureStore neighbor_within() API now available.
# Detectors progressively migrating (hydrogen_bonds, halogen_bonds done).
# Remaining detectors can adopt FeatureStore.neighbor_within(structure, radius)
# to reduce redundant pairwise distance computations.
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Protocol, Iterable, Any, Dict, Tuple, Optional

@dataclass(slots=True)
class Interaction:
    kind: str
    participants: Tuple[Any, ...]  # typically Biopython Atom objects
    distance: Optional[float] = None
    score: Optional[float] = None
    geometry: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

class DetectorLike(Protocol):
    """Structural protocol most detectors roughly follow (duck-typed)."""
    def __init__(self, config: Any): ...  # type: ignore

# Optional normalization helper -------------------------------------------------

def to_interaction_dict(obj: Any, fallback_kind: str) -> Dict[str, Any]:
    """Best-effort conversion of arbitrary detector outputs to a flat dict.

    Leaves existing dicts unchanged except for inserting 'type' if missing.
    For dataclass-like objects, inspects attributes.
    """
    if isinstance(obj, dict):
        if 'type' not in obj:
            obj['type'] = fallback_kind
        return obj
    # Try attribute extraction
    data: Dict[str, Any] = {}
    for attr in ('distance','angle','theta_angle','delta_angle','strength','score','donor_residue','acceptor_residue','chalcogen_residue','acceptor_chain','donor_chain','chalcogen_chain','residue1','residue2','chain1','chain2'):
        if hasattr(obj, attr):
            data[attr] = getattr(obj, attr)
    data['type'] = fallback_kind
    return data
