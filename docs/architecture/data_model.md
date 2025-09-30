# 2. Data Model & Schema Layer

## 2.1 Atom-Level Representation
Each atom record (conceptual) comprises:
| Field | Type | Purpose | Notes |
|-------|------|---------|-------|
| chain_id | str (compressed later) | Biological chain grouping | Encoded to small int for memory |
| resnum | int | Residue sequence number | Combined with insertion code when present |
| icode | str? | Insertion code | Optional (defaults blank) |
| resname | str | Residue three-letter code | Normalization ensures uppercase |
| atom_name | str | Atom name (e.g., CA, ND1) | Fixed-length ASCII subset |
| element | str | Element symbol | Derived from atom/pdb heuristics |
| x,y,z | float | Cartesian coordinates Å | Stored collectively as float32 arrays |
| occupancy | float | Crystallographic occupancy | First altLoc chosen if disambiguated |
| bfactor | float | Temperature factor | Retained for potential quality filters |
| flags | bitmask | Role classification | Derived (donor, acceptor, ring, charged, metal) |

## 2.2 Role Bitmask Encoding
16-bit layout (tentative):
```
bit 0  DONOR
bit 1  ACCEPTOR
bit 2  AROM_RING_ATOM
bit 3  CHARGED_POS
bit 4  CHARGED_NEG
bit 5  METAL_CENTER
bit 6  SULFUR_CONTAINING
bit 7  HALOGEN_CONTAINING
bit 8  POLAR_NEUTRAL
bit 9  HYDROPHOBIC
bits 10–15  Reserved
```
Rationale: constant-time role checks via bitwise operations.

## 2.3 Feature Store Objects
| Artifact | Shape / Type | Derived From | Used By |
|----------|--------------|-------------|---------|
| coords | (N,3) float32 | Atoms | All detectors |
| atom_flags | (N,) uint16 | Atoms | Fast role filtering |
| atom_index_by_name | dict[str, ndarray[int]] | Atoms | Domain subset retrieval |
| ring_centroids | (R,3) float32 | Aromatic residues | π–π, cation–π, CH–π |
| ring_normals | (R,3) float32 | Same | π alignment & offset |
| ring_members | list[list[int]] | Atom indices | Visual overlays, advanced filters |
| charged_pos_centroids | (P,3) float32 | ARG, LYS, etc. | Cation–π, ionic |
| charged_neg_centroids | (M,3) float32 | ASP, GLU, etc. | Ionic, salt bridge |
| donor_atoms | ndarray[int] | Flags | Hydrogen bonds |
| acceptor_atoms | ndarray[int] | Flags | Hydrogen bonds |
| sulfur_atoms | ndarray[int] | Element filter | Chalcogen/pnictogen/tetrel |
| halogen_atoms | ndarray[int] | Element filter | Halogen bonds |

## 2.4 Normalized Interaction Schema
A normalized record emphasizes structural neutrality and comparability across types.
```
{
  "interaction_type": str,
  "participants": [ {"chain":str,"resnum":int,"atom":str,"role":str}, ...],
  "metrics": { str: float, ... },
  "subtype": str?,
  "score": float?,
  "provenance_hash": str?,
  "flags": { "experimental": bool?, ... }
}
```

## 2.5 Canonical vs Normalized Tradeoffs
| Aspect | Canonical | Normalized |
|--------|----------|------------|
| Speed to produce | Faster | Slight overhead |
| Readability | Type-specific | Uniform |
| Downstream analytics | Requires join logic | Direct aggregation |
| Storage size | Smaller | Slightly larger (extra keys) |

## 2.6 Provenance Hash Inputs
Concatenated canonical form (stable ordering) of: detector list, parameter values, structure signature (chain/residue hash), version tags. SHA256 → hex digest truncated to 16–24 chars for display.

## 2.7 Integrity Invariants
| Invariant | Enforcement |
|----------|-------------|
| All participant indices valid bounds | Assertion / defensive checks |
| Distance metrics ≥ 0 | Constructed via norms |
| Angles within [0, 180] | Clamped arccos inputs |
| No duplicate participant pair ordering | Sorting + set membership or i<j constraint |
| Normalization stable ordering | Sort participants lexicographically |

## 2.8 Serialization Formats
| Format | Use Case | Notes |
|--------|----------|------|
| JSON lines | Streaming metrics | Append-only friendly |
| JSON (array) | Batch export | Human-readable |
| Columnar (future) | Large-scale analysis | Potential Arrow/Parquet integration |
| CSV | Simple spreadsheets | Flattened metrics only |

## 2.9 Potential Future Expansions
- Residue environment fingerprints (hydropathy radius profiles).
- Ring aromaticity scoring metadata.
- Confidence / quality flags (occupancy/bfactor threshold).
- Graph adjacency matrix export for network analysis.

Proceed to `detectors.md` for per-interaction algorithm specifics.
