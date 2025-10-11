# 4. Feature Extraction Pipeline

## 4.1 Stages Overview
| Stage | Input | Output | Purpose |
|-------|-------|--------|---------|
| Parse | PDB/CIF text | Atom objects | Raw atomic inventory |
| Sanitize | Atom objects | Filtered atom list | Remove altLoc duplicates, water (option configurable) |
| Classify | Filtered atoms | Flags bitmask | Role tagging |
| Group Derivations | Filtered atoms | Ring sets, charge groups | Higher-level semantics |
| Geometric Derivations | Ring sets, charge groups | Centroids, normals | Downstream geometry |
| Indexing | Atoms | Name→indices maps | Fast subset retrieval |
| Spatial Index | Coordinates | KD-tree | Neighbor queries |
| Packaging | All above | FeatureStore | Single access facade |

## 4.2 Parsing Details
- Primary: Biopython `PDBParser` or `MMCIFParser` depending on file extension.
- Error Tolerance: Non-fatal issues (unknown residue) logged; residue skipped if lacking coordinates.

## 4.3 Sanitization Rules
| Rule | Rationale |
|------|-----------|
| Keep first altLoc with highest occupancy | Avoid coordinate ambiguity |
| Drop hetero waters (HOH) if not requested | Performance gain, irrelevant to many detectors |
| Standardize residue & atom names (uppercase) | Stable hashing / provenance |

## 4.4 Role Classification Heuristics
| Role | Detection Logic |
|------|-----------------|
| Donor | Residue templates + presence of H or potential H site |
| Acceptor | Electronegative (O, N, S) with lone pair potential | 
| Aromatic Ring Atom | Membership in phenyl/imidazole/indole ring set | 
| Charged Positive | Side chain pattern (ARG, LYS, protonated HIS) | 
| Charged Negative | ASP, GLU side chain carboxylate | 
| Metal | Element in metal whitelist (Zn, Fe, Mg, etc.) | 

## 4.5 Aromatic Ring Derivation
1. Identify residues with known aromatic scaffolds.  
2. Extract ring atom coordinates (e.g., six for phenyl).  
3. Compute centroid = mean(atom coords).  
4. Normal vector = normalized cross product of two non-colinear edge vectors (averaged for stability).  
5. Store: centroid, normal, membership indices.

## 4.6 Charged Group Centroids
Group-specific atom selections (e.g., ARG: CZ, NH1, NH2).  
Centroid = average.  
Optional future weighting: partial charges from library.

## 4.7 KD-tree Construction
```
from scipy.spatial import cKDTree
kdtree = cKDTree(coords)  # coords is (N,3) float32
```
Reused across detectors; radius queries use pre-expanded maximum needed across all active detectors.

## 4.8 Index Maps
Examples: `atom_index_by_name['CA']` returns ndarray of indices of alpha carbons.  
Avoids repeated scanning.

## 4.9 FeatureStore Packaging Contract
Pseudo:
```
class FeatureStore:
    coords: ndarray
    atom_flags: ndarray
    ring_centroids: ndarray
    ring_normals: ndarray
    ... (accessor helpers)
```
Focus: lightweight container, no heavy logic; pure data + helper methods.

## 4.10 Validation Checks
| Check | Action on Fail |
|-------|----------------|
| Missing ring atoms count mismatch | Exclude ring; log warning |
| Zero-length normal (degenerate) | Recompute using alternate pair; else drop |
| Empty donor/acceptor sets | Detector skip pre-flag |

## 4.11 Performance Considerations
- Vectorized centroid & normal computation reduces Python overhead.
- Deduplicate repeated lookups (e.g., ring membership) with cached properties.
- Pre-filter atoms by element before bitmask application when large (branch prediction aid).

## 4.12 Future Enhancements
| Enhancement | Benefit |
|------------|---------|
| Lazy ring derivation (only if π detectors selected) | Save ~5–10% upfront time |
| Partial charge estimation cache | Better scoring for electrostatic interactions |
| Spatial hashing alternative | Lower build cost for extremely large N |

Proceed to `performance_engine.md`.
