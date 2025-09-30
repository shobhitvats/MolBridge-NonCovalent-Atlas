# 8. Caching & Persistence

## 8.1 Objectives
Reduce redundant computation, enable rapid re-analysis, and support reproducible parameter/application states.

## 8.2 Cache Layers
| Layer | Medium | Scope | TTL |
|-------|--------|-------|-----|
| Hot In-Memory | Python dict / LRU | Single process session | Session |
| Persistent Disk | diskcache | Cross session | Configurable / indefinite |
| Adaptive Threshold Store | JSON file | Detector state | Manual invalidation |

## 8.3 Cache Keys
```
structure_hash = SHA1(structure_bytes canonicalized)
feature_key = f"feat:{structure_hash}:v{FEATURE_VERSION}"
result_key = f"res:{structure_hash}:{detector_set_hash}:{param_hash}" 
```

## 8.4 Structure Hash Canonicalization
- Sort atoms by (chain, resnum, atom_name)
- Concatenate element + coordinates (rounded to fixed decimals)
- Hash bytes → reduces false cache misses due to file ordering differences.

## 8.5 Eviction Policy
Hot cache uses LRU with max size (configurable). Diskcache size-based eviction when byte quota exceeded (future: enforce programmatically).

## 8.6 Consistency & Staleness Prevention
Parameter hash incorporates all user-facing threshold values. Any change → new result key, preventing mismatched reuse.

## 8.7 Provenance Coupling
Normalization step attaches `provenance_hash`; mismatch detection possible if user attempts to merge incompatible result sets.

## 8.8 Failure Handling
| Failure | Response |
|---------|----------|
| Disk full | Warn + gracefully disable persistence |
| Corrupt entry | Delete key & recompute |
| Permission error | Fallback to in-memory only |

## 8.9 Export Persistence
Reports and exports optionally embed minimal metadata block including provenance and parameter snapshot for future matching.

## 8.10 Future Enhancements
| Idea | Benefit |
|------|--------|
| SQLite-backed feature store | Faster random access for huge datasets |
| Remote object store integration (S3) | Team sharing of caches |
| Automatic cache warm seeds | Pre-populate common PDB IDs |

Proceed to `normalization_and_provenance.md`.
