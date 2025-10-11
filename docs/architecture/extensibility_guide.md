# 11. Extensibility Guide

## 11.1 Philosophy
Extensibility emphasizes low friction: new interaction detectors plug in via a registry without core engine rewrites.

## 11.2 Extension Types
| Type | Examples | Mechanism |
|------|----------|----------|
| Detector | New interaction like cation–anion ring clamp | Registry decorator |
| Feature | Partial charge estimator | FeatureStore augmentation |
| Exporter | Custom graph format | Reporting plugin interface (future) |
| Visualization | Specialized plot | UI panel registering function |

## 11.3 Detector Contract (Formal)
```
@register_detector(name: str, version: str, category: str="standard")
def detect(store, params, context) -> (list[dict], dict):
    ...
```
`params` object resolves namespaced keys (e.g., `hbond.max_distance`).

## 11.4 Parameter Declaration Pattern
Central configuration file enumerates default scalar values and ranges. Custom detectors should extend this map to appear in UI automatically (planned dynamic discovery).

## 11.5 Checklist for New Detector
1. Define algorithm & minimal geometry metrics.  
2. Decide mandatory vs optional parameters.  
3. Implement detect function using FeatureStore (avoid re-parsing).  
4. Emit funnel metrics.  
5. Add parity test (if replacing existing approach) or golden baseline snapshot.  
6. Document criteria in scientific doc.  
7. Add to UI detectors list.  
8. Update normalization logic if schema requires new field.

## 11.6 Anti-Patterns
| Anti-Pattern | Why Bad | Correction |
|-------------|---------|-----------|
| Rebuilding KD-tree inside detector | Wastes CPU | Use shared tree from context |
| Manual path building to caches | Risk of collisions | Use cache API helpers |
| Embedding UI code in detector | Breaks separation | Keep detectors pure |
| Using global mutable state | Non-deterministic | Pass via params/context |

## 11.7 Example: Directional Interaction Template
```
vec = coords[b] - coords[a]
dist = np.linalg.norm(vec)
if not (min_d <= dist <= max_d): return
unit = vec / dist
angle = np.degrees(np.arccos(np.clip(np.dot(unit, axis_ref), -1, 1)))
if angle > angle_cut:
    return
# record assembly
```

## 11.8 Performance Considerations
- Pre-extract subsets (`subset = coords[candidate_indices]`).
- Avoid Python loops where a NumPy mask suffices.
- Batch computations (distance arrays) then apply boolean filters.

## 11.9 Testing Template
```
def test_new_detector_counts(sample_store):
    recs, m = detect(sample_store, params, context)
    assert m['accepted_pairs'] == EXPECTED
```

## 11.10 Versioning Strategy
Increment detector `version` when output semantics (criteria thresholds or record fields) change; normalization can incorporate version tag for clarity.

## 11.11 Future Extensibility Hooks
| Hook | Purpose |
|------|--------|
| Pre-normalization pipeline events | Insert custom transformations |
| Export plugin registry | Add output formats without editing core |
| Visualization registry | Lazy load heavy plotting libs when needed |

Proceed to `testing_and_quality.md`.
