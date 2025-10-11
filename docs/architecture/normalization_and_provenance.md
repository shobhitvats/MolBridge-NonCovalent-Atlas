# 9. Normalization & Provenance

## 9.1 Goals
Provide a consistent schema for downstream analytics and a cryptographic link between results and the exact computational context.

## 9.2 Normalization Steps
1. Receive heterogeneous canonical records per detector.
2. For each record build participant objects with explicit roles.
3. Sort participants deterministically.
4. Promote geometry metrics into `metrics` dict.
5. Attach provenance hash if enabled.
6. Optionally classify strength / subtype fields.

## 9.3 Provenance Hash Algorithm
Inputs (string joined with `|`):
- Detector names sorted
- Parameter JSON canonical dump (sorted keys)
- Structure signature (chain/resnum/atom ordering + coordinate rounding)
- Application version / git hash (if available)

Run SHA256 → hex digest; truncated for UI display. Full retained internally for collision safety.

## 9.4 Deterministic Ordering Guarantees
| Element | Strategy |
|---------|----------|
| Participants in an interaction | Sort by (chain, resnum, atom, role) |
| Detector iteration order | Pre-sorted registry keys |
| Metric key order (serialization) | JSON dumps with `sort_keys=True` |

## 9.5 Schema Evolution Strategy
| Change Type | Handling |
|------------|----------|
| Add optional field | Backwards compatible |
| Rename field | Dual write (old+new) across one version window |
| Remove field | Deprecation notice + major version bump |

## 9.6 Strength Classification Example (Hydrogen Bond)
```
if distance <= 2.9 and angle >= 155: strong
elif distance <= 3.2 and angle >= 140: moderate
else: weak
```
Applied post-normalization to ensure standardized metrics.

## 9.7 Provenance in Exports
Embedded block (YAML front-matter style in some formats):
```
provenance:
  hash: a1b2c3d4e5f6
  detectors: [hbond, pi_stacking]
  parameters:
    hbond.max_distance: 3.6
    pi_stacking.max_distance: 6.5
  generated: 2025-09-30T12:30:00Z
```

## 9.8 Verification Recipe
1. Extract provenance block from export A.  
2. Re-run analysis with same parameter set.  
3. Normalize results; compute new provenance hash.  
4. Hashes must match; if mismatch, identify drift cause (code version, parameter default changed, structure differences).

## 9.9 Future Enhancements
| Idea | Benefit |
|------|--------|
| Binary canonical serialization (MessagePack) | Faster hashing on large sets |
| Provenance diff CLI | Quick mismatch diagnostics |
| Git hash auto capture | Stronger reproducibility anchor |

Proceed to `api_and_ui_integration.md`.
