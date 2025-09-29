## MolBridge Optimization Backlog (Deferred Scope)

This document freezes the current (Sept 2025) state of identified performance / architecture improvement opportunities that are intentionally deferred. Treat this as a staging queue – pick the highest impact & lowest risk items first.

### Legend
| Tag | Meaning |
|-----|---------|
| [P] | Performance (speed / scalability) |
| [M] | Memory / footprint |
| [A] | Architecture / maintainability |
| [O] | Observability / metrics |
| [R] | Research / accuracy improvement |

Priority tiers: High (H), Medium (M), Low (L)

---
### 1. Detector Migrations (Remaining O(N²) Paths) [P][H]
Detectors still using legacy nested loops should be ported to the unified spatial pruning & vector geometry API:
* chalcogen_bonds
* pnictogen_bonds
* tetrel_bonds
* anion_pi_interactions
* n_pi_star_interactions
* london_dispersion
* sulfur_pi_interactions (verify fully migrated – partial vector geometry present)

Expected benefit: 1.5–8× speedups depending on candidate density.

Action Pattern:
1. Request coordinate slice & relevant feature arrays from `FeatureStore`.
2. Use `neighbor_within(radius, mask=<subset>)` or `neighbor_between(radius, group_a, group_b)`.
3. Vector filter geometry (angles / planarity) with NumPy arrays.
4. Materialize accepted interactions (possibly lazy – see item 3).

### 2. Multi-Radius Neighbor Cache Promotion to Shared Scope [P][M][M]
Current radius-pair cache lives per FeatureStore instance. For batch execution across structures in one session this is fine, but across processes we rebuild identical radii (e.g., 4.0Å, 5.5Å). Investigate:
* Shared read-only memory block containing compressed neighbor lists keyed by (structure_hash, radius).
* Memory map (.npz) with LRU eviction.

### 3. Lazy Interaction Materialization (Dict → Structured Array) [P][M][H]
Presently each accepted pair is turned immediately into a Python dict. Proposal:
* Maintain columnar arrays (int32 indices, float32 metrics, uint8 flags).
* Provide projection layer that yields dicts only when exporting / UI table requests.
* Add `MOLBRIDGE_LAZY_MATERIALIZATION=1` flag.

### 4. Pair Storage Representation Upgrade [P][M][H]
Neighbor pair lists (list[tuple[int,int]]) should become a contiguous `np.ndarray[int32]` shape (N,2). Reduces Python overhead & accelerates downstream filters.

### 5. Ring Normal Fast Path [P][M][M]
Replace covariance + eigen decomposition with direct cross-product normal (triangle fan) for 6-member rings; fallback to current method on pathological geometry / hetero aromatic cases.

### 6. Adaptive Threshold Model Enrichment [P][O][M]
Add new predictors: acceptance_ratio rolling mean, variance of candidate density, ring_count, charged_count (already partly present), maybe log atoms.
Use lightweight ridge regression to prevent coefficient explosion.

### 7. Cross-Process KD-tree Reuse [P][H]
Prototype using shared memory for the KD-tree raw coordinate block + serialized tree (SciPy cKDTree pickling speed?). If too heavy, implement coarse voxel grid reuse across processes.

### 8. Structured Logging Laziness [O][M]
Skip constructing detailed log dicts unless `settings.json_logging` or `verbose_detector_logs` true. Add guard macro-like helper.

### 9. Provenance Hash Incremental Mode [A][M]
Cache partial hashes per residue / chain; recompute diff only when structure subset changes.

### 10. OpenTelemetry Exporter Integration [O][M]
Replace stub with minimal OTLP HTTP exporter (config gated). Map detector spans & batch span.

### 11. GPU Feasibility Spike (Distance / Angle Kernels) [R][L]
CuPy prototype for large structures (>50k atoms) – measure crossover point.

### 12. Detector Plugin API [A][M]
Formal interface: registration metadata, required feature sets, optional vector path, fallback path, auto parameter schema injection into UI.

### 13. Memory Pressure Adaptive Behavior [P][M]
If RSS > threshold mid-run: temporarily disable vector geometry for low-yield detectors & shrink caches.

### 14. Unified Parameter Registry Finalization [A][M]
Move scattered defaults into single schema file with type, range, preset overrides, docstrings.

### 15. External Ligand Handling Expansion [R][M]
Add detection toggles for ligand-only vs. protein-ligand interface interactions.

---
## Recently Completed (September 2025 Freeze)
* Vector hydrogen bond angular batching
* Warm-start adaptive thresholds & env precedence
* Multi-radius LRU neighbor cache
* Logging verbosity gating (`verbose_detector_logs`)
* Fused normalization pass
* Adaptive model expansion (atoms + rings + charged)
* Shared memory feature distribution baseline

---
## Triage Guidance
| Impact | Typical Effort | Recommended Order |
|--------|----------------|-------------------|
| High (≥3× speed on hotspots) | Moderate | 1,4,3 |
| Medium (10–40% total improvement) | Moderate–High | 2,6,7 |
| Low / Strategic | Variable | 5,8,9,10 |

Start with detectors whose raw→candidate acceptance ratio <2% (high pruning waste) then move to lazy materialization to slash Python object overhead.

---
## Success Metrics (Define Before Implementation)
Record per change:
* Wall time delta (p95) across golden PDB set
* Memory (RSS) delta
* Detector acceptance ratio shifts
* Interaction count parity (must hold)
* CPU utilization profile (optional)

Maintain a CHANGELOG section referencing issue IDs for transparency.
