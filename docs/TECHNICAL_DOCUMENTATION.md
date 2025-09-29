## MolBridge Technical Documentation

Version Snapshot: September 2025

This document provides an engineering-centric view of the system: architecture, core modules, processing pipeline, performance infrastructure, configuration, and extensibility points.

### 1. High-Level Architecture
```
                +-------------------+
 User (UI/API)  |  Streamlit / REST |
                +---------+---------+
                          |
                          v
                 +------------------+           +------------------+
                 |  Batch / Unified |  <------> |  Feature Store    |
                 |  Processors      |           |  (per structure)  |
                 +---------+--------+           +------------------+
                           | (detector dispatch)
       +-------------------+-------------------+-------------------+
       v                   v                   v                   v
  Hydrogen Bonds     π-π / Cation-π       Ionic / Salt        Others (X...
  (vector path)      (vector path)        etc.

  (All detectors emit structured interaction records + funnel metrics)
```

### 2. Core Modules Overview
| Module | Responsibility |
|--------|----------------|
| `src/analysis/*` | Detector implementations, batch processors, feature extraction helpers |
| `src/utils/` | Settings, caching, logging config, geometry helpers, provenance, session mgmt |
| `src/performance/` | Auto-tuning, timing utilities, adaptive threshold logic |
| `src/reporting/` | Report generation, export serializers |
| `src/api/` | FastAPI endpoints bridging programmatic access |
| `app/ui/` | Streamlit panels & layout / performance profile controls |
| `src/tests/` | Parity, regression, metrics & config tests |

### 3. Processing Pipeline (Unified High Performance Path)
1. Input normalization (PDB fetch / parse via Biopython)
2. FeatureStore population:
   * Atomic coordinates array (float32 N×3)
   * Aromatic ring centroids & normals
   * Charged & acidic group centroids
   * Donor / acceptor participants
3. Spatial pruning: KD-tree neighbor queries (multi-radius LRU cache)
4. Detector candidate enumeration (vector / subset / cross-set APIs)
5. Geometry filtering (angles, planarity, centroid offsets)
6. Interaction acceptance & (optional) normalization
7. Funnel metrics instrumentation (raw→candidate→accepted)
8. Adaptive threshold feedback update
9. Aggregation & (optional) reporting/export

### 4. Adaptive Threshold System
Location: `utils/kdtree_thresholds.py`
* Maintains target candidate density per detector.
* Observes raw & candidate pair counts; adjusts thresholds (up or down) using heuristics plus optional warm-start scaling.
* Environment variables override at highest precedence.
* Optional persistence to JSON (`MOLBRIDGE_ADAPTIVE_CACHE_PATH`).

### 5. Parallel & Shared Memory Execution
* Process pool path enabled via profile or manual toggle (`MOLBRIDGE_USE_PROCESS_POOL`).
* Large static arrays (coords, ring centroids, normals, charged/acidic centroids) placed into POSIX shared memory blocks when `MOLBRIDGE_USE_SHM=1`.
* Fallback to threads for small workloads (auto profile).

### 6. Vector Geometry Acceleration
* Fast batched distance & angle math using NumPy for hydrogen bonds, π-π stacking, cation-π, hydrophobic (partial), etc.
* Future: extend to chalcogen/pnictogen/tetrel dispersion (see backlog).

### 7. Normalization & Columnar Data
* Fused single-pass creation of canonical + normalized representations when `MOLBRIDGE_ENABLE_NORMALIZATION=1`.
* Columnar storage optional (env flags) with direct serialization (`MOLBRIDGE_ENABLE_COLUMNAR`, `MOLBRIDGE_DIRECT_COLUMNAR_JSON`).

### 8. Logging & Metrics
* JSON mode: `MOLBRIDGE_JSON_LOGS=1`.
* Detector verbosity gating via `performance_mode` and `verbose_detector_logs` flag.
* Metrics funnel persisted optionally to `metrics.jsonl` and summarized by CLI.

### 9. Provenance
* Enabled via `MOLBRIDGE_ENABLE_PROVENANCE=1` to attach structural hash & parameter signature for reproducibility.

### 10. Configuration Layers
| Layer | Purpose |
|-------|---------|
| Environment Flags (`Settings`) | Fast runtime toggles / feature flags |
| AppConfig (`src/utils/config.py`) | User parameter & preset management |
| UI Session State | Transient interactive control & caching |

### 11. Extending a Detector (Checklist)
1. Define parameter schema additions in central registry (or module-level defaults if interim).
2. Acquire features from FeatureStore (avoid recomputation).
3. Use `neighbor_within` / `neighbor_between` for pruning.
4. Vectorize geometry filters.
5. Emit funnel metrics: raw_pairs, candidate_pairs, accepted_pairs.
6. Add unit test(s) (geometry parity, vector vs legacy if refactor).
7. Register in UI parameter table.

### 12. Error Handling & Resilience
* All acceleration flags are fail-safe; exceptions trigger fallback & log at debug.
* Cache read/write guarded; corrupt entries evicted.
* Adaptive system bounds threshold adjustments within reasonable literature ranges.

### 13. Security & Integrity
* Dependency pinning via `requirements.lock` optional.
* CI: interaction count & timing regression gates, vulnerability scans, Docker build verification.

### 14. Testing Strategy
* Golden snapshot parity tests
* Vector vs legacy geometry parity (π-π, H-bond, etc.)
* Cache hot layer behavior tests
* Threshold adaptation boundary tests

### 15. CLI Utilities
| CLI | Purpose |
|-----|---------|
| `cli_metrics` | Summarize metrics JSONL |
| `cli_golden` | Snapshot & compare golden baselines |
| (others) | Bench / perf regression scripts under `performance/` |

### 16. Deployment Notes
* Single entrypoint `server.py` for Streamlit (also instantiates API if configured).
* Dockerfile provides slim runtime environment.

### 17. Backlog Reference
See `OPTIMIZATION_BACKLOG.md` for deferred tasks.

---
## Appendix A: Environment Flag Index
Summarizes all recognized `MOLBRIDGE_` flags (may duplicate README for engineering convenience).

| Flag | Category | Effect |
|------|----------|--------|
| `MOLBRIDGE_ENABLE_VECTOR_GEOM` | Geometry | Enable vector paths |
| `MOLBRIDGE_USE_PROCESS_POOL` | Parallel | Use process pool |
| `MOLBRIDGE_USE_SHM` | Parallel | Shared memory arrays |
| `MOLBRIDGE_USE_NUMBA` | Geometry | Numba kernels |
| `MOLBRIDGE_USE_RUST` | Geometry | Prefer Rust extension |
| `MOLBRIDGE_TASK_GRAPH` | Feature | Precompute reusable features |
| `MOLBRIDGE_ENABLE_NORMALIZATION` | Data | Emit normalized records |
| `MOLBRIDGE_ENABLE_PROVENANCE` | Data | Include provenance block |
| `MOLBRIDGE_ENABLE_COLUMNAR` | Data | Columnar store mode |
| `MOLBRIDGE_DIRECT_COLUMNAR_JSON` | Data | Direct columnar JSON serialization |
| `MOLBRIDGE_VERBOSE_DETECTOR_LOGS` | Logging | Force verbose detector logs |
| `MOLBRIDGE_PERF_PROFILE` | Perf | Performance profile selection |
| `MOLBRIDGE_ADAPTIVE_CACHE_PATH` | Adaptive | Persist adaptive thresholds |

---
## Appendix B: Data Schemas (Simplified)
Hydrogen Bond (canonical):
```json
{
  "type": "hydrogen_bond",
  "donor_atom": "NE2",
  "acceptor_atom": "O",
  "distance": 3.12,
  "angle": 156.4,
  "chain_d": "A",
  "chain_a": "B"
}
```
Normalized (example fields – actual may include provenance & strength):
```json
{
  "interaction_type": "hydrogen_bond",
  "participants": [["A","45","NE2"],["B","66","O"]],
  "metrics": {"distance":3.12,"angle":156.4},
  "strength_class": "moderate"
}
```
