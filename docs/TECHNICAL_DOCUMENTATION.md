## MolBridge Technical Documentation (Beginner–Friendly Deep Dive)

Version: September 2025 (Expanded Edition)

> Goal of this document: explain “how MolBridge works” from the ground up so that a beginner who knows only basic Python can still follow, while providing links and depth for experienced engineers.

---
## Table of Contents
1. Introduction & Audience Levels  
2. What Problem Does MolBridge Solve?  
3. Core Concepts (Three-Minute Primer)  
4. Technology Stack & Key Libraries  
5. Annotated Directory Structure  
6. High-Level Architecture (Conceptual + Data Flow)  
7. Life of a Request / Analysis (Step‑by‑Step)  
8. Data Inputs (PDB IDs, Uploaded Files)  
9. Feature Extraction (What We Precompute & Why)  
10. Interaction Detectors (Design & Patterns)  
11. Geometry & Math Primer (Angles, Planes, Centroids)  
12. Spatial Acceleration (KD-Trees & Pruning)  
13. Vectorization & Performance Profiles  
14. Parallel Processing & Shared Memory  
15. Adaptive Threshold Algorithm (Explained Simply)  
16. Configuration System & Precedence  
17. Logging, Metrics & Funnel Instrumentation  
18. Provenance & Reproducibility  
19. Data Models & Schemas (Canonical + Normalized)  
20. UI Architecture (Streamlit Panels)  
21. REST API Architecture (FastAPI)  
22. CLI Tools (Usage Examples)  
23. Extending MolBridge (New Detector Tutorial)  
24. Testing Strategy (From Unit to Regression)  
25. Error Handling & Resilience Patterns  
26. Security & Dependency Hygiene  
27. Deployment Scenarios (Local, Docker, Cloud)  
28. Performance Tuning Cookbook  
29. Troubleshooting FAQ  
30. Glossary (Jargon → Plain Language)  
31. Learning Resources & References (Interactive Links)  
32. Version & Change Log  

### Architecture Suite Module Index (Deep Dives)
The monolithic document (this file) is complemented by a modular architecture suite for readers who prefer topic-focused deep dives. Each module is linked below (all paths relative to `docs/`):

| # | Module | Description |
|---|--------|-------------|
| 1 | [architecture/overview.md](architecture/overview.md) | Narrative system overview, principles, quality attributes |
| 2 | [architecture/data_model.md](architecture/data_model.md) | Atom schema, bitmasks, feature store artifacts, invariants |
| 3 | [architecture/detectors.md](architecture/detectors.md) | Detector taxonomy & per-interaction algorithms |
| 4 | [architecture/feature_pipeline.md](architecture/feature_pipeline.md) | Stepwise feature extraction stages & validations |
| 5 | [architecture/performance_engine.md](architecture/performance_engine.md) | Optimization pillars, profiling, vector & parallel strategies |
| 6 | [architecture/adaptive_thresholds.md](architecture/adaptive_thresholds.md) | Adaptive controller design, parameters, stability |
| 7 | [architecture/concurrency_and_parallelism.md](architecture/concurrency_and_parallelism.md) | Scheduling, shared memory model, determinism |
| 8 | [architecture/caching_and_persistence.md](architecture/caching_and_persistence.md) | Multi-layer caching, key design, invalidation |
| 9 | [architecture/normalization_and_provenance.md](architecture/normalization_and_provenance.md) | Canonical schema, provenance hash specification |
| 10 | [architecture/api_and_ui_integration.md](architecture/api_and_ui_integration.md) | FastAPI contract & Streamlit wiring |
| 11 | [architecture/extensibility_guide.md](architecture/extensibility_guide.md) | Extension contracts & patterns |
| 12 | [architecture/testing_and_quality.md](architecture/testing_and_quality.md) | Test matrix, coverage, failure forensics |
| 13 | [architecture/deployment_and_env.md](architecture/deployment_and_env.md) | Profiles, resource sizing, upgrade processes |
| 14 | [architecture/roadmap_and_decisions.md](architecture/roadmap_and_decisions.md) | ADRs, roadmap horizons, KPIs & revisit triggers |

Use these specialized documents for drilling deeper into a specific layer without scanning the entire technical narrative.

---
## 1. Introduction & Audience Levels
MolBridge is a toolkit + web application that finds noncovalent interactions inside protein structures (optionally with ligands). It supports researchers (structural biology, computational chemistry) and developers integrating data pipelines.

We intentionally layer explanations:
| Level | You Know… | Read These Sections First |
|-------|-----------|---------------------------|
| Beginner | Python basics, lists, functions | 2, 3, 4, 6, 7, 9 |
| Intermediate | NumPy, basic OOP, API calls | 6–15, 17–22 |
| Advanced | Performance tuning, concurrency | 13–15, 17–19, 22–28 |

---
## 2. What Problem Does MolBridge Solve?
Structural biologists want to understand how atoms “hold together” via subtle forces (hydrogen bonds, π–π, cation–π, halogen bonds, etc.). Doing this manually from 3D coordinates is slow and error‑prone. MolBridge automates:
* Parsing PDB/CIF structures
* Computing derived geometric features
* Applying literature‑backed criteria
* Reporting counts, hotspots, and visual overlays

---
## 3. Core Concepts (Three-Minute Primer)
| Concept | Plain Explanation | Why It Matters |
|---------|------------------|----------------|
| Interaction Detector | A Python function/class that filters atom pairs/triples to find one interaction type | Modular extensibility |
| Feature Store | Precomputed arrays (coordinates, ring centroids) | Avoid recomputation; speed |
| KD-Tree | Spatial index for “what’s near what?” queries | Reduces O(N²) brute force |
| Vector Path | Using NumPy arrays to do many calculations at once | Fast vs Python loops |
| Funnel Metrics | raw → candidate → accepted counts | Diagnose over/under filtering |
| Adaptive Threshold | Auto-tunes distance windows | Balances speed & completeness |
| Provenance | Hash of parameters + structure signature | Reproducibility/reporting |

---
## 4. Technology Stack & Key Libraries
| Domain | Library | Purpose | Link |
|--------|---------|---------|------|
| Core Language | Python 3.11 / 3.12 | Implementation | <a href="https://www.python.org/" target="_blank" rel="noopener">python.org</a> |
| Structural Parsing | Biopython | Read PDB/CIF & derive chains, residues | <a href="https://biopython.org/wiki/Documentation" target="_blank" rel="noopener">Biopython Docs</a> |
| Numerics | NumPy | Vectorized math | <a href="https://numpy.org/doc/" target="_blank" rel="noopener">NumPy Docs</a> |
| Statistics / Sci | SciPy | KDE density, math helpers | <a href="https://docs.scipy.org/doc/" target="_blank" rel="noopener">SciPy Docs</a> |
| DataFrames | pandas | Tabular aggregation & export | <a href="https://pandas.pydata.org/docs/" target="_blank" rel="noopener">pandas Docs</a> |
| Web UI | Streamlit | Interactive app interface | <a href="https://docs.streamlit.io/" target="_blank" rel="noopener">Streamlit Docs</a> |
| REST API | FastAPI | Programmatic HTTP endpoints | <a href="https://fastapi.tiangolo.com/" target="_blank" rel="noopener">FastAPI Docs</a> |
| Models | Pydantic | Input validation for API/config | <a href="https://docs.pydantic.dev/" target="_blank" rel="noopener">Pydantic</a> |
| Visualization | Plotly / py3Dmol | 2D plots & 3D molecular viewer | <a href="https://plotly.com/python/" target="_blank" rel="noopener">Plotly</a> / <a href="https://github.com/3dmol/3Dmol.js" target="_blank" rel="noopener">3Dmol.js</a> |
| Caching | diskcache | Persistent key/value storage | <a href="http://www.grantjenks.com/docs/diskcache/" target="_blank" rel="noopener">diskcache</a> |
| Benchmarking | pytest-benchmark | Performance regression detection | <a href="https://pytest-benchmark.readthedocs.io/" target="_blank" rel="noopener">Docs</a> |
| Packaging | PyPI / pyproject.toml | Distribution metadata | <a href="https://packaging.python.org/" target="_blank" rel="noopener">Packaging Guide</a> |

---
## 5. Annotated Directory Structure
```
src/
  analysis/        # Detectors, batch/unified processors, feature store helpers
  api/             # FastAPI endpoint definitions
  performance/     # Timing utilities, parallel helpers, adaptive threshold logic
  reporting/       # Report assembly (PDF, PPTX, Excel serializers)
  utils/           # Config, caching, logging, geometry helpers, provenance
  visualization/   # Plot construction (Ramachandran, networks, distributions)
app/ui/            # Streamlit panels, styling, interactive controls
docs/              # This technical & scientific documentation
scripts/           # Performance / regression maintenance scripts
tests/             # Pytest suites (parity, timing, cache, regression)
```

---
## 6. High-Level Architecture
```
User Action (UI / API) -> Validate Inputs -> Load/Fetch Structure -> Build Feature Store
  -> For each requested detector:
       1) Get candidate atom sets (spatial queries)
       2) Apply geometry filters (angles/distances)
       3) Classify / compute metrics
       4) Record funnel metrics
  -> Aggregate results -> (Optional) Normalize -> Persist / Cache -> Visualize / Export
```

ASCII Component Diagram:
```
            +---------------+                 +------------------+
  UI / API  |  Controller   |  ---- fetch --> |  Structure Loader |
            +-------+-------+                 +---------+--------+
                    |                                   |
                    v                                   v
            +---------------+                   +---------------+
            | Feature Store | <--- derive ----- |  Raw Atoms     |
            +-------+-------+                   +---------------+
                    |
          +---------+---------+
          |  Detector Registry |
          +----+--+--+--+-----+
               |  |  |  |
               v  v  v  v
        H-Bond  ππ  Cationπ  ... (each emits records + metrics)

```

---
## 7. Life of a Request / Analysis
1. User enters PDB ID(s) or uploads file.  
2. System fetches structure (remote if not cached).  
3. Biopython parses chains/residues/atoms.  
4. Feature Store builds arrays (coordinates, aromatic ring centroids, donors/acceptors, charged groups).  
5. For each detector: use KD-tree to get proximity-based candidate sets.  
6. Apply detector-specific rules (distance thresholds, angle windows, planarity).  
7. Keep accepted interactions and track counts at each funnel stage.  
8. Adaptive threshold module optionally fine-tunes pruning cutoff for next run.  
9. Results cached & rendered (tables, 3D overlay, plots).  
10. User exports PDF/CSV/Excel or inspects metrics.

---
## 8. Data Inputs
| Input Method | Description | Notes |
|--------------|-------------|-------|
| PDB ID(s) | e.g. “1CRN, 4HHB” fetched from RCSB | <a href="https://www.rcsb.org/" target="_blank" rel="noopener">RCSB PDB</a> |
| Local Upload | User-provided .pdb or .cif | Stored in temp cache |
| Text Paste | Raw coordinate block | Parsed similarly |

Missing atoms / alt locs: currently first conformer kept; future: occupancy weighting.

---
## 9. Feature Extraction Explained
| Feature | How | Why |
|---------|-----|-----|
| Aromatic Rings | Pattern residue side chains (PHE/TYR/TRP/HIS) -> ring atom coords -> centroid + normal (via cross product) | π–π, cation–π geometry |
| Charged Groups | Identify side chain groups (ARG, LYS, ASP, GLU, etc.) -> compute centroids | Ionic / salt / cation–π |
| Donor / Acceptor Atoms | Templates based on residue + atom naming | Hydrogen bond detection |
| Coordinate Array | N × 3 float32 | Shared foundation for distance calculations |
| KD-tree | SciPy / custom neighbor indexing | Fast neighbor queries |

---
## 10. Interaction Detectors
Each detector follows a loose template:
```
def detect(structure_features, params):
    raw_pairs = spatial_prune(features, params.distance_window)
    candidate_pairs = apply_simple_filters(raw_pairs)
    accepted = precise_geometry(candidate_pairs, params)
    return accepted, funnel_metrics
```
Design Goals:
* Minimize duplicate work (reuse Feature Store).
* Keep detector code declarative: “what qualifies” rather than “how to loop.”
* Emit standardized record dictionaries.

---
## 11. Geometry & Math Primer (Made Simple)
| Term | Plain Meaning | Implementation Hint |
|------|---------------|---------------------|
| Distance | Straight line between two 3D points | NumPy vector norm |
| Angle | Turn amount between two vectors | dot(u,v) / (|u||v|) -> arccos |
| Planarity | Whether two rings are aligned | angle between normals |
| Centroid | Average position of atoms in a group | mean(axis=0) |
| Offset (lateral) | Sideways displacement from normal axis | projection & subtraction |

Recommended Reading: <a href="https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps" target="_blank" rel="noopener">LibreTexts (molecular geometry)</a>.

---
## 12. Spatial Acceleration (KD-Trees)
Brute force pair checks scale as O(N²). We build a KD-tree so “find all atoms within R of atom X” becomes O(log N) + local points. See <a href="https://en.wikipedia.org/wiki/K-d_tree" target="_blank" rel="noopener">KD-tree (Wikipedia)</a>.

Simplified Pseudocode:
```
tree = build_kdtree(coords)
for atom i:
    neighbors = tree.query_radius(coords[i], r=max_cutoff)
    evaluate only those neighbors
```

---
## 13. Vectorization & Performance Profiles
Instead of Python loops:
```
dx = X[:,0] - X[:,0][:,None]  # broadcast differences
dist2 = dx*dx + ...            # compute matrix of squared distances
```
We choose faster vector paths when environment variable `MOLBRIDGE_ENABLE_VECTOR_GEOM=1` or profile = `auto/full`.

Profiles (`MOLBRIDGE_PERF_PROFILE`):
| Profile | Focus | Tradeoffs |
|---------|-------|-----------|
| auto | Balance | Heuristic toggles |
| full | Max speed | Higher RAM usage |
| minimal | Lowest overhead | Slow for huge structures |
| manual | User control | Expertise needed |

---
## 14. Parallel Processing & Shared Memory
Why processes? Python’s GIL blocks true multi-core CPU in threads for CPU-bound loops. Process pools avoid this. Shared memory (`MOLBRIDGE_USE_SHM=1`) prevents copying large arrays to each child.

Reference: <a href="https://docs.python.org/3/library/multiprocessing.html" target="_blank" rel="noopener">multiprocessing docs</a>.

---
## 15. Adaptive Threshold Algorithm (Friendly Explanation)
Goal: dynamic distance cutoffs to keep candidate counts “healthy” (not too huge, not too tiny).

Analogy: Like adjusting fishing net size—too wide (slow), too tight (miss fish). We monitor acceptance ratios and adjust ± a small step.

Pseudo Steps:
```
for each detector run:
    observe raw_pairs, candidate_pairs
    density = candidate_pairs / raw_pairs (if raw>0)
    if density << target_low: increase distance (loosen)
    elif density >> target_high: decrease distance (tighten)
    clamp distance within literature min/max
    persist if enabled
```

Target bands chosen empirically. Persistence path: `MOLBRIDGE_ADAPTIVE_CACHE_PATH`.

---
## 16. Configuration System & Precedence
1. Hard defaults (code)  
2. Preset (Literature / Conservative / Exploratory)  
3. UI user changes  
4. Environment variables  
5. Command-line / API overrides  

Highest wins. Config object (`AppConfig`) aggregates and caches.*

---
## 17. Logging, Metrics & Funnel Instrumentation
| Stage | Meaning |
|-------|---------|
| raw_pairs | After initial spatial prune |
| candidate_pairs | After coarse geometry filters |
| accepted_pairs | Final interactions |

Metrics serialized to `metrics.jsonl`. Summaries via `cli_metrics`.

Set `MOLBRIDGE_JSON_LOGS=1` for machine-readable logs. Logging library: <a href="https://loguru.readthedocs.io/" target="_blank" rel="noopener">Loguru</a>.

---
## 18. Provenance & Reproducibility
We hash:
* Structure (sorted residue/atom signature)
* Active parameter set
* Enabled detectors

So a published figure can be recreated. Enabled by `MOLBRIDGE_ENABLE_PROVENANCE=1`.

Further reading: <a href="https://reproducibilitychallenge.org/" target="_blank" rel="noopener">Reproducibility Principles</a>.

---
## 19. Data Models & Schemas
Canonical Interaction Record (example – hydrogen bond):
```json
{
  "type": "hydrogen_bond",
  "donor_chain": "A",
  "donor_resnum": 45,
  "donor_atom": "NE2",
  "acceptor_chain": "B",
  "acceptor_resnum": 66,
  "acceptor_atom": "O",
  "distance": 3.12,
  "angle": 156.4,
  "raw_score": 0.82
}
```
Normalized Form (structure-agnostic layout):
```json
{
  "interaction_type": "hydrogen_bond",
  "participants": [
    {"chain":"A","resnum":45,"atom":"NE2","role":"donor"},
    {"chain":"B","resnum":66,"atom":"O","role":"acceptor"}
  ],
  "metrics": {"distance":3.12,"angle":156.4},
  "strength_class": "moderate",
  "provenance_hash": "abc123..."
}
```

---
## 20. UI Architecture (Streamlit Panels)
Panels live under `app/ui/` and call service logic. State stored in `st.session_state`: caches (analysis results), flags (density on/off), and current structure id.

Streamlit Basics: <a href="https://docs.streamlit.io/library/get-started" target="_blank" rel="noopener">Get Started Guide</a>.

---
## 21. REST API Architecture (FastAPI)
FastAPI routers (in `src/api/`) expose endpoints like `/analyze` returning JSON interaction bundles. Pydantic models validate input.

Key Concept: asynchronous endpoints let long tasks run without blocking event loop (optionally with background workers).

---
## 22. CLI Tools (Examples)
| Command | Purpose | Example |
|---------|---------|---------|
| `python -m cli_golden snapshot` | Create baseline | `--output golden_baseline.json` |
| `python -m cli_golden compare` | Compare with baseline | `--baseline golden_baseline.json` |
| `python -m performance.microbench` | Micro benchmark run | `> bench.json` |
| `python -m cli_metrics summarize` | Summarize metrics | `--file metrics.jsonl` |

---
## 23. Extending MolBridge (Add a Detector) – Tutorial
Example: Detect “close CA–CA pseudo contacts.”
1. Create file `src/analysis/ca_contacts.py`.
2. Add function `detect_ca_contacts(feature_store, params)`.
3. Use KD-tree to find CA atoms within threshold (e.g., 6.0 Å).
4. Build records `{"type":"ca_contact", ...}`.
5. Register in a registry map (e.g., `src/analysis/registry.py`).
6. Add config defaults (distance cutoff).
7. Add unit test confirming count on a small known PDB.
8. Add to UI parameter list.

Reference patterns: open existing detectors (e.g., hydrogen bonds) and copy structure.

---
## 24. Testing Strategy
| Test Type | Purpose | Example |
|-----------|---------|---------|
| Unit | Function correctness | Ring centroid calculation |
| Parity | Vector vs legacy outputs | π–π stacking parity test |
| Regression | Golden baseline stability | `cli_golden compare` |
| Performance | Timing not degraded | Microbench gate |
| Integration | API + detector interplay | Analyze endpoint returns expected keys |

Resource: <a href="https://docs.pytest.org/en/stable/" target="_blank" rel="noopener">pytest docs</a>.

---
## 25. Error Handling & Resilience
Pattern: “Fail soft, log clearly.”
* Wrap optional acceleration in try/except → fallback.
* Validate user parameters (distance ranges) early.
* Evict corrupt cache entries gracefully.

---
## 26. Security & Dependency Hygiene
| Aspect | Practice |
|--------|----------|
| Dependencies | Pin versions when shipping releases (`pip freeze > requirements.lock`) |
| Vulnerabilities | Run `pip-audit` locally | <a href="https://github.com/pypa/pip-audit" target="_blank" rel="noopener">pip-audit</a> |
| Secrets | Use environment variables, never commit tokens |
| Supply Chain | Minimal direct transitive footprint |

---
## 27. Deployment Scenarios
| Scenario | Steps (Summary) | Notes |
|----------|-----------------|-------|
| Local Dev | `pip install -r requirements.txt` → `streamlit run server.py` | Quick iteration |
| Docker | `docker build -t molbridge .` → `docker run -p 8501:8501 molbridge` | Reproducible env |
| Streamlit Cloud | Connect repo, set entry to `server.py` | Auto-rebuild on push |
| Heroku / Cloud Run | Container or Procfile (already present) | Scale horizontally |

Docker Reference: <a href="https://docs.docker.com/get-started/" target="_blank" rel="noopener">Docker Get Started</a>.

---
## 28. Performance Tuning Cookbook
| Goal | Action | Flag / Tool |
|------|--------|-------------|
| Faster large batch | Enable vector + process pool | `MOLBRIDGE_PERF_PROFILE=full` |
| Lower memory | Disable process pool & SHM | `MOLBRIDGE_USE_PROCESS_POOL=0` |
| Stable counts | Lock adaptive thresholds | unset `MOLBRIDGE_ADAPTIVE_CACHE_PATH` |
| Investigate slowdown | Run microbench & compare | `performance.microbench` |
| Reduce noise in CI | Use minimal profile | `MOLBRIDGE_PERF_PROFILE=minimal` |

---
## 29. Troubleshooting FAQ
| Symptom | Likely Cause | Fix |
|---------|-------------|-----|
| No interactions found | Too strict cutoffs | Relax distance/angle sliders |
| Very slow run | Huge structure + all detectors | Use `auto` profile or disable least-needed detectors |
| Memory spike | Many structures + process pool | Reduce workers / disable SHM |
| “Angle NaN” | Colinear or degenerate vectors | Safe-guard present; counts unaffected |
| API 422 error | Bad JSON schema | Check FastAPI docs for model format |

---
## 30. Glossary
| Term | Meaning |
|------|---------|
| Detector | Module that finds one interaction type |
| Feature Store | Bundle of arrays precomputed once per structure |
| Funnel Metrics | raw → candidate → accepted counts |
| Normalization | Convert various record shapes into uniform schema |
| Provenance Hash | Fingerprint tying results to parameters + structure |
| Vectorization | Batch math on arrays |
| KD-tree | Spatial index for neighbor queries |

---
## 31. Learning Resources & References
* Structural Databases: <a href="https://www.rcsb.org/" target="_blank" rel="noopener">RCSB PDB</a>  
* Geometry Math: <a href="https://mathinsight.org/" target="_blank" rel="noopener">Math Insight</a>  
* Python Packaging: <a href="https://packaging.python.org/" target="_blank" rel="noopener">Guide</a>  
* Streamlit Tutorials: <a href="https://docs.streamlit.io/library/get-started" target="_blank" rel="noopener">Start</a>  
* FastAPI Course: <a href="https://fastapi.tiangolo.com/tutorial/" target="_blank" rel="noopener">Tutorial</a>  
* Biopython Tutorial: <a href="https://biopython.org/DIST/docs/tutorial/Tutorial.html" target="_blank" rel="noopener">PDF / HTML</a>  
* Performance Tips (NumPy): <a href="https://numpy.org/devdocs/user/quickstart.html" target="_blank" rel="noopener">Quickstart</a>  
* Reproducible Research: <a href="https://the-turing-way.netlify.app/" target="_blank" rel="noopener">The Turing Way</a>

---
## 32. Version & Change Log (Technical Doc)
| Version | Date | Summary |
|---------|------|---------|
| 1.0 | 2025-09-20 | Initial concise technical overview |
| 1.1 | 2025-09-30 | Massive beginner-friendly expansion; added TOC, tutorials, hyperlinks |

---
### Quick “You Got This” Recap
1. Load structure → build features.  
2. Each detector: prune, filter, accept.  
3. Metrics recorded; thresholds adapt.  
4. UI/API present normalized, exportable data.  

Happy exploring. Add a detector, tweak a profile, and iterate!

---
# Deep Architecture & Optimization Supplement (Advanced Reference)
This supplement drills into low-level design, data structures, optimization rationale, and extension internals. Use it when you need to: profile, extend core detectors, or evaluate architectural tradeoffs.

## A. Component Responsibility Matrix
| Component | File(s) / Area | Core Responsibilities | Key Interfaces Consumed | Key Interfaces Exposed |
|-----------|----------------|-----------------------|-------------------------|------------------------|
| Structure Loader | `utils/pdb_handler.py` | Fetch, parse, cache raw structure | Biopython PDB parser | FeatureStore builder inputs |
| FeatureStore | (implicit across `analysis/structure_index.py`, helpers) | Assemble typed arrays & semantic indices | Raw Biopython objects | Uniform arrays (coords, ring_centroids, donors) |
| Detector Registry | `analysis/registry.py` | Map detector name → callable & metadata | FeatureStore, config | Detector callable entrypoint |
| Batch Processor | `analysis/batch_processor.py` | Orchestrate multi-structure analysis | Registry, FeatureStore | Structured per-structure results |
| Unified Processor | `analysis/unified_batch.py` | Optimize cross-structure runs, reuse caches | KD-tree builder, detectors | Bulk result aggregation |
| Adaptive Thresholds | (internal utils) | Maintain dynamic distance/angle cutoffs | Funnel metrics, env | Updated thresholds state |
| Caching Layer | `utils/cache.py` | Multi-tier (mem/disk) retrieval & eviction | DiskCache library | Transparent get/put primitives |
| Provenance | `utils/provenance.py` | Hash configuration & structure signature | Config, raw structure | Provenance hash, metadata block |
| Normalization | `utils/normalization.py` | Transform raw detector records | Raw records | Canonical normalized entries |
| Reporting | `reporting/report_generator.py` | Build PDF/Excel/CSV bundles | Normalized data | Binary / textual artifacts |
| Visualization | `visualization/*` | Plot generation (Ramachandran, networks) | Normalized / raw metrics | Plot objects / HTML embeddings |
| API Layer | `api/endpoints.py` | Expose programmatic REST interface | Processors, models | JSON responses |
| UI Layer | `app/ui/*` | User interactions, parameter capture | Processors, config | Streamlit components |

## B. Core Data Structures
### 1. Coordinate Array
Type: `float32` contiguous NumPy array shape (N,3).  
Rationale: `float64` unnecessary precision for distance windows (typical cutoffs < 20 Å), halving memory footprint and improving cache locality.

### 2. Atom Metadata Table
Parallel arrays or structured dtype (implementation detail may evolve):
| Field | Type | Description |
|-------|------|-------------|
| chain_id | uint8 / object | Encoded chain (mapping for memory efficiency) |
| resnum | int32 | Residue sequence number |
| atom_name | fixed-length string | Atom label (e.g., "CA") |
| element | uint8 or category | Periodic symbol index |
| flags | bitmask | Donor / acceptor / aromatic / charged roles |

Bitmask Example (uint16): bits 0..5 encode roles (DONOR=1, ACCEPTOR=2, AROM_RING=4, CHARGED_POS=8, CHARGED_NEG=16, METAL=32). Future roles can extend higher bits.

### 3. Ring Centroids & Normals
Arrays shape (R,3) each. Built once; each ring also holds index list of constituent atom row indices for optional recalculation or highlight rendering.

### 4. Charged Group Centroids
Pre-averaged positions for side chain heavy atoms contributing to charge center (e.g., guanidinium group in ARG). Reduces multi-atom evaluation to single vector math.

### 5. KD-tree Structure
Implementation: SciPy `cKDTree` or a lightweight custom wrapper.  
Stored: reference to coordinate array (no copy).  
Cache Key: `hash(coords_shape, float_precision, structure_id_version)`.

### 6. Funnel Metrics Record
Dictionary per detector: `{ "raw_pairs": int, "candidate_pairs": int, "accepted_pairs": int, "threshold_snapshot": {...} }`.

## C. Detector Algorithm Patterns
We support three canonical detector archetypes:
| Archetype | Examples | Loop Primitive | Vectorization Strategy |
|-----------|----------|----------------|------------------------|
| Pairwise Distance | Hydrogen bonds (initial), hydrophobic contacts | KD-tree neighbor list → pair filter | Broadcast distance arrays or incremental filtered loops |
| Pairwise Planar | π–π stacking, cation–π | Pre-filter by centroids within R, then angle between normals | Vector dot normalization → boolean mask |
| Multi-Stage (Directional) | Halogen bonds, chalcogen, pnictogen | Step1: distance prune; Step2: angle + σ-hole alignment | Two-phase mask refinement |

### Example: π–π Stacking Detailed Steps
1. Input ring centroids C (R×3) and normals N (R×3).  
2. KD-tree radius search with `R_max` (e.g., 6.5 Å) produces candidate index pairs P.  
3. Vector compute: center distance `d = ||C_i - C_j||`.  
4. Normal angle `θ = arccos( clamp( dot(N_i, N_j) ) )`.  
5. Offset vector `v = (C_j - C_i) - proj_{N_i}(C_j - C_i)` to evaluate lateral displacement.  
6. Accept if: `d_min ≤ d ≤ d_max`, `θ ≤ θ_parallel_cut` (parallel) OR within skew band, and lateral offset ≤ L_max.  
7. Classify subtype (parallel vs T-shaped) using angle + offset heuristics.  
8. Append record, update metrics.

### Example: Hydrogen Bond Advanced Filters
1. Donor-H (D-H) and acceptor (A) atom identification via residue templates.  
2. Donor heavy atom vs acceptor distance prune (≤ 3.6 Å).  
3. Compute D-H-A angle (if H exact position estimated via idealized geometry when not present).  
4. Compute H-A-Base angle for linearity discrimination.  
5. Score = weighted function of distance + angle alignment.  
6. Strength classification mapping (tight thresholds → strong).

## D. Mathematical Derivations (Selected)
### 1. Angle Between Ring Normals
Given normals n1, n2:  
`θ = arccos( max(-1, min(1, (n1 · n2) / (||n1|| ||n2||))) )`  
Optimization: normals pre-normalized → skip denominator.

### 2. Lateral Offset for π–π
`v = c2 - c1`  
`parallel_component = (v · n1) * n1`  
`lateral = v - parallel_component`  
`offset = ||lateral||`

### 3. Projected Hydrogen Position (If Missing)
For donor heavy atom D and attached hetero atom X approximate hydrogen position H' along normalized vector to acceptor candidate A:  
`u = (A - D) / ||A - D||`  
`H' = D + d_DH * u` where `d_DH` is ideal bond length constant table.

### 4. Dihedral (n→π*) Style Interactions
Computed using four atoms (A,B,C,D).  
Vectors: `b1 = B-A`, `b2 = C-B`, `b3 = D-C`.  
Normals: `n1 = b1 × b2`, `n2 = b2 × b3`.  
Angle: `φ = atan2( ||b2|| * (b1 · n2), n1 · n2 )`.

## E. Complexity Analysis
| Stage | Complexity (N atoms, R rings) | Notes |
|-------|-------------------------------|-------|
| Feature extraction | O(N) | Linear pass over atoms |
| KD-tree build | O(N log N) | SciPy implementation |
| Neighbor queries | O(k log N) | k = number of query centers (≈N or subsets) |
| π–π detection | O(P) | P = pruned pairs (typically << N²) |
| Hydrogen bonds | O(P + refinements) | H present subset smaller |

Goal: keep P scalable by tight but safe pruning windows.

## F. Memory Footprint Anatomy
| Component | Approx Size Formula | Example (N=50k atoms) |
|-----------|---------------------|-----------------------|
| Coordinates | 12N bytes (float32×3) | ≈ 1.8 MB |
| Atom metadata | ~ (chain/resnum/flags) 12–24N bytes | ≈ 1.2 MB |
| Rings (centroids+normals) | 24R bytes | ~ small (kB) |
| KD-tree | Implementation specific (~ 2–4× coords) | ≈ 7 MB |
| Cache (disk) | Variable (serialized JSON / binary) | Configurable |

## G. Optimization Catalog
| Category | Technique | Justification | Tradeoff |
|----------|-----------|---------------|----------|
| Algorithmic | KD-tree pruning | Avoid O(N²) | Build overhead |
| Algorithmic | Two-stage filtering (distance then angle) | Early discard | Slight complexity |
| Data Layout | float32 coords | Cache-friendly | Precision reduced (acceptable) |
| Data Layout | Bitmask roles | Compact, vectorizable | Requires mask ops |
| Vectorization | NumPy broadcasting | Speed vs Python loops | Temporary memory spikes |
| Parallelism | Process pool | Multi-core speed | IPC overhead for small N |
| Parallelism | Shared memory arrays | Zero-copy across workers | Complexity in setup |
| Adaptivity | Dynamic thresholds | Balanced workload | Warm-up oscillations |
| Caching | DiskCache layering | Avoid recomputation | Stale risk (hashed keys mitigate) |
| Logging | Structured JSON optional | Tooling integration | Slight overhead |
| Normalization | Single pass merge | Avoid multiple traversals | More complex record builder |

## H. Adaptive Threshold Parameters
| Parameter | Meaning | Typical Value |
|-----------|---------|---------------|
| target_density_low | Lower candidate density bound | 0.05–0.08 |
| target_density_high | Upper density bound | 0.20–0.25 |
| step_size | Distance adjustment Å | 0.1 |
| angle_step | Angle window adjustment (deg) | 2–3 |
| min_distance_window | Hard lower clamp | Literature min |
| max_distance_window | Hard upper clamp | Literature max |

Stability Tactic: Exponential smoothing can be applied (α ≈ 0.3) to observed density before decision (future enhancement).

## I. Concurrency Model Detail
1. Main process builds FeatureStore (shared arrays optional).  
2. Work unit partition: list of detectors × structures.  
3. Dispatcher enqueues tasks in process pool if `MOLBRIDGE_USE_PROCESS_POOL=1`.  
4. Each worker references shared memory blocks (names derived from structure hash).  
5. Results (accepted interactions + metrics) returned via queue; merging ensures deterministic ordering by detector list index.

Deadlock Avoidance: large result serialization avoided by returning primitive lists; normalization delayed until aggregation.

## J. Caching Strategy Internal
| Layer | Key Construction | Expiration | Notes |
|-------|------------------|------------|-------|
| In-Memory LRU | `(structure_hash, feature_version)` | Session lifetime | Fast ephemeral |
| DiskCache Persistent | Same key + detection profile | Configurable / manual | Survives restarts |
| Adaptive Thresholds | Detector name | Manual invalidation | JSON snapshot |

Invalidation: structure hash changes if file contents differ (content-based hashing vs filename).

## K. Error Taxonomy
| Category | Examples | Handling |
|----------|----------|----------|
| Input | Missing chain, malformed PDB | Raise user-facing warning, skip structure |
| Geometry | NaN angle (degenerate vectors) | Log debug, discard pair |
| Config | Unsupported flag combo | Fallback to safe defaults |
| Parallel | Worker crash | Retry sequentially; disable pool for session |
| Cache | Corrupt serialization | Evict entry, rebuild |

## L. Extension Guidelines (Formal Contract)
Detector Callable Signature:
```
def detect(feature_store, params, context) -> tuple[list[dict], dict]:
  """Return (records, funnel_metrics)."""
```
Required funnel metrics keys: `raw_pairs`, `candidate_pairs`, `accepted_pairs` (ints). Additional keys permitted.

Record Requirements:
| Field | Required | Description |
|-------|----------|-------------|
| type | Yes | Interaction type slug |
| participants / atom fields | Yes | Either unified list or explicit donor/acceptor | 
| metrics | Yes | Numeric measurements (distance, angles) |
| subtype | Optional | Interaction subtype classification |
| score | Optional | Ranking metric |

## M. Configuration Resolution Algorithm
Pseudo:
```
cfg = defaults()
apply(preset)
apply(user_ui_changes)
apply(env_overrides)
apply(cli_or_api_overrides)
validate(cfg)
freeze(cfg)
```
Validation ensures numeric ranges, mutual exclusion (e.g., cannot enable both experimental Rust & incompatible Numba kernel simultaneously if present).

## N. Metrics Persistence Flow
1. Detector run produces metrics dict.  
2. Aggregator enriches with timestamp & structure id.  
3. Append JSON line to `metrics.jsonl` if enabled.  
4. CLI summarizer groups by detector, computes averages, ranges, std dev, anomaly outliers (IQR rule prospective).

## O. Testing Coverage Mapping (Illustrative)
| Test File | Coverage Focus |
|-----------|---------------|
| `test_hbond_subtypes.py` | Hydrogen bond classification & angles |
| `test_cation_pi_vector_parity.py` | Vector vs reference parity (cation–π) |
| `test_pipi_vector_parity.py` | π–π stacking geometry parity |
| `test_cache_hot_layer.py` | LRU + disk caching interactions |
| `test_registry_and_cache_metrics.py` | Detector registry integrity + metrics funnel |
| `test_key_resolution_and_normalization.py` | Config + normalization pathways |
| `test_timing_utils.py` | Performance timing instrumentation baseline |

Potential Additions: Negative tests for malformed PDB, adaptive threshold oscillation stability, shared memory race resilience.

## P. Scaling Behavior & Limits
| Dimension | Current Practical Limit | Notes |
|----------|-------------------------|-------|
| Atoms per structure | ~150k (depends on RAM) | KD-tree memory grows; consider chunking |
| Structures per batch | Dozens (moderate) | Parallel overhead beyond certain point |
| Detectors simultaneously | All current + future planned | Incremental; watch memory footprint |

Future: On-demand ring/charge features only if relevant detectors selected.

## Q. Future Optimization Roadmap (Aspirational)
| Idea | Expected Benefit | Complexity |
|------|------------------|-----------|
| SIMD via NumPy ufunc rewrite | 1.2–1.5× certain vector ops | Medium |
| Numba JIT for directional detectors | 1.3–2× speed on loops | Medium |
| Rust extension (PyO3) for KD pruning | Lower Python overhead | High |
| GPU distance pruning (CuPy) | Large N acceleration | High (data transfer) |
| Adaptive smoothing / PID control | Stable thresholds faster | Low |
| Ring normal pre-clustering | Fewer pair checks for π–π | Medium |

## R. Observability Enhancements (Planned)
| Feature | Description |
|---------|-------------|
| Prometheus Exporter | Expose metrics over /metrics endpoint |
| Latency Heatmap | Per detector timing distribution |
| Config Snapshot Artifact | Auto-save config + provenance with export |

## S. Design Tradeoffs Summary
| Decision | Rationale | Alternative | Revisit Criteria |
|----------|-----------|------------|------------------|
| Streamlit UI | Rapid development, research audience | React front-end | If complex collaborative features needed |
| Python core | Ecosystem, readability | C++ core | If CPU becomes bottleneck irreducibly |
| KD-tree approach | Standard, reliable | Uniform grid spatial hashing | If KD-tree build cost dominates small N |
| Dynamic thresholds | Balanced performance | Fixed literature constants | If reproducibility priority > performance |
| float32 coords | Memory & speed | float64 | If precision complaints emerge |

## T. Decision Log (Excerpt)
| Date | Decision | Outcome |
|------|----------|---------|
| 2025-07 | Adopt adaptive thresholds | Reduced outlier runtimes 18% (internal benchmark) |
| 2025-08 | Introduce funnel metrics | Faster mis-tuning diagnosis |
| 2025-09 | Square Ramachandran plot | Visualization consistency improved |

## U. Extension Example (Concrete Code Sketch)
```
from analysis.registry import register_detector

@register_detector(name="ca_contacts", version="1.0", category="experimental")
def detect_ca_contacts(store, params, context):
  coords = store.coords  # (N,3)
  ca_idx = store.atom_name_index['CA']
  pairs = context.kdtree.radius_pairs(coords[ca_idx], params.max_distance)
  raw = len(pairs)
  # Filter self, duplicates, chain equality optional
  filtered = [p for p in pairs if p[0] < p[1]]
  records = []
  for i,j in filtered:
    d = context.fast_distance(coords[ca_idx][i], coords[ca_idx][j])
    if d < params.min_distance:  # avoid steric artifact
      continue
    records.append({
      "type": "ca_contact",
      "participants": store.participant_tuple_from_indices(ca_idx[i], ca_idx[j]),
      "metrics": {"distance": float(d)},
    })
  metrics = {"raw_pairs": raw, "candidate_pairs": len(filtered), "accepted_pairs": len(records)}
  return records, metrics
```

## V. Profiling Workflow Example
1. Enable JSON logs + full profile: `MOLBRIDGE_JSON_LOGS=1 MOLBRIDGE_PERF_PROFILE=full`.  
2. Run analysis on representative large structure set.  
3. Extract timing blocks from logs (grep `TIMING`).  
4. Aggregate average per-detector time; identify top 2 hotspots.  
5. Introduce micro-optimization (e.g., preallocate array) → re-run → compare diff.

## W. Maintainability Metrics (Qualitative)
| Aspect | Current State | Risk Level |
|--------|---------------|-----------|
| Detector Cohesion | High (clear boundaries) | Low |
| Cross-Cutting Concerns (logging, metrics) | Encapsulated | Low |
| Configuration Drift | Moderate (flags vs UI) | Medium |
| Test Coverage Breadth | Good (core detectors) | Medium (advanced edge cases) |
| Complexity Concentration | Adaptive & vector paths | Medium |

Mitigation: auto-generate flag matrix docs from introspection in future.

## X. Reproducibility Recipe (Exact Reconstruction)
1. Record provenance hash from exported report.  
2. Check out commit matching stored git hash (planned inclusion).  
3. Set environment flags from captured config block.  
4. Re-run batch with identical structure input set.  
5. Compare normalized JSON via stable sort → diff must be empty except permissible floating tolerance (±1e-5 Å).

## Y. Frequently Overlooked Details
| Item | Why It Matters |
|------|----------------|
| Normalization order stability | Ensures consistent downstream diffs |
| Angle clamp in arccos | Avoids NaNs due to rounding | 
| Bitmask usage | Fast role checks (branch reduction) |
| Early bail on empty ring list | Saves constructing unused arrays |
| Consolidated export pipeline | Prevents format drift across PDF/CSV |

## Z. Meta: How to Read the Source Efficiently
1. Start at `analysis/registry.py` to see available detectors.  
2. Open one detector file and trace its imports into `utils/geom.py` / feature helpers.  
3. Inspect `utils/cache.py` for retrieval logic.  
4. Follow `normalization.py` to learn final record shape.  
5. Use tests as executable specifications.

---
End of Deep Supplement. If something is still unclear, add a “Documentation Gap” issue with a pointer to the code region.
