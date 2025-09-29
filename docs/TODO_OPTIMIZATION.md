Optimization & Refactor Roadmap
================================

Legend: [x] done / [~] partial / [ ] pending

Architecture
------------
[x] Split batch processor monolith
[x] Add detector registration scaffold
[x] Decorator adoption for all detectors
[x] Backward key compatibility tests

Geometry
--------
 [x] Vector hydrogen bonds (vector path + instrumentation; hydrogen angle batch vectorization complete; benchmark suite extended)
[x] KD-tree pruning layer (unified helper + env tuning + adaptive runtime threshold tuning; hydrogen bonds + ionic + hydrophobic instrumented; shared memory serialization for participants added)
[x] Ionic centroid prefilter (charged vs acidic group centroids; prunes atom-level candidate pairs)
[x] Salt bridge detection (centroid-based positive vs acidic groups)
[x] Shared aromatic ring & centroid feature extraction (FeatureStore + reuse across π-π, cation-π, sulfur-π)
[x] Vectorized angle filtering generalization (halogen & pnictogen integrated; chalcogen θ & δ now vectorized; tetrel θ1/θ2 vectorized)

Parallelism
-----------
[x] Refactored high performance orchestration
 [x] Process pool for CPU-heavy detectors (experimental flag, coords + ring + charged + acidic + HBond donor/acceptor SHM prototype + consumption path)
 [x] Task graph precompute stage (shared feature extraction & wiring; worker consumption via shared meta for all major features)

Caching
-------
 [x] Feature-level ephemeral caching (rings + charged (LYS/ARG/HIS) + acidic (ASP/GLU) centers; hydrogen bond donors/acceptors; centroid prefilter feeding ionic)
[x] Bundle result caching (aggregate) with manifest
[x] Parameter subset hashing per detector (stable subset hash utility)

Observability
-------------
[x] Env-gated system tuning
 [x] Structured logging & metrics export (batch-complete + per-detector streaming JSON + metrics channel + batch metrics emission; enriched CLI with per-detector p95 & adaptive threshold stats)
 [x] Actual per-detector timing in parallel path (thread + process variants captured)

Testing & Benchmarking
----------------------
 [x] pytest-benchmark suite (HBond vector vs legacy benchmark added; extended curated structure list driving thresholds)
 [x] Golden dataset for regression (snapshot utility + compare function; curated list file and CI integration)
 [x] Performance CI threshold (workflow enhanced: artifact reuse + metrics export; strict gating with non-zero exit on regression)

Packaging
--------
[x] Optional extras for performance & profiling
[x] Requirements lock generation doc (guidance added; integrate into CI pending)

Future
------
 [x] Numba / Cython acceleration (optional numba distance kernel gated by MOLBRIDGE_USE_NUMBA)
 [x] Rust geometry kernel prototype (initial pairwise_sq_dists via PyO3; integration fallback wrapper)
 [x] Distributed execution option (Ray/Dask) scaffold (DistributedProcessor facade; baseline coverage added)

This document evolves with implementation progress.
