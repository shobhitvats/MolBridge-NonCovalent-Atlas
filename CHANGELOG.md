# Changelog

## 1.1.0 (2025-09-28)
### Added
- Unified funnel metrics across all detectors (raw_pairs, candidate_pairs, accepted_pairs, acceptance_ratio, candidate_density).
- Adaptive KD-tree thresholds extended to chalcogen & tetrel detectors with persistence (`MOLBRIDGE_ADAPTIVE_CACHE_PATH`).
- Microbenchmark harness (`performance.microbench`) and regression gate (`performance.bench_regression_gate`).
- CI integration for microbench p95 + acceptance ratio gating.
- Shared pruning heuristic `should_flag_kdtree` + boundary unit tests.
- Persistence tests, failure-path test, bench gate exit-code test, microbench schema test.
- Dependency lock file `requirements.lock` for reproducibility.
- Quick Performance Validation section in README.

### Changed
- Refactored microbench to align with registry tuple format; normalized detector keys.
- Bumped project version to 1.1.0 and updated license metadata to MIT.

### Fixed
- Microbench baseline JSON malformed duplication removed; now a single valid object.

### Known Issues
- Funnel metrics currently not populated in synthetic microbench output (instrumentation key mapping may be extended later).
- Some detectors (anion-π) raise attribute errors in synthetic runs due to missing feature precomputation; tolerated for timing envelope but future improvement planned.

## 1.0.1
(See README Version History for earlier release notes.)
