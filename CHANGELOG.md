# Changelog

All notable changes to this project will be documented in this file. This project follows a semantic-ish versioning approach (major.minor.patch) where detector additions & default parameter changes usually increment the minor version; fixes and documentation-only updates increment the patch.

## [Unreleased]
- (placeholder for pending changes)

## [1.2.0] - 2025-09-30
### User Experience, Overlay & Productivity Expansion
- Command Palette (Ctrl/Cmd+K) with fuzzy search (rapidfuzz) for presets, profiles, layout snapshot ops, and UI toggles.
- Scenario Profiles (`templates/scenario_profiles.yaml`) providing task-focused interaction sets & linked parameter presets.
- Layout Snapshot save / restore UI (sidebar) + palette actions for fast environment switching.
- Multi-Structure Ramachandran overlay with structure provenance in hover tooltips & stratified performance safeguards.
- Interaction Network enhancements: edge filtering, collapse by interaction type, label visibility toggle, aggregated edge stats.
- Diff improvements: coordinate-aware hashing, filtering, diff summary counters, CSV export.
- Persistent detector progress stream with live status pill & timestamped events (rotating log).
- Progressive Ramachandran density rendering (auto-disable on large residue sets) + high-contrast density colorbar toggle.
- Virtualized interaction table (AgGrid fallback) for large interaction lists.
- Parameter modified indicators + one-click revert-to-preset; autosave of session metadata.
- Unified dark theme only plus contrast & classic Ramachandran color toggles.

## [1.1.0] - 2025-07-??
### Performance & Observability Expansion
- Unified funnel metrics across all detectors.
- Adaptive KD-tree threshold persistence (opt-in via `MOLBRIDGE_ADAPTIVE_CACHE_PATH`).
- Microbenchmark harness enhancements: memory RSS delta, funnel aggregation totals, per-detector warnings, deterministic seeding.
- Baseline management script `scripts/update_microbench_baseline.py` with embedded metadata.
- OpenTelemetry export stub (`MOLBRIDGE_OTEL_EXPORT=1`).
- CI hardening: security scan, Docker smoke test, microbench regression gate integration.
- Documentation updates & MIT headers for regression scripts.
- Detector robustness fix (anion–π attribute alias).

## [1.0.1] - 2025-05-??
- New detectors: Cation–π, Salt Bridges, Sulfur–π, Metal Coordination.
- Hydrogen Bond Subtype extension.
- Added new configurable cutoffs (cation_pi, salt_bridge, sulfur_pi, metal_coordination...).
- Feature toggle: `enable_hbond_subtypes`.
- Extended reporting & exports for new analytics.

## [1.0.0] - 2025-03-??
- Initial release with 11 interaction types, Streamlit interface, REST API, reporting system, 3D visualization, batch processing.

---

### Format Guidance
- Add newest version at the top.
- Use present tense for summaries ("Add", "Fix") and group by categories if large (Added / Changed / Fixed / Removed).

### Compare Links (add after tagging releases)
<!-- Example (update repository URLs once tags exist) -->
<!-- [Unreleased]: https://github.com/shobhitvats/Protein-Interaction-Analysis-Server/compare/v1.2.0...HEAD -->
<!-- [1.2.0]: https://github.com/shobhitvats/Protein-Interaction-Analysis-Server/releases/tag/v1.2.0 -->
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
