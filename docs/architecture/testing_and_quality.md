# 12. Testing & Quality

## 12.1 Quality Objectives
Ensure scientific correctness, performance stability, and reproducibility across releases.

## 12.2 Test Layers
| Layer | Purpose | Examples |
|-------|---------|----------|
| Unit | Validate small functions | Angle clamp, centroid calc |
| Component | Detector logic correctness | Hydrogen bond criteria |
| Parity | Vector vs legacy output equivalence | π–π stacking tests |
| Regression | Detect unintended changes | Golden snapshot comparisons |
| Performance | Timing stability | Microbench harness |
| Integration | Cross-layer orchestration | API analyze endpoint |

## 12.3 Golden Snapshot Flow
1. Run standardized structure set; capture normalized output.  
2. Store canonical JSON as baseline.  
3. Future runs diff against baseline; differences flagged for review.

## 12.4 Metrics Validation
Check monotonic funnel: `raw_pairs >= candidate_pairs >= accepted_pairs`. Any violation triggers failure.

## 12.5 Performance Guard Strategy
Store previous mean timings; allow tolerance window (e.g., +15%). Exceeding triggers investigation.

## 12.6 Failure Forensics Toolkit
| Symptom | Investigation Path |
|---------|--------------------|
| Interaction count drop | Compare adaptive threshold values first |
| Spike in raw pairs | KD-tree radius or mis-set parameter |
| NaN angles | Look for zero-length normals or degenerate vectors |
| Memory bloat | Large intermediate arrays retained; check vectorization chunking |

## 12.7 Test Data Design
Use small curated PDB fragments that each intentionally exercise a distinct edge case (e.g., ring adjacency, mixed charge cluster, metal center geometry).

## 12.8 Mocking & Isolation
- Mock network fetch for remote PDB retrieval.  
- Provide in-memory FeatureStore fixtures.  
- Isolate adaptive thresholds by resetting state before tests.

## 12.9 Coverage Goals
| Area | Target % | Rationale |
|------|----------|-----------|
| Core geometry functions | 95 | Critical correctness |
| Detector decision branches | 90 | Avoid silent logic drift |
| Normalization | 90 | Schema integrity |
| Adaptive thresholds | 85 | Stability under extremes |

## 12.10 Continuous Quality (Future)
- Optional reinstated CI pipeline running: lint, unit, parity, performance smoke.
- Artifact upload of metrics JSON for trend dashboards.

Proceed to `deployment_and_env.md`.
