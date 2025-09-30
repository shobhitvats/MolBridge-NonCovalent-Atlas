# 14. Roadmap & Architectural Decisions

## 14.1 Decision Record Format
Each ADR (Architectural Decision Record) captures context, decision, consequences, and revisit triggers.

## 14.2 Selected ADRs
### ADR-001: Use KD-tree for Spatial Pruning
- Context: Need scalable neighbor search.  
- Decision: SciPy cKDTree chosen.  
- Consequences: O(N log N) build; good for repeated queries.  
- Revisit: If build cost dominates for extremely sparse detectors.

### ADR-002: Float32 Coordinates
- Context: Memory + speed tradeoff.  
- Decision: Adopt float32 for coords, keep float64 only where numerically sensitive (rare).  
- Consequences: Half memory, slight precision loss acceptable.  
- Revisit: If precision error reported in edge scientific cases.

### ADR-003: Adaptive Threshold Controller
- Context: Workload variability across structures.  
- Decision: Proportional controller with clamped step.  
- Consequences: More stable runtimes; minor reproducibility variance.  
- Revisit: Need strict reproducibility or oscillation observed.

### ADR-004: Streamlit for UI
- Context: Rapid prototyping & domain scientist accessibility.  
- Decision: Streamlit chosen over custom React front-end.  
- Consequences: Faster iteration; limited deep UI customization.  
- Revisit: If collaborative multi-user editing demanded.

### ADR-005: Registry-Based Detector Discovery
- Context: Ease of extension.  
- Decision: Decorator registration pattern.  
- Consequences: Low coupling; dynamic listing.  
- Revisit: If plugin distribution externalization required.

## 14.3 Near-Term Roadmap
Scoring: Impact (1–5), Effort (1–5, lower = easier), Priority = Impact * (1 / Effort) heuristic.

| Item | Benefit | Status | Impact | Effort | Priority |
|------|---------|--------|--------|--------|----------|
| Provenance diff CLI | Faster mismatch debugging | Planned | 4 | 2 | 2.0 |
| Export plugin registry | Custom formats | Planned | 3 | 3 | 1.0 |
| GPU distance kernel experiment | Performance spike for large N | Research | 5 | 4 | 1.25 |
| Partial charge integration | Improved scoring fidelity | Backlog | 4 | 5 | 0.8 |

## 14.4 Mid-Term Roadmap
| Item | Benefit | Status | Impact | Effort | Priority |
|------|---------|--------|--------|--------|----------|
| PID adaptive refinement | Stability & faster convergence | Planned | 4 | 3 | 1.33 |
| Async prefetch pipeline | Reduced idle CPU time | Planned | 4 | 4 | 1.0 |
| Parquet columnar export | Large-scale analytics | Planned | 5 | 3 | 1.67 |
| Visualization registry | Reduce initial import cost | Planned | 3 | 3 | 1.0 |

## 14.5 Long-Term Vision
| Theme | Direction |
|-------|-----------|
| Distributed Execution | Multi-node analysis of massive structure sets |
| ML Integration | Learn predictive filters to pre-prune candidates |
| Knowledge Graph | Persist interactions into graph DB for query |
| Interactive Workflows | Notebook integration & API synergy |

## 14.6 Revisit Triggers
| Area | Trigger |
|------|--------|
| KD-tree | Frequent build cost >50% runtime |
| Adaptive thresholds | Unstable oscillation or scientific pushback |
| float32 precision | Reported false negatives due to rounding |
| Streamlit UI | Multi-user concurrency scaling issues |

## 14.7 Deprecation Policy (Draft)
- Announce in docs one minor version before removal.
- Provide compatibility shim where feasible.

## 14.8 Measuring Success
| KPI | Target |
|-----|-------|
| Median analysis latency (standard set) | < 3s |
| Detector parity failure rate | 0 per release |
| Performance regression incidents | < 1 per quarter |
| Documentation gap issues opened | Trending downward |

## 14.9 Feedback Loop
Encourage issue templates: Performance Regression, Detector Accuracy, Documentation Gap, Feature Request—each with reproducibility checklist.

End of architecture suite.
