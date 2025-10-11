# 1. System Overview

## 1.1 Mission Statement
MolBridge automates identification, classification, and reporting of noncovalent interactions in macromolecular structures to accelerate exploratory structural biology and structure‑based design.

## 1.2 Architectural Style
Hybrid layered + component registry pattern:
- Presentation (Streamlit UI, optional REST API)
- Orchestration (Batch / Unified processors, detector registry)
- Domain Logic (Detectors, feature builders, adaptive threshold engine)
- Data Access & Services (Structure loader, caching, provenance, normalization)
- Cross‑Cutting Concerns (Logging, metrics, configuration)

## 1.3 Core Principles
| Principle | Description | Example Enforcement |
|----------|-------------|---------------------|
| Idempotence | Same inputs → same outputs (unless adaptive tuning explicitly enabled) | Provenance hash & normalization order |
| Progressive Optimization | Start correct → add vectorization/parallelism behind flags | Vector path gated by `MOLBRIDGE_ENABLE_VECTOR_GEOM` |
| Observability | Every stage emits measurable metadata | Funnel metrics & timing decorators |
| Extensibility | Add detectors without touching core orchestrator | Detector registry + contract |
| Reproducibility | Ability to reconstruct published result | Config snapshot + provenance hash |

## 1.4 High-Level Data Flow
```
Input (PDB IDs / files) → Structure Loader → Feature Extraction → Detector Dispatch Loop
  → (Per detector: Prune → Filter → Classify → Record) → Aggregation → Normalization → Reporting/Visualization/Export
```

## 1.5 Interaction Classes Supported
Hydrogen bonds, π–π stacking, cation–π, anion–π, CH–π, sulfur–π, sulfur-mediated chalcogen/pnictogen/tetrel bonds, halogen bonds, ionic interactions, salt bridges, hydrophobic contacts, London dispersion approximations, n→π* interactions, metal coordination.

## 1.6 Evolution Overview
Milestones:
| Phase | Focus | Result |
|-------|-------|--------|
| 0.x Prototype | Correctness & breadth | Baseline detector set |
| 0.9 Pre‑Release | Performance & normalization | Vector paths + canonical schema |
| 1.0 Foundation | Reproducibility & metrics | Provenance, funnel metrics |
| 1.1 Expansion | Deep documentation | Architecture suite | 

## 1.7 Non-Goals
| Non-Goal | Rationale |
|----------|-----------|
| Full MD trajectory analysis | Scope limited to static structures |
| Quantum mechanical scoring | Out of performance envelope |
| In-browser 3D editing | Focus remains analysis not modeling |

## 1.8 Constraints & Assumptions
- Python 3.11+ baseline (pattern matching, performance gains).
- Typical structure sizes ≤ 150k atoms single structure (RAM-based heuristics tuned here).
- Sufficient disk space for caching (configurable path).
- Network access required for remote PDB fetch unless cached.

## 1.9 Quality Attribute Priorities
| Attribute | Priority (1–5) | Notes |
|-----------|----------------|-------|
| Correctness | 5 | Scientific validity is foundational |
| Maintainability | 4 | Detector modularity & docs reduce friction |
| Performance | 4 | Large structures actionable interactively |
| Reproducibility | 4 | Provenance + config capture |
| Extensibility | 5 | New interaction paradigms anticipated |
| Usability | 3 | Research audience tolerant of complexity |
| Security | 3 | Local/offline use; basic hygiene maintained |

## 1.10 Reading Guide
Proceed next to `data_model.md` for concrete representations before diving into detectors.
