# 7. Concurrency & Parallelism

## 7.1 Goals
Exploit multi-core CPUs to reduce wall-clock time while avoiding excessive overhead for small tasks and preserving deterministic aggregation order.

## 7.2 Execution Modes
| Mode | Trigger | Characteristics |
|------|--------|-----------------|
| Sequential | Small N_atoms or few detectors | Lowest overhead |
| Parallel (process pool) | Many detectors or large structures | Multi-core speedup |
| Hybrid (future) | Mix heavy & light detectors | Work-stealing scheduler |

## 7.3 Work Unit Definition
Tuple: `(structure_id, detector_name, parameters_snapshot)`.

## 7.4 Scheduling Strategy
1. Pre-scan detectors to estimate cost (heuristic: expected pair counts, past timings).  
2. Sort heavier tasks earlier to keep workers busy (reduces tail latency).  
3. Submit tasks to pool with bounded queue size to prevent memory blow-up.

## 7.5 Shared Memory Layout
| Segment | Contents | Naming Convention |
|---------|----------|-------------------|
| COORDS | Raw coordinate float32 array | `mb_{struct_hash}_coords` |
| FLAGS | uint16 flags array | `mb_{struct_hash}_flags` |
| RINGS | Centroids + normals concatenated | `mb_{struct_hash}_rings` |

Workers attach via `SharedMemory(name=segment_name)`; construct NumPy view.

## 7.6 Result Aggregation
Results returned as small Python objects (lists of dict + metrics). Aggregator merges preserving configured detector order. Order invariants maintained for reproducibility.

## 7.7 Failure Handling
| Failure | Detection | Response |
|---------|-----------|----------|
| Worker crash | Future exception | Log; rerun sequentially; optionally disable pool |
| Shared memory leak | Segment count mismatch on teardown | Force unlink all matched prefix |
| Timeout (future) | Duration > threshold | Cancel & fallback |

## 7.8 Contention & Overhead Analysis
| Source | Impact | Mitigation |
|--------|--------|-----------|
| Process spawn | Startup latency | Warm pool / persistent processes (future) |
| Large result pickling | CPU overhead | Delay normalization; keep records simple |
| GIL in heavy Python loops | Limits thread scaling | Prefer processes; vectorize loops |

## 7.9 Determinism Strategies
- Stable ordering of detector names.
- Sorting pair indices before record emission.
- Avoid random seeding (no stochastic algorithms currently).

## 7.10 Future Parallel Enhancements
| Idea | Benefit |
|------|--------|
| Async I/O for network fetch stage | Overlap fetch + compute |
| Heterogeneous task cost model | Better load balancing |
| Thread + SIMD combination | Fine-grained CPU utilization |

Proceed to `caching_and_persistence.md`.
