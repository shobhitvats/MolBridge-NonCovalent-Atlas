# 5. Performance Engine

## 5.1 Objectives
Deliver near real-time interaction analysis for medium/large proteins while preserving correctness, reproducibility, and maintainability.

## 5.2 Performance Pillars
| Pillar | Techniques |
|--------|-----------|
| Reduce Work | KD-tree pruning, role pre-filtering, lazy derivations |
| Do Work Faster | Vectorization, parallel processes, shared memory |
| Avoid Repeating Work | Multi-layer caching, adaptive thresholds |
| Observe & Tune | Funnel metrics, timing logs, microbench harness |

## 5.3 Profiling Methodology
1. Baseline run with representative structure set.  
2. Collect timing spans (decorators or context managers).  
3. Identify top 2–3 hotspots (typically KD queries or geometry classification).  
4. Apply optimization; re-run microbench; record delta.  
5. Guard with parity tests (ensure outputs unchanged within tolerance).

## 5.4 Timing Instrumentation Pattern
Simplified:
```
with timing("pi_stacking:geometry"):
    compute_geometry()
```
Outputs structured log line with duration ms.

## 5.5 KD-tree Radius Strategy
Compute global maximum needed radius across enabled detectors to build once. Example:
```
max_radius = max(detector.max_prune_radius for detector in active_detectors)
```
This prevents rebuilding tree for differing radii while potentially retrieving slightly more candidates (tradeoff accepted).

## 5.6 Vectorization Strategies
| Pattern | Example | Risk |
|---------|---------|------|
| Broadcasting differences | Distance matrix slabs | High temporary memory |
| Masked indexing | Filtered candidate arrays | Low |
| Pre-normalized normals | Skip per-pair normalization | None (pre-check) |
| Chunked vector ops | Large ring sets | Slight overhead |

## 5.7 Parallel Processing Heuristic
| Condition | Action |
|-----------|--------|
| N_atoms < threshold_small | Keep sequential |
| Detectors > 3 & N_atoms >= threshold_small | Enable process pool |
| Large coordinate arrays | Enable shared memory |

Thresholds empirically derived; environment overrides permitted.

## 5.8 Shared Memory Implementation Sketch
```
from multiprocessing import shared_memory
shm = shared_memory.SharedMemory(create=True, size=coords.nbytes)
backing = np.ndarray(coords.shape, dtype=coords.dtype, buffer=shm.buf)
backing[:] = coords
```
Workers attach by name; no copy performed.

## 5.9 Adaptive Threshold Feedback Loop
Monitors ratio candidate/raw. Out-of-band adjustments occur after detector completes to avoid in-loop instability.

## 5.10 Memory Optimization Tactics
| Tactic | Description |
|--------|-------------|
| float32 coordinates | Halve memory vs float64 |
| On-demand arrays | Build only if active detector requires them |
| Release references | Drop intermediate large temporaries post-use |

## 5.11 Benchmark Dimensions
| Dimension | Metric |
|----------|--------|
| Throughput | interactions/sec |
| Latency | time per structure |
| Memory Peak | MB resident |
| Candidate Efficiency | candidate/raw ratio |
| CPU Utilization | core usage % |

## 5.12 Representative Baseline (Illustrative Template)
| Structure Set | Atoms | Detectors | Latency (s) | Peak Mem (MB) |
|---------------|-------|-----------|-------------|---------------|
| Set A | 25k avg | 8 | 2.3 | 180 |
| Set B | 60k avg | 10 | 5.9 | 410 |

(Actual numbers depend on hardware; fill in via benchmarks.)

## 5.13 Microbenchmark Harness
Focuses on isolated geometry kernel segments (e.g., ring pair offset computation). Helps attribute regressions to specific functions.

## 5.14 Regression Guard Philosophy
A performance regression is only accepted if accompanied by documented correctness gain or feature addition; otherwise blocked in review.

## 5.15 Future Performance Ideas
- Persistent process workers (avoid startup cost)
- Ring centroid spatial grid hashing for π detectors
- Partial GPU offload for distance slabs (CuPy)
- Pre-fetch & streaming analysis for multi-structure batches

Proceed to `adaptive_thresholds.md`.
