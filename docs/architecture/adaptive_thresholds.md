# 6. Adaptive Thresholds Engine

## 6.1 Purpose
Dynamically tune pruning windows (distances, sometimes angles) to balance speed and sensitivity without manual retuning for varying structure sizes and densities.

## 6.2 Control Theory Analogy
Acts like a coarse proportional controller:
```
error = observed_density - target_density
adjustment = -k_p * error
```
Where observed_density = candidate_pairs / max(raw_pairs,1). Clamped within safety bounds.

## 6.3 Parameters
| Name | Description | Default |
|------|-------------|---------|
| target_low | Lower acceptable density | 0.05 |
| target_high | Upper acceptable density | 0.20 |
| k_p | Proportional scalar | 0.5 |
| max_step | Max distance change / iteration (Å) | 0.15 |
| angle_k_p | Angle adjustment proportional scalar | 0.3 |

## 6.4 State Persistence
Optionally stored JSON keyed by detector name:
```
{
  "pi_stacking": {"distance": 6.0, "updated": "2025-09-30T12:00:00Z"},
  "hbond": {"distance": 3.4}
}
```
Enabled via `MOLBRIDGE_ADAPTIVE_CACHE_PATH`.

## 6.5 Update Algorithm Pseudocode
```
obs = candidate_pairs / max(raw_pairs, 1)
if obs < target_low:
    distance += min(max_step, k_p*(target_low - obs) * base_window)
elif obs > target_high:
    distance -= min(max_step, k_p*(obs - target_high) * base_window)
clamp(distance, min_window, max_window)
```
Angles (if adaptive) similar with degree clamps.

## 6.6 Stability Considerations
| Issue | Mitigation |
|-------|-----------|
| Oscillation | Use asymmetric deadband (target_low, target_high) |
| Overshoot | Cap max_step |
| Drift from literature | Clamp to physical min/max |
| Multi-detector interference | Independent state per detector |

## 6.7 Telemetry Fields
Recorded in metrics snapshot under `threshold_snapshot` key:
```
{
  "distance": 6.1,
  "prev_distance": 6.0,
  "obs_density": 0.07,
  "action": "+0.1"
}
```

## 6.8 Failure Modes
| Mode | Cause | Fallback |
|------|------|----------|
| JSON corruption | External edit | Reinitialize defaults |
| Path unwritable | Permissions | Disable persistence; warn |
| Extreme density zero | Pathological prune | Increase distance aggressively (double until min raw pairs) |

## 6.9 When to Disable
- Benchmarking reproducibility runs.
- Scientific validation requiring fixed criteria.

## 6.10 Future Enhancements
| Idea | Benefit |
|------|--------|
| Exponential moving average smoothing | Reduce jitter |
| PID (add integral term) | Better steady-state accuracy |
| Detector coupling awareness | Coordinated global target budget |

Proceed to `concurrency_and_parallelism.md`.
