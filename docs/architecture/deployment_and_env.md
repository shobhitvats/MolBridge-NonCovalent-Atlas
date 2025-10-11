# 13. Deployment & Environment

## 13.1 Environment Profiles
| Profile | Description | Flags Influenced |
|---------|-------------|------------------|
| minimal | Baseline correctness, lowest memory | Disables process pool, normalization |
| auto | Balanced defaults | Enables adaptive moderate thresholds |
| full | Maximum performance | Vector + parallel + shared memory |
| manual | User-specified | Honors explicit environment variables |

## 13.2 Environment Variable Precedence
Explicit env > profile default > code hard default.

## 13.3 Deployment Targets
| Target | Method | Notes |
|--------|-------|-------|
| Local Dev | `streamlit run server.py` | Immediate feedback |
| Docker | Build image, run container | Reproducible; pinned dependencies |
| Cloud PaaS | Use Docker/Procfile | Horizontal scaling with multiple workers |
| Research Cluster | Virtualenv + potentially SLURM job | Batch analyses |

## 13.4 Container Topology
Single container runs both UI and API to simplify environment; can be split for scaling (UI front container + API service) if request volume rises.

## 13.5 Configuration Secrets
Currently minimal (no credentials stored). Future API keys should use environment variables injected at runtime (never committed).

## 13.6 Observability (Future)
Add optional Prometheus endpoint for metrics export; logs structured JSON for ingestion by ELK stack if enabled.

## 13.7 Resource Sizing Guidelines
| Vector | Guideline |
|--------|----------|
| CPU | ≥ 4 cores recommended for parallel profile |
| RAM | ~1 GB per 100k atoms analysis headroom |
| Disk | Cache size scaling with #structures (monitor) |

## 13.8 Failure Domains
| Domain | Failure | Mitigation |
|--------|--------|-----------|
| Network | PDB fetch timeout | Retry/backoff or offline cache |
| Disk | Cache full | Evict LRU / alert |
| CPU Saturation | Too many workers | Cap worker count to logical cores |

## 13.9 Upgrade Process
1. Pull new code.  
2. Rebuild container or reinstall dependencies.  
3. Invalidate adaptive thresholds if criteria changed.  
4. Re-run golden snapshots to verify parity.

## 13.10 Disaster Recovery
- Backup caches optional; recomputation deterministic if inputs preserved.
- Store exported reports + provenance for rehydration.

Proceed to `roadmap_and_decisions.md`.
