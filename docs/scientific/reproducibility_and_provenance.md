# Reproducibility & Provenance (Scientific Framing)

## Motivation
Ensuring that derived interaction landscapes can be recreated underpins credible comparative studies and publications.

## Components Captured
| Component | Rationale |
|-----------|-----------|
| Structure signature | Guarantees atom ordering & coordinates unchanged |
| Detector list | Confirms analytical scope |
| Parameter set | Ensures thresholds identical |
| Application version | Ties results to implementation state |
| (Optional) Adaptive state disabled | Stabilizes thresholds |

## Recommended Reproducibility Mode
Set environment: `MOLBRIDGE_ENABLE_PROVENANCE=1` and disable adaptive thresholds for final publication analysis. Archive exported JSON + provenance block.
Sample provenance JSON snippet:
```json
{
	"provenance_hash": "a1b2c3d4e5f6",
	"detectors": ["hbond", "pi_stacking"],
	"parameters": {"hbond.max_distance": 3.5, "pi_stacking.max_distance": 5.5},
	"structure_signature": "sha256:9f...",
	"version": "1.1.0",
	"generated_utc": "2025-09-30T13:05:22Z"
}
```

## Audit Procedure
1. Extract provenance hash from original dataset.
2. Re-run with identical parameters; obtain new hash.
3. Compare interaction counts & random spot-check geometry metrics.
4. Document any intentional differences (e.g., updated detector algorithm) in changelog.

## Handling Evolving Criteria
When criteria updated (e.g., refined halogen angle), dual reporting of counts under old vs new thresholds during transition window improves continuity.
Consider archiving a “criteria manifest” file enumerating each detector → parameter → reference DOI. This enables transparent peer review and rapid replication audits.

Further reading: reproducible computational chemistry practices (<a href="https://doi.org/10.1021/acs.jcim.0c01102" target="_blank" rel="noopener noreferrer">Best Practices JCIM 2021</a>).
