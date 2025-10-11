# Security Policy

## Supported Versions
Active development branch: `main`. Security-related fixes will be applied to the latest release branch when applicable.

## Reporting a Vulnerability
If you discover a potential vulnerability (dependency CVE, unsafe file parsing, injection vector):
1. Do **not** open a public issue with exploit details.
2. Email the maintainer (placeholder: add contact) or open a minimal issue titled `SECURITY DISCLOSURE` indicating you have a report.
3. Provide reproduction steps privately when contacted.

A response target of 5 business days applies for triage acknowledgment.

## Handling
- Triage: validate impact & scope
- Mitigation: patch or apply configuration workaround
- Disclosure: coordinated public disclosure after fix release (unless actively exploited—then accelerated)

## Best Practices Implemented
- Dependency pinning via `requirements.txt`
- Optional vulnerability scan (historically via `pip-audit` in CI)
- Separation of untrusted user input (PDB text) from execution logic (no dynamic code eval)

## Hardening Roadmap
| Area | Planned Improvement |
|------|---------------------|
| Dependency Scanning | Re-introduce automated scan GitHub Action (manual for now) |
| File Validation | Enhanced PDB/CIF sanitization & size limits |
| Rate Limiting | Optional middleware for API job submission |
| Sandbox | Evaluate restricted execution for future plugin system |

## Responsible Disclosure
We value coordinated disclosure. Please allow reasonable time to address issues before publicizing details.
