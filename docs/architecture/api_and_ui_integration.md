# 10. API & UI Integration

## 10.1 Dual Interface Strategy
- Streamlit UI for interactive exploration.
- FastAPI for programmatic / automated workflows.

## 10.2 Responsibilities Split
| Layer | Responsibility |
|-------|---------------|
| UI Panel | Input capture, parameter adjustments, visualization rendering |
| API Endpoint | Validation, triggering analysis, returning serialized results |
| Core Engine | Stateless computation based on provided parameters |

## 10.3 FastAPI Contract (Illustrative)
```
POST /analyze
{
  "structures": ["1CRN"],
  "detectors": ["hbond", "pi_stacking"],
  "parameters": {"hbond.max_distance": 3.6}
}
→ 200
{
  "results": [... normalized records ...],
  "metrics": {...},
  "provenance": {...}
}
```

## 10.4 Validation via Pydantic
- Ensures detectors exist.
- Ensures parameters within legal ranges.
- Coerces numeric types; rejects unknown fields (strict mode advisable).

## 10.5 Streamlit Parameter Binding
Widget state updates feed a config assembly function producing runtime parameter map used by engine.

## 10.6 Visualization Binding
UI obtains normalized interaction list → groups by type → passes to plotting utilities (e.g., network graph, Ramachandran overlay filters).

## 10.7 Error Surface Strategy
User errors surfaced with gentle warnings (Streamlit `st.warning`), internal errors logged at debug with truncated stack traces to user.

## 10.8 API Versioning Plan
Prefix future breaking changes under `/v2` path; maintain `/v1` stable for deprecation window.

## 10.9 Authentication (Future)
Token-based (Bearer) header for multi-user or cloud deployment; omitted locally for frictionless research usage.

## 10.10 Response Size Management
| Technique | Use Case |
|----------|----------|
| Pagination | Very large result sets |
| Compression (gzip) | Automated pipelines |
| Field filtering | Clients only needing metrics |

Proceed to `extensibility_guide.md`.
