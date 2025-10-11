# Contributing to MolBridge

Thank you for your interest in improving MolBridge! This document outlines the end-to-end workflow for contributing code, documentation, or scientific criteria.

## 1. Code of Conduct
All participants are expected to uphold respectful, inclusive collaboration. (A `CODE_OF_CONDUCT.md` will be added; until then, follow the principles of the Contributor Covenant.)

## 2. Issue First (Strongly Recommended)
Open an issue (or GitHub Discussion) before large changes:
- Type: Bug / Feature / Performance / Documentation / Scientific Criteria / Refactor
- Include: Motivation, minimal reproduction (if bug), acceptance criteria, performance implications, references (for criteria changes)

## 3. Development Environment
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
# (Optional) dev extras if defined later
```

## 4. Branching & Commits
| Aspect | Guideline |
|--------|-----------|
| Branch naming | `feature/<slug>` or `fix/<slug>` or `docs/<slug>` |
| Commit style | Imperative present tense ("Add detector", "Refactor cache layer") |
| Granularity | One logical change per commit; tests + code together |

## 5. Adding a New Detector (Checklist)
1. Draft geometric / physicochemical rationale with citations.
2. Implement detector in `src/analysis/`.
3. Register via existing registry decorator.
4. Add parameters + defaults (config + README parameter glossary update).
5. Add tests:
   - Unit geometry edge cases
   - Parity test (if alternate path exists)
   - Performance sanity (counts within expected range on curated PDB sample)
6. Update `docs/scientific/criteria_breakdown.md` + rationale file.
7. Add to interaction overview tables (README condensed + scientific summary).
8. Run full test suite `pytest -v`.
9. Include benchmark notes (optional) + any adaptive threshold considerations.

## 6. Modifying Criteria / Defaults
Provide:
- Before vs after values
- Literature citations (DOI preferred)
- Rationale (precision vs recall impact)
- If backward incompatible, note in PR description (will trigger CHANGELOG entry)

## 7. Performance Contributions
If optimizing:
- Capture baseline timings (microbench + representative real PDBs)
- Report percent change (mean + p95)
- Confirm no change in interaction counts unless explicitly intended
- Flag new environment variables in docs

## 8. Documentation Improvements
- Keep top-level README concise; deep detail belongs in `docs/`
- Add new architecture deep dives under `docs/architecture/`
- Add new scientific elaborations under `docs/scientific/`
- Use tables for parameter matrices; keep lines < 120 chars

## 9. Testing
```bash
pytest -v
```
Optional focused runs:
```bash
pytest tests/test_hbond_subtypes.py -k strong
```
Add tests near existing patterns; avoid creating duplicate utility helpers without consolidation.

## 10. Style & Lint (Future Placeholder)
A formatting / linting config may be introduced (e.g., Ruff, Black). Until formalized:
- Follow existing code style conventions in edited modules
- Keep imports grouped (stdlib, third-party, local)

## 11. Opening a Pull Request
Checklist for your PR description:
- Summary of change
- Linked issue number(s)
- Test evidence (paste microbench delta if perf change)
- Documentation updated? (Y/N + files)
- Backward compatibility notes
- Screenshots (UI changes)

## 12. Review Expectations
- At least one maintainer approval required
- Address review comments with follow-up commits (avoid force-push unless requested)
- Squash commits on merge if many granular WIP commits

## 13. Release & Versioning
Semantic-ish versioning (major.minor.patch). Detector additions or parameter default changes generally bump minor; regression fixes bump patch.

### Dependency Lock Workflow
- Edit top-level specs only in `requirements.txt`.
- Regenerate the frozen snapshot when intentionally updating: 
   ```bash
   bash scripts/update_lock.sh
   ```
- Commit the resulting `requirements.lock` in the same PR with a short note (e.g., "Update lock after adding X").
- Do not hand-edit `requirements.lock`.

## 14. Security / Disclosure
For potential security issues (dependency vulnerability, unsafe file handling) contact the maintainer privately first (placeholder: add email or SECURITY.md).

## 15. Attribution
Add yourself to `PROJECT_SUMMARY.txt` or forthcoming `AUTHORS.md` once a PR is merged.

Thanks for contributing! 🚀
