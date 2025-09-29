# Reference & Asset Migration Guide

This document tracks the plan to relocate large binary reference materials (PDF, PPTX) out of the core application tree for a cleaner repository and faster installs.

## Current Locations

```
Protein-Interaction-Analysis-Server/IUPAC/*.pdf
Protein-Interaction-Analysis-Server/Other References/*.pdf
protein_interaction_presentation_1_structures.pptx (root)
```

## Rationale
- Large binaries bloat source distribution and slow CI/CD.
- References change infrequently and are not executed as code.
- Clear separation improves license auditing and packaging.

## Target Structure
```
docs/
  references/
    iupac/
    literature/
    presentations/
```

## Proposed Moves
| Source | Destination | Notes |
|--------|-------------|-------|
| IUPAC/*.pdf | docs/references/iupac/ | Keep original filenames |
| Other References/*.pdf | docs/references/literature/ | Group by topic if desired |
| protein_interaction_presentation_1_structures.pptx | docs/references/presentations/ | Rename with date if multiple |

## Git LFS (Optional)
If these documents must remain in repo, consider enabling Git LFS:
1. Install git-lfs locally.
2. `git lfs install`
3. `git lfs track "*.pdf" "*.pptx"`
4. Commit updated `.gitattributes`.

## Implementation Steps (Manual Outside This Automation)
1. Create target folders:
   - `mkdir -p docs/references/{iupac,literature,presentations}`
2. Move files:
   - `git mv IUPAC/*.pdf docs/references/iupac/`
   - `git mv "Other References"/*.pdf docs/references/literature/`
   - `git mv protein_interaction_presentation_1_structures.pptx docs/references/presentations/`
3. Update any code or README links (search for filenames).
4. Commit changes with message: "ref(docs): relocate reference PDFs & presentation".

## Post-Migration Cleanup
- Verify no stale empty directories remain (`IUPAC/`, `Other References/`).
- Add a README in each new directory if contextual notes are useful.

## Future Enhancements
- Add citation metadata (BibTeX) in `docs/references/metadata/`.
- Provide a summarized bibliography section in main README.

---
Maintainer Note: This guide is generated as part of automated cleanup; adjust as needed for internal policies.
