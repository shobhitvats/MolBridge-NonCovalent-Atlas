# Interpretation Guidelines

## Practical Analysis Workflow
1. Run default criteria to establish baseline interaction landscape.
2. Examine funnel metrics for each interaction; unexpectedly low acceptance may indicate over-strict filters.
3. Adjust a single parameter axis (e.g., distance) at a time to understand sensitivity.
4. Compare normalized counts across homologous structures or mutant vs wild-type.
5. Focus on interaction clusters near active or binding sites for functional inference.

## Enrichment Analysis
- Define region of interest (active site) vs background (rest of structure).
- Compute per-interaction-type density (count per 100 residues) in region vs background.
- Enrichment ratio = region_density / background_density; values >1.5 may indicate functional specialization.
Formula:
```
ER_type = ( Count_type_region / Residues_region ) / ( Count_type_background / Residues_background )
```
Confidence estimation (approximate): treat counts as Poisson; compute 95% CI via sqrt(n) and propagate ratio (delta method) for preliminary significance screen (<a href="https://doi.org/10.1038/nmeth.2613" target="_blank" rel="noopener noreferrer">Basic statistical guidance</a>).

## Network Perspective
Construct graph where nodes are residues and edges interactions; analyze:
- Degree distribution (hubs may stabilize folds)
- Betweenness (bridging residues in communication pathways)
- Community detection for interaction modules.
Optional tools: apply Louvain or Leiden algorithms on residue interaction graph (external libraries). Weighted edges by interaction type multiplicity or strength classification.

## Structural Quality Checks
| Indicator | Potential Issue | Action |
|-----------|-----------------|-------|
| Many long hydrogen bonds | Low resolution or generous distance cutoff | Tighten distance |
| Sparse hydrophobic contacts in core | Misparsed chains or missing residues | Verify input structure |
| Excess cation–π vs literature norms | Parameter offset too wide | Reduce lateral offset cutoff |

## Comparative Mutational Assessment
Track lost or gained interactions after mutation modeling to hypothesize stability or affinity effects.
ΔInteractionScore (prototype): sum(weight_per_type * delta_count). Calibrate weights empirically (e.g., hydrogen bond 1.0, cation–π 1.2, π–π 0.9, hydrophobic 0.3). Not an energy; heuristic prioritization aid.

## Reporting Best Practices
- Include provenance hash and parameter table in supplementary materials.
- Provide both raw and normalized counts for transparency.
- Use consistent color coding across figures (e.g., hydrogen bonds always blue).

## Limitations Reminder
Computational criteria do not guarantee energetic favorability; experimental corroboration (mutagenesis, crystallography at higher resolution) recommended for critical claims.
Add parallel validation with MD or QM cluster calculations for ambiguous or borderline high-impact interactions (<a href="https://doi.org/10.1021/cr200106v" target="_blank" rel="noopener noreferrer">Desiraju 2011</a> for hydrogen bond conceptual refinements).
