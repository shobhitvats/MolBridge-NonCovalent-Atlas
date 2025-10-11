# Scientific Overview Module

This suite complements `SCIENTIFIC_DOCUMENTATION.md` by providing deep, focused analyses. Each module is self-contained and cross-links primary, high-impact literature (DOI hyperlinks open in new tabs). Use this overview as a navigational hub.

## Module Index & Scope Matrix
| File | Scope Focus | Key Additions | Back-Link |
|------|-------------|---------------|-----------|
| `interaction_rationale.md` | Biological & physicochemical roles | Energetic ranges, occurrence notes | Main doc table |
| `criteria_breakdown.md` | Geometric + energetic justification | Literature sourcing per parameter | Threshold section |
| `edge_cases.md` | Ambiguities & classification pitfalls | Resolution strategies & references | Caveats |
| `interpretation_guidelines.md` | Analytical workflows, statistics | Enrichment math, network metrics | Data interpretation |
| `reproducibility_and_provenance.md` | Provenance science | JSON schema & audit recipe | Provenance section |
| `future_directions.md` | Roadmap (scientific) | ML/energy integration refs | Future enhancements |

## Reading Path Recommendations
Beginner: Rationale → Criteria → Interpretation.
Advanced analysis: Criteria → Edge Cases → Reproducibility → Future Directions.

## Global Energetic Reference (Approximate Ranges)
| Interaction | Typical Stabilization (kcal/mol)* | Notes |
|-------------|----------------------------------|-------|
| Hydrogen bond | 1–5 (context-dependent) | Strong short linear <3 Å often >3 kcal/mol (<a href="https://doi.org/10.1021/bi00396a001" target="_blank" rel="noopener noreferrer">Cleland & Kreevoy 1994</a>) |
| Halogen bond | 1–5 | Heavier halogens (I>Br>Cl) stronger (<a href="https://doi.org/10.1039/c3cp00054k" target="_blank" rel="noopener noreferrer">Politzer 2013</a>) |
| Chalcogen / Pnictogen / Tetrel | 0.5–3 | Strength rises with polarizability (<a href="https://doi.org/10.1515/pac-2020-1002" target="_blank" rel="noopener noreferrer">IUPAC 2023</a>) |
| π–π stacking | 1–3 per pair | Aromatic substitution modulates (<a href="https://doi.org/10.1002/anie.200390319" target="_blank" rel="noopener noreferrer">Meyer 2003</a>) |
| Cation–π | 2–6 | High quadrupole rings (Trp>Phe>Tyr) (<a href="https://doi.org/10.1073/pnas.96.17.9459" target="_blank" rel="noopener noreferrer">Dougherty 1999</a>) |
| Anion–π | 0.5–2 | Sensitive to ring electronics (<a href="https://doi.org/10.1039/c2cs35251h" target="_blank" rel="noopener noreferrer">Gale 2013</a>) |
| CH–π | 0.2–1 | Cumulative networks significant (<a href="https://doi.org/10.1039/C3CP53729G" target="_blank" rel="noopener noreferrer">Nishio 2014</a>) |
| n→π* | ~0.5–1 | Subtle local stabilization (<a href="https://doi.org/10.1038/nchembio.157" target="_blank" rel="noopener noreferrer">Choudhary 2009</a>) |
| Hydrophobic (aggregate) | Contextual (surface burial) | Entropic + dispersion (<a href="https://doi.org/10.1016/0968-0004(90)90045-I" target="_blank" rel="noopener noreferrer">Dill 1990</a>) |
| London dispersion (pair) | 0.1–0.5 | Highly additive (<a href="https://doi.org/10.1002/anie.201500057" target="_blank" rel="noopener noreferrer">Wagner 2015</a>) |
| Sulfur–π | 0.5–2 | Polarizable S plays role (<a href="https://doi.org/10.1002/anie.201411510" target="_blank" rel="noopener noreferrer">Iwaoka 2015</a>) |
| Salt bridge | 1–8 (desolvation context) | Environment dependent (<a href="https://doi.org/10.1002/prot.10372" target="_blank" rel="noopener noreferrer">Kumar 2002</a>) |
| Metal coordination | Variable (can be >10) | Not purely noncovalent |

*Approximate; actual free energy depends on environment, desolvation, cooperativity.

Return to main: see `../SCIENTIFIC_DOCUMENTATION.md`.
