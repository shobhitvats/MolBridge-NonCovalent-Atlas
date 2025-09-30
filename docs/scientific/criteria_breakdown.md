# Criteria Breakdown (Geometric & Energetic Justification)

## Philosophy
Criteria reflect consensus ranges balancing sensitivity (detect true interactions) and specificity (avoid noise). Adjustable bands support exploratory analysis while literature defaults anchor reproducibility. Each section below provides: core parameters, exploration ranges, rationale, and key peer-reviewed sources (hyperlinked DOIs open in new tabs).

## Hydrogen Bonds
Distance <3.5 Å captures the majority of canonical hydrogen bonds; <3.0 Å often indicates strong linear bonds. Angular constraints ensure orbital alignment.

| Parameter | Default | Exploration Range | Rationale | Key Sources |
|-----------|---------|-------------------|-----------|-------------|
| Heavy atom distance | 3.5 Å | 2.6–3.6 | Avoid excessive noise beyond 3.6 | IUPAC 2011; McDonald & Thornton 1994 |
| D–H···A angle | ≥120° | 100–180° | Directionality & overlap | Baker & Hubbard 1984 |
| A–Base angle (optional) | 90–140° | 80–160° | Lone pair geometry approximation | Empirical backbone surveys |
| Strong subclass | ≤2.9 & ≥155° | Tunable | Identify potential catalytic or LBHB | Cleland & Kreevoy 1994 |

## Halogen Bonds
σ-hole driven directionality; heavier halogens (I > Br > Cl >> F) show increased electropositive cap enabling selective design.

| Parameter | Default | Range | Rationale | Sources |
|-----------|---------|-------|-----------|---------|
| X···A distance | ΣvdW + 0.2 Å | ΣvdW+0.0–0.4 | Specific interaction vs close packing | IUPAC 2013 |
| C–X···A angle | ≥150° | 140–180° | Enforces σ-hole alignment | Politzer 2013 |
| Strong subclass | ≥165° | 160–180° | Higher linearity correlates with strength | Cavallo 2016 |

## Chalcogen / Pnictogen / Tetrel Bonds
Analogous anisotropic electrostatics; polarizability increases with atomic number, broadening viable distance slightly beyond hydrogen bonds.

| Parameter | Default | Range | Rationale | Source |
|-----------|---------|-------|-----------|--------|
| Donor–acceptor distance | ΣvdW + 0.5 | +0.3–0.7 | Polarizability tolerance | IUPAC 2023 |
| Axis angle | 115–155° | 110–170° | Captures anisotropic positive potential | IUPAC 2023 |
| Out-of-plane deviation | ≤50° | 30–60° | Maintain axis focus | Survey analyses |

## π–π Stacking
Dispersion + electrostatic quadrupole complementarity; classification by distance, angle, and lateral offset.

| Parameter | Parallel | Offset Parallel | T-shaped | Rationale | Sources |
|-----------|---------|-----------------|---------|-----------|---------|
| Centroid distance | 3.3–4.1 (strong), classify ≤5.5 | ≤6.0 | ≤6.0 | Favorable dispersion zone | Hunter & Sanders 1990; Meyer 2003 |
| Interplanar angle | ≤30° | ≤30° | 60–120° | Orientation family separation | Hunter & Sanders 1990 |
| Lateral offset | <1.5 Å | 1.5–3.0 Å | N/A | Avoid direct electron cloud clash | PDB surveys |

## Cation–π
Electrostatic attraction augmented by aromatic quadrupole; tryptophan rings often dominate.

| Parameter | Default | Range | Rationale | Sources |
|-----------|---------|-------|-----------|---------|
| Centroid distance | 3–6 Å | 2.5–6.5 | Balance attraction/desolvation | Dougherty 1999 |
| Lateral offset | ≤1.5 Å | 0.5–2.0 | Central π density contact | Dougherty 2013 |
| Axial alignment (optional) | cosθ ≥0.6 | 0.5–0.9 | Distinguish axial vs peripheral | Model sets |

## Anion–π
Requires an electron-deficient aromatic face; approach guided by quadrupole axis alignment.

| Parameter | Default | Range | Rationale | Sources |
|-----------|---------|-------|-----------|---------|
| Anion–centroid distance | ≤5.0 Å | 3.5–5.2 | Optimal near 4.3 | Gale 2013 |
| Approach angle (vs normal) | ≤30° | 0–35° | Maximize charge–quadrupole complement | Gale 2013 |
| Electronics filter | optional | substituent aware | Reduce false positives | Review metrics |

## CH–π
Weak donor; cumulative contributions matter in hydrophobic clustering.

| Parameter | Default | Range | Rationale | Sources |
|-----------|---------|-------|-----------|---------|
| C–centroid distance | 2.0–4.8 Å | 2.0–5.2 | Extended tail allowed | Nishio 2014 |
| C–H···centroid angle | ≥90° | 80–180° | Approach orientation | Nishio 2014 |
| Perpendicular offset | ≤2.5 Å | 2.0–3.0 | Ring face engagement | Empirical |

## n→π*
Subtle orbital interaction aligning with Bürgi–Dunitz trajectory.

| Parameter | Default | Range | Rationale | Sources |
|-----------|---------|-------|-----------|---------|
| O→C=O distance | 2.9–3.6 Å | 2.8–3.8 | Overlap zone | Choudhary 2009 |
| Approach angle | 95–140° | 90–145° | Bürgi–Dunitz approximation | Choudhary 2009 |

## Hydrophobic Contacts
Nonpolar heavy atom proximity indicative of packing optimization.

| Parameter | Default | Range | Rationale | Sources |
|-----------|---------|-------|-----------|---------|
| Heavy atom distance | ≤5.0 Å | 4.5–6.0 | Empirical shell | Dill 1990 |
| Element filter | C/S nonpolar | adjustable | Remove polar noise | Heuristic |

## London Dispersion (Heuristic)
Heuristic cluster detection within van der Waals shell to approximate additive dispersion stabilization.

| Parameter | Default | Range | Rationale | Sources |
|-----------|---------|-------|-----------|---------|
| Shell inner radius | 3.3 Å | 3.0–3.5 | Avoid steric region | vdW minima |
| Shell outer radius | 5.0 Å | 4.8–5.5 | Include decaying zone | Wagner 2015 |
| Cluster min size | 3 atoms | 2–5 | Emphasize cumulative effect | Heuristic |

## Sulfur–π
Polarizable sulfur expands permissible distance while maintaining orientation constraints.

| Parameter | Default | Range | Rationale | Sources |
|-----------|---------|-------|-----------|---------|
| S–centroid distance | ≤6.0 Å | 5.0–6.2 | Extended reach | Iwaoka 2015 |
| Perpendicular offset | ≤3.0 Å | 2.0–3.5 | Maintain engagement | Empirical |
| Normal alignment (optional) | cosθ ≥0.3 | 0–0.5 | Orientation quality | Surveys |

## Metal Coordination
Metal–ligand distances and geometry vary by element; generic defaults provide envelope detection, with optional geometry RMSD classification.

| Parameter | Generic Default | Metal-Specific Example | Rationale | Sources |
|-----------|-----------------|-----------------------|-----------|---------|
| Primary shell radius | 2.6 Å | Zn–N/O 2.0–2.3 | Coordinate bond length envelope | Harding 2004 |
| Secondary shell | 3.0 Å | Ca 2.8–3.2 | Extended electrostatic sphere | Zheng 2014 |
| Coordination number window | 4–6 | Zn (4), Mg (6) | Identify geometry class | Harding 2004 |
| Geometry RMSD tolerance | 0.25 Å | 0.20–0.35 | Allow mild distortion | Empirical |

## Notes on Parameter Adjustments
Adjust one dimension at a time; track funnel metrics to detect over-broadening. Record provenance hash before/after exploratory tuning for audit trails.

## Derived Metrics (Future Extensions)
| Metric | Description | Potential Use |
|--------|-------------|---------------|
| Interaction network clustering coefficient | Graph property of residue interaction graph | Detect cooperative motifs |
| Per-residue interaction energy proxy | Weighted sum of interactions | Mutation prioritization |
| Motif recurrence frequency | Frequency of patterns (e.g., stacked triads) | Fold comparison |

## Reference Cross-Mapping Table (Abbrev → Full)
| Abbrev | Full Citation (DOI linked) |
|--------|---------------------------|
| IUPAC 2011 | <a href="https://doi.org/10.1351/PAC-REP-10-01-01" target="_blank" rel="noopener noreferrer">IUPAC Hydrogen Bonding Rec. 2011</a> |
| McDonald & Thornton 1994 | <a href="https://doi.org/10.1021/bi00627a012" target="_blank" rel="noopener noreferrer">Hydrogen bonding patterns</a> |
| Baker & Hubbard 1984 | <a href="https://doi.org/10.1098/rspb.1984.0023" target="_blank" rel="noopener noreferrer">Hydrogen bond geometry</a> |
| Cleland & Kreevoy 1994 | <a href="https://doi.org/10.1021/bi00396a001" target="_blank" rel="noopener noreferrer">Low-barrier hydrogen bonds</a> |
| IUPAC 2013 | <a href="https://doi.org/10.1351/PAC-REC-12-05-10" target="_blank" rel="noopener noreferrer">Halogen bonding</a> |
| Politzer 2013 | <a href="https://doi.org/10.1039/c3cp00054k" target="_blank" rel="noopener noreferrer">σ-hole perspective</a> |
| Cavallo 2016 | <a href="https://doi.org/10.1107/S2052252516012636" target="_blank" rel="noopener noreferrer">Halogen bond survey</a> |
| IUPAC 2023 | <a href="https://doi.org/10.1515/pac-2020-1002" target="_blank" rel="noopener noreferrer">σ-hole classification</a> |
| Meyer 2003 | <a href="https://doi.org/10.1002/anie.200390319" target="_blank" rel="noopener noreferrer">π–π stacking energetics</a> |
| Hunter & Sanders 1990 | <a href="https://doi.org/10.1021/ja00167a008" target="_blank" rel="noopener noreferrer">Aromatic interactions model</a> |
| Dougherty 1999 | <a href="https://doi.org/10.1073/pnas.96.17.9459" target="_blank" rel="noopener noreferrer">Cation–π in proteins</a> |
| Dougherty 2013 | <a href="https://doi.org/10.1021/ar300265y" target="_blank" rel="noopener noreferrer">Cation–π perspective</a> |
| Gale 2013 | <a href="https://doi.org/10.1039/c2cs35251h" target="_blank" rel="noopener noreferrer">Anion–π review</a> |
| Nishio 2014 | <a href="https://doi.org/10.1039/C3CP53729G" target="_blank" rel="noopener noreferrer">CH–π interactions</a> |
| Choudhary 2009 | <a href="https://doi.org/10.1038/nchembio.157" target="_blank" rel="noopener noreferrer">n→π* context</a> |
| Dill 1990 | <a href="https://doi.org/10.1016/0968-0004(90)90045-I" target="_blank" rel="noopener noreferrer">Hydrophobic effect</a> |
| Wagner 2015 | <a href="https://doi.org/10.1002/anie.201500057" target="_blank" rel="noopener noreferrer">London dispersion</a> |
| Iwaoka 2015 | <a href="https://doi.org/10.1002/anie.201411510" target="_blank" rel="noopener noreferrer">Sulfur mediated interactions</a> |
| Harding 2004 | <a href="https://doi.org/10.1107/S0907444903025839" target="_blank" rel="noopener noreferrer">Metal geometry</a> |
| Zheng 2014 | <a href="https://doi.org/10.1038/nchem.1541" target="_blank" rel="noopener noreferrer">Metal-binding analysis</a> |
