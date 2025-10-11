## MolBridge Scientific Documentation

> Version: 2025-09 (synchronized with Technical Documentation)

### Modular Scientific Suite Index
For focused deep dives, the scientific content has been decomposed into topic modules (paths relative to `docs/`):

| # | Module | Focus |
|---|--------|-------|
| 1 | [scientific/interaction_rationale.md](scientific/interaction_rationale.md) | Biological & physicochemical roles of each interaction family |
| 2 | [scientific/criteria_breakdown.md](scientific/criteria_breakdown.md) | Justification of geometric thresholds & energetic context |
| 3 | [scientific/edge_cases.md](scientific/edge_cases.md) | Ambiguities, borderline scenarios, classification caveats |
| 4 | [scientific/interpretation_guidelines.md](scientific/interpretation_guidelines.md) | Practical usage, enrichment analysis, network perspective |
| 5 | [scientific/reproducibility_and_provenance.md](scientific/reproducibility_and_provenance.md) | Provenance framing, audit procedures |
| 6 | [scientific/future_directions.md](scientific/future_directions.md) | Planned scientific enhancements & research directions |

Use these modules when you need additional nuance without scanning the entire consolidated document below.

### Purpose
MolBridge provides reproducible, literature-grounded identification of noncovalent interactions across protein (and protein–ligand) structures for structural biology, computational chemistry, and design studies.

### Interaction Classes & Rationale
| Class | Biophysical Role | Core Criteria (Simplified) | Ref |
|-------|------------------|----------------------------|-----|
| Hydrogen Bond | Stabilizes secondary & tertiary structure; specificity | Donor–acceptor heavy atom ≤3.5 Å; D–H···A ≥120° | <a href="https://doi.org/10.1351/PAC-REP-10-01-01" target="_blank" rel="noopener noreferrer">1</a> |
| Halogen Bond | Directional σ-hole interaction; ligand optimization | X···A ≤ Σ(vdW) (+0.2 Å slack); C–X···A ≥150° | <a href="https://doi.org/10.1351/PAC-REC-12-05-10" target="_blank" rel="noopener noreferrer">2</a> |
| Chalcogen / Pnictogen / Tetrel | Emerging σ-hole families; packing & recognition | Donor element–acceptor ≤(ΣvdW+0.5); axis angle window | <a href="https://doi.org/10.1515/pac-2020-1002" target="_blank" rel="noopener noreferrer">3</a>, <a href="https://doi.org/10.1017/qrd.2023.3" target="_blank" rel="noopener noreferrer">4</a> |
| π-π Stacking | Aromatic stabilization, interface packing | Centroid ≤5.5 Å; interplanar ≤30° (parallel) / 60–120° (T) | <a href="https://doi.org/10.1021/ja00167a008" target="_blank" rel="noopener noreferrer">7</a> |
| Cation–π | Electrostatic + quadrupole stabilization | Centroid distance 3–6 Å; lateral offset ≤1.5 Å | <a href="https://doi.org/10.1073/pnas.96.17.9459" target="_blank" rel="noopener noreferrer">6</a> |
| Anion–π | Charge–quadrupole modulation | Anion–centroid ≤5.0 Å; approach vs ring normal ≤30° | <a href="https://doi.org/10.1021/acsomega.3c00205" target="_blank" rel="noopener noreferrer">5</a> |
| n→π* | Local backbone stabilization | O(i)→C(i+1) 2.9–3.6 Å; approach 95–140° | <a href="https://doi.org/10.1021/acsomega.3c00205" target="_blank" rel="noopener noreferrer">5</a> |
| Salt Bridge | Strong electrostatic pair | Heavy atom ≤4.0 Å; sidechain centroid ≤6.0 Å | <a href="https://doi.org/10.1021/acsomega.3c00205" target="_blank" rel="noopener noreferrer">5</a> |
| Ionic (general) | Broader electrostatic landscape | Charge centroid distance cutoff | <a href="https://doi.org/10.1021/acsomega.3c00205" target="_blank" rel="noopener noreferrer">5</a> |
| Hydrophobic | Packing & desolvation | C···C ≤5.0 Å (explore up to 5.5–6.0) | <a href="https://doi.org/10.1021/acsomega.3c00205" target="_blank" rel="noopener noreferrer">5</a> |
| London Dispersion | Weak cumulative packing | vdW shell band 3.3–5.0 Å cluster | <a href="https://doi.org/10.1021/acsomega.3c00205" target="_blank" rel="noopener noreferrer">5</a> |
| Sulfur–π | Polarizable sulfur orientation | S–centroid ≤6.0 Å; offset ≤3.0 Å | <a href="https://doi.org/10.1021/acsomega.3c00205" target="_blank" rel="noopener noreferrer">5</a> |
| Metal Coordination | Catalysis / structure | Metal–ligand primary ≤2.6 Å (element specific) | <a href="https://doi.org/10.1107/S0907444903025839" target="_blank" rel="noopener noreferrer">Harding 2004</a>, <a href="https://doi.org/10.1038/nchem.1541" target="_blank" rel="noopener noreferrer">Zheng 2014</a> |

<sub>Clicking a reference number opens the source in a new browser tab.</sub>

### Additional Context & Foundational Literature
The table above lists a core curated reference per class (or cluster). For deeper study and historical evolution of criteria, below are supplementary highly cited sources (retaining original numbered references unchanged):

* Hydrogen Bonds (historical & geometric refinements): <a href="https://doi.org/10.1098/rspb.1984.0023" target="_blank" rel="noopener noreferrer">Baker & Hubbard 1984</a>, <a href="https://doi.org/10.1021/bi00627a012" target="_blank" rel="noopener noreferrer">McDonald & Thornton 1994</a>, <a href="https://doi.org/10.1021/cr200106v" target="_blank" rel="noopener noreferrer">Desiraju 2011 Review</a>.
* Halogen / σ-Hole Interactions: <a href="https://doi.org/10.1039/c3cp00054k" target="_blank" rel="noopener noreferrer">Politzer & Murray 2013</a>, <a href="https://doi.org/10.1107/S2052252516012636" target="_blank" rel="noopener noreferrer">Cavallo et al. 2016</a>.
* Chalcogen / Pnictogen / Tetrel: <a href="https://doi.org/10.1021/acs.chemrev.5b00527" target="_blank" rel="noopener noreferrer">Mahmudov et al. Chem. Rev. 2017</a>.
* π–π Stacking Foundations: <a href="https://doi.org/10.1002/anie.200390319" target="_blank" rel="noopener noreferrer">Meyer et al. 2003</a>, <a href="https://doi.org/10.1039/C2SC00644C" target="_blank" rel="noopener noreferrer">Martinez & Iverson 2012</a>.
* Cation–π: (classic already cited) plus <a href="https://doi.org/10.1021/ar300265y" target="_blank" rel="noopener noreferrer">Dougherty 2013 Perspective</a>.
* Anion–π: <a href="https://doi.org/10.1021/cr400273w" target="_blank" rel="noopener noreferrer">Garcia-Raso & co. (Comprehensive Review)</a>, <a href="https://doi.org/10.1039/c2cs35251h" target="_blank" rel="noopener noreferrer">Gale et al. 2013</a>.
* CH–π: <a href="https://doi.org/10.1039/C3CP53729G" target="_blank" rel="noopener noreferrer">Nishio 2014 Review</a>.
* n→π*: <a href="https://doi.org/10.1038/nchembio.157" target="_blank" rel="noopener noreferrer">Choudhary et al. 2009</a>, <a href="https://doi.org/10.1038/nchembio.124" target="_blank" rel="noopener noreferrer">Raines Group Insights</a>.
* Hydrophobic Effect: <a href="https://doi.org/10.1002/bip.360010501" target="_blank" rel="noopener noreferrer">Kauzmann 1959</a>, <a href="https://doi.org/10.1016/0968-0004(90)90045-I" target="_blank" rel="noopener noreferrer">Dill 1990</a>.
* London Dispersion in Biology: <a href="https://doi.org/10.1002/anie.201500057" target="_blank" rel="noopener noreferrer">Wagner & Schreiner 2015</a>.
* Sulfur–π & Sulfur-mediated: <a href="https://doi.org/10.1002/anie.201411510" target="_blank" rel="noopener noreferrer">Iwaoka & Takemoto (survey)</a>.
* Salt Bridges & Electrostatics: <a href="https://doi.org/10.1002/prot.10372" target="_blank" rel="noopener noreferrer">Kumar & Nussinov 2002</a>.
* Metal Coordination Geometry: <a href="https://doi.org/10.1107/S0907444903025839" target="_blank" rel="noopener noreferrer">Harding 2004</a>, <a href="https://doi.org/10.1038/nchem.1541" target="_blank" rel="noopener noreferrer">Zheng et al. 2014</a>.


### Parameterization Sources
Parameters anchored to peer-reviewed ranges (see README references) with adjustable sliders permitting controlled exploration. Each detection module exposes defaults aligned to frequently cited midpoints.

### Hydrogen Bond Subtypes
Subtype classification distinguishes structural roles (backbone vs sidechain, water mediation). This enables secondary analyses: fold stabilization mapping, interface enrichment studies.
Additional considerations: Bifurcated hydrogen bonds and low-barrier (short, highly linear) bonds can signal catalytic residues (<a href="https://doi.org/10.1021/bi00396a001" target="_blank" rel="noopener noreferrer">Cleland & Kreevoy 1994</a>). Backbone-mediated networks propagate stability across β-sheets.

### Provenance & Reproducibility
Optional provenance hashing captures:
* Structure digest (sorted residue + atom signatures)
* Parameter signature
* Detector set
Use this in manuscripts to assert exact analytic conditions.

### Data Interpretation Guidance
* High hydrogen bond counts with low average angle may indicate relaxed angle cutoff or low-resolution structure.
* Elevated hydrophobic contact density can correlate with compact interface or artifactually high B-factor smoothing.
* Cation–π enrichment near active sites often signals ligand or catalytic stabilization.
* Outlier metal coordination distances should be verified for alternate conformations / occupancy issues.

### Caveats
* Geometry-only heuristics (no QM correction) – borderline interactions may require manual review.
* Missing hydrogens may shift certain donor classifications; current method infers from heavy atom templates.
* Alternate location indicators: first conformer retained; future enhancement may weight occupancy.

### Recommended Workflow For Publication-Grade Results
1. Use Literature Default preset.
2. Inspect hydrogen bond subtype distribution; adjust angle if inconsistent with expected secondary structure content.
3. Validate metal coordination distances vs known ion-specific radii.
4. Export provenance-enabled normalized dataset.
5. Archive configuration file + commit SHA.

### Statistical Summaries
MolBridge provides acceptance ratios and funnel metrics enabling quality control: abrupt drops in acceptance can imply overly strict angular filters or coordinate anomalies.
Energy perspective (approximate, context-dependent typical magnitudes): Hydrogen bonds ~1–5 kcal/mol; cation–π similar order; π–π stacking 1–3 kcal/mol per pair; individual CH–π and dispersion contacts <1 kcal/mol but cumulatively significant. (Representative aggregate studies: <a href="https://doi.org/10.1039/C3CP00054K" target="_blank" rel="noopener noreferrer">Politzer Review</a>, <a href="https://doi.org/10.1002/anie.200390319" target="_blank" rel="noopener noreferrer">Meyer 2003</a>.)

### Future Scientific Enhancements
* Ligand-specific interaction enrichment analysis
* Interface vs core differential interaction profiling
* Cooperative network motifs (triads, ladders)
* Energetic scoring overlays (empirical potentials)
* Empirical energy calibration via knowledge-based potentials (<a href="https://doi.org/10.1002/prot.20033" target="_blank" rel="noopener noreferrer">Zhou & Zhou 2002</a>)
* Integration with high-accuracy structure prediction frameworks (<a href="https://doi.org/10.1038/s41586-021-03819-2" target="_blank" rel="noopener noreferrer">AlphaFold2</a>) for comparative evolutionary interaction conservation.

---

## Detailed Detection Criteria (Full Specification)

This section mirrors and expands the in-application Info tab with explicit default thresholds, recommended exploration ranges, and literature anchors. All angles are in degrees (°) and distances in Å unless otherwise noted. Where vdW radii sums are referenced, a small tolerance (+0.2–0.3 Å) is often permitted to accommodate crystallographic uncertainty and thermal motion.

### 1. Hydrogen Bonds
**Defaults:** Donor–Acceptor heavy atom distance ≤ 3.5; D–H···A angle ≥ 120.  
**Exploration Range:** Distance 2.6–3.6; Angle 100–180 (tighten ≥150 for high specificity).  
**Notes:** Heavy-atom proxy geometry used if hydrogens absent. Low-barrier variants (short, near-linear) are not separately subclassed yet.  
**Reference:** [<a href="https://doi.org/10.1351/PAC-REP-10-01-01" target="_blank" rel="noopener noreferrer">1</a>]

### 2. Halogen Bonds
**Defaults:** X···A ≤ Σ(vdW radii) + 0.2; C–X···A ≥ 150 (soft lower bound); preferred ~165–180.  
**Directionality:** σ-hole alignment enforces near-linear C–X···A axis.  
**Reference:** [<a href="https://doi.org/10.1351/PAC-REC-12-05-10" target="_blank" rel="noopener noreferrer">2</a>]

### 3. Chalcogen Bonds (S / Se / Te)
**Defaults:** S/Se/Te···A ≤ 4.0 (typical strong 3.2–3.8); C–S···A (or substituent axis) 115–155; out-of-plane deviation |φ| ≤ 50.  
**Rationale:** σ-hole opposite covalent bond(s) yields anisotropic positive electrostatic potential. Tighter angle windows increase precision.  
**Reference:** [<a href="https://doi.org/10.1515/pac-2020-1002" target="_blank" rel="noopener noreferrer">3</a>], [<a href="https://doi.org/10.1017/qrd.2023.3" target="_blank" rel="noopener noreferrer">4</a>]

### 4. Pnictogen Bonds (N / P / As)
**Defaults:** Donor pnictogen–acceptor ≤ Σ(vdW) + 0.4; alignment ≥ 150.  
**Reference:** [<a href="https://doi.org/10.1515/pac-2020-1002" target="_blank" rel="noopener noreferrer">3</a>]

### 5. Tetrel Bonds (C / Si / Ge)
**Defaults:** C/Si/Ge···A ≤ Σ(vdW) + 0.5 (protein context: carbon-centered cases rare); approach (σ-hole axis) ≥ 160 for stringent classification (≥150 relaxed).  
**Reference:** [<a href="https://doi.org/10.1515/pac-2020-1002" target="_blank" rel="noopener noreferrer">3</a>]

### 6. π–π Stacking
**Parallel (Face-to-Face):** Centroid distance 3.3–4.1 strongly favored; classify ≤5.5. Interplanar angle ≤ 30.  
**T-shaped / Edge-to-Face:** Centroid distance ≤ 6.0; ring normal angle ~60–120 (center above ring rim).  
**Reference:** [<a href="https://doi.org/10.1021/ja00167a008" target="_blank" rel="noopener noreferrer">7</a>]

### 7. C–H···π Interactions
**Defaults:** C (donor carbon) projected to ring centroid 2.0–4.8; C–H···centroid angle ≥ 90; perpendicular offset ≤ 2.5.  
**Reference:** [<a href="https://doi.org/10.1021/acsomega.3c00205" target="_blank" rel="noopener noreferrer">5</a>]

### 8. Cation–π Interactions
**Defaults:** Cation (sidechain charged group centroid) to aromatic centroid 3.0–6.0 (optimum 4.0–5.0); lateral offset (distance from ring normal projection) ≤ 1.5.  
**Note:** Metal cations (e.g., Zn) can be optionally included depending on parameter flags.  
**Reference:** [<a href="https://doi.org/10.1073/pnas.96.17.9459" target="_blank" rel="noopener noreferrer">6</a>]

### 9. Anion–π Interactions
**Defaults:** Anion center to ring centroid ≤ 5.0 (enrichment <4.3); angle between anion→centroid vector and ring normal ≤ 30 (approach along positive quadrupole axis).  
**Reference:** [<a href="https://doi.org/10.1021/acsomega.3c00205" target="_blank" rel="noopener noreferrer">5</a>]

### 10. n→π* Interactions
**Defaults:** Carbonyl O(i) → C(i+1) distance 2.9–3.6 (cutoff ≤3.5 default); approach (proxy for lone pair to π* orientation) 95–140.  
**Reference:** [<a href="https://doi.org/10.1021/acsomega.3c00205" target="_blank" rel="noopener noreferrer">5</a>]

### 11. Hydrophobic Contacts
**Defaults:** Nonpolar sidechain heavy atom (sp3/sp2 carbon) pair distance ≤ 5.0 (exploratory expansion 5.5–6.0).  
**Note:** Counts scale with surface burial; interpret alongside SASA if available.

### 12. London Dispersion (Heuristic)
**Approach:** Identify nonpolar atom pairs in vdW window shell 3.3–5.0 and cluster to emphasize cumulative packing rather than isolated geometry.  
**Note:** Not an energy calculation; heuristic signal only.

### 13. Sulfur–π Interactions
**Defaults:** Sulfur (CYS SG / MET SD) to aromatic centroid ≤ 6.0 (strong <5.0); perpendicular offset ≤ 3.0; projection onto ring normal near center enhances classification confidence.  
**Reference:** [<a href="https://doi.org/10.1021/acsomega.3c00205" target="_blank" rel="noopener noreferrer">5</a>]

### 14. Salt Bridges
**Defaults:** Opposite formal charge sidechain heavy atom pair ≤ 4.0 (e.g., NH₁/₂/₃ to carboxylate O); centroid (positive group vs carboxylate) ≤ 6.0.  
**Note:** Histidine considered protonated only under selected protonation assumptions; adjustable by parameter flags.  
**Reference:** [<a href="https://doi.org/10.1021/acsomega.3c00205" target="_blank" rel="noopener noreferrer">5</a>]

### 15. Ionic (General Electrostatics)
**Defaults:** Charged group centroids ≤ 6.5 (looser) or ≤ 5.5 (stricter).  
**Use Case:** Broader charge network mapping distinct from strict salt bridges.  
**Reference:** [<a href="https://doi.org/10.1021/acsomega.3c00205" target="_blank" rel="noopener noreferrer">5</a>]

### 16. Metal Coordination
**Defaults (Generic):** Primary coordination shell ≤ 2.6 (Zn/N/O donors often 2.0–2.3; Mg 2.0–2.2; Ca 2.3–2.6); secondary (outer) shell exploration ≤ 3.0.  
**Note:** Element-specific tuning advisable; current defaults represent conservative upper bounds.  
**Reference:** (specialized ion-specific literature to be added)

---

## References
Below are numbered primary sources; each citation is hyperlinked to the DOI for direct access.

1. Hydrogen bonding definitions – IUPAC Recommendations 2011. <a href="https://doi.org/10.1351/PAC-REP-10-01-01" target="_blank" rel="noopener noreferrer">https://doi.org/10.1351/PAC-REP-10-01-01</a>  
2. Halogen bonding definition – IUPAC Recommendations 2013. <a href="https://doi.org/10.1351/PAC-REC-12-05-10" target="_blank" rel="noopener noreferrer">https://doi.org/10.1351/PAC-REC-12-05-10</a>  
3. σ-Hole interactions (pnictogen, chalcogen, tetrel) – IUPAC perspective 2023. <a href="https://doi.org/10.1515/pac-2020-1002" target="_blank" rel="noopener noreferrer">https://doi.org/10.1515/pac-2020-1002</a>  
4. Chalcogen bonding in protein structural context (review). <a href="https://doi.org/10.1017/qrd.2023.3" target="_blank" rel="noopener noreferrer">https://doi.org/10.1017/qrd.2023.3</a>  
5. Survey of unconventional noncovalent interactions in proteins (covers C–H···π, anion–π, sulfur–π, n→π*, salt bridge context, dispersion heuristics). <a href="https://doi.org/10.1021/acsomega.3c00205" target="_blank" rel="noopener noreferrer">https://doi.org/10.1021/acsomega.3c00205</a>  
6. Protein cation–π interaction paradigms. <a href="https://doi.org/10.1073/pnas.96.17.9459" target="_blank" rel="noopener noreferrer">https://doi.org/10.1073/pnas.96.17.9459</a>  
7. Aromatic π–π stacking conceptual framework (Hunter–Sanders model). <a href="https://doi.org/10.1021/ja00167a008" target="_blank" rel="noopener noreferrer">https://doi.org/10.1021/ja00167a008</a>  

Additional specialized literature (anion–π focused treatments, refined n→π* energetics, metal-specific coordination radii) can be appended upon request; current selection anchors defaults surfaced in the UI.

---

## Mapping to Implementation
Each detector applies these thresholds in a vectorized geometry pipeline. Adjustable parameters (distance, angle windows, offsets) are exposed via settings or environment flags; provenance hashing records the active numeric set to ensure reproducibility.

---

## Change Log (Scientific Documentation)
| Version | Date | Update |
|---------|------|--------|
| 1.1.0 | 2025-09-29 | Added full criteria & reference summary; aligned with UI Info tab. |

