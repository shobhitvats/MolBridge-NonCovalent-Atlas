## MolBridge Scientific Documentation

### Purpose
MolBridge provides reproducible, literature-grounded identification of noncovalent interactions across protein (and protein–ligand) structures for structural biology, computational chemistry, and design studies.

### Interaction Classes & Rationale
| Class | Biophysical Role | Core Criteria (Simplified) |
|-------|------------------|----------------------------|
| Hydrogen Bond | Stabilizes secondary & tertiary structure; specificity | Donor–acceptor heavy atom ≤3.5 Å; D–H···A ≥120° |
| Halogen Bond | Directional σ-hole interaction; ligand optimization | X···A ≤ Σ(vdW) (+0.2 Å slack); C–X···A ≥150° |
| Chalcogen / Pnictogen / Tetrel | Emerging σ-hole families; packing & recognition | Donor element–acceptor ≤(ΣvdW+0.5); axis angle window |
| π-π Stacking | Aromatic stabilization, interface packing | Centroid ≤5.5 Å; interplanar ≤30° (parallel) / 60–120° (T) |
| Cation–π | Electrostatic + quadrupole stabilization | Centroid distance 3–6 Å; lateral offset ≤1.5 Å |
| Anion–π | Charge–quadrupole modulation | Anion–centroid ≤5.0 Å; approach vs ring normal ≤30° |
| n→π* | Local backbone stabilization | O(i)→C(i+1) 2.9–3.6 Å; approach 95–140° |
| Salt Bridge | Strong electrostatic pair | Heavy atom ≤4.0 Å; sidechain centroid ≤6.0 Å |
| Ionic (general) | Broader electrostatic landscape | Charge centroid distance cutoff |
| Hydrophobic | Packing & desolvation | C···C ≤5.0 Å (explore up to 5.5–6.0) |
| London Dispersion | Weak cumulative packing | vdW shell band 3.3–5.0 Å cluster |
| Sulfur–π | Polarizable sulfur orientation | S–centroid ≤6.0 Å; offset ≤3.0 Å |
| Metal Coordination | Catalysis / structure | Metal–ligand primary ≤2.6 Å (element specific) |

### Parameterization Sources
Parameters anchored to peer-reviewed ranges (see README references) with adjustable sliders permitting controlled exploration. Each detection module exposes defaults aligned to frequently cited midpoints.

### Hydrogen Bond Subtypes
Subtype classification distinguishes structural roles (backbone vs sidechain, water mediation). This enables secondary analyses: fold stabilization mapping, interface enrichment studies.

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

### Future Scientific Enhancements
* Ligand-specific interaction enrichment analysis
* Interface vs core differential interaction profiling
* Cooperative network motifs (triads, ladders)
* Energetic scoring overlays (empirical potentials)

---

## Detailed Detection Criteria (Full Specification)

This section mirrors and expands the in-application Info tab with explicit default thresholds, recommended exploration ranges, and literature anchors. All angles are in degrees (°) and distances in Å unless otherwise noted. Where vdW radii sums are referenced, a small tolerance (+0.2–0.3 Å) is often permitted to accommodate crystallographic uncertainty and thermal motion.

### 1. Hydrogen Bonds
**Defaults:** Donor–Acceptor heavy atom distance ≤ 3.5; D–H···A angle ≥ 120.  
**Exploration Range:** Distance 2.6–3.6; Angle 100–180 (tighten ≥150 for high specificity).  
**Notes:** Heavy-atom proxy geometry used if hydrogens absent. Low-barrier variants (short, near-linear) are not separately subclassed yet.  
**Reference:** IUPAC Recommendations 2011 (DOI: 10.1351/PAC-REP-10-01-01).

### 2. Halogen Bonds
**Defaults:** X···A ≤ Σ(vdW radii) + 0.2; C–X···A ≥ 150 (soft lower bound); preferred ~165–180.  
**Directionality:** σ-hole alignment enforces near-linear C–X···A axis.  
**Reference:** IUPAC Recommendations 2013 (DOI: 10.1351/PAC-REC-12-05-10).

### 3. Chalcogen Bonds (S / Se / Te)
**Defaults:** S/Se/Te···A ≤ 4.0 (typical strong 3.2–3.8); C–S···A (or substituent axis) 115–155; out-of-plane deviation |φ| ≤ 50.  
**Rationale:** σ-hole opposite covalent bond(s) yields anisotropic positive electrostatic potential. Tighter angle windows increase precision.  
**Reference:** Review of unconventional noncovalent interactions (DOI: 10.1017/qrd.2023.3).

### 4. Pnictogen Bonds (N / P / As)
**Defaults:** Donor pnictogen–acceptor ≤ Σ(vdW) + 0.4; alignment ≥ 150.  
**Reference:** IUPAC 2023 σ-hole families (DOI: 10.1515/pac-2020-1002).

### 5. Tetrel Bonds (C / Si / Ge)
**Defaults:** C/Si/Ge···A ≤ Σ(vdW) + 0.5 (protein context: carbon-centered cases rare); approach (σ-hole axis) ≥ 160 for stringent classification (≥150 relaxed).  
**Reference:** IUPAC 2023 (DOI: 10.1515/pac-2020-1002).

### 6. π–π Stacking
**Parallel (Face-to-Face):** Centroid distance 3.3–4.1 strongly favored; classify ≤5.5. Interplanar angle ≤ 30.  
**T-shaped / Edge-to-Face:** Centroid distance ≤ 6.0; ring normal angle ~60–120 (center above ring rim).  
**Reference (aromatic interactions framework):** Hunter & Sanders 1990 (DOI: 10.1021/ja00167a008).

### 7. C–H···π Interactions
**Defaults:** C (donor carbon) projected to ring centroid 2.0–4.8; C–H···centroid angle ≥ 90; perpendicular offset ≤ 2.5.  
**Reference:** Incorporated in broad unconventional interaction surveys (DOI: 10.1021/acsomega.3c00205).

### 8. Cation–π Interactions
**Defaults:** Cation (sidechain charged group centroid) to aromatic centroid 3.0–6.0 (optimum 4.0–5.0); lateral offset (distance from ring normal projection) ≤ 1.5.  
**Note:** Metal cations (e.g., Zn) can be optionally included depending on parameter flags.  
**Reference (protein cation–π paradigms):** Gallivan & Dougherty 1999 (DOI: 10.1073/pnas.96.17.9459).

### 9. Anion–π Interactions
**Defaults:** Anion center to ring centroid ≤ 5.0 (enrichment <4.3); angle between anion→centroid vector and ring normal ≤ 30 (approach along positive quadrupole axis).  
**Reference:** Captured within surveys of unconventional interactions (DOI: 10.1021/acsomega.3c00205).

### 10. n→π* Interactions
**Defaults:** Carbonyl O(i) → C(i+1) distance 2.9–3.6 (cutoff ≤3.5 default); approach (proxy for lone pair to π* orientation) 95–140.  
**Reference:** Backbone stabilization discussions (survey context DOI: 10.1021/acsomega.3c00205). Additional specialized literature may be integrated in future revisions.

### 11. Hydrophobic Contacts
**Defaults:** Nonpolar sidechain heavy atom (sp3/sp2 carbon) pair distance ≤ 5.0 (exploratory expansion 5.5–6.0).  
**Note:** Counts scale with surface burial; interpret alongside SASA if available.

### 12. London Dispersion (Heuristic)
**Approach:** Identify nonpolar atom pairs in vdW window shell 3.3–5.0 and cluster to emphasize cumulative packing rather than isolated geometry.  
**Note:** Not an energy calculation; heuristic signal only.

### 13. Sulfur–π Interactions
**Defaults:** Sulfur (CYS SG / MET SD) to aromatic centroid ≤ 6.0 (strong <5.0); perpendicular offset ≤ 3.0; projection onto ring normal near center enhances classification confidence.  
**Reference:** Unconventional interaction survey (DOI: 10.1021/acsomega.3c00205).

### 14. Salt Bridges
**Defaults:** Opposite formal charge sidechain heavy atom pair ≤ 4.0 (e.g., NH₁/₂/₃ to carboxylate O); centroid (positive group vs carboxylate) ≤ 6.0.  
**Note:** Histidine considered protonated only under selected protonation assumptions; adjustable by parameter flags.

### 15. Ionic (General Electrostatics)
**Defaults:** Charged group centroids ≤ 6.5 (looser) or ≤ 5.5 (stricter).  
**Use Case:** Broader charge network mapping distinct from strict salt bridges.

### 16. Metal Coordination
**Defaults (Generic):** Primary coordination shell ≤ 2.6 (Zn/N/O donors often 2.0–2.3; Mg 2.0–2.2; Ca 2.3–2.6); secondary (outer) shell exploration ≤ 3.0.  
**Note:** Element-specific tuning advisable; current defaults represent conservative upper bounds.

---

## Reference Summary
| Topic | Citation / DOI |
|-------|----------------|
| Hydrogen Bonds | 10.1351/PAC-REP-10-01-01 |
| Halogen Bonds | 10.1351/PAC-REC-12-05-10 |
| σ-Hole (Pnictogen / Chalcogen / Tetrel) | 10.1515/pac-2020-1002 |
| Chalcogen in Proteins Review | 10.1017/qrd.2023.3 |
| Unconventional Noncovalent Interactions (broad survey) | 10.1021/acsomega.3c00205 |
| Cation–π in Proteins | 10.1073/pnas.96.17.9459 |
| Aromatic / π–π Conceptual Framework | 10.1021/ja00167a008 |

Additional specialized literature (anion–π, n→π* deep dives, metal-specific coordination radii) can be appended upon request; current selection anchors defaults already surfaced in the UI.

---

## Mapping to Implementation
Each detector applies these thresholds in a vectorized geometry pipeline. Adjustable parameters (distance, angle windows, offsets) are exposed via settings or environment flags; provenance hashing records the active numeric set to ensure reproducibility.

---

## Change Log (Scientific Documentation)
| Version | Date | Update |
|---------|------|--------|
| 1.1.0 | 2025-09-29 | Added full criteria & reference summary; aligned with UI Info tab. |

