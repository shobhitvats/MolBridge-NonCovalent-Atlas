## MolBridge Scientific Documentation

### Purpose
MolBridge provides reproducible, literature-grounded identification of noncovalent interactions across protein (and protein–ligand) structures for structural biology, computational chemistry, and design studies.

### Interaction Classes & Rationale
| Class | Biophysical Role | Core Criteria (Simplified) |
|-------|------------------|----------------------------|
| Hydrogen Bond | Stabilizes secondary & tertiary structure; specificity | Donor–acceptor heavy atom distance; D–H···A angle |
| Halogen Bond | Directional σ-hole interaction; ligand optimization | X···A distance + C–X···A angle |
| Chalcogen / Pnictogen / Tetrel | Emerging σ-hole families influencing packing & recognition | Element-specific distance + alignment angle (and optional dihedral / out-of-plane) |
| π-π Stacking | Aromatic stabilization, interface packing | Centroid distance + inter-planar angle |
| Cation–π | Electrostatic + quadrupole stabilization | Cation centroid–ring centroid distance + perpendicular offset |
| Anion–π | Charge–quadrupole modulation | Anion–centroid distance + ring normal alignment |
| n→π* | Local backbone / sidechain stabilization | Carbonyl O→C distance + approach angle |
| Salt Bridge | Strong electrostatic pair | Closest heavy atom distance + centroid distance heuristic |
| Ionic (general) | Broader electrostatic landscape | Charge centroid distance only |
| Hydrophobic | Packing & desolvation | C–C distance window |
| London Dispersion | Weak cumulative packing | vdW window band |
| Sulfur–π | Polarizable sulfur orientation | S–centroid distance + perpendicular offset |
| Metal Coordination | Catalysis / structural integrity | Metal–ligand atom distances (primary + secondary shell) |

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
