# 3. Detector Taxonomy & Deep Dives

## 3.1 Taxonomy Overview
| Category | Detectors | Shared Geometric Theme |
|----------|-----------|------------------------|
| Directional Hydrogen-like | Hydrogen bonds, chalcogen, halogen, pnictogen, tetrel | Distance + angular linearity / σ-hole alignment |
| Aromatic & π Systems | π–π stacking, cation–π, anion–π, CH–π, n→π* | Planarity, ring normals, centroid offsets |
| Electrostatic | Ionic interactions, salt bridges | Charge center distances & optional orientation |
| Dispersion / Van der Waals | Hydrophobic contacts, London dispersion proxy | Proximity & element class filters |
| Coordination | Metal coordination | Sphere geometry + ligand set ranking |

## 3.2 Common Detector Interface
```
(records, metrics) = detect(feature_store, params, context)
```
- records: list[dict] canonical entries
- metrics: funnel counts + optional extras

## 3.3 Hydrogen Bonds
**Criteria Components**:
1. Donor–Acceptor heavy atom distance: typical window 2.5–3.6 Å (adaptive center ~3.2).  
2. D–H–A angle > 120° (liberal) or > 140° (strict).  
3. If hydrogen absent, idealized projection along donor-acceptor vector.
4. Optional second angle A–Base (acceptor lone pair direction) when geometry known.

**Algorithm**:
1. Subset donors/acceptors using bitmasks.  
2. KD-tree distance prune with max_cutoff.  
3. Vector compute distances; mask outside min/max.  
4. Estimate hydrogen coordinates (if needed).  
5. Compute angles; apply threshold mask.  
6. Score = f(distance, angle) → classification strong/moderate/weak.  
7. Emit records.

**Edge Cases**: Symmetry duplicates (i<j ordering), ambiguous donors (e.g., ARG side chain) resolved by including all hydrogens if explicit.

## 3.4 π–π Stacking
**Subtypes**: Parallel face-to-face, Offset parallel, T-shaped (edge-to-face).  
**Key Metrics**:
- Centroid distance d
- Normal angle θ
- Lateral offset o (projection difference)

Rules (illustrative):
| Subtype | d Range (Å) | θ Range (deg) | Offset |
|---------|-------------|---------------|--------|
| Parallel | 3.5–5.5 | 0–20 | o < 1.5 |
| Offset Parallel | 3.5–6.0 | 0–25 | 1.5 ≤ o ≤ 3.0 |
| T-shaped | 4.0–6.0 | 60–120 | N/A (angle dominated) |

Algorithm mirrors earlier description with subtype classification branch.

## 3.5 Cation–π Interactions
Use positive charge centroids vs ring centroids.  
Distance cutoff ~6.0 Å.  
Orientation heuristic: vector from ring centroid to cation should align with ring normal within angle tolerance for “axial” classification.

## 3.6 Anion–π Interactions
Similar to cation–π but with negative centroids. Consider electrostatic alignment: ring quadrupole axis alignment (future improvement: quadrupole tensor approximation).

## 3.7 CH–π Interactions
Hydrogen (on aliphatic carbon) projecting toward aromatic ring face.  
Simplified: distance from hydrogen (or carbon proxy) to ring centroid and angle of C–H vector vs ring normal.

## 3.8 n→π* Interactions
Backbone carbonyl oxygen lone pair donating into adjacent carbonyl π* orbital.  
Geometric signature: O_i → C=O_{i+1} distance ~3.0–3.8 Å and specific Bürgi–Dunitz angle (~107°) approximation.

## 3.9 Halogen / Chalcogen / Pnictogen / Tetrel Bonds
Unified pattern: electrophilic “σ-hole” region aligned with acceptor (often lone pair donor).  
Steps:
1. Identify heavy atom (e.g., Cl, Br, S, Se) with bound substituents.  
2. Estimate σ-hole axis along extension of covalent bond.  
3. Distance prune to potential acceptors.  
4. Angle (heavy–σ-hole–acceptor) within narrow cone.  
5. Secondary angle (acceptor base) optional.  
6. Strength classification via normalized directional alignment score.

## 3.10 Ionic Interactions & Salt Bridges
Positive vs negative charge centroids distance cutoffs.  
Salt bridge if distance ≤ threshold (e.g., 4.0 Å) and side chains appropriately oriented (orientation optional for base classification).

## 3.11 Hydrophobic Contacts
Residue-level contact network approximated by counting nonpolar heavy atom pairs within a window (e.g., 4.5 Å) excluding polar element flags.

## 3.12 London Dispersion Proxy
Not a rigorous physics model; counts close hydrophobic atom clusters. Potential future direct Lennard-Jones energy estimation.

## 3.13 Metal Coordination
Central metal atom (flagged) + ligand atoms (O, N, S) within coordination sphere (radius table by metal type).  
Optionally compute coordination geometry deviation (octahedral, tetrahedral) by angle set RMSD vs ideal templates.

## 3.14 Funnel Metrics Across Detectors (Comparative Insight)
| Detector | Typical Raw→Candidate Reduction | Heaviest Geometry Phase |
|----------|---------------------------------|--------------------------|
| Hydrogen bond | 10×–25× | Angle calculations |
| π–π stacking | 50×–200× | Offset + angle classification |
| Cation–π | 30×–80× | Normal alignment |
| Halogen bond | 15×–40× | σ-hole vector estimation |
| Metal coordination | 5×–15× | Geometry template fitting |

Efficiency improves by tightening initial prune radii while monitoring missed reference interactions (parity tests).

Proceed to `feature_pipeline.md`.
