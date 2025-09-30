# Edge Cases & Ambiguities

## Alternate Conformations
Multiple occupancy states may lead to ambiguous geometry. Current strategy: retain first highest-occupancy altLoc; future weighting could refine borderline interactions.
Recommended reading: <a href="https://doi.org/10.1107/S0907444909042073" target="_blank" rel="noopener noreferrer">Validation of alternate conformations</a>.

## Missing Hydrogens
Hydrogen positions often absent in X-ray structures; donors approximated by projecting along heavy atom-bond vector using ideal bond lengths.
Error source mitigation: limit strong subclass marking when hydrogens inferred.

## Protonation States
Histidine tautomer/protonation ambiguous at moderate resolution. Optional flag (future) to treat HIS as protonated only if forming plausible H-bond network.
Reference context: <a href="https://doi.org/10.1002/prot.20331" target="_blank" rel="noopener noreferrer">Histidine protonation surveys</a>.

## Metal Coordination Variability
Coordination numbers and geometries vary (e.g., Zn tetrahedral vs octahedral). Mis-assignment risk when ligand identities ambiguous.
Validation tip: compare observed angles vs ideal geometry templates (Harding 2004; Zheng 2014).

## Boundary Distances
Interactions at upper distance threshold edges may be weak; classification labels indicate relative confidence (strong/moderate/weak when applicable).
Suggestion: flag counts near boundary for sensitivity analysis passes.

## π–π Classification Overlap
Configurations can straddle parallel vs offset boundaries; algorithm assigns subtype based on priority order (parallel>offset>edge) to maintain deterministic labeling.
Future improvement: continuous scoring that interpolates between archetypes (Meyer 2003).

## σ-hole Direction Estimation Error
Approximation of σ-hole axis uses primary substituent vector; multi-substituted environments may skew actual electrostatic direction.
Potential fix: compute electrostatic potential map sampling (future performance permitting).

## Backbone Disorder
High B-factor or unresolved loops reduce reliability of inferred n→π* interactions due to coordinate uncertainty.
Mitigation: require B-factor percentile filter for strong classification.

## Carboxylate Symmetry
Equivalent oxygen atoms in ASP/GLU may produce duplicate geometry candidates; deduplicated via sorted pair constraint.
Corner case: bridging bidentate interactions may appear as two near-parallel entries—hand-labeled as single site in future normalization.

## Charged Group Averaging
Using centroid of multi-atom group may slightly misrepresent actual electrostatic potential center; tradeoff favors performance.
Alternative (future): assign weighted center using empirical partial charges.
