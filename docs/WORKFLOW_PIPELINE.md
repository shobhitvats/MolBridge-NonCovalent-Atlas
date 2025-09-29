## MolBridge Workflow & Logic Pipeline

### Overview Diagram (Textual)
```
User Input (PDB IDs / Upload) -> Parse & Validate -> Feature Extraction (Task Graph) ->
Spatial Pruning (KD-tree) -> Candidate Generation (vector / subset / cross-set) ->
Geometry Filtering -> Interaction Acceptance -> (Optional) Normalization & Provenance ->
Funnel Metrics + Adaptive Threshold Feedback -> Aggregation -> Visualization / Reporting / API Output
```

### Detailed Block Descriptions
1. Input Acquisition: PDB fetch or file upload; Biopython parsing.
2. Validation: Detect missing atoms, alt loc pruning, chain integrity.
3. Feature Extraction: Aromatic rings, charged groups, donors/acceptors, centroids.
4. Spatial Pruning: Reusable KD-tree + multi-radius neighbor cache supplies coarse candidate pairs.
5. Candidate Generation: Detector-specific subset or cross-set selection.
6. Geometry Filtering: Vectorized angle, distance, planarity, offset checks.
7. Acceptance: Interaction record created (or deferred if lazy materialization planned).
8. Normalization (flagged): Convert to unified schema + optional strength classification.
9. Provenance (flagged): Hash structure & parameter signature.
10. Metrics: Log raw/candidate/accepted counts + timings.
11. Adaptive Feedback: Adjust search radius thresholds toward density target.
12. Aggregation: Collect per-detector arrays for UI + exports.
13. Presentation: 3D viewer overlays, tables, plots, network graphs, reports.

### Data Flow
```
coords -> KD-tree -> neighbor pairs -> geometry filter -> accepted interactions
         ^                |                  |                 |
         |                |                  v                 v
    FeatureStore   detector selection   funnel metrics     normalized rows
```

### Failure / Fallback Paths
* KD-tree build failure → brute-force distance checks (logged debug)
* Vector path error → legacy loop path fallback
* Rust extension missing → Numba → pure NumPy

### Performance Hooks
* Auto profile chooses process pool & shared memory for large workloads.
* Adaptive thresholds reduce over-broad radii.
* Columnar storage reduces Python object overhead.

### Export Targets
| Target | Mechanism |
|--------|-----------|
| CSV | Serializer converts normalized rows |
| PDF | Report generator + templated figures |
| Excel | Multi-sheet workbook exporter |
| JSON | Raw or normalized payload |

### Miro AI Prompt (See Below) will generate a visual diagram.

---
## Miro AI Diagram Prompt
Use the following prompt in Miro AI Diagram to generate the workflow:

```
Create a detailed system data flow & processing pipeline diagram titled "MolBridge Interaction Analysis Pipeline" with the following blocks and directional edges:

Blocks (arrange left-to-right logical sequence):
1. User Input
2. Parsing & Validation
3. Feature Extraction (Task Graph)
4. Spatial Pruning (KD-tree + Multi-Radius Cache)
5. Candidate Generation (Subset / Cross-Set)
6. Geometry Filtering (Vector Math)
7. Interaction Acceptance
8. Normalization & Provenance (Optional)
9. Funnel Metrics Logging
10. Adaptive Threshold Feedback Loop (back edge to Spatial Pruning)
11. Aggregation & Storage
12. Visualization Layer (3D Viewer, Tables, Plots)
13. Reporting & Exports (PDF/CSV/Excel/JSON)
14. REST API Interface

Edges:
* Linear flow from 1 -> 13
* A side edge from 13 to 14 (API alternative output)
* Feedback loop from Adaptive Threshold Feedback to Spatial Pruning
* Feature Extraction feeds both Spatial Pruning and Candidate Generation
* Aggregation & Storage feeds Visualization and Reporting

Annotations:
* Mark performance accelerators: vector geometry, process pool + SHM, adaptive thresholds, columnar storage.
* Indicate fallback paths: Vector Geometry -> Legacy Loop, KD-tree -> Brute Force.
* Include a legend for optional blocks (dashed borders): Normalization & Provenance.

Style:
* Use color coding: Blue (Core Processing), Green (Performance), Orange (Output), Grey (Optional).
* Include arrow labels where helpful: "pairs", "filtered pairs", "interactions".
```
