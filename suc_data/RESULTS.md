# 16 Plate SomaScan-inspired calibration: Results

**Data:** 6,137 protein groups; 80 analyzed samples (16 LM-QC, 32 QC, 32 study), 16 plates, 87/97 HELP used for sample normalization. Raw data contain substantial missingness (38.31% among analyzed samples).

## Median per-PG CV of the two QC groups

| Calibration reference | Evaluation type | Raw CV | Final stage-4 CV | Finding |
|---|---|---:|---:|---|
| LM-QC | LM-QC calibrator in-sample | 16.81% | ≈0% | Overfit; invalid validation |
| LM-QC | QC (different QC type) | 41.43% | 75.22% | Worsens |
| QC | LM-QC (different QC type) | 16.81% | 30.35% | Worsens |
| QC | QC calibrator in-sample | 41.43% | 26.96% | Optimistic training-set result |
| QC | QC 2-fold held-out (one QC per plate per fold) | 41.43% | 81.83% | Worsens |

Raw LM-QC: 5,912 proteins with ≥3 detections across its 16 samples. Raw QC: 4,886 proteins with ≥3 detections across its 32 samples. All CVs on linear intensity scale, sd/mean ×100, sample sd (ddof=1); no missing-value imputation.

## Controls and sensitivity checks

QC held-out 2-fold cross-validation, with models recalibrated in every fold:

| Calibration steps | QC raw CV | Held-out QC CV |
|---|---:|---:|
| No HELP, plate + per-protein calibration | 41.43% | 49.51% |
| HELP, plate + per-protein calibration | 41.43% | 81.83% |
| HELP with >2.5x factors ignored, plate + per-protein | 41.43% | 50.24% |

`P7-73` (HELP factor 4.89x) and `P9-37` (10.89x) are extreme relative to all other analyzed samples. Blind HELP-based multiplication is harmful for these samples in particular. See `sensitivity_CV_summary.tsv` for alternate filters and detection limits, and `QC_intraplate_pair_CV.tsv` for within-plate stability.

## Interpretation

1. The apparent 0% LM-QC calibration CV is mathematically expected when calibrating a single LM-QC on each plate to the per-protein reference; it is **not evidence of technical success**.
2. Using the QC group as the calibrator improves QC CV *on the same QC samples*, but worsens the **held-out QC**. Avoid reporting the in-sample CV as performance.
3. Calibrator transfer between `LM-QC` and `QC` is unsatisfactory; if these are different biological materials, cross-type transfer is not necessarily an appropriate ground truth. The like-for-like held-out QC result is also negative.
4. This is a practical SomaScan-inspired method, not the complete SomaScan pipeline: lacks its 12 specific hybridization controls and global healthy reference (ANML); LM-QC provides no within-plate calibrator replicates.
5. The calibration should **not** replace the PACS correction process for these samples without independent technical replication and QC type/reference verification.

The detailed calculations, per-protein effect estimates and stagewise matrices are in `results/`.