# suc_data — SomaScan-inspired PACS normalization (16-plate dataset)

**Status: Exploratory comparison completed; no correction method is validated for deployment.**

This is a transparent **SomaScan-inspired** calibration experiment, **not** a reproduction of SomaLogic's proprietary SomaScan normalization. It uses HELP heavy-peptide signals as an approximate sample-level technical reference and either `LM-QC` or `QC` as an experiment-specific plate calibrator.

## Data and inputs

- Source metadata: `metadata(1).xlsx` (`sample_id`, `group`, `condition`, `batch`, `Plate`)
- Protein quantities: `Plate1-16_QC_report.pg_matrix.tsv`
- HELP quantification: `Accurate_Quantification_Report_v4_苏创.xlsx`, `Peptide_Long`, field `Accurate_AUC`
- Observed: 16 plates, 127 listed sample IDs, 6,137 protein groups.
- Analyzed: `LM-QC` 16 (one/plate), `QC` 32 (two/plate), study `sample` 32 (two/plate), total 80. Other groups (`Internal-QC` 32, `Blank` 15) excluded.
- HELP 97 sequences/precursors, 87 satisfying selection rule (at least 75% of included samples with `Peak_Detection_QC_Flag=PASS`, `Peak_Scan_Points>=8`, `N_Common_Transitions>=2`, positive `Accurate_AUC`).
- Protein values are raw linear intensities. Missing values are left missing (no imputation). Analyzed PG missingness: 38.31%.

## Model steps

1. **HELP sample normalization:** for each peptide and plate, find the median `log2(Accurate_AUC)` over its 5 non-blank, non-Internal-QC samples. Each sample receives the *negative median log2 deviation across the 87 HELP peptides*, then the sample shifts are re-centered within each plate. This is not SomaScan's 12 hybridization-control readout method; HELP peptides report different technical stages.
2. **Intraplate calibrator median normalization:** calibrator replicates within each plate are aligned by a median ratio across protein groups, for *calibrator samples only*. With `LM-QC`, **one calibrator per plate means this step is the identity**. With `QC`, 2 calibrators per plate allow this stage.
3. **Plate scaling:** for each protein, median of calibrators in the plate is compared with a reference median across the 16 plates; the across-protein median of those log2 differences yields one plate factor.
4. **Per-protein plate calibration:** after step 3, each protein gets its own residual plate factor, using a reference from calibrators over at least 8 plates. Plate-protein combinations without a reference retain the plate-level result.
5. **ANML:** **not run**. Requires the appropriate external healthy-population reference distribution; naive across-study median normalization can erase genuine biological changes.
6. **QC assessment:** comparison of protein CV on the original linear intensity scale (sample SD, ddof=1, minimum 3 detections per protein); separate within-calibrator, cross-QC-type and *QC 2-fold out-of-fold* assessments.

## Main findings

All numbers below are **median per-protein QC CV**, calculated from detected values, **not** a median of sample CVs.

| Calibrator | Evaluation | Raw CV | Final CV | Interpretation |
|---|---|---:|---:|---|
| LM-QC | Same LM-QC used in calibration | 16.81% | ~0% | Severe in-sample overfit (one LM-QC per plate; not valid evidence) |
| LM-QC | Independent QC group | 41.43% | 75.22% | Worse |
| QC | Separate LM-QC group | 16.81% | 30.35% | Worse |
| QC | Same QC used for calibration | 41.43% | 26.96% | Appears better **in-sample** only |
| QC | QC held out one/plate, 2-fold OOF | 41.43% | 81.83% | Worse on held-out replicates |

**Sensitivity to HELP normalization:** For QC 2-fold out-of-fold evaluation, a calibration-only run without HELP sample normalization produced 49.51% CV (vs 41.43% raw); a HELP-driven run with factors above 2.5x *disabled* produced 50.24%; full HELP correction produced 81.83%. Thus the tested correction is **not supported** for general use on this dataset.

**Anomalous HELP factors:** `P7-73` and `P9-37` require positive multipliers of 4.89x and 10.89x, respectively, under the sample-level rule. They have much lower median heavy HELP signals than other samples on the same plates, but cannot automatically be treated as uniformly low protein signals; aggressive correction increased variability. See `results/sample_HELP_factors.tsv`.

**Important limitations:** `LM-QC` and `QC` might not be identical reference material; cross-QC-type evaluation measures transfer to a different QC type, not necessarily expected biological equivalence. QC out-of-fold is the stronger like-for-like check, but only two QC replicates per plate means calibrating from one and testing on one is noisy. Any per-protein plate correction is unable to reduce *within-plate* CV for independent samples, since every sample on the plate receives the same per-protein shift. Only step 1 or a genuinely sample-specific method can change within-plate precision.

## Reproduction

Run from the local folder containing the three input files (absolute paths supported), with `numpy` installed:

```bash
python suc_data/run_soma_like.py \
  --metadata 'metadata(1).xlsx' \
  --pg 'Plate1-16_QC_report.pg_matrix.tsv' \
  --help-report 'Accurate_Quantification_Report_v4_苏创.xlsx' \
  --out suc_data/results
python suc_data/sensitivity_checks.py \
  --metadata 'metadata(1).xlsx' \
  --pg 'Plate1-16_QC_report.pg_matrix.tsv' \
  --help-report 'Accurate_Quantification_Report_v4_苏创.xlsx' \
  --out suc_data/results
```

**Do not use the corrected PG matrices for downstream biological association or fold-change analyses without additional validation.**

## Outputs

- `results/QC_CV_summary.tsv`: overall CV by calibrator, evaluation set, and stage
- `results/QC_CV_per_protein.tsv`: per-protein raw and corrected CVs
- `results/sensitivity_CV_summary.tsv`: no HELP, full HELP, gated/capped HELP, detection thresholds and QC cross-validation
- `results/QC_intraplate_pair_CV.tsv`: within-plate QC replicate CV
- `results/sample_HELP_factors.tsv`, `results/HELP_coverage.tsv`: HELP normalization factors and peptide coverage
- `results/plate_factors_*.tsv`, `results/per_protein_plate_factors_*.tsv`, `results/calibrator_intraplate_factors_*.tsv`: transparent correction factors
- `results/pg_matrix_*`: linear-intensity matrices after individual stages; QC OOF matrix contains **only the held-out QC columns** (all other intensity cells missing)
- `results/sample_manifest.tsv`, `results/run_manifest.json`: audit metadata

## References

- SomaLogic official technical note: https://somalogic.com/wp-content/uploads/2022/07/SL00000442_Rev4_2021-07_-SomaScan-v4.0-and-v4.1-Data-Standardization.pdf
- SomaScan normalization and quality control methods: https://pubmed.ncbi.nlm.nih.gov/40601147/