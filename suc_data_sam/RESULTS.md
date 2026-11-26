# Results — dual-QC PACS and HELP latent-factor evaluation

## Executive summary

Both requested PACS runs completed successfully. **LM-QC performed materially better as the pooled-QC definition in this dataset**, even though it provides only one QC injection per plate. The separate `QC` class contains two injections per plate, but its HELP and protein signals are substantially less stable.

A QC-only HELP latent-factor model contains useful information about sample-level global PG shifts, including a strong plate-associated structure in the `QC` class. However, applying the predicted global or protein-wise shifts did not improve held-out protein CV. The prototype is therefore retained as a diagnostic result, not a final correction method.

## 1. Input audit

| Item | Result |
|---|---:|
| Protein groups | 12,028 |
| Protein matrix sample columns | 196 |
| Metadata sample IDs | 127 |
| Metadata groups | 16 LM-QC; 32 QC; 32 sample; 32 Internal-QC; 15 Blank |
| HELP long-table records | 12,319 |
| Unique HELP peptide sequences | 97 |
| Duplicate Protein.Group rows | 0 |
| Duplicate HELP Sample + Peptide rows | 0 |
| Metadata IDs missing from HELP | 0 |
| Metadata IDs missing from protein matrix | 0 |
| Protein-matrix sample columns not represented in metadata | 69 |

All 127 metadata IDs matched both the HELP report and protein matrix. The 69 protein-matrix columns without metadata matches were not analyzed. One all-zero protein column, `P1-61`, belongs to Internal-QC and was excluded under both requested strategies.

The metadata CSV had one stray non-UTF8 byte in an unnamed trailing column. The unnamed columns were dropped; the five named metadata fields and all sample assignments were retained.

## 2. QC mapping and sample counts

| Run | Pooled QC | Study samples | Excluded Internal-QC | Excluded other QC | Excluded Blank |
|---|---:|---:|---:|---:|---:|
| LM-QC | 16 | 32 | 32 | 32 QC | 15 |
| QC | 32 | 32 | 32 | 16 LM-QC | 15 |

There are 16 plates. LM-QC has one injection per plate; QC has two injections per plate.

## 3. PACS results

Protein CV is calculated on linear intensities using sample standard deviation (`ddof=1`) and at least three matched observations per protein. The same protein/observation positions are used before and after correction.

| Pooled QC | QC n | Proteins with CV | Raw median CV | PACS median CV | Relative median reduction | Proteins with lower CV |
|---|---:|---:|---:|---:|---:|---:|
| LM-QC | 16 | 5,912 | 16.81% | 14.95% | 11.05% | 82.32% |
| QC | 32 | 4,886 | 41.43% | 40.17% | 3.03% | 56.04% |

The LM-QC run also increased the fraction of proteins with CV ≤15% from 40.93% to 50.29%. In the QC run, that fraction changed from 0.63% to 0.47%; the median moved slightly, but the low-CV tail did not improve.

### HELP selection by PACS

| Run | Corrected group | Selected HELP | HELP–PG Pearson r | RMSE | Selected HELP QC CV |
|---|---|---|---:|---:|---:|
| LM-QC | pooled QC | EAQLPVIENK | 0.765 | 0.102 | 10.55% |
| LM-QC | study samples | HVVPNEVVVQR | -0.136 | 0.343 | 8.02% |
| QC | pooled QC | IVIEYVDR | -0.733 | 1.082 | 25.57% |
| QC | study samples | HVVPNEVVVQR | -0.136 | 0.343 | 25.73% |

The HELP extraction step retained 97 peptides in both runs. For LM-QC, 84 were complete across all retained samples and the default scan-point threshold of 8 was accepted. For QC, the threshold was relaxed from 8 to 7, yielding 97 complete peptides, but the extraction remained unstable: the lowest pooled-QC HELP CV was 24.65%, above the script's 15% acceptance criterion.

The PACS QC-variability warning was also triggered in both runs:
- LM-QC: median raw protein CV 16.84% (slightly above the 15% warning threshold).
- QC: median raw protein CV 41.45%.

### Study-sample corrections are identical in both runs

All 32 study samples received exactly the same PACS offsets in the LM-QC and QC runs. Both runs selected `HVVPNEVVVQR` for the study-sample group, with Pearson r = -0.136 and RMSE = 0.343. This means the choice of pooled-QC class changed the QC branch and its QC evaluation, but did not change the 32 study-sample corrections in this dataset.

The study-sample median CV changed from 64.57% to 63.77% in both runs. Because all 32 study samples are assigned to one `condition=sample` group, this CV includes biological heterogeneity and is not a pure technical precision metric.

## 4. HELP latent-factor diagnostics

The new `HELP_GMF/` directory in this PACS run is a PCA-based latent-factor diagnostic; it is distinct from the earlier Gaussian factor pilot implemented in `run_help_gmf.py`. The exploratory prototype used only the selected pooled-QC class to fit a `StandardScaler + PCA` model on complete HELP features. No missing values were imputed. The fitted factors were then projected to the retained study samples.

| Metric | LM-QC | QC |
|---|---:|---:|
| Selected QC samples | 16 | 32 |
| HELP features complete across all retained samples | 84 | 97 |
| PC1 variance explained | 77.65% | 92.73% |
| HELP global shift vs PG global shift (QC Pearson r) | +0.660 | -0.825 |
| Grouped 8-fold CV R², 1 factor | 0.349 | 0.668 |
| Grouped CV RMSE, 1 factor | 0.103 | 0.303 |
| Baseline mean-predictor RMSE | 0.136 | 0.553 |
| RMSE improvement vs baseline | 24.62% | 45.19% |

The cross-validation split was grouped by plate, so the two QC injections from the same plate were never split between training and test folds.

### Plate structure

In the QC run, the plate-associated variance fractions (`η²`) for HELP factors were:
- PC1: 0.419
- PC2: 0.856
- PC3: 0.936

There are two QC injections per plate, so these are descriptive estimates from two replicates per plate. In the LM-QC run there is only one QC per plate; a within-plate variance estimate is not available, and a plate `η²` based on one observation per plate is not interpretable.

For QC, the median standardized HELP distance between plates was 2.08 times the median within-plate paired-QC distance. This supports a plate-associated HELP structure, but does not by itself establish that all of that structure is technical.

The HELP-to-PG global-shift correlations have opposite signs in the two QC types (+0.66 for LM-QC; -0.83 for QC). This is a warning against treating the two QC types as interchangeable or fitting one shared HELP-to-PG mapping across them.

Two QC injections are particularly influential in the QC-class run:
- `P9-37`: HELP global shift = -3.275 log2 (observed HELP level ≈0.103× the QC median; inverse correction multiplier ≈9.68×), while the PG global median shift is +2.138 log2 (≈4.40×).
- `P7-73`: HELP global shift = -2.272 log2 (observed HELP level ≈0.207×; inverse correction multiplier ≈4.83×), while the PG global median shift is +1.626 log2 (≈3.09×).

These are QC samples and are also extreme on the HELP PCA scores. They help explain the inverse HELP–PG relationship, but the pattern does not establish a causal loading, spike-in, or instrument error. See `comparison/QC_HELP_global_shift_outliers.tsv` and the existing `DATA_QUALITY.md` for related diagnostics.

## 5. Does the HELP factor model improve protein CV?

No, not in the held-out tests run here.

| QC strategy | Method | Held-out / evaluated median protein CV |
|---|---|---:|
| LM-QC | Raw | 16.81% |
| LM-QC | PACS | 14.95% |
| LM-QC | HELP-GMF global 1-factor OOF | 15.74% |
| LM-QC | HELP-GMF protein-wise 1-factor OOF | 15.44% |
| QC | Raw | 41.43% |
| QC | PACS | 40.17% |
| QC | HELP-GMF global 1-factor OOF | 42.11% |
| QC | HELP-GMF protein-wise 1-factor OOF | 44.51% |
| QC | QC-derived per-protein plate correction, paired-QC 2-fold OOF | 50.35% |

The global HELP factors can predict a sample-level PG median shift in grouped cross-validation, but subtracting that predicted global shift does not reduce protein-level CV to the PACS level. The protein-wise regression prototype also failed to improve held-out CV. A simple per-plate per-protein correction estimated from one QC replicate and evaluated on the paired QC replicate made CV worse.

This distinction matters: **predicting a global shift is not the same as correcting protein-specific technical variation**. The present HELP factor model is not ready to be used as a normalization matrix.

## 6. Important implementation limitation in current PACS

The supplied `PACS_norm.py` groups samples by `(group, condition)`. It checks that batch information is present and uses batch/plate for visualization, but batch/plate is not part of the correction-offset calculation. Therefore these two runs are faithful PACS runs, but **they are not an explicit inter-plate calibration model**.

The HELP latent-factor diagnostics reveal plate-associated structure, especially in the QC group, but the attempted OOF correction did not pass the protein-CV validation.

## 7. Recommendation

1. Use LM-QC as the provisional pooled-QC definition for this dataset's PACS baseline; it is clearly more stable in the observed results.
2. Do not combine LM-QC and QC into one pooled-QC group. Their HELP-to-PG relationships have opposite signs and very different CV distributions.
3. Do not deploy the current HELP-GMF candidate correction. Keep its latent scores as diagnostic features only.
4. To build a plate-aware HELP model, add at least two independent LM-QC injections per plate (or otherwise create a held-out same-material QC replicate), train plate/technical factors on QC only, and evaluate on held-out QC before applying the correction to study samples.
5. Before biological fold-change or association analysis, test whether the correction preserves known biological contrasts. The current study samples all share one metadata condition, so biological-signal preservation cannot be assessed from this test alone.

## 8. Validation

- Both saved PACS runner scripts reproduced the final protein matrices exactly; maximum absolute numeric difference was 0.0.
- The HELP latent-factor diagnostics script completed for both QC strategies.
- Raw input files were not copied into this public folder.
- The QC-run final protein matrix and the full PACS per-protein CV detail tables are split into ordered parts in GitHub because of the upload transport size limit; each parts folder includes reassembly instructions.
- The exploratory factor model is not an exact omicsGFM implementation.
