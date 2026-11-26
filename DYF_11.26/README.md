# DYF_11.26 — PACS normalization benchmark

**Exploratory only; no method validated for production.** 6,675 PG rows, 22 analyzed samples: 15 QC (5 per plate × 3), 3 HC, 4 S. Three internal calibrators and three blanks excluded (absent in intensity tables). 6,120 observed PGs, 3,930 complete PGs, 17.77% missing. 97 HELP, 76 complete, 86 eligible for GMF. Minimum QC HELP CV 4.09%; HELP median multipliers 0.8474–1.0947.

## QC results

| Method | Raw CV % | Full fit CV % | Held-out QC CV % | Held-out plate CV % |
|---|---:|---:|---:|---:|
| Raw | 11.8219 | 11.8219 | 11.8219 | 11.8219 |
| Current_PACS | 11.8219 | 11.8213 | — | — |
| HELP_median | 11.8219 | 12.5679 | 12.6458 | 12.9923 |
| Soma_like | 11.8219 | 10.3659 | 12.0970 | — |
| Gaussian_GMF | 11.8219 | 11.1720 | — | — |
| HELP_GMF | 11.8219 | 10.1118 | 12.5073 | 13.9816 |
| HELP_GMF_plate | 11.8219 | 9.6609 | 12.5247 | — |
| Gated_0.25 | 11.8219 | 11.4810 | 11.7507 | 11.8588 |
| Gated_0.5 | 11.8219 | 11.1538 | 11.7636 | 11.9520 |
| Gated_1 | 11.8219 | 10.8267 | 11.8961 | 12.3795 |

CV uses linear intensity, sample SD (ddof=1), >=3 matched observations per PG, 5,965 PGs, no imputation. Full fit is not independent validation. Five balanced folds hold out one QC per plate (12 training, 3 test per fold); leave-one-plate-out holds out 5 (10 training). HELP extraction threshold recalculated with training QC; unlabeled HELP measurements across all samples still used for completeness screening. M1 is group-transductive; M4 uses all protein data; neither has an independent CV here. M3 and M6 plate corrections require same-plate reference; not tested on unseen plates. M6 rank 2 with only 3 plates is particularly prone to overfitting. Gated 25% has only 0.071 pp OOF gain, not independently confirmed.

## Methods

- M0 Raw: match IDs, exclude unused groups, log2 positive intensities, preserve missing observations.
- M1 Current PACS: supplied PACS_norm.py; HELP PASS, >=2 transitions, scan points 8 then 7; complete HELP; QC CV ranking; group/condition-wise PG median; RMSE matching; adaptive offset. QC TLYSSSPR, HC TFTLLDPK, S STVEIFK.
- M2 Multi-HELP median: QC median reference per HELP, median log2 HELP deviation per sample; subtract from all PGs.
- M3 Soma-like: per-PG QC global and local plate medians; subtract local-minus-global. No HELP or ANML. Not proprietary SomaScan.
- M4 Gaussian GMF: one-hot QC/HC/S condition and plate, ridge plate penalty 2, rank-3 residual, three iterations; remove plate term only. Not official omicsGMF. Condition labels preserved.
- M5 HELP-GMF: HELP with >=80% QC coverage and finite CV; QC median reference; global median deviation plus 2 shape PCs trained on QC; standardize using QC; PG-wise ridge penalty 2, >=4 QC observations; subtract predicted technical term.
- M6 HELP-GMF + plate: M5 residual QC medians; global plate median; rank-2 shared residual; subtract 0.5 global + 0.5 shared + n/(n+4) remainder.
- Gated HELP-GMF: skip samples with max absolute standardized factor >2.5; otherwise subtract 25%, 50% or 100% of technical term.
- Optional imputation: disjoint 5%/5% masked study entries, select rank 0/2/5 by validation RMSE; test on 1815 entries. Only PGs with >=50% coverage eligible; original observations unchanged. RMSE 2.2295 raw, 2.2281 HELP-GMF, 2.2314 M6. Synthetic 1.5x shift recovery is not real biological validation.

## Reproduction

Install requirements, place original supplied PACS_norm.py beside run_dyf.py and run_help_gmf.py (identical model core copied from suc_data_sam). Run `python run_dyf.py /path/to/inputs results`. Input names: pg_matrix(2).tsv, metadata(1).xlsx, Accurate_Quantification_Report_v4_10ppm_no_msconvert_remove_8ScanPoints_transitionAgreement.xlsx.

This repository contains compact numerical result tables, a comparison chart and reproduction script. Full observed/imputed matrices, model snapshot, per-protein tables, original PACS script and PNG figures are provided in the separate downloadable archive; the GitHub connection cannot directly transfer the local binary archive. No raw input files republished.
