# PACS HELP-GMF pilot (16 plates)

**Exploratory, not validated for production.** Two separate fits: LM-QC (16) and QC (32), each evaluated against the other QC type and 32 study samples. Internal-QC and blank excluded. 12,028 PGs; 80 retained samples; 97 HELP peptides. Raw input files: Plate1-16_QC_Hela_merged.pg_matrix.tsv, Accurate_Quantification_Report_v4_苏创.csv, metadata.csv. Metadata contains an invalid UTF-8 byte in an unused column; use GB18030 fallback. Original PACS baseline uses supplied PACS_norm.py, with other QC mapped to a separate study condition (transductive). Gaussian GMF is a Python approximation, not official omicsGMF. Soma-like is direct log2 protein-wise reference calibration, not proprietary SomaScan.

## Results

See results/summary.csv, results/validation.csv, results/gated_sensitivity.csv, results/imputation_benchmark.csv and results/synthetic_shift.csv. CV is linear sample SD/mean x100, minimum 3 observations per PG; same observed support across methods. Study CV is descriptive, not a technical precision target.

- LM-QC reference: HELP-GMF plate training CV 8.93% vs raw 16.81%; alternate QC 44.56% vs raw 41.43%. LM-QC leave-one-plate-out HELP-GMF 15.64%. Soma-like yields ~0% training CV (one calibrator per plate), not real generalization.
- QC reference: HELP-GMF plate training CV 29.41% vs raw 41.43%; alternate LM-QC 23.92% vs raw 16.81%. Same-type QC 2-fold held-out well 37 vs 73: raw 41.43%, HELP-GMF 43.01%, plate variant 43.23%, Soma-like 50.35%, HELP median 71.58%. No validated improvement.
- Masked study observation RMSE (5,341 test entries): raw low-rank 1.0765 log2, anchored plate models 1.0940/1.0912. No improvement. Fixed corrections preserve synthetic 1.5x shift by construction, not proof of biological accuracy.
- HELP warning: LM-QC best complete HELP CV 8.02%, QC 24.65%. P9-37 and P7-73 have median HELP multipliers ~9.2x and ~4.9x. Conservative gate and shrinkage still do not establish benefit.
- Recommendation: investigate HELP spike-in/LC-MS response mismatch and QC material identity; compare within same material with independent controls; do not deploy this model.

## Reproduction

Place the original supplied PACS_norm.py beside the scripts and the three inputs in one folder. Install numpy, pandas, scipy, scikit-learn, matplotlib, openpyxl. Run:
```bash
python prepare.py results /path/to/input
python finish_benchmark.py results /path/to/input LM-QC
python finish_benchmark.py results /path/to/input QC
python sensitivity.py results /path/to/input
```

Model: log2 observed-only ridge with global HELP median and two shape PCs, coefficients trained on selected QC only, require >=4 observations/protein; plate residual rank 2, shrink common component by 0.5 and residual by n/(n+4). HELP extraction uses existing PACS quality rules (8 then 7 scan points); threshold fitted across all selected QC, so OOF is not fully nested. HELP-GMF plate applies calibration only to plates with a QC. Checkpoints store fitted sample scores and corrections, not transferable production model. Rank 0/2/5 selected using disjoint 5% validation/5% test observed study values, only >=50% overall detected PGs eligible for optional imputation; original observations unchanged. No claims about true MNAR values. Cross-QC types may be different materials.

**Delivery status (original HELP-GMF pilot):** The aggregate result tables and reproduction scripts from that pilot are committed here; its compressed model snapshots and per-protein matrices remain separate from the new PACS run. The new dual-QC PACS outputs are documented in `RESULTS.md` and committed under `LM-QC/`, `QC/`, `HELP_GMF/`, `comparison/`, and `scripts/`. Raw input files are not published.


## New PACS dual-QC baseline (2026-11-26)

A separate run was performed on the uploaded `Plate1-16_QC_Hela_merged.pg_matrix.tsv` matrix:

- `LM-QC/`: only LM-QC is pooled QC; QC, Internal-QC, and Blank are excluded.
- `QC/`: only QC is pooled QC; LM-QC, Internal-QC, and Blank are excluded.
- `RESULTS.md` contains the full comparison, HELP diagnostics, cross-validation results, and interpretation.
- `scripts/run_two_qc_pacs.py` reproduces the two PACS runs; `scripts/help_gmf_diagnostics.py` reproduces the PCA-based HELP factor diagnostics.

The `HELP_GMF/` directory created for this PACS run contains a PCA-based diagnostic prototype. It is **not** the earlier Gaussian factor model in `run_help_gmf.py`, and it is not an exact implementation of omicsGFM.

Important: the supplied PACS code checks batch/plate metadata but does not use plate in its correction-offset calculation. Therefore the PACS outputs are baselines, not an inter-plate calibration model. The QC-run final matrix and the two full PACS per-protein CV tables are split into ordered parts; each parts directory includes instructions for reassembly.

