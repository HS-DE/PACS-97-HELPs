# Data quality and additional diagnostics

- 12,028 input PG rows, 6,137 detected in at least one of the 80 retained samples, 713 complete across all 80. 302,865 observed cells; 68.525% missing over full 12,028 × 80 matrix (includes all-empty PG rows).
- LM-QC: 5,541–5,666 PG detected/sample (median 5,582.5). QC: 1,341–4,155 (median 3,459). Study: 1,592–5,103 (median 3,394.5).
- LM-QC selected 8 scan-point threshold, 84 complete HELP, min linear HELP CV 8.02%. QC required relaxation to 7 scan points, 97 complete HELP, min HELP CV still 24.65%; existing PACS flagged unstable HELP.
- P9-37 and P7-73 median HELP multipliers under QC reference 9.234 and 4.876; detected PG counts 1,341 and 1,689. Their detected-only log2 PG medians are 15.56 and 15.04 (other QC median sample medians 13.42), which may reflect detection censoring and not true loading. Do not infer causal spike error from this.
- QC type may denote different sample materials or stages. Cross-type transfer is not a like-for-like test. Within-QC two-fold holdout and LM-QC leave-one-plate-out are stronger evidence, although scan-point threshold selection was not fully nested.
- No production model recommended; no experimental biological ground truth or reliable uncertainty calibration.

Run `python plot_embeddings.py results /path/to/input` after the four commands in README to produce common-Raw-axis PCA and HELP factor score CSVs. PCA uses only 713 PGs complete across all 80 samples; not an unbiased biological embedding.
