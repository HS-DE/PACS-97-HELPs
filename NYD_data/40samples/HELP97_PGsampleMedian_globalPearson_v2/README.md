# NYD: one global Pearson per Heavy HELP (v2)

Each HELP provides a 40-observation Heavy log2 vector. Each sample provides
one original log2 PG median, computed across its own detected PGs.
Pearson across these TWO 40-value vectors yields ONE coefficient per HELP.
It is NOT per-protein Pearson and NOT per-sample Pearson.

- fi_i = median_all_samples(log2 Heavy_j) - log2 Heavy_ji.
- r_j > 0: ALL proteins in sample i receive +fi_i.
- r_j < 0: ALL proteins in sample i receive -fi_i.
- Undefined/zero r: no correction.

## Main full40 results

- 40/40 detected Heavy HELP: **74 / 97**.
- Positive Pearson HELP: 12; negative: 62; zero or undefined: 0.
- Better overall clinical median PG CV: 0 / 74.

| HELP | Pearson r | sign | Raw median CV | Corrected median CV | change pp |
|---|---:|---:|---:|---:|---:|
| ALDFAVGEYNK | -0.3782 | -1 | 25.18% | 27.65% | +2.47 |
| YWGVASFLQK | -0.5749 | -1 | 25.18% | 28.10% | +2.92 |
| VTVLGQPK | -0.6834 | -1 | 25.18% | 28.15% | +2.97 |
| LLIYWASTR | -0.1538 | -1 | 25.18% | 28.33% | +3.15 |
| STDYGIFQINSR | -0.4571 | -1 | 25.18% | 28.40% | +3.22 |
| IVIEYVDR | +0.2052 | +1 | 25.18% | 28.90% | +3.72 |
| SLDFTELDVAAEK | -0.1143 | -1 | 25.18% | 28.92% | +3.74 |
| AADDTWEPFASGK | -0.4418 | -1 | 25.18% | 28.92% | +3.74 |

## Historical exclusion sensitivity (39 samples)

- Full detection 39/39: 91 / 97.
- Better overall median PG CV: 0 / 91.

## File layout

- full40/01_global_Pearson_HELP_ranking_and_CV.csv
- full40/02_every_HELP_sample_factors.csv.gz
- full40/00_summary_all_HELP.png
- full40/per_HELP/SEQUENCE/corrected_PG_log2.tsv.gz
- full40/per_HELP/SEQUENCE/01_Raw_vs_Corrected_boxplots.png
- full40/per_HELP/SEQUENCE/02_Pearson_Heavy_vs_PGmedian_and_factors.png
- full40/per_HELP/SEQUENCE/03_Raw_vs_Corrected_common_PCA.png
- exclude_T003_bad_sample_39/01_global_Pearson_HELP_ranking_and_CV.csv
- exclude_T003_bad_sample_39/02_every_HELP_sample_factors.csv.gz

## Caveats

No identical-composition pooled QC in NYD. CV over clinical samples includes biology.
Pearson derived using clinical PG sample medians is not an independent technical reference.
A single per-sample additive log2 shift leaves sample-to-sample protein-profile
Pearson correlations invariant; PCA, sample median and CV are the useful diagnostics.
Previously completed per-protein Pearson v1 remains unchanged.
