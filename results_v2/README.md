# PACS Trial 2, new uploaded protein matrix (stepwise)

Previously generated results/ has been archived at results_v1/.
The new results/ is an independent rerun on the newly uploaded pg_matrix.tsv.
The root pg_matrix.tsv and older archive files are not changed.

## Sample screening

- New protein rows: 6677.
- Original sample columns: 75.
- Kept: 75 (Internal-QC 6, pooled QC 15, HC 27, S 27).
- Removed sample columns: 0.
- No imputation, filtering of protein rows, normalization or transformation of input raw data.

## Definitions (log2 except CV)

For each Study sample with yeast peptide log2 intensities y1 and y2:
- factor1a = median_study((y1+y2)/2) - (y1+y2)/2.
- factor1b = mean(median_study(y1), median_study(y2)) - (y1+y2)/2.
- factor2 = median_of_pooled_QC_sample_medians(log2 HELP L/H) - sample_median(log2 HELP L/H).
- Raw: original protein intensity.
- After factor1: log2 corrected protein = log2 raw protein + factor1.
- After factor1+factor2: log2 corrected protein = log2 raw + factor1 + factor2.
- Only HC and S are corrected. Internal-QC and pooled QC are entirely unchanged.
- CV uses unlogged positive intensity, SD(ddof=1)/arithmetic mean*100%,
  detection >= max(3,ceil(0.7*group n)), minimum mean 1e-8.
- PCA and boxplots use log2 intensities. PCA-only per-protein median imputation;
  each view selects >=70% detected proteins, same set across stages, centered not z-scored.
- PCA is separately fitted at each stage, so PC coordinates are not aligned between stages.

## References

- Yeast ref a: 23.6562336357.
- Yeast ref b: 23.6525507562.
- Common factor1 b minus a: -0.00368287950429.
- Pooled QC HELP reference log2: -5.69001011213.

## Per-group median protein CV (%)

| Method | Stage | Internal-QC | pooled QC | HC | S |
| --- | --- | ---: | ---: | ---: | ---: |
| a | Raw | 63.67% | 10.40% | 23.06% | 28.83% |
| a | Factor1 | 63.67% | 10.40% | 39.25% | 45.73% |
| a | Factor1 + Factor2 | 63.67% | 10.40% | 37.65% | 57.06% |
| b | Raw | 63.67% | 10.40% | 23.06% | 28.83% |
| b | Factor1 | 63.67% | 10.40% | 39.25% | 45.73% |
| b | Factor1 + Factor2 | 63.67% | 10.40% | 37.65% | 57.06% |

## Joint figures: Raw / Factor1 / Factor1+Factor2

Two yeast reference methods x two plotting views = 4 figure sets.
Each set contains six PNG charts (24 PNGs total), no SVG.

### Method a / all_groups

![Sample boxplots](method_a/figures/all_groups/01_sample_boxplots_3stages.png)

![Group sample-median boxplots](method_a/figures/all_groups/02_group_boxes_3stages.png)

![PCA](method_a/figures/all_groups/03_pca_3stages.png)

![Per-protein CV violin/box/scatter](method_a/figures/all_groups/04_cv_3stages.png)

![Median CV changes](method_a/figures/all_groups/05_cv_median_comparison.png)

![Factors](method_a/figures/all_groups/06_factor_components.png)

### Method a / study_only

![Sample boxplots](method_a/figures/study_only/01_sample_boxplots_3stages.png)

![Group sample-median boxplots](method_a/figures/study_only/02_group_boxes_3stages.png)

![PCA](method_a/figures/study_only/03_pca_3stages.png)

![Per-protein CV violin/box/scatter](method_a/figures/study_only/04_cv_3stages.png)

![Median CV changes](method_a/figures/study_only/05_cv_median_comparison.png)

![Factors](method_a/figures/study_only/06_factor_components.png)

### Method b / all_groups

![Sample boxplots](method_b/figures/all_groups/01_sample_boxplots_3stages.png)

![Group sample-median boxplots](method_b/figures/all_groups/02_group_boxes_3stages.png)

![PCA](method_b/figures/all_groups/03_pca_3stages.png)

![Per-protein CV violin/box/scatter](method_b/figures/all_groups/04_cv_3stages.png)

![Median CV changes](method_b/figures/all_groups/05_cv_median_comparison.png)

![Factors](method_b/figures/all_groups/06_factor_components.png)

### Method b / study_only

![Sample boxplots](method_b/figures/study_only/01_sample_boxplots_3stages.png)

![Group sample-median boxplots](method_b/figures/study_only/02_group_boxes_3stages.png)

![PCA](method_b/figures/study_only/03_pca_3stages.png)

![Per-protein CV violin/box/scatter](method_b/figures/study_only/04_cv_3stages.png)

![Median CV changes](method_b/figures/study_only/05_cv_median_comparison.png)

![Factors](method_b/figures/study_only/06_factor_components.png)

## PCA explained variance

| Method | View | Stage | Proteins | PC1 | PC2 |
| --- | --- | --- | ---: | ---: | ---: |
| a | all_groups | raw | 5510 | 25.81% | 15.43% |
| a | all_groups | factor1 | 5510 | 34.82% | 13.28% |
| a | all_groups | factor1_factor2 | 5510 | 41.18% | 12.70% |
| a | study_only | raw | 5596 | 27.36% | 16.48% |
| a | study_only | factor1 | 5596 | 43.84% | 14.92% |
| a | study_only | factor1_factor2 | 5596 | 52.91% | 12.56% |
| b | all_groups | raw | 5510 | 25.81% | 15.43% |
| b | all_groups | factor1 | 5510 | 34.84% | 13.28% |
| b | all_groups | factor1_factor2 | 5510 | 41.18% | 12.70% |
| b | study_only | raw | 5596 | 27.36% | 16.48% |
| b | study_only | factor1 | 5596 | 43.84% | 14.92% |
| b | study_only | factor1_factor2 | 5596 | 52.91% | 12.56% |

## Data products

- pg_matrix_filtered.tsv: screened source, unmodified protein intensities.
- sample_factors.csv: per-sample a/b factor1, factor2, and HELP ratio statistics.
- yeast_factor_details.csv: peptide-level log2 intensities and both reference methods.
- all_stage_cv_summary.csv and all_stage_pca_summary.csv: method x stage x group/view summaries.
- method_a/ and method_b/: linear and log2 factor1/factor1+factor2 corrected matrices;
  stage_offsets.csv, per_protein_cv_stages.csv, group_cv_stages.csv,
  pca_scores.csv, pca_summary.csv, sample_log2_medians.csv and PNG figures.

## Caution

CV for both QC groups must remain unchanged by design, and cannot validate
normalization effectiveness. HC/S CV includes biological variability.
Methods a/b differ by a common log2 offset and are expected to have
identical within-group CV; their all-groups PCA can differ because QC remains fixed.
PCA explained variance increasing does not imply better normalization.
