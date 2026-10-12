# PACS Trial 2 - yeast spike-in + log2 HELP offsets

This repository experiment is exploratory. It does NOT select a final PACS method.
Only study samples (HC/S; 54) are corrected. All 6 Internal-QC (LM) and
15 pooled QC intensities are unchanged across all four correction scenarios.
Trial 1A / 1B and all original input files remain untouched.

## Design and formulas (log2 scale)

Let y1_s = log2(AADALLLK intensity), y2_s = log2(VNQIGTLSESIK intensity).
Both peptides are required to be finite and positive in all 54 study samples.
Define mean_ye_s = (y1_s + y2_s)/2.

**Method a:** ref_ye_a = median_study(mean_ye_s); factor1_a_s = ref_ye_a - mean_ye_s.
**Method b:** ref_ye_b = (median_study(y1_s)+median_study(y2_s))/2;
factor1_b_s = mean(ref_ye_b-y1_s, ref_ye_b-y2_s).
ref_ye_a = 23.6562336357, ref_ye_b = 23.6525507562,
and factor1_b - factor1_a = -0.00368287950429 for **every** study sample.
Because both yeast peptides are detected for every study sample, the two
factors differ by a common additive constant: their within-group protein
CV and centered PCA are expected to match; the absolute intensity level differs.

For every sample, mR_s = median over 97 HELP **log2(L/H)** ratios (positive,
finite L/H only). ref_mR = median(mR_s for 15 pooled QC).
ref_mR = -5.69001011213 in log2 ratio units.
factor2_s = ref_mR - mR_s, applied only to study samples in combined scenarios.

Study sample protein log2 intensity is corrected as:
- factor1 only: log2(I') = log2(I) + factor1_a or factor1_b.
- factor1+factor2: log2(I') = log2(I) + factor1_a/b + factor2.
- Both QC types: no applied offset (0), so original linear intensities remain exact.

Each scenario exports corrected protein intensities on log2 and back-
transformed linear scales. The linear scale is only for ordinary CV and
downstream raw-scale comparisons; it is NOT subjected to a second correction.

## Four correction scenarios x two chart views

There are FOUR correction scenarios (a F1, a F1+F2, b F1, b F1+F2).
Each generates an all-groups version (75 samples) and a study-only
version (54 HC/S samples): **8 figure collections, each 5 PNG images**.
No SVG files are generated.

### a_factor1

| Group | n | Eligible proteins | Median CV original | Median CV corrected | Median paired Δ (pp) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Internal-QC | 6 | 672 | 56.51% | 56.51% | +0.00 |
| pooled QC | 15 | 5966 | 10.33% | 10.33% | +0.00 |
| HC | 27 | 5581 | 23.35% | 39.23% | +16.84 |
| S | 27 | 5822 | 28.71% | 45.84% | +17.85 |

PCA feature counts/variance:

- all_groups: n=75, proteins=5585; before PC1/PC2 22.3%/16.0%; after 34.8%/12.2%.
- study_only: n=54, proteins=5662; before PC1/PC2 27.2%/17.0%; after 43.5%/15.1%.

**all_groups:**
![Sample protein distributions](a_factor1/figures/all_groups/01_sample_boxplots.png)

![Group boxplots](a_factor1/figures/all_groups/02_group_boxplots.png)

![PCA](a_factor1/figures/all_groups/03_pca.png)

![CV violin box jitter](a_factor1/figures/all_groups/04_cv_violin_box.png)

![Applied log2 sample offsets](a_factor1/figures/all_groups/05_applied_offsets.png)

**study_only:**
![Sample protein distributions](a_factor1/figures/study_only/01_sample_boxplots.png)

![Group boxplots](a_factor1/figures/study_only/02_group_boxplots.png)

![PCA](a_factor1/figures/study_only/03_pca.png)

![CV violin box jitter](a_factor1/figures/study_only/04_cv_violin_box.png)

![Applied log2 sample offsets](a_factor1/figures/study_only/05_applied_offsets.png)

### a_factor1_plus_factor2

| Group | n | Eligible proteins | Median CV original | Median CV corrected | Median paired Δ (pp) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Internal-QC | 6 | 672 | 56.51% | 56.51% | +0.00 |
| pooled QC | 15 | 5966 | 10.33% | 10.33% | +0.00 |
| HC | 27 | 5581 | 23.35% | 37.70% | +16.04 |
| S | 27 | 5822 | 28.71% | 57.34% | +26.81 |

PCA feature counts/variance:

- all_groups: n=75, proteins=5585; before PC1/PC2 22.3%/16.0%; after 42.1%/11.1%.
- study_only: n=54, proteins=5662; before PC1/PC2 27.2%/17.0%; after 52.5%/12.7%.

**all_groups:**
![Sample protein distributions](a_factor1_plus_factor2/figures/all_groups/01_sample_boxplots.png)

![Group boxplots](a_factor1_plus_factor2/figures/all_groups/02_group_boxplots.png)

![PCA](a_factor1_plus_factor2/figures/all_groups/03_pca.png)

![CV violin box jitter](a_factor1_plus_factor2/figures/all_groups/04_cv_violin_box.png)

![Applied log2 sample offsets](a_factor1_plus_factor2/figures/all_groups/05_applied_offsets.png)

**study_only:**
![Sample protein distributions](a_factor1_plus_factor2/figures/study_only/01_sample_boxplots.png)

![Group boxplots](a_factor1_plus_factor2/figures/study_only/02_group_boxplots.png)

![PCA](a_factor1_plus_factor2/figures/study_only/03_pca.png)

![CV violin box jitter](a_factor1_plus_factor2/figures/study_only/04_cv_violin_box.png)

![Applied log2 sample offsets](a_factor1_plus_factor2/figures/study_only/05_applied_offsets.png)

### b_factor1

| Group | n | Eligible proteins | Median CV original | Median CV corrected | Median paired Δ (pp) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Internal-QC | 6 | 672 | 56.51% | 56.51% | +0.00 |
| pooled QC | 15 | 5966 | 10.33% | 10.33% | +0.00 |
| HC | 27 | 5581 | 23.35% | 39.23% | +16.84 |
| S | 27 | 5822 | 28.71% | 45.84% | +17.85 |

PCA feature counts/variance:

- all_groups: n=75, proteins=5585; before PC1/PC2 22.3%/16.0%; after 34.9%/12.1%.
- study_only: n=54, proteins=5662; before PC1/PC2 27.2%/17.0%; after 43.5%/15.1%.

**all_groups:**
![Sample protein distributions](b_factor1/figures/all_groups/01_sample_boxplots.png)

![Group boxplots](b_factor1/figures/all_groups/02_group_boxplots.png)

![PCA](b_factor1/figures/all_groups/03_pca.png)

![CV violin box jitter](b_factor1/figures/all_groups/04_cv_violin_box.png)

![Applied log2 sample offsets](b_factor1/figures/all_groups/05_applied_offsets.png)

**study_only:**
![Sample protein distributions](b_factor1/figures/study_only/01_sample_boxplots.png)

![Group boxplots](b_factor1/figures/study_only/02_group_boxplots.png)

![PCA](b_factor1/figures/study_only/03_pca.png)

![CV violin box jitter](b_factor1/figures/study_only/04_cv_violin_box.png)

![Applied log2 sample offsets](b_factor1/figures/study_only/05_applied_offsets.png)

### b_factor1_plus_factor2

| Group | n | Eligible proteins | Median CV original | Median CV corrected | Median paired Δ (pp) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Internal-QC | 6 | 672 | 56.51% | 56.51% | +0.00 |
| pooled QC | 15 | 5966 | 10.33% | 10.33% | +0.00 |
| HC | 27 | 5581 | 23.35% | 37.70% | +16.04 |
| S | 27 | 5822 | 28.71% | 57.34% | +26.81 |

PCA feature counts/variance:

- all_groups: n=75, proteins=5585; before PC1/PC2 22.3%/16.0%; after 42.1%/11.1%.
- study_only: n=54, proteins=5662; before PC1/PC2 27.2%/17.0%; after 52.5%/12.7%.

**all_groups:**
![Sample protein distributions](b_factor1_plus_factor2/figures/all_groups/01_sample_boxplots.png)

![Group boxplots](b_factor1_plus_factor2/figures/all_groups/02_group_boxplots.png)

![PCA](b_factor1_plus_factor2/figures/all_groups/03_pca.png)

![CV violin box jitter](b_factor1_plus_factor2/figures/all_groups/04_cv_violin_box.png)

![Applied log2 sample offsets](b_factor1_plus_factor2/figures/all_groups/05_applied_offsets.png)

**study_only:**
![Sample protein distributions](b_factor1_plus_factor2/figures/study_only/01_sample_boxplots.png)

![Group boxplots](b_factor1_plus_factor2/figures/study_only/02_group_boxplots.png)

![PCA](b_factor1_plus_factor2/figures/study_only/03_pca.png)

![CV violin box jitter](b_factor1_plus_factor2/figures/study_only/04_cv_violin_box.png)

![Applied log2 sample offsets](b_factor1_plus_factor2/figures/study_only/05_applied_offsets.png)

## Numeric results

- yeast_factor_details.csv: both peptide log2 intensities, per-sample means,
  peptide medians, method-a and method-b reference levels/factor1 values.
- help_factor_details.csv: per-sample HELP log2 median, pooled-QC ref_mR,
  candidate factor2 and number of observed HELP ratios.
- factor_summary.csv: per-sample candidate factor1a/b, factor2, and actual
  offset for each of the four scenarios; zero offsets for both QC types.
- comparisons.csv: median protein CV per group per scenario (paired proteins).
- <scenario>/applied_factors.csv: log2 and linear applied multipliers.
- <scenario>/pg_matrix_corrected_log2.tsv: log2 corrected protein matrix.
- <scenario>/pg_matrix_corrected_linear.tsv: back-transformed matrix for CV.
- <scenario>/per_protein_cv.csv: per-protein/group original and corrected CV.
- <scenario>/group_cv_summary.csv: CV median comparison.
- <scenario>/pca_scores.csv and pca_summary.csv: per-view PCA results.
- <scenario>/sample_log2_medians.csv: values underlying group boxplots.

## Statistical definitions and safeguards

- All offsets, yeast factors, HELP reference/factor2, boxplots and PCA
  calculations are on log2 intensities/ratios.
- CV (%) = sample SD(ddof=1) / arithmetic mean of **linear** positive intensities ×100;
  require detected n >= max(3, ceil(70% of group n)) and mean >= 1e-8.
- Missing/nonpositive protein observations are not imputed for CV and
  positive detections are unchanged by the multiplicative correction.
- PCA: choose >=70%-detected proteins within each view; use identical
  protein set before and after in that view, median-impute missing log2
  values for PCA only, center but do not z-scale features.
- all_groups PCA and study_only PCA are **separately fitted** and do not
  have directly comparable PC coordinate axes.
- Both QC groups retain original linear intensity values. LM duplicate
  pairs are intentionally preserved; nominal n=6 is not independent n=6.
- Biological variation influences HC and S CV as well as technical factors.
- This design cannot establish biological validity of either factor without
  further evaluation (QC stability, yeast spike-in, HC/S biology).
