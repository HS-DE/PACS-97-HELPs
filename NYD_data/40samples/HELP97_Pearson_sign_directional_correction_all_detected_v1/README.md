# NYD: Every fully detected Heavy HELP with protein-specific Pearson-sign correction

## Protocol

Primary: all **40 original NYD samples**; each Heavy peptide must have 40/40 detected.
Separately rerun after excluding known T003 outlier MAP26092326E01 (39 samples);
do not silently mix candidate peptide sets, Pearson references, or CV results.

For one peptide j and sample i:

- h_ji = log2 Heavy intensity from the single peptide measurement.
- ref_j = median over samples of h_ji.
- fi = ref_j - h_ji (signed; distance not absolute value).
- For each protein g, Pearson r_gj is computed across sample axis
  between its original log2 PG signal and log2 Heavy signal j, pairwise available.
- Positive Pearson: add fi to original log2 PG_gi.
- Negative Pearson: subtract fi from original log2 PG_gi.
- Undefined or exact zero: keep original PG_gi unchanged (no invented direction).

IMPORTANT: a sample only has one Heavy reading for one peptide, so
a per-sample Pearson is undefined. Signs attach to protein groups,
not to individual samples. Missing original protein values stay missing.

## Results

- Primary fully detected: 74/97 Heavy HELP.
- Sensitivity 39 fully detected: 91/97 Heavy HELP.
- Primary HELP variants improving overall within-cohort median PG CV: 0/74.

## Main ranking (full 40)

| HELP | raw median CV | corrected median CV | change (percentage points) | positive proteins | negative proteins |
|---|---:|---:|---:|---:|---:|
| LLIYWASTR | 25.18 | 26.05 | +0.87 | 2327 | 1825 |
| ALDFAVGEYNK | 25.18 | 26.72 | +1.54 | 2221 | 1931 |
| YWGVASFLQK | 25.18 | 26.76 | +1.58 | 2193 | 1959 |
| SLDFTELDVAAEK | 25.18 | 26.81 | +1.63 | 2372 | 1780 |
| LSITGTYDLK | 25.18 | 26.97 | +1.79 | 2283 | 1869 |
| IVIEYVDR | 25.18 | 27.25 | +2.07 | 2190 | 1962 |
| LVVLPFPK | 25.18 | 27.28 | +2.10 | 2057 | 2095 |
| DGAGDVAFVK | 25.18 | 27.28 | +2.10 | 2106 | 2046 |

## Outputs

- full40/01_HELP_all_candidates_metrics_sorted_by_deltaCV.csv
- full40/02_every_HELP_per_sample_factors.csv.gz
- full40/03_every_HELP_per_protein_Pearson_signs.csv.gz
- full40/per_HELP/<sequence>/corrected_PG_log2.tsv.gz
- full40/per_HELP/<sequence>/01_Raw_vs_corrected_boxplots.png
- full40/per_HELP/<sequence>/02_Heavy_factor_and_PG_sample_medians.png
- full40/per_HELP/<sequence>/03_Raw_vs_corrected_common_PCA.png
- full40/00_overview_deltaCV_all_HELP.png
- exclude_T003_bad_sample_39/01_HELP_all_candidates_metrics_sorted_by_deltaCV.csv
- exclude_T003_bad_sample_39/02_every_HELP_per_sample_factors.csv.gz
- exclude_T003_bad_sample_39/03_every_HELP_per_protein_Pearson_signs.csv.gz

## Limitations

- This method chooses correction direction from raw protein-Heavy correlations
  **on the same clinical samples**; that can encode biology and is not
  independent evidence of removal of technical batch effects.
- CV over 40 clinical individuals is not replicate-only technical CV.
- Corrected sample-to-sample Pearson correlations need not remain unchanged:
  signs vary by protein so transformation is NOT one constant per sample.
- No QC or calibration reference beyond clinical samples is available.
- No source files or earlier results overwritten.
