# PACS five-HELP reruns — same existing methods, only HELP set reduced

Previous results/ archived intact as **results_v4/**; all previous results_vN/ retained.
The five peptides were selected explicitly by the user, not algorithmically optimized on these results.
This is a controlled re-run of the existing three types of experiments,
NOT a new batch-correction model or batch-effect audit.

## Exactly five HELP sequences, all from ratio_LH sheet

1. FSVVYAK
2. STGSWSTLK
3. VTGWGNLK
4. VVEESELAR
5. VYFAGFPR

The original 97-row 3-sheet workbook is unchanged. Exactly these five sequence
rows are selected from the existing 97, for every 75-sample HELP calculation.
Only valid finite positive ratios are included in each sample median.
The 75 × 5 ratio grid contains 0 missing/nonpositive values.
Valid HELP counts per sample range 5
to 5; none has zero usable HELP ratios.
All five-HELP raw-ratio pooled-QC reference: 0.0185423054179.
All five-HELP log2-ratio pooled-QC reference: -5.75303556026.

## Input provenance

- pg_matrix.tsv 6675 protein groups: sha256 112d59739481b2fbf1d5252d0fc754c4a4494af60c77ff050af153fe423012be.
- Full HELP workbook used only as row source: LH_ratio_tables_2026.10.13—ALL-HELP.xlsx: sha256 5523ffca836fd7244610322f781084619541482835f108a92cdb350a6e3a4e53.
- Yeast peptide source, unchanged: sha256 4813ebd8387ee704e55a2c9b2ae5ed8cebd86d3e3db8c4792d320d0e3aaccee8.
- Yeast protein source, unchanged: sha256 08122d5a044689a7863cb85768d96797566d976a76925ecd3a71bdedcef09c45.
- Metadata: 75 samples, HC=27, S=27, Internal-QC=6, pooled QC=15.

## Three retained algorithms re-executed

| Folder | Exact inherited procedure | QC treatment |
| --- | --- | --- |
| [trial1b](trial1b/README.md) | HELP-only median(raw L/H), pooled-QC reference; two PCA views | Internal-QC unchanged, pooled QC corrected |
| [trial2_peptide](trial2_peptide/README.md) | Yeast peptide Factor1 methods a/b (unchanged) + median(log2 HELP L/H) Factor2; Raw/F1/F1+F2 | Both QC groups unchanged |
| [trial3_eno1](trial3_eno1/README.md) | EN01 protein Factor1: higher row sum vs linear sum of both (unchanged) + median(log2 HELP L/H) Factor2 | Both QC groups unchanged |

HELP-only correction uses median(raw L/H), whereas two-step Trial2/3
methods use median(log2 L/H); these are distinct calculations.

## Trial 2 and EN01 median protein CV (%)

CV on linear protein intensities, HC/S include biological variation;
unchanged QC CV in Trial2/3 cannot independently validate normalization.

| Trial | Method | Stage | HC | S |
| --- | --- | --- | ---: | ---: |
| trial2_peptide | a | raw | 23.04% | 28.82% |
| trial2_peptide | a | factor1 | 39.25% | 45.72% |
| trial2_peptide | a | factor1_factor2 | 41.87% | 56.01% |
| trial2_peptide | b | raw | 23.04% | 28.82% |
| trial2_peptide | b | factor1 | 39.25% | 45.72% |
| trial2_peptide | b | factor1_factor2 | 41.87% | 56.01% |
| trial3_eno1 | single_max_rowsum | raw | 23.04% | 28.82% |
| trial3_eno1 | single_max_rowsum | factor1 | 49.10% | 68.95% |
| trial3_eno1 | single_max_rowsum | factor1_factor2 | 54.57% | 85.70% |
| trial3_eno1 | sum_two | raw | 23.04% | 28.82% |
| trial3_eno1 | sum_two | factor1 | 49.01% | 68.60% |
| trial3_eno1 | sum_two | factor1_factor2 | 54.55% | 85.51% |

## 5-vs-97 historical context — interpretation limit

See [five_vs_97_historical_cv.csv](five_vs_97_historical_cv.csv).
Historical Trial2 97-HELP data are read only from results_v3/;
historical EN01 97-HELP data are read only from results_v4/.
**Comparison design:** Prior Trial2 in results_v3 contains 6677 proteins,
whereas the current data contain 6675 after removing two yeast protein rows.
Thus Trial2 changes are descriptive rather than a strictly isolated 5-vs-97 test.
Prior EN01 Trial3 in results_v4 and current Trial3 both use the same 6675
protein matrix, so Trial3 is directly comparable for the HELP-subset change.
Factor1-only results should not change within a fixed protein set, because
the HELP subset affects only Factor2. Do not claim causal batch reduction.
Raw, Factor1, and Factor1+Factor2 are shown together for both Trial2/3.

## Files and reproducibility

- [five_help_sample_diagnostics.csv](five_help_sample_diagnostics.csv): 5 requested
  ratios, valid-count, raw/log2 medians, and QC/reference per 75 samples.
- [five_vs_97_historical_cv.csv](five_vs_97_historical_cv.csv): prior-run comparison
  with clear protein-row caveat.
- [run_manifest.json](run_manifest.json): sources, selected sequences, methods.
- Each trial directory has its complete matrices, factors, CVs, PCA and original plots.
- Replay only via scripts/pacs_help5_reruns.py; the historical scripts remain
  configured for all 97 as a preserved historical baseline.
