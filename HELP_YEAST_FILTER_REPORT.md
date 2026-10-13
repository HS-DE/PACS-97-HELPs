# Metadata filtering of 2026-10-13 HELP and yeast workbooks

Using all_sample_metadata.xlsx and all_QC_metadata.xlsx, retain 75 samples:
27 HC + 27 S + 6 Internal-QC (LM) + 15 pooled QC (QC1).
Blank/unlisted sample columns are removed.

Only sample columns are removed; no log2 transformation, normalization,
imputation, ratio recomputation or peptide filtering is performed.
Every retained sample/peptide value was compared to the original before
modification; peptide sequences and three HELP sheets are preserved.
Two yeast peptides, AADALLLK and VNQIGTLSESIK, remain independent.
Original versions remain recoverable from the repository's Git history.

## LH_ratio_tables_2026.10.13—ALL-HELP.xlsx

- Worksheets: L, H, ratio_LH
- Peptide rows: 97 (unchanged)
- Sample columns: 109 -> 75
- Removed columns: 34
- Status: filtered
- Previous SHA256: 9c3a493577fab8b5ef988b248df650e60696905a770da0f3346c7651adc39e4a
- Current SHA256: 5523ffca836fd7244610322f781084619541482835f108a92cdb350a6e3a4e53
- Previous size: 323819 bytes
- Current size: 204954 bytes

### Removed sample IDs

- k_26_9_2_1
- k_26_9_2_10
- k_26_9_2_101
- k_26_9_2_103
- k_26_9_2_12
- k_26_9_2_15
- k_26_9_2_17
- k_26_9_2_19
- k_26_9_2_20
- k_26_9_2_25
- k_26_9_2_27
- k_26_9_2_28
- k_26_9_2_29
- k_26_9_2_4
- k_26_9_2_40
- k_26_9_2_41
- k_26_9_2_42
- k_26_9_2_45
- k_26_9_2_5
- k_26_9_2_50
- k_26_9_2_53
- k_26_9_2_56
- k_26_9_2_6
- k_26_9_2_64
- k_26_9_2_65
- k_26_9_2_67
- k_26_9_2_70
- k_26_9_2_71
- k_26_9_2_72
- k_26_9_2_73
- k_26_9_2_79
- k_26_9_2_80
- k_26_9_2_81
- k_26_9_2_89

## sample-QC-InternalQC-HELP_yeast.xlsx

- Worksheets: L
- Peptide rows: 2 (unchanged)
- Sample columns: 75 -> 75
- Removed columns: 0
- Status: already filtered
- Previous SHA256: 4813ebd8387ee704e55a2c9b2ae5ed8cebd86d3e3db8c4792d320d0e3aaccee8
- Current SHA256: 4813ebd8387ee704e55a2c9b2ae5ed8cebd86d3e3db8c4792d320d0e3aaccee8
- Previous size: 8036 bytes
- Current size: 8036 bytes

### Removed sample IDs

- None

## Important note

Older PACS trial outputs were NOT recomputed with these updated files.
Previous scripts still pointing to the deleted 2026.9.22 HELP filename
must be updated before the next normalization analysis.
