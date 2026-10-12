# Latest pg_matrix.tsv metadata-based refilter

This run filters ONLY pg_matrix.tsv, using the existing current
all_sample_metadata.xlsx and all_QC_metadata.xlsx sample IDs.
No other input files or prior trial result files were modified by this script.

- Status: **filtered and rewritten**.
- Protein.Group rows: **6677** before and after.
- Original TSV columns: **114** (including five protein annotations).
- Retained TSV columns: **80** = five protein annotations + 75 samples.
- Original named sample columns: **109**.
- Original unnamed sample columns: **0**.
- Removed metadata-unlisted named columns: **34**.
- Retained groups: **Internal-QC (6), pooled QC (15), HC (27), S (27)**.
- Blank QC and other metadata-unlisted sample columns are excluded.
- Sample order follows the newly uploaded matrix, without reordering.
- All protein annotations and retained intensity strings are unchanged.
- No normalization, log transform, or protein-row filtering is performed.
- Uploaded matrix SHA-256: 8f2af2cf40ef2be1e96c983f148fbf0c5bc065900ed1b797c5daf2625067f169.
- Filtered matrix SHA-256: 1890ae9651b67e8ec944f8a72f76579e94b9660b1ff37fe0111e9d3ec1f787c9.
- Original size: 4901253 bytes; filtered size: 3896244 bytes.

## Excluded named sample columns

- k_26_9_2_1
- k_26_9_2_4
- k_26_9_2_5
- k_26_9_2_6
- k_26_9_2_10
- k_26_9_2_12
- k_26_9_2_15
- k_26_9_2_17
- k_26_9_2_19
- k_26_9_2_20
- k_26_9_2_25
- k_26_9_2_27
- k_26_9_2_28
- k_26_9_2_29
- k_26_9_2_40
- k_26_9_2_41
- k_26_9_2_42
- k_26_9_2_45
- k_26_9_2_50
- k_26_9_2_53
- k_26_9_2_56
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
- k_26_9_2_101
- k_26_9_2_103

Note: downstream Trial 1/2 results generated from an earlier matrix
should be considered historical until explicitly recomputed.
