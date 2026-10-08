# 97-HELP metadata-based sample filtering

Retain all 54 study samples (27 HC + 27 S), 6 Internal-QC, and 15 QC1 pooled plasma; exclude 6 Blank QC records.
All three HELP sheets retain 97 peptides x 75 sample columns.
pg_matrix.tsv retains 6613 protein rows, 5 annotation columns + 75 sample columns.
Original HELP sample columns: 109; removed 34. Original pg_matrix: 87 named sample columns plus 22 unnamed trailing columns; removed 12 named and 22 unnamed.
This operation only deletes unused columns and Blank metadata rows; it does NOT alter retained intensity or ratio values.

## Removed HELP IDs

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

## Removed pg_matrix IDs

- k_26_9_2_1
- k_26_9_2_4
- k_26_9_2_27
- k_26_9_2_40
- k_26_9_2_41
- k_26_9_2_42
- k_26_9_2_71
- k_26_9_2_80
- k_26_9_2_81
- k_26_9_2_89
- k_26_9_2_101
- k_26_9_2_103
