# ENO1_YEAST separation and EN01.xlsx metadata filtering

Data preparation only; no correction, statistical analysis,
normalization, protein intensity conversion or new trial results.

## Metadata selection

- Study: HC 27, S 27
- Internal-QC (LM): 6
- pooled QC (QC1): 15
- Total approved sample IDs: 75; Blank and other IDs removed

## pg_matrix.tsv

- Protein rows: 6677 → 6675.
- Removed exactly two yeast-enolase rows, ENO1_YEAST and ENO2_YEAST.
- All other protein rows, sample columns, and original TSV values preserved.
- Removed Protein.Group identifiers and original line numbers:
  - line 881, Protein.Group: P00924, Protein.Names: ENO1_YEAST
  - line 882, Protein.Group: P00925, Protein.Names: ENO2_YEAST
- Input SHA-256: 1890ae9651b67e8ec944f8a72f76579e94b9660b1ff37fe0111e9d3ec1f787c9
- Output SHA-256: 112d59739481b2fbf1d5252d0fc754c4a4494af60c77ff050af153fe423012be

## EN01.xlsx

- Sheet: Sheet1
- Protein rows: 2; source IDs P00924_ENO1_YEAST and P00925_ENO1_YEAST; retained verbatim.
- Note: pg_matrix.tsv calls P00925 ENO2_YEAST, while EN01.xlsx calls the same accession P00925_ENO1_YEAST. Names were NOT rewritten.
- Sample columns: 109 → 75.
- Removed sample columns: 34.
- Preserved the protein identifier column and numeric value strings
  of all 75 retained samples for both yeast-enolase protein rows.
- Input SHA-256: 3f8247c2be21a8a20d770ec464611a62ab9b7f631f44ec5a3134d9e272e6e151
- Output SHA-256: 08122d5a044689a7863cb85768d96797566d976a76925ecd3a71bdedcef09c45

### Removed sample IDs

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

## Other inputs and existing results

- sample-QC-InternalQC-HELP_yeast.xlsx was not modified.
- Other metadata and HELP files were not modified.
- results/, results_v1/ and results_v2/ were not modified.
- Prior result summaries describe the earlier protein matrix.
- Original source versions remain available in Git history.
