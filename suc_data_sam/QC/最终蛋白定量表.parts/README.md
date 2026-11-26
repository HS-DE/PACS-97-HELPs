# QC final matrix parts

The QC-run final protein matrix is stored as four sequential TSV parts because the single text blob exceeded the upload transport limit.

Rebuild the original matrix by concatenating the files in order:

```bash
cat part1.tsv part2.tsv part3.tsv part4.tsv > ../最终蛋白定量表.tsv
```

- `part1.tsv` contains the column header and first data rows.
- `part2.tsv`–`part4.tsv` continue the data rows and do not repeat the header.
- The original row order and values are preserved.
