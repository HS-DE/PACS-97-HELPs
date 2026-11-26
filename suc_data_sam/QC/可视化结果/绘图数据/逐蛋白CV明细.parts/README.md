# Rebuild the original PACS per-protein CV detail table

The file is split into four sequential TSV parts to fit the GitHub upload transport. Concatenate in order:

```bash
cat part1.tsv part2.tsv part3.tsv part4.tsv > ../逐蛋白CV明细.tsv
```

Part 1 contains the header; parts 2–4 continue the rows without repeating the header.
