# suc_data_sam — PACS dual-QC and HELP latent-factor evaluation

This folder contains the analysis of the uploaded **`Plate1-16_QC_Hela_merged.pg_matrix.tsv`** matrix with the supplied HELP CSV and metadata CSV. It is a new run and is kept separate from the older `suc_data/` experiment, which used a different protein matrix.

## Scope

Two independent PACS runs were completed:

- **LM-QC run:** only `group=LM-QC` is mapped to pooled `qc`.
- **QC run:** only `group=QC` is mapped to pooled `qc`.

In both runs:
- `group=sample` remains biological/study samples.
- `Internal-QC` and the non-selected QC class are mapped to `internal-calibrator` and excluded by the PACS input logic.
- `Blank` remains `blank` and is excluded.
- No missing protein values are imputed.

The original input files are not copied into this public folder. Derived normalized matrices and analysis outputs are included.

## Key result

On this matrix, **LM-QC is the more stable pooled-QC definition**: its median protein CV changed from 16.81% to 14.95% after PACS. The QC-class run started at 41.43% and changed to 40.17%, while PACS also emitted HELP-instability and high-QC-variability warnings.

The exploratory HELP latent-factor prototype can predict QC global PG shifts in plate-grouped cross-validation, but applying its predicted shifts did **not** improve protein-level QC CV over PACS. It is therefore diagnostic only and is not promoted as a final normalization.

## Folder map

```text
suc_data_sam/
├── README.md
├── RESULTS.md
├── run_manifest.json
├── scripts/
│   ├── PACS_norm.py
│   ├── run_two_qc_pacs.py
│   ├── help_gmf_diagnostics.py
│   └── requirements.txt
├── LM-QC/
│   ├── 最终蛋白定量表.tsv
│   ├── 蛋白鉴定量统计柱状图.png
│   └── 可视化结果/
├── QC/
│   ├── 最终蛋白定量表.tsv
│   ├── 蛋白鉴定量统计柱状图.png
│   └── 可视化结果/
├── HELP_GMF/
│   ├── LM-QC/
│   └── QC/
└── comparison/
    ├── *.tsv
    └── figures/
```

Each PACS result directory contains the final protein matrix (the QC-run matrix is split into four sequential parts under `QC/最终蛋白定量表.parts/`), protein-count chart, before/after boxplots, PCA, CV plots, QC CV histogram, and plotting data.

The full 32k-row PACS `逐蛋白CV明细.tsv` tables are also split into four sequential parts under each run's `可视化结果/绘图数据/逐蛋白CV明细.parts/`. Each parts folder has a README with the concatenation command. Part 1 contains the header; later parts do not repeat it.

The source PACS script writes two QC CV plotting-data files in GBK. UTF-8 copies ending in `_utf8.tsv` are included for reliable viewing on GitHub.

## Reproduce

Install dependencies:

```bash
python -m pip install -r scripts/requirements.txt
```

Run PACS twice:

```bash
python scripts/run_two_qc_pacs.py \
  --protein-table /path/to/Plate1-16_QC_Hela_merged.pg_matrix.tsv \
  --help-table /path/to/Accurate_Quantification_Report_v4_苏创.csv \
  --metadata /path/to/metadata.csv \
  --out-dir ./suc_data_sam
```

Run the exploratory HELP latent-factor diagnostics:

```bash
python scripts/help_gmf_diagnostics.py \
  --protein-table /path/to/Plate1-16_QC_Hela_merged.pg_matrix.tsv \
  --help-table /path/to/Accurate_Quantification_Report_v4_苏创.csv \
  --metadata /path/to/metadata.csv \
  --out-dir ./suc_data_sam/HELP_GMF_reproduced
```

## Important interpretation notes

1. The PACS script groups by `group` and `condition`. The `batch`/`plate` field is checked and used for plotting, but is **not** used to estimate plate-specific offsets in the current PACS correction.
2. The HELP latent-factor analysis is a transparent PCA-based exploratory prototype inspired by factor models. It is **not** an exact implementation or reproduction of omicsGFM.
3. The QC factor model is trained on the selected QC class only. Its grouped cross-validation keeps all QC replicates from a plate in the same fold.
4. The 32 study samples all have `condition=sample`; their CV includes biological heterogeneity and must not be interpreted as pure technical CV.
5. Do not use the exploratory HELP-GMF candidate shifts as a production normalization. The held-out QC protein CV did not improve over PACS.
