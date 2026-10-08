# 酵母 spike-in 蛋白 QC｜Yeast spike-in protein QC

此目录结果由 [`scripts/yeast_protein_qc.py`](../../scripts/yeast_protein_qc.py) 从更新后的 [`pg_matrix.tsv`](../../pg_matrix.tsv)、`all_sample_metadata.xlsx` 和 `all_QC_metadata.xlsx` 生成。**无插补、无归一化，CV 按原始强度计算。**

## 样本、统计单位与总体检出

- 75 个 metadata 合法样本：Internal-QC 6，pooled QC（QC1）15，HC 27，S 27；未发现多余样本列和重复样本 ID。
- 6,830 条蛋白组定量记录中，**217 条** `Protein.Names` 含有 `YEAST`。
- 217 条酵母注释记录中，**216 条**在至少一个保留样本内检出（有限且大于 0 的强度）。
- 仅 `FADH_YEAST`（`Protein.Group=P32771`）在全部 75 个样本中均未检出。
- 统计单位是 `pg_matrix.tsv` 的一行 **Protein.Group**；带分号的 `Protein.Names` 条目按一个定量蛋白组计数，**不能解释为逐个同源蛋白均单独被鉴定**。此数据中 217 条记录对应 217 个不同的 `Protein.Names` 字符串。

## 四组结果

| 组别 | 样本 n | ≥1 个组内样本检出的蛋白组 | 每样本平均检出 | 可计算单蛋白 CV 数量 | 单蛋白 CV 中位数 (%) | 单样本检出数量的 CV (%) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Internal-QC | 6 | 24 | 13.50 | 10 | 68.34 | 31.69 |
| pooled QC | 15 | 210 | 147.93 | 128 | 18.38 | 4.98 |
| HC | 27 | 212 | 154.56 | 134 | 49.34 | 10.33 |
| S | 27 | 216 | 149.11 | 125 | 66.52 | 10.93 |

**特别说明：** Internal-QC 的酵母检出覆盖率极低，且 `LM_2_1/LM_2_2`、`LM_3_2/LM_3_3` 为已确认的**有意重复**。本表依 metadata 保留名义 n=6 计算描述性指标，但不能把它视为 6 个独立测量重复。HC/S 组的 CV 也可能同时包含真实生物差异，不应直接等同于纯技术误差。

## 检出与 CV 计算口径

1. 筛选 `Protein.Names` **包含 YEAST** 的行（不区分大小写）。
2. 单样本检出：强度为**有限、正数**。NA、空白、0、负数、Inf 均按未检出处理。
3. **组内检出率** = 组内检出样本数 / 该组 metadata 样本总数 × 100%。检测率分母包括未检出样本。
4. **单蛋白组内 CV (%)** = `sample SD (ddof=1) / mean × 100`，仅使用该组内已检出的原始强度，**不先取 log2**。
5. 参照所附 `CV作图(1).R` 默认设置，CV 计算条件为：有效点 `n_detected ≥ max(3, ceil(0.70 × n_group))` 且 `mean ≥ 1e-8`。阈值分别为 Internal-QC 5/6、pooled QC 11/15、HC 19/27、S 19/27。不达标的 CV 记为缺失，**不以 0 代替**。
6. 每样本检出数量的 CV 使用同组所有样本的每样本检出种类数量（包括 0），同样使用样本 SD/均值 ×100%。
7. 不同组内的检出总数可能重叠，不能将四组蛋白数量相加得到全体蛋白数。

## 可视化

四组颜色在所有图中一致：**蓝色 Internal-QC、橙色 pooled QC、绿色 HC、紫色 S**。

### 一、各组检出酵母蛋白数（柱状图）

![各组检出数](figures/01_group_detected_counts.png)

### 二、每种酵母蛋白在四组中的检出率（分组横向柱状图）

按 30 条蛋白组分页，包含所有 217 条记录：

[第 1 页](figures/02_per_protein_detection_rates_01.png) · [第 2 页](figures/02_per_protein_detection_rates_02.png) · [第 3 页](figures/02_per_protein_detection_rates_03.png) · [第 4 页](figures/02_per_protein_detection_rates_04.png) · [第 5 页](figures/02_per_protein_detection_rates_05.png) · [第 6 页](figures/02_per_protein_detection_rates_06.png) · [第 7 页](figures/02_per_protein_detection_rates_07.png) · [第 8 页](figures/02_per_protein_detection_rates_08.png)

### 三、每种酵母蛋白在四组中的 CV（小提琴图 + 箱线图 + 散点）

![CV 分布](figures/03_group_cv_violin_box.png)

具体每条蛋白组的四组 CV 以散点图分页展示（CV 不满足计算门槛则不显示点）：

[第 1 页](figures/04_per_protein_cv_01.png) · [第 2 页](figures/04_per_protein_cv_02.png) · [第 3 页](figures/04_per_protein_cv_03.png) · [第 4 页](figures/04_per_protein_cv_04.png) · [第 5 页](figures/04_per_protein_cv_05.png) · [第 6 页](figures/04_per_protein_cv_06.png) · [第 7 页](figures/04_per_protein_cv_07.png) · [第 8 页](figures/04_per_protein_cv_08.png)

### 四、补充图（非 CV 均用柱状图）

![每样本平均检出数量](figures/05_group_mean_identification.png)

![各样本检出数量](figures/06_per_sample_identification.png)

![可计算 CV 的蛋白组数量](figures/07_cv_eligible_counts.png)

所有图片同时包含 **SVG 矢量版**，可在 `figures/` 中下载用于论文或汇报。

## 完整统计文件与复现

- [`group_summary.csv`](group_summary.csv)：四组总体检出数、每样本均值、中位数、可估计 CV 数及 CV 汇总。
- [`per_protein.csv`](per_protein.csv)：**217 行**酵母蛋白组；四组各自的检出数、检出率百分比和 CV 百分比。
- [`per_sample.csv`](per_sample.csv)：**75 行**每个样本检出的酵母蛋白组数。
- [`scripts/yeast_protein_qc.py`](../../scripts/yeast_protein_qc.py)：完整可复现的 Python 分析及作图脚本。
- [GitHub Actions workflow](../../.github/workflows/yeast-protein-qc.yml)：数据或脚本更新后自动重算并提交结果。

本目录输出是 raw intensity 的描述性质量检查，并非 PACS 归一化结果。
