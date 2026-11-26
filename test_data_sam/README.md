# 蛋白质组学数据说明与处理指南

Proteomics Data Description and Processing Guide

包含 LC-MS 蛋白质组学分析所需的三类数据：HELP 重标肽定量报告、样本信息表和蛋白组定量矩阵。它们分别提供外加技术内标信号、实验设计信息和内源性蛋白组定量信号，可用于数据质量控制、技术变异分析、归一化（normalization）、批次效应校正（batch correction）、缺失值分析、统计分析及机器学习建模。本说明不依赖特定分析软件或校正算法，可用于 Python、R 或其他蛋白质组学分析流程。

## 1. 文件概览

| 文件名模式 | 内容 | 数据结构 | 主要用途 |
| --- | --- | --- | --- |
| `Accurate_Quantification_Report*` | HELP 重标肽定量及质控报告 | 长表：样本 × 肽段记录 | 技术内标质量评估、技术偏移估计 |
| `sampleinfo*` / `metadata*` | 样本信息与实验设计 | 每行一个样本 | 样本匹配、QC 识别、分组与批次分析 |
| `*_pg_matrix.*` / `pg_matrix*` | DIA-NN 蛋白组定量矩阵 | 蛋白组 × 样本宽表 | 蛋白定量、归一化、统计及建模 |

文件名及后缀可能随项目变化，读取时应以实际表头、工作表名称及数据内容为准。常见格式包括 XLSX、CSV 和 TSV。

## 2. `Accurate_Quantification_Report*`：HELP 重标肽定量报告

### 2.1 数据是什么？

HELP（本项目使用的稳定同位素标记重肽 panel）为外加技术内标。本项目 panel 包含 97 条目标肽段，但每个样本实际通过质量筛选的肽段数量可能不同。报告记录其色谱峰积分面积（AUC）、峰检测质量、保留时间及碎片离子等信息。

HELP 的作用是提供独立于内源性蛋白丰度的技术参考。在加样量固定、定量响应可靠的条件下，其信号变化可用于评估 LC-MS 响应、上机量、保留时间漂移、离子抑制等技术因素，但不能单独确定这些因素的具体来源。如果 HELP 在消化后加入，它不能直接监控加入之前的提取、富集或消化损失。HELP 加样误差本身也可能改变信号，不应直接将所有 HELP 波动等同于内源性蛋白的技术偏移。

### 2.2 数据结构及关键字段

Excel 报告优先读取 `Peptide_Long` 工作表；CSV/TSV 通常对应同结构的长表。每行是一个样本中的一个目标肽段定量记录。

| 字段 | 含义 |
| --- | --- |
| `Sample` | 样本名称或路径 |
| `Peptide` | 肽段氨基酸序列 |
| `Peptide_Key`（如有） | 肽段定量目标标识，用于区分不同目标形式 |
| `Protein` | 目标肽段关联的蛋白注释，不代表内源性蛋白丰度 |
| `Precursor_mz` / `Precursor_charge` | 前体离子质荷比 / 电荷数 |
| `Accurate_AUC` | HELP 色谱峰积分面积，主要定量字段 |
| `Peak_Detection_QC_Flag` | 峰检测质量标记，例如 `PASS` |
| `Peak_Scan_Points` | 色谱峰内扫描点数 |
| `N_Common_Transitions` | 共同碎片离子 transitions 数量 |
| `Reference_RT_min` / `Apex_RT_min` | 参考保留时间 / 实测峰顶保留时间（分钟） |
| `AUC_Time_Unit` | 积分时间单位 |
| `QC_Flag` / `Failure_Reason` | 其他质控标记 / 失败原因（如有） |

某些 Excel 报告还包含 `Peptide_Wide`、`Peptide_QC`、`Transition_Long`、`Sample_QC`、`Run_Config` 等工作表，可用于进一步诊断。不同版本字段不一定完全一致；不能将 `Detected_AUC` 或 `QC_Flag` 无条件视为 `Accurate_AUC` 或 `Peak_Detection_QC_Flag` 的等价字段。

### 2.3 如何处理？

1. 读取 `Peptide_Long`（Excel）或相同结构的 CSV/TSV 长表，检查必需字段及数值格式。
2. 统一 `Sample` 的命名，清理路径、首尾空格及末尾 `.raw` 后缀，与 metadata 的 `sample_id` 对齐。保留原始名称与标准化 ID 的映射。
3. 确认肽段标识。如果一个 `Peptide` 序列对应多个电荷态或定量目标，应使用 `Peptide_Key` 等能够唯一识别目标的字段，不能盲目合并重复记录。
4. 筛选有效定量记录。本项目可采用以下现有质量筛选规则作为起点：
   - `Peak_Detection_QC_Flag == PASS`；
   - `N_Common_Transitions >= 2`；
   - `Peak_Scan_Points >= 8`；
   - `Accurate_AUC` 为有限且大于 0 的数值。
   - 如严格筛选后可用内标不足，可预先定义放宽至 7 个扫描点的敏感性分析，并记录原因及影响。阈值是当前报告的经验性规则，不是所有 LC-MS 实验的通用标准。
5. 检查同一个样本与同一个肽段目标的重复记录。确认重复来源后再处理，不应默认求和或取平均。
6. 将长表转换为 HELP 宽矩阵，行是肽段，列是样本，数值是 `Accurate_AUC`。未检出、未通过筛选或无有效强度的位置保留为 NA，而非 0。
7. 统计每条 HELP 的检出率、每个样本的有效 HELP 数量，以及每类同质 QC 中的 HELP CV。可进一步检查 RT 偏移、峰形和极端 AUC。
8. 根据分析方法选取 HELP：单内标方法可选择覆盖率高且重复性好的肽段；多内标方法可使用符合预设覆盖率和质控标准的肽段；矩阵分解方法可在训练 QC 中估计多个技术因子。无需所有方法都强制使用全样本完整检出的肽段。
9. 按需要保留线性矩阵和 log2 矩阵；log2 仅对正强度计算。不同积分单位、HELP 配方或加样量的报告应确认可比性后再合并。

示意（非真实强度）：

| Peptide | Sample_01 | Sample_02 | QC_01 |
| --- | ---: | ---: | ---: |
| HELP_A | 100000 | 105000 | 98000 |
| HELP_B | 250000 | NA | 245000 |

建议输出：`help_matrix.tsv`、`help_quality.tsv`、`help_log2_matrix.tsv`、筛选参数及排除记录。

## 3. `sampleinfo*` / `metadata*`：样本信息表

### 3.1 数据是什么？

该表记录每个样本的唯一身份、样本类型、生物学分组和实验批次。它是连接 HELP 报告与蛋白定量矩阵的索引，也是确定校正参考、验证样本、统计设计和批次变量的依据。

| 字段 | 含义 |
| --- | --- |
| `sample_id` | 唯一样本 ID |
| `group` | 样本类型，如 QC、sample、internal-calibrator、blank |
| `condition` | 生物学分组，如 HC、S、疾病组、处理组 |
| `batch` | 实验批次 |
| `Plate` | 实验板号，与 batch 不一定相同 |
| 其他可选字段 | 仪器、进样顺序、制备日期、受试者 ID、重复关系等（若有） |

### 3.2 样本类型

| 类型 | 处理原则 |
| --- | --- |
| `sample` | 正式研究样本；保留其真实 condition，不应仅为了降低组内 CV 而强制对齐蛋白分布 |
| `QC` | 用于技术重复性评估；若为相同 pooled material，可用于建立参考或验证校正 |
| `LM-QC` 等其他 QC 类型 | 需要确认材料及实验阶段；可分别用于建模和评价，不默认与其他 QC 合并 |
| `Internal-QC` / `internal-calibrator` | 根据实验目的决定是否用于校准或监控；不是所有流程都需要纳入 |
| `Blank` | 用于背景、污染或 carryover 检查，通常不纳入生物学分析或常规 pooled QC CV |

### 3.3 如何处理？

1. 读取样本表，统一表头的大小写及首尾空格，例如将 `Plate` 统一为 `plate`。
2. 将 `sample_id` 作为字符串读取，保留前导零，检查是否存在空值或重复 ID。
3. 标准化 group 标签（如 `QC` 和 `qc`）；保存原始标签。不同 QC 类型必须保留其区别。
4. 根据实验目的定义分析样本、校正参考样本、独立验证样本及排除样本，记录选择原因。
5. 检查纳入样本的 condition、batch、plate 等建模所需字段。若无板号，不应自行假设所有 batch 都等于 plate。
6. 统计各 group、condition、batch 和 plate 的样本数量，并检查生物学分组与批次是否混杂。如果某个疾病组仅存在于某块板，通常无法单靠数据可靠地区分疾病效应与板间效应。
7. 核对与 HELP 报告及蛋白矩阵的匹配情况，输出缺失、多余及冲突的样本列表，不静默丢弃样本。
8. 按蛋白矩阵的样本列顺序重排 metadata，或者明确指定统一样本顺序。
9. 若做模型验证，提前划分训练 QC 和验证 QC；可使用同板留出重复样本或整板留出，避免用验证蛋白强度估计待评价的校正参数。

建议输出：`metadata_clean.tsv`、`sample_manifest.tsv`、`sample_matching_report.tsv`、`group_batch_counts.tsv` 及验证划分表（如适用）。

## 4. `*_pg_matrix.*`：DIA-NN 蛋白组定量矩阵

### 4.1 数据是什么？

该文件是 DIA-NN 输出的 protein group 定量结果，每行对应一个 protein group，每个样本列记录该蛋白组在对应样本中的定量强度。

| 字段 | 含义 |
| --- | --- |
| `Protein.Group` | 蛋白组标识，用于唯一定位矩阵行 |
| `Protein.Ids` | 蛋白 accession / ID |
| `Protein.Names` | 蛋白名称 |
| `Genes` | 基因名称 |
| `First.Protein.Description` | 蛋白描述 |
| 其他样本定量列 | 各 LC-MS 样本的蛋白组强度，列名可能带有原始文件路径 |

Protein group 并不一定对应单个唯一蛋白；不要直接以 `Genes` 替代 `Protein.Group` 作为唯一键。当前输入示例为线性定量强度，但其他项目应检查 DIA-NN 的定量及 normalization 设置，确认是否已经归一化、转换尺度或进行过其他处理，避免重复校正。

### 4.2 如何处理？

1. 按实际分隔符读取文件，检查表头、重复列名和 `Protein.Group` 唯一性。
2. 分离蛋白注释列与样本定量列。通过 metadata 的标准化样本 ID 识别定量列，而不是假设所有非 `Protein.Group` 列都是数值。
3. 清理样本列名中的路径、空格及末尾 `.raw`，与 metadata 一一匹配。检查清理后是否发生名称冲突。
4. 依据实验设计选择所需样本，按统一样本顺序建立 protein group × sample 矩阵。
5. 将定量列转换为数值；在线性正强度数据中，空值、NA、0、负值、无穷大及非法字符串应标记为缺失或异常，并记录数量。若实际文件采用其他缺失编码或已处于 log2 尺度，应按其定义处理，不能机械地删除 0 或负数。
6. 统计每个样本检出的蛋白数、每个蛋白的检出率和整体缺失率；保留未经蛋白过滤的分析矩阵。
7. 根据下游任务定义蛋白筛选标准，例如某组内检出率达到 70% 或至少有若干有效重复。不同任务可使用不同的筛选范围，不能只依据所有样本的总体检出率删除疾病或组别特异蛋白。
8. 对有效线性强度计算 log2，保留缺失值及原始线性矩阵。多数回归、PCA 和矩阵分解可在 log2 尺度上建模。
9. 在校正前完成基线质量分析：QC CV、蛋白检出数、样本强度分布、相关性、PCA、批次或板间差异。
10. 将校正后的实测矩阵、插补矩阵及原始矩阵分别保存；保留蛋白注释及原始 ID。

建议输出：`protein_annotation.tsv`、`protein_matrix_linear.tsv`、`protein_matrix_log2.tsv`、`protein_detection.tsv` 及 `observed_mask.tsv`。

## 5. 三类数据的匹配与整合

```text
Accurate_Quantification_Report  --> HELP long table --> HELP matrix (H x N)
sampleinfo / metadata          --> cleaned sample metadata (N rows)
pg_matrix                      --> protein matrix (P x N)
                                      |
                               Sample ID matching
                                      |
                        Aligned matrices + sample metadata
                                      |
                       Baseline QC / normalization / modeling
```

- H：有效 HELP 定量目标数；N：纳入分析的样本数；P：蛋白组数。
- 两种矩阵保持相同的样本列顺序，metadata 按此顺序排列。算法若要求 sample × feature 格式，可转置后使用。
- HELP 与蛋白矩阵之间通过样本 ID 关联，不需要按 HELP 的 `Protein` 与 DIA-NN 的 `Protein.Group` 强制一一对应。
- 统一命名示例：`D:\\data\\Sample_01.raw` 与 `C:/data/Sample_01` 均可映射为 `Sample_01`，前提是该映射不会造成冲突。样本 ID 中的其他字符或大小写不要任意改变。
- 如果某个样本缺少 HELP 定量，应标记为不能执行依赖 HELP 的方法；是否保留在不依赖 HELP 的分析中取决于任务。不同方法比较时应注明样本集合是否相同。

## 6. 推荐的数据分析与校正流程（算法无关）

### Step 1 — 原始数据质量控制

按第 2–5 节完成文件读取、数据清洗和匹配。分别检查 HELP 稳定性和蛋白定量稳定性；对不同 QC 类型分别统计。不能仅凭 HELP CV 较低就认定所有蛋白均可采用相同校正因子。

### Step 2 — 建立 Raw baseline

使用未进行本轮校正的蛋白矩阵作为共同对照，记录 QC CV、检出率、缺失率、PCA、样本相关性以及板内/板间变异。记录输入数据是否已进行上游 normalization。

### Step 3 — 选择并实施 normalization / batch correction

以下方法可以独立比较，也可以按实验设计组合；并非要求顺序执行全部方法。

| 方法 | 基本处理过程 | 关键限制 |
| --- | --- | --- |
| Raw | 不校正，仅完成清洗与匹配 | 作为 baseline |
| Endogenous median normalization | 估计各样本的蛋白中位数偏移并校正 | 假设多数蛋白整体丰度稳定；可能消除真实全局生物变化 |
| Single-HELP normalization | 选择可靠的 HELP，估计样本偏移并校正 | 单条肽段可能无法代表其他蛋白 |
| Multi-HELP normalization | 每条 HELP 建立参考，计算样本偏移后进行稳健聚合 | 需确认 HELP 与蛋白技术响应一致，防止过度校正 |
| QC-based plate scaling | 比较各板同质 QC 与共同参考，估计每板整体校正因子 | 不能校正蛋白特异性偏移 |
| Protein-specific plate calibration | 对每个蛋白估计板内 QC 与共同参考的差异，应用于该板 | QC 少或材料不一致时可能过拟合；缺乏参考时不强行校正 |
| Regression / batch correction | 使用 batch、plate 或技术指标为协变量，估计并移除技术项 | 需考虑生物分组与技术变量混杂 |
| Matrix factorization / GMF | 利用观测掩码、协变量及低秩结构联合建模 | 未知潜在因子不一定是技术因素；不可全部当作噪声移除 |
| HELP-guided model | 从 HELP 提取技术因子，利用 QC 或控制实验估计蛋白响应，必要时加入正则化与板间项 | 需要独立验证，不能把训练 QC 对齐当成准确性提升 |

**示例 A：Multi-HELP sample correction（log2 尺度）**

1. 用同质 reference QC 计算每条 HELP 的参考值 `R[h] = median(H[q,h])`。
2. 计算样本相对参考的偏移 `d[i,h] = H[i,h] - R[h]`。
3. 稳健聚合有效偏移 `t[i] = median_h(d[i,h])`。
4. 若技术响应假设成立，可按 `Y_corrected[i,j] = Y[i,j] - t[i]` 校正；也可以估计不同蛋白的响应系数，并设置质量门控或缩减校正幅度。

**示例 B：QC-based protein-specific plate correction（log2 尺度）**

1. 对相同 QC 材料计算每个蛋白的 global median `G[j]`。
2. 计算每板每个蛋白的 local QC median `L[p,j]`。
3. 估计偏移 `b[p,j] = L[p,j] - G[j]`。
4. 在同板技术偏移可迁移至目标样本的假设下，使用 `Y_corrected[i,j] = Y[i,j] - b[p,j]`。若有效重复不足，可使用共享参数、shrinkage 或保留原值。

上述为可比较的算法示例，不是默认推荐所有数据执行完整校正。两种校正叠加时应避免重复移除同一技术效应。

### Step 4 — 缺失值处理（可选）

校正与插补应分开评价。先保留校正后实测矩阵的 NA，再根据后续用途决定是否预测缺失值。可比较 KNN、低秩矩阵分解、适合低丰度缺失的统计模型等；不要默认用 0、全局最小值或任意常数填充全部缺失位置。

如果执行插补，记录蛋白覆盖率阈值、预测方法、参数及填充位置，禁止将真实观测值静默替换为模型估计值。随机遮蔽观测值可用于测试预测误差，但不能证明对真实非随机缺失（MNAR）的插补同样准确。

### Step 5 — 校正效果评价

| 指标 | 计算/评价原则 |
| --- | --- |
| QC Median Protein CV | 使用线性强度；每个蛋白分别计算 QC CV，再取所有有效蛋白的中位数 |
| CV <10% / <20% 蛋白比例 | 使用校正前后相同的有效蛋白及观测位置 |
| 板内、板间 CV | 分开评估同质 QC 的板内重复性和跨板一致性 |
| PCA / correlation | 检查批次聚集及生物学结构；比较时尽量使用共同蛋白集合和一致的变换 |
| Correction magnitude | 检查每个样本、每个蛋白的校正量、极端值及未校正位置 |
| Independent QC validation | 验证 QC 不参与其校正参数拟合；分别评价同板留出和跨板泛化能力 |
| Biological preservation | 用已知混合比例、真实 fold change、技术重复或独立实验验证生物变化是否保留 |
| Imputation RMSE / MAE | 对独立遮蔽的实测值计算预测误差，不与未插补 CV 混合统计 |

统一 CV 口径建议：`CV (%) = 100 × sample_SD / mean`，标准差采用 `ddof=1`，每个蛋白至少 3 个有效 QC 观测。记录有效样本及蛋白数量，校正前后使用相同的有效观测，不使用插补值计算实测 QC CV。数学上，一个板内所有样本的同一蛋白若乘以相同的正校正因子，该蛋白板内 CV 不会改变。

参与建模的 QC CV 下降可能只是拟合结果，不能证明校正方法提高了独立样本的定量准确性；正式研究样本组内 CV 也包含真实生物差异，不宜直接作为越低越好的优化目标。

## 7. 推荐输出目录及可重复性

```text
project/
├── README.md
├── input/                     # 原始三类输入文件，不覆盖
├── processed/
│   ├── metadata_clean.tsv
│   ├── sample_matching_report.tsv
│   ├── help_matrix.tsv
│   ├── help_quality.tsv
│   ├── protein_annotation.tsv
│   ├── protein_matrix_linear.tsv
│   ├── protein_matrix_log2.tsv
│   └── observed_mask.tsv
├── results/
│   ├── normalization_summary.tsv
│   ├── qc_cv_per_protein.tsv
│   ├── validation_summary.tsv
│   ├── correction_factors.tsv
│   ├── normalized_observed_matrix.tsv
│   ├── normalized_imputed_matrix.tsv  # 可选
│   ├── imputed_mask.tsv              # 可选
│   └── figures/
├── scripts/                   # 分析代码、依赖及版本
└── run_config.json             # 参数、样本选择、随机种子、输入版本
```

上述为建议的数据组织方式，并非固定软件的输出要求。多种方法比较时，建议按方法名建立子目录，并记录参考 QC 类型、筛选阈值、校正公式、模型参数、训练/验证划分、输入文件版本及缺失值处理规则。所有校正应能追溯至原始观测和对应校正量。

## 8. 说明与适用范围

本指南的数据字段依据本项目实际 HELP 报告、metadata 和 DIA-NN 矩阵格式整理。筛选阈值沿用当前数据处理经验作为示例；其他分析方法、输出目录及验证流程为通用建议，需根据实验设计、仪器、内标加入阶段和数据质量调整。没有任何单一 normalization 方法适合所有数据集，应通过独立技术验证及生物学信号保留测试选择方法。
