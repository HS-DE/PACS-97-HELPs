#!/usr/bin/env python3
"""PACS：HELP 蛋白强度校正及自动生成校正前后可视化。"""

import argparse
import csv
import hashlib
from datetime import datetime, timedelta, timezone
from pathlib import Path
import re
import sys
from tempfile import TemporaryDirectory
from textwrap import fill

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy.stats import pearsonr
from sklearn.impute import KNNImputer
from sklearn.decomposition import PCA

EPS = 1e-8
THRESHOLDS = (70, 65, 60, 55, 50)
TARGET_PROTEINS = 2800
TARGET_SAMPLE_PERCENT = 70
HELP_COLUMNS = (
    "Sample", "Peptide", "Peak_Scan_Points", "Peak_Detection_QC_Flag",
    "Accurate_AUC", "N_Common_Transitions",
)
HELP_WARNING = (
    "HELP体系检出不稳定，建议核查实验是否存在异常，建议人工检查是否存在异常检出的HELP"
    "（推荐用箱线图、PCA等方式进行细察），本次矫正结果可信度可能较低"
)
COLORS = {"QC": "#356D96", "HC": "#C28B3B", "S": "#89638F", "T2DM": "#6A8A4E",
          "T002": "#C28B3B", "T003": "#89638F", "T004": "#6A8A4E"}


STAGES = (
    ("raw", "Raw", "校正前", "", "#78858E"),
    ("corrected", "Corrected", "校正后", "", "#356D96"),
)


MARKERS = ("o", "s", "^", "D", "v", "P", "X", "<", ">", "h")


def clean_sample_id(value):
    """跨平台提取 sample_id：兼容 Windows/Linux 路径及全角符号，移除末尾 .raw（不区分大小写）。"""
    basename = re.split(r"[/\\／＼:：]", str(value).strip().strip("\"'"))[-1]
    return re.sub(r"\.raw$", "", basename, flags=re.IGNORECASE).strip()


def validate_header(header, filename):
    if any(value is None or not str(value).strip() for value in header):
        raise ValueError(f"表头存在空列名：{filename}")
    names = [str(value) for value in header]
    if len(names) != len(set(names)):
        raise ValueError(f"表头存在重复列名：{filename}")


def read_table(path):
    """使用第一行表头，忽略无列名的全空列和全空行；Excel 读取第一张表。"""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".xlsx":
        book = load_workbook(path, read_only=True, data_only=True)
        try:
            iterator = book.worksheets[0].iter_rows(
                min_row=1, max_row=1, values_only=True
            )
            try:
                header = list(next(iterator, []))
            finally:
                iterator.close()
            validate_header(
                [name for name in header if name is not None and str(name).strip()],
                path.name,
            )
        finally:
            book.close()
        table = pd.read_excel(path, sheet_name=0, dtype=object, keep_default_na=False)
    elif suffix in (".csv", ".tsv", ".txt"):
        delimiter = "," if suffix == ".csv" else "\t"
        with path.open(encoding="utf-8-sig", newline="") as stream:
            header = next(csv.reader(stream, delimiter=delimiter), [])
        validate_header(
            [name for name in header if name is not None and str(name).strip()],
            path.name,
        )
        table = pd.read_csv(
            path,
            sep=delimiter,
            dtype=object,
            keep_default_na=False,
            encoding="utf-8-sig",
        )
    else:
        raise ValueError(f"蛋白定量表和 metadata 支持 xlsx、csv、tsv、txt：{path.name}")
    kept = []
    for index in range(table.shape[1]):
        name = header[index]
        if name is None or not str(name).strip():
            values = table.iloc[:, index]
            empty = values.isna() | values.astype(str).str.strip().eq("")
            if not empty.all():
                raise ValueError(f"第 {index + 1} 列没有列名但包含数据：{path.name}")
        else:
            kept.append(index)
    table = table.iloc[:, kept].copy()
    table.columns = [str(header[index]) for index in kept]
    empty_rows = table.apply(
        lambda values: values.isna() | values.astype(str).str.strip().eq("")
    ).all(axis=1)
    table = table.loc[~empty_rows].reset_index(drop=True)
    if table.empty or not table.columns.is_unique:
        raise ValueError(f"表格为空或表头重复：{path.name}")
    return table


def validate_batch_info(meta):
    if "batch" not in meta.columns:
        raise ValueError("请在metadata中填写完整Batch信息")
    labels = meta.batch.astype(str).str.strip().str.casefold()
    if (meta.batch.isna() | labels.isin(["", "na", "n/a", "nan", "null", "none", "<na>", "#n/a", "#na"])).any():
        raise ValueError("请在metadata中填写完整Batch信息")


def validate_qc_samples(meta):
    validate_batch_info(meta)
    if not meta.group.eq("qc").any():
        raise ValueError("metadata 中必须包含 pooled QC 样本")


def load_metadata(path, plate_col=None, min_samples=2):
    """列名忽略大小写；保留原始样本记录，condition 用于生物学分组。"""
    source = read_table(path)
    columns = source.columns.str.strip().str.casefold()
    duplicates = columns[columns.duplicated()].unique().tolist()
    if duplicates:
        raise ValueError(f"metadata 列名忽略大小写并去除首尾空白后存在重复：{duplicates}")
    source.columns = columns
    if plate_col is not None:
        plate_col = plate_col.strip().casefold()
    required = ["sample_id", "group", "condition"] + (
        [plate_col] if plate_col is not None else []
    )
    missing = [name for name in required if name not in source.columns]
    if missing:
        raise ValueError(f"metadata 缺少必需列：{missing}")
    sample_ids = source.sample_id.astype(str).str.strip()
    if sample_ids.eq("").any() or not sample_ids.is_unique:
        raise ValueError("metadata 的 sample_id 不得为空或重复")
    types = source.group.astype(str).str.strip().str.casefold()
    types = types.str.replace(r"^internal[-_\s]+calibrator$", "internal-calibrator", regex=True)
    allowed = {"qc", "internal-calibrator", "blank", "sample"}
    invalid = sorted(set(types) - allowed)
    if invalid:
        raise ValueError(f"metadata group 仅允许 qc、internal-calibrator、blank、sample；不支持：{invalid}")
    retained = types.isin(["qc", "sample"])
    excluded_ids = set(sample_ids[~retained])
    # 先剔除不用的样本，只要求保留的 QC/sample 填写完整 batch。
    validate_batch_info(source.loc[retained])
    rows = []
    for index in source.index[retained]:
        row = source.loc[index]
        group = types.loc[index]
        value = row["condition"]
        condition = "" if pd.isna(value) else str(value).strip()
        condition = "QC" if group == "qc" else (condition or "study_sample")
        batch = str(row["batch"]).strip()
        record = {
            "sample_id": sample_ids.loc[index],
            "group": group,
            "condition": condition,
            "batch": batch,
        }
        if plate_col is not None:
            record["plate"] = row[plate_col]
        rows.append(record)
    meta = pd.DataFrame(rows)
    if meta.empty:
        raise ValueError("metadata 中没有可分析的 QC 或 sample 组")
    if min_samples < 2:
        raise ValueError("HELP 校正所需每组最低样本数不得少于 2")
    validate_qc_samples(meta)
    return meta, excluded_ids


def read_help_report(path):
    """读取 Peptide_Long；CSV/TSV 输入须为同结构长表。"""
    path = Path(path)
    if path.suffix.lower() == ".xlsx":
        try:
            source = pd.read_excel(
                path, sheet_name="Peptide_Long", dtype=object, keep_default_na=False
            )
        except ValueError as exc:
            raise ValueError("HELP数据输入有误") from exc
    elif path.suffix.lower() in (".csv", ".tsv", ".txt"):
        source = read_table(path)
    else:
        raise ValueError("HELP数据输入有误")
    if not set(HELP_COLUMNS).issubset(source.columns):
        raise ValueError("HELP数据输入有误")
    source = source[list(HELP_COLUMNS)].copy()
    for column in ("Sample", "Peptide"):
        if (source[column].isna() | source[column].astype(str).str.strip().eq("")).any():
            raise ValueError("HELP数据输入有误")
    source["Sample"] = source.Sample.map(clean_sample_id)
    source["Peptide"] = source.Peptide.astype(str).str.strip()
    if source.Sample.eq("").any():
        raise ValueError("HELP数据输入有误")
    for column in ("Peak_Scan_Points", "Accurate_AUC", "N_Common_Transitions"):
        source[column] = pd.to_numeric(source[column], errors="coerce")
    return source


def extract_help_table(report, min_scan_points=8):
    """HELP表提取：PASS、共同 transition≥2、扫描点≥阈值，然后转宽表。"""
    keep = (
        report.Peak_Detection_QC_Flag.eq("PASS")
        & report.N_Common_Transitions.ge(2)
        & report.Peak_Scan_Points.ge(min_scan_points)
    )
    selected = report.loc[keep, ["Sample", "Peptide", "Accurate_AUC"]]
    if selected.duplicated(["Sample", "Peptide"]).any():
        raise ValueError("HELP数据输入有误：同一样本和肽段存在重复记录")
    heavy = selected.pivot(index="Peptide", columns="Sample", values="Accurate_AUC")
    # 保留筛选前的样本列：未通过质量筛选的定量位置仍为缺失，不补值。
    heavy = heavy.reindex(columns=report.Sample.drop_duplicates())
    # 宽表左上角留空，不再执行旧版的转置操作。
    heavy.index.name = None
    heavy.columns.name = None
    return heavy


def extract_help_report(report, meta):
    """按可用于校正的全样本检出 HELP 数量及 pooled QC CV，最多从8放宽至7。"""
    missing = sorted(set(meta.sample_id) - set(report.Sample))
    if missing:
        raise ValueError(f"HELP 缺少需要保留的 QC/sample 样本：{missing}")
    active = report.loc[report.Sample.isin(meta.sample_id)]
    trials = []
    for threshold in (8, 7):
        heavy = extract_help_table(active, threshold).reindex(columns=meta.sample_id)
        heavy.index.name = heavy.columns.name = None
        logged = positive_log(heavy.to_numpy(float))
        complete = np.isfinite(logged).all(axis=1)
        qc = meta.group.eq("qc").to_numpy()
        cv, _ = linear_cv(logged[complete][:, qc], min_observations=2)
        finite_cv = cv[np.isfinite(cv)]
        acceptable = int(complete.sum()) >= 5 and bool((finite_cv <= 15).any())
        trials.append({
            "min_scan_points": threshold,
            "n_extracted_HELPs": len(heavy),
            "n_complete_HELPs": int(complete.sum()),
            "n_finite_QC_CV": int(finite_cv.size),
            "min_QC_CV_percent": float(finite_cv.min()) if finite_cv.size else np.nan,
            "acceptable": acceptable,
        })
        if acceptable:
            break
    unstable = not trials[-1]["acceptable"]
    if unstable:
        print(HELP_WARNING, file=sys.stderr)
    if heavy.empty:
        raise ValueError("HELP表提取后没有HELP，无法进行矫正")
    heavy.attrs["HELP_extraction"] = {"min_scan_points": threshold, "unstable": unstable, "trials": trials}
    return heavy


def load_inputs(
    protein_path,
    help_path,
    metadata_path,
    plate_col=None,
    min_samples=2,
):
    meta, excluded_ids = load_metadata(metadata_path, plate_col, min_samples)
    report = read_help_report(help_path)
    linear_heavy = extract_help_report(report, meta)
    protein = read_table(protein_path)
    if "Protein.Group" not in protein.columns:
        raise ValueError("蛋白定量表缺少 Protein.Group 列")
    protein_ids = protein["Protein.Group"].astype(str).str.strip()
    if not protein_ids.ne("").all() or not protein_ids.is_unique:
        raise ValueError("Protein.Group 必须非空且每行唯一")
    metadata_ids = set(meta.sample_id)
    sample_ids = metadata_ids | excluded_ids | set(report.Sample)
    sample_columns, rename = [], {}
    discarded = []
    for name in protein.columns:
        if name == "Protein.Group":
            continue
        sample_id = clean_sample_id(name)
        if sample_id in metadata_ids:
            sample_columns.append(sample_id)
            rename[name] = sample_id
        elif sample_id in sample_ids or re.search(r"[/\\／＼]|\.raw", name, flags=re.IGNORECASE):
            discarded.append(name)
    if len(sample_columns) != len(set(sample_columns)):
        raise ValueError("蛋白表的样本列在清理路径和 .raw 后重复")
    missing = sorted(metadata_ids - set(sample_columns))
    if missing:
        raise ValueError(f"以下 metadata 样本在蛋白表中没有匹配列：{missing}")
    if not sample_columns:
        raise ValueError("蛋白表中没有匹配 metadata 的样本列")
    protein = protein.drop(columns=discarded).rename(columns=rename)
    if not protein.columns.is_unique:
        raise ValueError("清理后的样本列名与其他列名重复")
    meta = meta.set_index("sample_id").loc[sample_columns].reset_index()
    required = metadata_ids
    missing_heavy = sorted(required - set(linear_heavy.columns))
    if missing_heavy:
        raise ValueError(f"HELP 缺少需要保留的 QC/sample 样本：{missing_heavy}")
    raw = positive_log(
        protein[sample_columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    )
    heavy = positive_log(linear_heavy.reindex(columns=sample_columns).to_numpy(float))
    invalid = ~np.isfinite(raw).any(axis=0)
    if invalid.any():
        raise ValueError("存在未检出任何有效正强度蛋白的分析样本")
    return protein, meta, raw, heavy, linear_heavy.index.tolist(), sample_columns


def positive_log(values):
    values = np.asarray(values, float)
    result = np.full(values.shape, np.nan)
    valid = np.isfinite(values) & (values > 0)
    result[valid] = np.log2(values[valid])
    return result


def factors(p, h, x):
    p, h = np.asarray(p, float), np.asarray(h, float)
    if (
        p.shape != h.shape
        or not np.isfinite(p).all()
        or not np.isfinite(h).all()
        or not np.isfinite(x)
        or not 1 <= x <= 10
    ):
        raise ValueError("Finite offsets and group x in [1,10] are required")
    zero = (abs(p) <= EPS) | (abs(h) <= EPS)
    same = (np.sign(p) == np.sign(h)) & ~zero
    small = same & (abs(h) <= abs(p))
    large = same & ~small
    opposite = ~zero & ~same
    f = np.zeros(len(p))
    branch = np.full(len(p), "zero", dtype=object)
    f[small] = abs(h[small])
    branch[small & (p > 0)] = "a"
    branch[small & (p < 0)] = "b"
    f[large] = abs(p[large]) / x
    branch[large & (p > 0)] = "c"
    branch[large & (p < 0)] = "d"
    f[opposite] = np.minimum(abs(h[opposite]), abs(p[opposite])) / x
    branch[opposite] = "e"
    shift = -np.sign(p) * f
    assert ((f >= 0) & (f <= abs(p) + 1e-12)).all()
    assert (abs(p + shift) <= abs(p) + 1e-12).all()
    return f, shift, branch


def choose_help(metrics):
    """依次检查 RMSE 前五条；HELP_j 为全部候选中 pooled QC CV 最小者。"""
    if not metrics:
        return None, "no_eligible_HELP", None
    by_rmse = sorted(metrics, key=lambda m: (m["RMSE"], m["HELP"]))
    by_cv = sorted(
        (m for m in metrics if np.isfinite(m["QC_CV_percent"])),
        key=lambda m: (m["QC_CV_percent"], m["HELP"]),
    )
    if not by_cv:
        # 例如仅有一个 QC：CV 无法估计。提取阶段已警告，仍按 RMSE 选择。
        return by_rmse[0], "QC_CV_unavailable", None
    reference = by_cv[0]
    cv_top5 = {m["HELP"] for m in by_cv[:5]}
    for candidate in by_rmse[:5]:
        cv = candidate["QC_CV_percent"]
        if candidate["HELP"] in cv_top5:
            return candidate, "QC_CV_top5", reference
        if np.isfinite(cv) and cv <= 10:
            return candidate, "QC_CV_le10", reference
        if np.isfinite(cv) and cv - reference["QC_CV_percent"] < 2:
            return candidate, "QC_CV_gap_lt2", reference
    return reference, "fallback_min_QC_CV", reference


def correct_group(raw, heavy, sequences, qc_cv=None):
    if not np.isfinite(raw).any(axis=0).all():
        raise ValueError("Sample has no detected PG")
    med = np.nanmedian(raw, axis=0)
    ref = float(np.median(med))
    p = med - ref
    candidates = []
    for j, (seq, h) in enumerate(zip(sequences, heavy)):
        if np.isfinite(h).all():
            delta = h - np.median(h)
            rmse = float(np.sqrt(np.mean((delta - p) ** 2)))
            pearson = (
                float(pearsonr(h, med).statistic)
                if np.std(h) >= 1e-12 and np.std(med) >= 1e-12 else np.nan
            )
            candidates.append({
                "HELP": seq, "index": j, "RMSE": rmse, "Pearson": pearson,
                "QC_CV_percent": float(qc_cv.get(seq, np.nan)) if qc_cv is not None else np.nan,
            })
    info = {
        "PG_ref": ref,
        "P_delta": p,
        "PG_median_raw": med,
        "H_delta": np.full(len(med), np.nan),
        "HELP": "",
        "Pearson": np.nan,
        "RMSE": np.nan,
        "x": np.nan,
        "status": "no_eligible_HELP",
        "HELP_metrics": candidates,
        "QC_CV_percent": np.nan,
        "HELP_j": "",
        "HELP_j_QC_CV_percent": np.nan,
    }
    shift = np.zeros(len(med))
    chosen, reason, reference = choose_help(candidates)
    info["selection_reason"] = reason
    if reference is not None:
        info.update(HELP_j=reference["HELP"], HELP_j_QC_CV_percent=reference["QC_CV_percent"])
    if chosen is not None:
        seq, j = chosen["HELP"], chosen["index"]
        h = heavy[j] - np.median(heavy[j])
        info.update(
            HELP=seq,
            H_delta=h,
            Pearson=chosen["Pearson"],
            RMSE=chosen["RMSE"],
            QC_CV_percent=chosen["QC_CV_percent"],
            status="no_valid_ratio_for_x",
        )
        valid = (abs(p) > EPS) & (abs(h) > EPS) & np.isfinite(p) & np.isfinite(h)
        if valid.any():
            ratio = np.maximum(abs(p[valid]), abs(h[valid])) / np.minimum(
                abs(p[valid]), abs(h[valid])
            )
            if not np.isfinite(ratio).all():
                raise ValueError("Offset ratio overflow")
            x = min(10.0, float(np.median(ratio)))
            _, shift, _ = factors(p, h, x)
            info.update(x=x, status="ok")
    after = raw + shift[None, :]
    np.testing.assert_array_equal(np.isfinite(raw), np.isfinite(after))
    info.update(
        PG_median_corrected=np.nanmedian(after, axis=0), corrected_offset=p + shift
    )
    return after, info


def linear_cv(logged, min_observations=3):
    if min_observations < 2:
        raise ValueError("Sample CV needs at least 2 observations")
    observed = np.isfinite(logged)
    x = np.exp2(logged)
    assert np.isfinite(x[observed]).all() and (x[observed] > 0).all()
    n = observed.sum(axis=1)
    mean = np.divide(
        np.where(observed, x, 0).sum(axis=1), n, out=np.zeros(len(n)), where=n > 0
    )
    ss = np.where(observed, x - mean[:, None], 0) ** 2
    var = np.divide(ss.sum(axis=1), n - 1, out=np.zeros(len(n)), where=n > 1)
    cv = np.divide(
        100 * np.sqrt(var), mean, out=np.full(len(n), np.nan), where=mean > 0
    )
    cv[n < min_observations] = np.nan
    return cv, n


def qc_cv_top5(heavy, sequences, meta):
    """合并全部 pooled QC，返回所有候选的线性 CV 及 CV 前五条。"""
    assert heavy.shape == (len(sequences), len(meta))
    qc = meta.group.astype(str).str.casefold().eq("qc").to_numpy()
    if not qc.any():
        return None
    values = heavy[:, qc]
    cv, _ = linear_cv(values, min_observations=2)
    complete = np.isfinite(values).all(axis=1)
    ranked = sorted(
        (float(cv[j]), str(seq))
        for j, seq in enumerate(sequences)
        if complete[j] and np.isfinite(cv[j])
    )
    return {
        "n_QC": int(qc.sum()),
        "HELPs": [seq for _, seq in ranked[:5]],
        "CV_percent": {seq: value for value, seq in ranked},
    }


def detected(values):
    """检出：有限、非零数值；空值、NA、非数值在入口转换为 NaN。"""
    values = np.asarray(values, dtype=float)
    return np.isfinite(values) & (values != 0)


def target_met(counts, target_proteins=TARGET_PROTEINS):
    counts = np.asarray(counts)
    if counts.ndim != 1 or not counts.size:
        raise ValueError("鉴定量统计必须包含至少一个样本")
    return bool(
        np.count_nonzero(counts > target_proteins) * 100
        >= TARGET_SAMPLE_PERCENT * counts.size
    )


def partial_impute(values, target_proteins=TARGET_PROTEINS):
    """70→65→60→55→50%，严格大于阈值；原始观测和前轮插补值均不改写。"""
    original = np.asarray(values, dtype=float)
    if original.ndim != 2 or not all(original.shape):
        raise ValueError("部分插补需要非空的蛋白×样本矩阵")
    if target_proteins < 0:
        raise ValueError("蛋白鉴定量目标不得为负数")
    result = original.copy()
    baseline = detected(original)
    before_counts = baseline.sum(axis=0)
    report = {
        "before_counts": before_counts,
        "after_counts": before_counts.copy(),
        "stages": [],
        "threshold": None,
        "target_reached": target_met(before_counts, target_proteins),
        "stop_reason": "already_satisfied",
    }
    if report["target_reached"]:
        return result, report
    # 检出率基于插补前真实观测固定计算；分母为当前类别的样本数。
    detected_per_protein = baseline.sum(axis=1)
    for threshold in THRESHOLDS:
        eligible = detected_per_protein * 100 > threshold * original.shape[1]
        subset = result[eligible].copy()
        missing = ~detected(subset)
        added = 0
        if missing.any():
            log_values = np.full(subset.shape, np.nan)
            positive = np.isfinite(subset) & (subset > 0)
            log_values[positive] = np.log2(subset[positive])
            available_columns = np.isfinite(log_values).any(axis=0)
            # 行为蛋白、列为样本，直接由 sklearn 完成距离、邻居选择和填补。
            estimated = KNNImputer(
                n_neighbors=10,
                weights="uniform",
                metric="nan_euclidean",
                keep_empty_features=True,
            ).fit_transform(log_values)
            # sklearn 保留全缺失特征时会填 0；此处仍保留缺失，避免还原成强度 1。
            estimated[:, ~available_columns] = np.nan
            fill = missing & np.isfinite(estimated)
            with np.errstate(over="ignore", under="ignore"):
                linear = np.exp2(estimated[fill])
            if not np.isfinite(linear).all() or not (linear > 0).all():
                raise ValueError("KNN 插补值转换为线性强度时溢出或下溢")
            subset[fill] = linear
            result[eligible] = subset
            added = int(fill.sum())
        counts = detected(result).sum(axis=0)
        report["stages"].append(
            {
                "threshold_percent": threshold,
                "eligible_proteins": int(eligible.sum()),
                "new_values": added,
                "samples_above_target": int((counts > target_proteins).sum()),
            }
        )
        report.update(
            after_counts=counts,
            threshold=threshold,
            target_reached=target_met(counts, target_proteins),
        )
        if report["target_reached"]:
            report["stop_reason"] = "target_reached"
            break
    else:
        report["stop_reason"] = "threshold_floor_reached"
    np.testing.assert_array_equal(result[baseline], original[baseline])
    assert (report["after_counts"] >= before_counts).all()
    return result, report


def partial_impute_by_type(values, meta, target_proteins=TARGET_PROTEINS):
    """QC 和 sample 分别判断、筛选和插补，再按原样本列序合并。"""
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or not all(values.shape) or values.shape[1] != len(meta):
        raise ValueError("部分插补矩阵与 metadata 的样本数量不一致")
    types = meta.group.astype(str).str.casefold()
    result = values.copy()
    before_counts = detected(values).sum(axis=0)
    after_counts = before_counts.copy()
    reports = {}
    for name, mask in (
        ("QC", types.eq("qc").to_numpy()),
        ("sample", types.eq("sample").to_numpy()),
    ):
        if not mask.any():
            continue
        result[:, mask], reports[name] = partial_impute(values[:, mask], target_proteins)
        after_counts[mask] = reports[name]["after_counts"]
    return result, {
        "before_counts": before_counts,
        "after_counts": after_counts,
        "by_type": reports,
    }


def qc_variability(raw, meta):
    """校正前合并 pooled QC 计算蛋白线性 CV；原有质量警报继续保留。"""
    qc = meta.group.eq("qc").to_numpy()
    cv, _ = linear_cv(raw[:, qc], min_observations=2)
    valid = cv[np.isfinite(cv)]
    median = float(np.median(valid)) if valid.size else np.nan
    if median > 15:
        print(
            "QC 样本差异过大，建议检测实验过程是否出现问题，本次矫正结果可信度可能较低。",
            file=sys.stderr,
        )
    return {
        "n_QC": int(qc.sum()),
        "n_proteins_with_CV": int(valid.size),
        "median_CV_percent": median,
    }


def correct_data(meta, raw, heavy, sequences, min_samples=2):
    """统一视为一个分析批次；QC 合并筛选 HELP，再按 group/condition 校正。"""
    if min_samples < 2:
        raise ValueError("HELP 校正所需每组最低样本数不得少于 2")
    validate_qc_samples(meta)
    corrected = raw.copy()
    # batch 标签仅作输入完整性检查，不参与筛选、分组或偏移计算。
    complete = np.isfinite(heavy).all(axis=1)
    complete_sequences = [seq for seq, keep in zip(sequences, complete) if keep]
    screen = qc_cv_top5(heavy[complete], complete_sequences, meta)
    models = {
        "groups": {},
        "QC_CV_top5": screen,
        "globally_complete_HELPs": complete_sequences,
        "QC_variability": qc_variability(raw, meta),
    }
    for key, members in meta.groupby(["group", "condition"], sort=False):
        indices = members.sort_values("sample_id").index.to_numpy()
        if len(indices) < min_samples:
            models["groups"][key] = {
                "status": "too_few_samples_for_HELP", "n_samples": len(indices), "HELP": "",
            }
            continue
        corrected[:, indices], models["groups"][key] = correct_group(
            raw[:, indices], heavy[complete][:, indices], complete_sequences, screen["CV_percent"]
        )
    np.testing.assert_array_equal(np.isfinite(corrected), np.isfinite(raw))
    return corrected, models


def corrected_quantification(protein, sample_columns, raw, corrected):
    after = protein.copy()
    for j, name in enumerate(sample_columns):
        numeric = pd.to_numeric(protein[name], errors="coerce").to_numpy(float)
        finite = np.isfinite(numeric)
        after.loc[finite, name] = numeric[finite]
        changed = np.isfinite(raw[:, j]) & (corrected[:, j] != raw[:, j])
        if changed.any():
            shifts = corrected[changed, j] - raw[changed, j]
            np.testing.assert_allclose(shifts, shifts[0], atol=1e-12, rtol=0)
            values = numeric[changed] * np.exp2(shifts[0])
            if not np.isfinite(values).all() or not (values > 0).all():
                raise ValueError(f"校正后的线性强度溢出或下溢：{name}")
            after.loc[changed, name] = values
    return after


def plot_counts(meta, report, path):
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
    after = np.asarray(report["after_counts"])
    if len(meta) != len(after) or not np.isfinite(after).all() or (after < 0).any():
        raise ValueError("鉴定量柱状图的样本数或鉴定量数值不正确")
    qc = meta.group.eq("qc").to_numpy()
    order = np.r_[np.flatnonzero(qc), np.flatnonzero(~qc)]
    meta = meta.iloc[order].reset_index(drop=True)
    after = after[order]
    groups = meta.condition.drop_duplicates().tolist()
    colors = {
        group: COLORS.get(group, matplotlib.colormaps["tab20"](i % 20))
        for i, group in enumerate(groups)
    }
    fills = [colors[group] for group in meta.condition]
    positions = np.arange(len(meta))
    fig, ax = plt.subplots(figsize=(max(12, 0.28 * len(meta) + 4), 9))
    ax.bar(positions, after, width=0.75, color=fills)
    upper = max(1, int(after.max())) * 1.13
    for x, count in zip(positions, after):
        ax.text(
            x,
            count + upper * 0.009,
            f"{int(count)}",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax.set_xticks(positions, meta.sample_id, rotation=90, fontsize=8)
    ax.set(
        xlim=(-0.7, len(meta) - 0.3),
        ylim=(0, upper),
        xlabel="Sample",
        ylabel="Quantified proteins (N)",
    )
    ax.grid(axis="y", alpha=0.18)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    handles = [Patch(facecolor=colors[group], label=group) for group in groups]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.99),
        ncol=min(6, len(handles)),
        frameon=False,
    )
    bottom = max(
        0.13, min(0.45, 0.07 + 0.007 * max(len(str(value)) for value in meta.sample_id))
    )
    fig.subplots_adjust(
        top=0.93 - 0.03 * ((len(groups) - 1) // 6), bottom=bottom, left=0.06, right=0.99
    )
    fig.savefig(path, dpi=180, bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)


def normalize_data(
    protein_path,
    help_path,
    metadata_path,
    output_dir,
    plate_col=None,
    min_samples=2,
):
    output_dir = Path(output_dir).expanduser().resolve()
    protein, meta, raw, heavy, sequences, sample_columns = load_inputs(
        protein_path,
        help_path,
        metadata_path,
        plate_col,
        min_samples,
    )
    corrected, _ = correct_data(meta, raw, heavy, sequences, min_samples)
    final = corrected_quantification(protein, sample_columns, raw, corrected)
    numeric = (
        final[sample_columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    )
    # 当前版本固定跳过插补，仅统计校正后已有的有效观测。
    report = {"after_counts": detected(numeric).sum(axis=0)}
    output_dir.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".PACS_", dir=output_dir) as temporary:
        staging = Path(temporary)
        # 仅在导出时限制样本强度的小数位，保留注释和原始缺失标记。
        for name in sample_columns:
            values = pd.to_numeric(final[name], errors="coerce")
            finite = np.isfinite(values)
            final.loc[finite, name] = values.loc[finite].map("{:.2f}".format)
        final.to_csv(
            staging / "最终蛋白定量表.tsv",
            sep="\t",
            index=False,
            encoding="utf-8",
            lineterminator="\n",
        )
        plot_counts(meta, report, staging / "蛋白鉴定量统计柱状图.png")
        if len(list(staging.iterdir())) != 2:
            raise RuntimeError("输出文件数量不正确，应为最终 TSV 表和最后一张柱状图")
        for name in ("最终蛋白定量表.tsv", "蛋白鉴定量统计柱状图.png"):
            (staging / name).replace(output_dir / name)
    return output_dir


def top_legends(fig, groups, max_columns=3):
    """按实际宽度换行，并增加图高，为图例、子图标题分别留出空间。"""
    groups = [(handles, title) for handles, title in groups if handles]
    if not groups:
        return
    width, height = fig.get_size_inches()
    renderer = fig.canvas.get_renderer()
    legends = []
    available_width = (width / len(groups) - .3) * fig.dpi
    for i, (handles, title) in enumerate(groups):
        for columns in range(min(max_columns, len(handles)), 0, -1):
            legend = fig.legend(
                handles=handles, labels=[fill(str(handle.get_label()), width=28) for handle in handles],
                title=title, loc="upper center", bbox_to_anchor=((i + .5) / len(groups), 1.),
                ncol=columns, frameon=False, borderaxespad=0, columnspacing=1.2, handletextpad=.5,
            )
            if legend.get_window_extent(renderer).width <= available_width or columns == 1:
                break
            legend.remove()
        legends.append(legend)
    legend_height = max(legend.get_window_extent(renderer).height for legend in legends) / fig.dpi
    new_height = height + legend_height + .45
    bottom = fig.subplotpars.bottom * height / new_height
    fig.set_size_inches(width, new_height)
    for i, legend in enumerate(legends):
        legend.set_bbox_to_anchor(((i + .5) / len(groups), 1. - .12 / new_height),
                                  transform=fig.transFigure)
    # 保持底部标签空间和绘图区高度，标题与图例之间额外留白。
    fig.subplots_adjust(top=(height - .22) / new_height, bottom=bottom)


def save(fig, folder, filename):
    fig.savefig(folder / filename, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def style(axis):
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", color="#D9DDE1", linewidth=.6, alpha=.7)
    axis.set_axisbelow(True)


def read_matrix(path, meta, proteins=None):
    table = read_table(path)
    if "Protein.Group" not in table.columns:
        raise ValueError(f"蛋白定量表缺少 Protein.Group 列：{path.name}")
    protein_ids = table["Protein.Group"].astype(str).str.strip()
    if not protein_ids.ne("").all() or not protein_ids.is_unique:
        raise ValueError(f"Protein.Group 必须非空且每行唯一：{path.name}")
    names = set(meta.sample_id)
    table = table.rename(columns={
        c: clean_sample_id(c) for c in table.columns if clean_sample_id(c) in names
    })
    if not table.columns.is_unique or not table["Protein.Group"].is_unique:
        raise ValueError(f"蛋白或样本 ID 重复：{path}")
    missing = names - set(table.columns)
    if missing:
        raise ValueError(f"样本未匹配：{path.name}: {sorted(missing)}")
    table = table.set_index("Protein.Group")
    if proteins is not None:
        if set(table.index) != set(proteins):
            raise ValueError(f"校正前后蛋白集合不一致：{path}")
        table = table.loc[proteins]
    values = table[meta.sample_id].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    return table.index.to_numpy(), values


def boxplots(logged, meta, scopes, folder):
    rows = []
    for number, (scope, indices) in enumerate(scopes.items(), 1):
        selected = meta.iloc[indices]
        group_colors = [COLORS.get(c, "#78858E") for c in selected.condition]
        fig, axes = plt.subplots(len(STAGES), 1, figsize=(max(12, .43 * len(indices) + 3), 4 * len(STAGES)),
                                 sharex=True, sharey=True)
        for axis, (stage, label, _, _, _) in zip(axes, STAGES):
            vectors = [logged[stage][np.isfinite(logged[stage][:, i]), i] for i in indices]
            if not all(len(v) for v in vectors):
                raise ValueError("箱线图存在无有效蛋白的样本")
            boxes = axis.boxplot(
                vectors, positions=np.arange(len(indices)), widths=.65, patch_artist=True,
                medianprops={"color": "#252525", "linewidth": 1.1},
                flierprops={"marker": ".", "markersize": 1.3, "alpha": .15,
                            "markerfacecolor": "#65717D", "markeredgecolor": "#65717D",
                            "markeredgewidth": 0},
                whiskerprops={"color": "#68717A", "linewidth": .7},
                capprops={"color": "#68717A", "linewidth": .7},
            )
            for patch, color in zip(boxes["boxes"], group_colors):
                patch.set(facecolor=color, edgecolor="#46505A", linewidth=.7, alpha=.8)
            axis.set(title=label, ylabel="Protein intensity (log2)")
            style(axis)
            blocks = (selected.group + "/" + selected.batch).to_numpy()
            boundaries = np.r_[0, np.flatnonzero(blocks[1:] != blocks[:-1]) + 1, len(indices)]
            for left, right in zip(boundaries[:-1], boundaries[1:]):
                if left:
                    axis.axvline(left - .5, color="#AEB6BD", linestyle="--", linewidth=.7)
                axis.text((left + right - 1) / 2, .985, blocks[left].replace("qc", "QC"),
                          transform=axis.get_xaxis_transform(), ha="center", va="top", fontsize=9)
            for i, vector in zip(indices, vectors):
                row = meta.iloc[i].to_dict()
                q = np.percentile(vector, [0, 25, 50, 75, 100])
                rows.append({"stage": stage, **row, "n_proteins": len(vector),
                             **dict(zip(("min_log2", "Q1_log2", "median_log2", "Q3_log2", "max_log2"), q))})
        axes[-1].set_xticks(np.arange(len(indices)), selected.sample_id,
                           rotation=65, ha="right", fontsize=9)
        axes[-1].set_xlabel("Sample")
        axes[-1].set_xlim(-.7, len(indices) - .3)
        fig.subplots_adjust(bottom=.18, hspace=.27, left=.07, right=.99)
        top_legends(fig, [([Patch(facecolor=COLORS.get(c, "#78858E"), label=c)
                           for c in selected.condition.drop_duplicates()], None)], max_columns=4)
        save(fig, folder, f"{number:02d}_{scope}_蛋白强度箱线图.png")
    # 同一样本在全部/QC/sample 视图中重复出现，仅保存一份统计。
    return pd.DataFrame(rows).drop_duplicates(["stage", "sample_id"])


def pca_plots(logged, meta, proteins, scopes, folder, shape_col="batch", shape_title="Batch"):
    rows, features, usage = [], [], []
    batches = meta[shape_col].drop_duplicates().tolist()
    markers = {b: MARKERS[i % len(MARKERS)] for i, b in enumerate(batches)}
    for number, (scope, indices) in enumerate(scopes.items(), 4):
        complete = np.logical_and.reduce([
            np.isfinite(logged[s][:, indices]).all(axis=1) for s, *_ in STAGES
        ])
        informative = np.std(logged["raw"][complete][:, indices], axis=1) > 1e-12
        chosen = np.flatnonzero(complete)[informative]
        status = {"scope": scope, "n_samples": len(indices), "n_complete_proteins": int(complete.sum()),
                  "n_used_proteins": len(chosen), "status": "ok"}
        usage.append(status)
        if len(indices) < 3 or len(chosen) < 2:
            status["status"] = "insufficient_samples_or_complete_varying_proteins"
            fig, axes = plt.subplots(1, len(STAGES), figsize=(6 * len(STAGES), 6.5))
            for axis, (_, label, _, _, _) in zip(axes, STAGES):
                axis.set_title(label)
                axis.text(.5, .5, "PCA unavailable\nNeed at least 3 samples and\n2 varying complete proteins",
                          ha="center", va="center", transform=axis.transAxes)
                axis.set_axis_off()
            save(fig, folder, f"{number:02d}_{scope}_PCA.png")
            continue
        # 两个结果共用 Raw 拟合的均值和主轴，不补值、不做逐蛋白方差标准化。
        model = PCA(n_components=2, svd_solver="full").fit(logged["raw"][chosen][:, indices].T)
        coords = {s: model.transform(logged[s][chosen][:, indices].T) for s, *_ in STAGES}
        all_coords = np.concatenate(list(coords.values()))
        low, high = all_coords.min(axis=0), all_coords.max(axis=0)
        margin = np.maximum(high - low, 1.) * .13
        selected = meta.iloc[indices].reset_index(drop=True)
        fig, axes = plt.subplots(1, len(STAGES), figsize=(6 * len(STAGES), 6.5), sharex=True, sharey=True)
        for axis, (stage, label, _, _, _) in zip(axes, STAGES):
            for condition in selected.condition.drop_duplicates():
                for batch in batches:
                    mask = (selected.condition.eq(condition) & selected[shape_col].eq(batch)).to_numpy()
                    axis.scatter(coords[stage][mask, 0], coords[stage][mask, 1], s=66,
                                 marker=markers[batch], color=COLORS.get(condition, "#78858E"),
                                 edgecolors="white", linewidths=.55, alpha=.9)
            axis.set(title=label, xlim=(low[0] - margin[0], high[0] + margin[0]),
                     ylim=(low[1] - margin[1], high[1] + margin[1]),
                     xlabel=f"PC1 ({100 * model.explained_variance_ratio_[0]:.1f}% of Raw variance)")
            style(axis)
            axis.set_aspect("equal", adjustable="box")
            for j, row in selected.iterrows():
                rows.append({"scope": scope, "stage": stage, **row.to_dict(),
                             "PC1": coords[stage][j, 0], "PC2": coords[stage][j, 1],
                             "n_complete_proteins": int(complete.sum()), "n_used_proteins": len(chosen),
                             "raw_PC1_variance_fraction": model.explained_variance_ratio_[0],
                             "raw_PC2_variance_fraction": model.explained_variance_ratio_[1]})
        axes[0].set_ylabel(f"PC2 ({100 * model.explained_variance_ratio_[1]:.1f}% of Raw variance)")
        condition_handles = [Line2D([], [], marker="o", linestyle="", color=COLORS.get(c, "#78858E"), label=c)
                             for c in selected.condition.drop_duplicates()]
        batch_handles = [Line2D([], [], marker=markers[b], linestyle="", color="#505961", label=b)
                         for b in batches if selected[shape_col].eq(b).any()]
        fig.subplots_adjust(bottom=.16, left=.07, right=.99, wspace=.16)
        top_legends(fig, [(condition_handles, "Condition (color)"),
                          (batch_handles, f"{shape_title} (shape)")])
        save(fig, folder, f"{number:02d}_{scope}_PCA.png")
        features.extend({"scope": scope, "Protein.Group": proteins[i]} for i in chosen)
    coordinate_columns = ["scope", "stage", *meta.columns, "PC1", "PC2", "n_complete_proteins", "n_used_proteins",
                          "raw_PC1_variance_fraction", "raw_PC2_variance_fraction"]
    return (pd.DataFrame(rows, columns=coordinate_columns),
            pd.DataFrame(features, columns=["scope", "Protein.Group"]), pd.DataFrame(usage))


def cv_statistics(linear, logged, meta, proteins):
    qc = meta.group.eq("qc").to_numpy()
    groups = {"QC": ("All QC", qc)}
    if (~qc).any():
        groups["sample"] = ("All study samples", ~qc)
    groups.update({f"QC_{b}": (f"QC in {b}", qc & meta.batch.eq(b).to_numpy())
                   for b in meta.batch.unique() if (qc & meta.batch.eq(b).to_numpy()).any()})
    groups.update({f"sample_{c}": (f"Study {c}", (~qc) & meta.condition.eq(c).to_numpy())
                   for c in meta.loc[~qc, "condition"].unique()})
    matched = np.logical_and.reduce([np.isfinite(logged[s]) for s, *_ in STAGES])
    rows, summaries = [], []
    for group, (label, mask) in groups.items():
        for stage, _, chinese, _, _ in STAGES:
            values = np.where(matched[:, mask], logged[stage][:, mask], np.nan)
            cv, n = linear_cv(values, min_observations=3)
            valid = np.isfinite(cv)
            # 独立地从输入线性矩阵直接复核每一个蛋白的 CV。
            direct = pd.DataFrame(np.where(matched[:, mask], linear[stage][:, mask], np.nan))
            independent = (100 * direct.std(axis=1, ddof=1) / direct.mean(axis=1)).to_numpy()
            np.testing.assert_allclose(cv[valid], independent[valid], rtol=1e-10, atol=1e-10)
            distribution = cv[valid]
            for i in np.flatnonzero(valid):
                rows.append({"scope": group, "stage": stage, "method": chinese,
                             "Protein.Group": proteins[i], "n_matched_observations": int(n[i]),
                             "CV_percent": cv[i]})
            summaries.append({"scope": group, "label": label, "stage": stage, "method": chinese,
                              "n_samples": int(mask.sum()), "n_proteins": int(valid.sum()),
                              "median_CV_percent": float(np.median(distribution)) if distribution.size else np.nan,
                              "Q1_CV_percent": float(np.percentile(distribution, 25)) if distribution.size else np.nan,
                              "Q3_CV_percent": float(np.percentile(distribution, 75)) if distribution.size else np.nan,
                              "max_CV_percent": float(distribution.max()) if distribution.size else np.nan})
    return (pd.DataFrame(rows, columns=["scope", "stage", "method", "Protein.Group", "n_matched_observations", "CV_percent"]),
            pd.DataFrame(summaries))


def cv_plots(cv, summary, meta, folder):
    available = set(summary.scope)
    qc_groups = ["QC"] + [f"QC_{b}" for b in meta.batch.unique() if f"QC_{b}" in available]
    conditions = [f"sample_{c}" for c in meta.loc[meta.group.eq("sample"), "condition"].unique()]
    views = {"全部样本": ["QC"] + conditions, "QC": qc_groups}
    if "sample" in available:
        views["研究样本"] = ["sample"] + conditions
    for number, (scope, groups) in enumerate(views.items(), 7):
        fig, axes = plt.subplots(1, len(groups), figsize=(max(12, 4.3 * len(groups)), 6.8),
                                 sharey=True, squeeze=False)
        upper = float(cv.loc[cv.scope.isin(groups), "CV_percent"].max()) * 1.13
        if not np.isfinite(upper) or upper <= 0:
            upper = 1.
        for axis, group in zip(axes[0], groups):
            if not (cv.scope.eq(group)).any():
                axis.text(.5, .5, "CV unavailable\nNeed at least 3 paired\nobservations per protein",
                          ha="center", va="center", transform=axis.transAxes)
            for position, (stage, _, _, _, color) in enumerate(STAGES, 1):
                values = cv.loc[cv.scope.eq(group) & cv.stage.eq(stage), "CV_percent"].to_numpy()
                if not len(values):
                    continue
                if len(values) > 1 and np.ptp(values) > 0:
                    violin = axis.violinplot([values], positions=[position], widths=.75,
                                             showextrema=False, points=300, bw_method=.2)
                    for body in violin["bodies"]:
                        body.set(facecolor=color, edgecolor=color, alpha=.65, linewidth=.8)
                axis.boxplot([values], positions=[position], widths=.17, showfliers=False,
                             patch_artist=True, boxprops={"facecolor": "white", "linewidth": .8},
                             medianprops={"color": "#202020", "linewidth": 1.2},
                             whiskerprops={"linewidth": .8}, capprops={"linewidth": .8})
                axis.text(position, .94, f"{np.median(values):.2f}%", ha="center", va="top",
                          transform=axis.get_xaxis_transform(), color=color, fontweight="bold", fontsize=11)
            row = summary.loc[summary.scope.eq(group)].iloc[0]
            axis.set(title=f"{row.label}\n{row.n_samples} samples, {row.n_proteins:,} proteins",
                     xlim=(.45, len(STAGES) + .55), ylim=(0, upper))
            axis.set_xticks(np.arange(1, len(STAGES) + 1), [s[0].capitalize() for s in STAGES])
            style(axis)
        axes[0, 0].set_ylabel("Protein CV (%)")
        fig.text(.5, .02, "Numbers: median CV. Linear intensities; matched observed entries; at least 3 measurements per protein.",
                 ha="center", fontsize=10, color="#505961")
        fig.subplots_adjust(left=.06, right=.99, bottom=.13, top=.88, wspace=.18)
        save(fig, folder, f"{number:02d}_{scope}_CV提琴图.png")


def qc_cv_histogram(cv, folder):
    """使用本轮计算的 pooled QC 蛋白 CV 绘图，配套 TSV 使用 GBK。"""
    folder = Path(folder)
    data_dir = folder / "绘图数据"
    data_dir.mkdir(parents=True, exist_ok=True)
    qc = cv.loc[cv.scope.eq("QC") & cv.stage.isin(["raw", "corrected"])].copy()
    values = qc.pivot(index="Protein.Group", columns="stage", values="CV_percent")
    values = values.reindex(columns=["raw", "corrected"])
    if not np.isfinite(values.to_numpy(float)).all() or (values.to_numpy(float) < 0).any():
        raise ValueError("QC CV 前后必须使用相同蛋白集合，且 CV 必须为有限非负数")

    # 左闭右开：[0,1)、……、[39,40)，最后一栏包含所有 CV >=40% 的蛋白。
    edges = np.r_[np.arange(41, dtype=float), np.inf]
    counts = {stage: np.histogram(values[stage].to_numpy(float), bins=edges)[0]
              for stage in ("raw", "corrected")}
    if any(int(count.sum()) != len(values) for count in counts.values()):
        raise ValueError("QC CV 分箱计数总和与输入蛋白数量不一致")
    qc.to_csv(data_dir / "QC蛋白CV绘图数据.tsv", sep="\t", index=False, encoding="gbk")
    pd.DataFrame({
        "CV区间": [f"[{i}%, {i+1}%)" for i in range(40)] + ["[40%, +inf)"],
        "CV下界百分比": np.arange(41),
        "CV上界百分比": list(range(1, 41)) + [np.inf],
        "校正前蛋白数量": counts["raw"],
        "校正后蛋白数量": counts["corrected"],
    }).to_csv(data_dir / "QC蛋白CV分箱统计.tsv", sep="\t", index=False, encoding="gbk")

    with plt.rc_context({
        "font.family": "DejaVu Sans", "font.size": 11,
        "axes.labelcolor": "#354052", "xtick.color": "#536075",
        "ytick.color": "#536075", "pdf.fonttype": 42, "ps.fonttype": 42,
    }):
        fig, axis = plt.subplots(figsize=(10.4, 4.4), facecolor="white")
        fig.subplots_adjust(left=.085, right=.985, bottom=.19, top=.84)
        x = np.arange(41) + .5
        before_color, before_edge, after_color = "#DCE3EE", "#AAB8CC", "#5790E5"
        axis.bar(x, counts["raw"], width=.96, color=before_color,
                 edgecolor=before_edge, linewidth=.6, zorder=2)
        axis.bar(x, counts["corrected"], width=.64, color=after_color,
                 edgecolor=after_color, linewidth=.3, zorder=3)
        fig.legend(
            handles=[Patch(facecolor=before_color, edgecolor=before_edge, label="Before"),
                     Patch(facecolor=after_color, edgecolor=after_color, label="After")],
            loc="upper center", bbox_to_anchor=(.53, .97), ncol=2, frameon=False,
            handlelength=1.3, handleheight=.85, columnspacing=2.4,
        )
        axis.set_xlim(0, 41)
        upper = max(int(count.max()) for count in counts.values())
        axis.set_ylim(0, 1.12 * upper if upper else 1)
        if values.empty:
            axis.text(.5, .5, "CV unavailable\nNo proteins with at least 3 matched QC observations",
                      ha="center", va="center", transform=axis.transAxes)
        axis.set_xticks([*range(0, 40, 5), 40.5],
                        [*[str(i) for i in range(0, 40, 5)], "40+"])
        axis.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
        axis.set_xlabel("Protein CV (%) · last bin ≥40%", labelpad=11)
        axis.set_ylabel("Proteins", labelpad=10)
        axis.set_axisbelow(True)
        axis.grid(axis="y", color="#E8ECF1", linewidth=.8)
        axis.tick_params(axis="both", which="both", length=3, color="#9AA5B5")
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)
        for spine in ("bottom", "left"):
            axis.spines[spine].set_color("#A7B1BE")
            axis.spines[spine].set_linewidth(.8)
        for suffix in ("png", "pdf"):
            fig.savefig(folder / f"10_QC_蛋白CV分布图.{suffix}", dpi=300,
                        bbox_inches="tight", pad_inches=.12, facecolor="white")
        plt.close(fig)


def configure_style(meta):
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    preferred = ("Microsoft YaHei", "Noto Sans CJK SC", "SimHei", "WenQuanYi Micro Hei", "DejaVu Sans")
    family = [name for name in preferred if name in available_fonts]
    plt.rcParams.update({"font.family": family or ["DejaVu Sans"], "font.size": 11,
                         "figure.facecolor": "white", "axes.labelcolor": "#252A30", "axes.unicode_minus": False})
    unknown = [value for value in meta.condition.drop_duplicates() if value not in COLORS]
    for i, value in enumerate(unknown):
        COLORS[value] = matplotlib.colors.to_hex(matplotlib.colormaps["tab20"](i % 20))


def generate_visualizations(raw_path, corrected_path, metadata_path, output_dir=None, plate_col=None):
    """只读取三个输入表；在暂存目录生成全部图表后写入指定目录。"""
    raw_path, corrected_path, metadata_path = [Path(p).expanduser().resolve()
                                               for p in (raw_path, corrected_path, metadata_path)]
    stamp = datetime.now(timezone(timedelta(hours=8))).strftime("%Y%m%d_%H%M%S_%f")
    folder = (Path(output_dir).expanduser() if output_dir is not None
              else corrected_path.parent / f"校正前后可视化_{stamp}").resolve()
    inputs = (raw_path, corrected_path, metadata_path)
    hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    if plate_col is None:
        plate_col = "plate" if any(str(c).strip().casefold() == "plate" for c in read_table(metadata_path).columns) else None
    has_plate = plate_col is not None
    meta, excluded = load_metadata(metadata_path, plate_col=plate_col)
    meta = meta.assign(_kind=meta.group.ne("qc").astype(int),
                       _batch_order=pd.Categorical(meta.batch, categories=meta.batch.drop_duplicates(), ordered=True))
    meta = meta.sort_values(["_kind", "_batch_order"], kind="stable").drop(columns=["_kind", "_batch_order"]).reset_index(drop=True)
    meta["plot_order"] = np.arange(1, len(meta) + 1)
    if has_plate:
        meta["plate"] = meta.plate.fillna("").astype(str).str.strip().replace("", "Unspecified")
    proteins, before = read_matrix(raw_path, meta)
    _, after = read_matrix(corrected_path, meta, proteins)
    linear = {"raw": before, "corrected": after}
    logged = {stage: positive_log(values) for stage, values in linear.items()}
    for stage, values in logged.items():
        if not np.isfinite(values).any(axis=0).all():
            raise ValueError(f"{stage} 存在未检出有效正强度蛋白的样本")
    qc = meta.group.eq("qc").to_numpy()
    scopes = {"全部样本": np.arange(len(meta)), "QC": np.flatnonzero(qc)}
    if (~qc).any():
        scopes["研究样本"] = np.flatnonzero(~qc)
    configure_style(meta)
    folder.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".PACS_plot_", dir=folder) as temporary:
        staging = Path(temporary)
        data_dir = staging / "绘图数据"
        data_dir.mkdir()
        boxes = boxplots(logged, meta, scopes, staging)
        coordinates, features, pca_usage = pca_plots(
            logged, meta, proteins, scopes, staging,
            shape_col="batch", shape_title="Batch",
        )
        cvs, summary = cv_statistics(linear, logged, meta, proteins)
        cv_plots(cvs, summary, meta, staging)
        qc_cv_histogram(cvs, staging)
        for name, frame in (("样本强度统计.tsv", boxes), ("PCA坐标.tsv", coordinates),
                            ("PCA共同蛋白集合.tsv", features), ("PCA使用情况.tsv", pca_usage),
                            ("逐蛋白CV明细.tsv", cvs), ("CV统计汇总.tsv", summary)):
            frame.to_csv(data_dir / name, sep="\t", index=False, encoding="utf-8-sig")
        lines = ["# 校正前后绘图说明", "", f"校正前：{raw_path}", f"校正后：{corrected_path}",
                 f"metadata：{metadata_path}", "",
                 f"保留 {len(meta)} 个样本（QC {int(qc.sum())}，研究样本 {int((~qc).sum())}），排除 {len(excluded)} 个 Blank/Internal-calibrator；共 {len(proteins)} 个蛋白组。", "",
                 "- 箱线图：有效正强度转 log2，每个样本一个箱；上方 Raw，下方 Corrected，显示离群点。",
                 "- PCA：每个视图分别固定前后共同全检出、且 Raw 中有变化的蛋白；仅用 Raw 拟合均值和主轴，再投影校正后数据。仅均值中心化，不做蛋白方差标准化，不额外插补。轴百分比为 Raw 的解释方差。",
                 "- PCA 颜色表示 condition，点形按 batch 区分。",
                 "- CV：从线性强度计算样本标准差/均值×100%，使用前后共同有效的观测位置；组内每个蛋白至少3个有效观测。前后使用同一蛋白集合和观测位置，提琴展示完整范围，顶部数字是中位 CV。",
                 "- QC 图另列出各 batch 的 QC CV。研究样本合并 CV 和生物学组 CV 包含生物差异，不能解释为纯技术误差。",
                 "- QC CV 分布图：从原始蛋白表与本轮 PACS 校正后的最终蛋白表重新计算 pooled QC 的蛋白 CV，复用上述 CV 口径；所有 QC 合并统计。每栏1个百分点，最后一栏合并 CV≥40%。",
                 "- QC蛋白CV绘图数据.tsv 和 QC蛋白CV分箱统计.tsv 使用 GBK 编码，位于绘图数据目录；其余绘图数据表保持原编码。CV_percent 已为百分数，例如12.5表示12.5%。",
                 "- QC 在前、研究样本在后；批次按 metadata 中首次出现的顺序排列，同批内保留 metadata 顺序。",
                 "- PCA 的样本或共同蛋白不足时，图中标注不可计算原因；CV 有效观测不足时标注不可计算，不填零、不虚构结果。",
                 "",
                 "## 蛋白中位 CV（%）", "", "| 分组 | 校正前 | 校正后 |", "|---|---:|---:|"]
        for group in summary.scope.drop_duplicates():
            part = summary.loc[summary.scope.eq(group)].set_index("stage")
            values = [part.loc[s, "median_CV_percent"] for s in ("raw", "corrected")]
            text = [f"{v:.4f}" if np.isfinite(v) else "不可计算" for v in values]
            lines.append(f"| {group} | {text[0]} | {text[1]} |")
        lines.extend(["", "所有可计算 CV 均从线性输入独立复核。输入文件 SHA256 在绘图前后逐一核对，原始表未改写。"])
        (staging / "绘图说明.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
        if len(list(staging.glob("*.png"))) != 3 * len(scopes) + 1:
            raise RuntimeError("绘图数量与分析视图不一致")
        artifacts = [p for p in staging.rglob("*") if p.is_file()]
        if any((folder / p.relative_to(staging)).resolve() in inputs for p in artifacts):
            raise ValueError("输出文件与输入文件重名，请使用其他输出目录")
        if any(hashlib.sha256(p.read_bytes()).hexdigest() != before for p, before in hashes.items()):
            raise ValueError("绘图期间输入文件发生变化，请保存输入文件后重新运行")
        for artifact in artifacts:
            destination = folder / artifact.relative_to(staging)
            destination.parent.mkdir(parents=True, exist_ok=True)
            artifact.replace(destination)
    return folder


def run(
    protein_path,
    help_path,
    metadata_path,
    output_dir,
    plate_col=None,
    min_samples=2,
):
    """先完成数据处理，再用原始表和最终表绘图；全部成功后发布结果。"""
    inputs = tuple(Path(path).expanduser().resolve() for path in (protein_path, help_path, metadata_path))
    output_dir = Path(output_dir).expanduser().resolve()
    if any((output_dir / name).resolve() in inputs
           for name in ("最终蛋白定量表.tsv", "蛋白鉴定量统计柱状图.png")):
        raise ValueError("输出文件与输入文件重名，请使用其他输出目录")
    hashes = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}
    output_dir.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".PACS_all_", dir=output_dir) as temporary:
        staging = Path(temporary)
        normalize_data(*inputs, staging, plate_col, min_samples)
        final_table = staging / "最终蛋白定量表.tsv"
        plots = generate_visualizations(inputs[0], final_table, inputs[2], staging / "可视化结果", plate_col=plate_col)
        # 绘图说明引用最终交付路径，避免保留暂存目录的位置。
        explanation = plots / "绘图说明.md"
        explanation.write_text(
            explanation.read_text(encoding="utf-8").replace(str(final_table.resolve()), str(output_dir / final_table.name)),
            encoding="utf-8",
        )
        artifacts = [path for path in staging.rglob("*") if path.is_file()]
        if any((output_dir / path.relative_to(staging)).resolve() in inputs for path in artifacts):
            raise ValueError("输出文件与输入文件重名，请使用其他输出目录")
        if any(hashlib.sha256(path.read_bytes()).hexdigest() != hashes[path] for path in inputs):
            raise ValueError("处理期间输入文件发生变化，请保存输入文件后重新运行")
        for artifact in artifacts:
            destination = output_dir / artifact.relative_to(staging)
            destination.parent.mkdir(parents=True, exist_ok=True)
            artifact.replace(destination)
    return output_dir


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--protein-table",
        type=Path,
        required=True,
        help="DIA-NN 蛋白线性强度表，文件名任意",
    )
    parser.add_argument(
        "--help-table",
        "--heavy-table",
        dest="help_table",
        type=Path,
        required=True,
        help="HELP 定量报告，xlsx 必须包含 Peptide_Long；csv/tsv 须为同结构长表",
    )
    parser.add_argument(
        "--metadata", type=Path, required=True, help="含 sample_id、group、condition、batch 的单份 metadata，列名不区分大小写；须包含 pooled QC"
    )
    # 兼容既有命令中同名的固定列参数；新命令无需指定。
    parser.add_argument("--group-col", type=lambda value: value.strip().casefold(), default="group", choices=("group",), help=argparse.SUPPRESS)
    parser.add_argument("--sample-id-col", type=lambda value: value.strip().casefold(), default="sample_id", choices=("sample_id",), help=argparse.SUPPRESS)
    parser.add_argument("--plate-col", help="可选的 plate 列名；本版本不参与计算")
    parser.add_argument("--output-dir", type=Path, help="结果目录，自动创建；已存在时写入并覆盖同名结果")
    parser.add_argument(
        "--skip-imputation",
        action="store_true",
        help=argparse.SUPPRESS,  # 仅兼容旧命令，不再控制流程；当前固定跳过插补。
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        choices=(2, 3),
        default=2,
        help="每个分组进行 HELP 校正的最低样本数，默认 2；不足时保留该组原始强度",
    )
    args = parser.parse_args()
    output = (
        args.output_dir
        or args.protein_table.parent / f"PACS结果_{datetime.now():%Y%m%d_%H%M%S_%f}"
    )
    try:
        run(
            args.protein_table.expanduser(),
            args.help_table.expanduser(),
            args.metadata.expanduser(),
            output,
            args.plate_col,
            args.min_samples,
        )
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(2, f"错误：{exc}\n")


if __name__ == "__main__":
    main()
