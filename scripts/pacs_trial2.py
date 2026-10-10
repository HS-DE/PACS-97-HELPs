#!/usr/bin/env python3
"""PACS Trial 2 - two yeast-factor methods, log2 HELP factor, additive corrections.

Study samples only (HC/S, n=54) receive correction. 6 Internal-QC and 15 pooled
QC samples are strictly unchanged on the linear protein intensity scale.

Method a: ref_ye_a = median_study(mean(log2 yeast peptide 1, log2 yeast peptide 2))
Method b: ref_ye_b = mean(median_study(log2 peptide 1), median_study(log2 peptide 2))
factor1_s = ref_ye - mean(log2 yeast peptide intensities for that study sample)
factor2_s = median_pooledQC(median_97HELP(log2 L/H)) - median_97HELP(log2 L/H)_s

For method a and b independently, compare two corrections on study samples:
  (i) log2(I_corrected) = log2(I_raw) + factor1
 (ii) log2(I_corrected) = log2(I_raw) + factor1 + factor2

Each correction gets two visualization views: four groups, HC/S only.
Protein CV is on unlogged positive intensities with the reference R thresholds.
PCA-only missing-value median imputation; CV has no imputation. PNG only.
Input files / previous trials are never modified.
"""
import math
import statistics
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from filter_samples import ci, readbook, tab
from pacs_trial1 import ANN, PG, ROOT, HELP, cv_by_feature, finite_positive
from yeast_protein_qc import COLORS, GROUPS, metadata

OUT = ROOT / "results" / "pacs_trial2"
YEAST = ROOT / "sample-QC-InternalQC-HELP_yeast.xlsx"
VIEWS = {"all_groups": GROUPS, "study_only": ["HC", "S"]}
COLOR = dict(zip(GROUPS, COLORS))
STUDY_GROUPS = {"HC", "S"}
PEPTIDES = ["AADALLLK", "VNQIGTLSESIK"]
SCENARIOS = (
    ("a_factor1", "a", False),
    ("a_factor1_plus_factor2", "a", True),
    ("b_factor1", "b", False),
    ("b_factor1_plus_factor2", "b", True),
)
PCA_MIN_DETECTION = .70
CV_VISUAL_CAP = 200.


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def save_plot(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=210, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def values_from_workbook(path, sheet_index, expected_peptides, expected_samples):
    strings, roots = readbook(path)
    if sheet_index >= len(roots):
        raise ValueError(f"Missing sheet index {sheet_index} for {path}")
    rows = tab(roots[sheet_index], strings)
    names = {sid: col for col, sid in rows[0].items() if col != "A"}
    if len(names) != len(expected_samples) or set(names) != set(expected_samples):
        raise ValueError(f"Sample IDs do not match metadata for {path}")
    observed = [r.get("A", "") for r in rows[1:]]
    if len(observed) != len(set(observed)) or set(observed) != set(expected_peptides):
        raise ValueError("Incorrect or duplicated peptide identities in " + str(path))
    return {r["A"]: {sid: r.get(col, "") for sid, col in names.items()}
            for r in rows[1:]}


def log_positive(value):
    try:
        n = float(value)
        return math.log2(n) if math.isfinite(n) and n > 0 else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def compute_yeast_factors(samples, labels):
    yeast = values_from_workbook(YEAST, 0, PEPTIDES, samples)
    study_ids = [sid for sid in samples if labels[sid] in STUDY_GROUPS]
    if len(study_ids) != 54:
        raise ValueError("Trial 2 must have exactly 54 HC/S study samples")
    logged = {seq: {sid: log_positive(yeast[seq][sid]) for sid in study_ids}
              for seq in PEPTIDES}
    for seq in PEPTIDES:
        missing = [sid for sid, v in logged[seq].items() if not math.isfinite(v)]
        if missing:
            raise ValueError(f"Invalid spike-in intensity for {seq}: {missing}")
    # A and B are intentionally computed separately without substituting one.
    study_means = {sid: statistics.mean(logged[seq][sid] for seq in PEPTIDES)
                   for sid in study_ids}
    ref_a = float(np.median([study_means[sid] for sid in study_ids]))
    peptide_medians = {seq: float(np.median([logged[seq][sid] for sid in study_ids]))
                       for seq in PEPTIDES}
    ref_b = float(statistics.mean(peptide_medians.values()))
    factor_a = {sid: ref_a - study_means[sid] for sid in study_ids}
    factor_b = {sid: float(statistics.mean(ref_b - logged[seq][sid]
                                          for seq in PEPTIDES)) for sid in study_ids}
    # When both log2 peptides are detected, factor differences are an additive constant.
    expected_constant = ref_b - ref_a
    if not np.allclose(
        np.asarray([factor_b[s] - factor_a[s] for s in study_ids]),
        expected_constant, rtol=0, atol=1e-10
    ):
        raise AssertionError("Yeast a/b factor difference should be constant")
    table = []
    for sid in study_ids:
        table.append({
            "sample_id": sid, "group": labels[sid],
            "log2_AADALLLK": logged[PEPTIDES[0]][sid],
            "log2_VNQIGTLSESIK": logged[PEPTIDES[1]][sid],
            "mean_log2_yeast": study_means[sid],
            "ref_ye_method_a": ref_a,
            "ref_ye_method_b": ref_b,
            "median_log2_AADALLLK": peptide_medians[PEPTIDES[0]],
            "median_log2_VNQIGTLSESIK": peptide_medians[PEPTIDES[1]],
            "factor1_a": factor_a[sid],
            "factor1_b": factor_b[sid],
            "factor1_b_minus_a": factor_b[sid] - factor_a[sid],
        })
    write_csv(OUT / "yeast_factor_details.csv", table)
    return factor_a, factor_b, ref_a, ref_b, expected_constant


def compute_help_factor2(samples, labels):
    import zipfile
    from xml.etree import ElementTree as ET
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    with zipfile.ZipFile(HELP) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        sheet_names = [x.get("name") for x in wb.find(ns+"sheets")]
    if sheet_names != ["L", "H", "ratio_LH"]:
        raise ValueError(f"Unexpected HELP sheets: {sheet_names}")
    ss, sheets = readbook(HELP)
    entries = tab(sheets[2], ss)
    if len(entries) != 98 or len(set(r.get("A", "") for r in entries[1:])) != 97:
        raise ValueError("Expected exactly 97 unique HELP peptides")
    colmap = {sid: col for col, sid in entries[0].items() if col != "A"}
    if len(colmap) != 75 or set(colmap) != set(samples):
        raise ValueError("HELP sample ids not matched")
    medians, counts = {}, {}
    for sid in samples:
        vs = [log_positive(row.get(colmap[sid], "")) for row in entries[1:]]
        vs = [v for v in vs if math.isfinite(v)]
        if not vs:
            raise ValueError(f"No positive HELP L/H values for {sid}")
        medians[sid] = float(np.median(vs))  # log2 transform first
        counts[sid] = len(vs)
    qc_ids = [sid for sid in samples if labels[sid] == "pooled QC"]
    if len(qc_ids) != 15:
        raise ValueError("Expected 15 pooled QC samples")
    ref = float(np.median([medians[sid] for sid in qc_ids]))
    return {sid: ref - medians[sid] for sid in samples}, medians, counts, ref


def log2_matrix(raw):
    return np.log2(np.where(finite_positive(raw), raw, np.nan))


def export_protein_matrices(folder, pg, samples, raw, corrected_log2, corrected_linear, labels):
    annotation = pg.loc[:, ANN]
    frame_log = annotation.copy()
    frame_linear = annotation.copy()
    source = pg[samples]
    for j, sid in enumerate(samples):
        log_values = corrected_log2[:, j]
        frame_log[sid] = [format(x, ".15g") if math.isfinite(x) else ""
                          for x in log_values]
        if labels[sid] not in STUDY_GROUPS:
            # Exact source text; no roundtrip through float for either QC group.
            frame_linear[sid] = source[sid]
            if not frame_linear[sid].equals(source[sid]):
                raise AssertionError("Uncorrected QC source values unexpectedly changed")
        else:
            original_text = source[sid].to_numpy()
            frame_linear[sid] = [
                format(x, ".16g") if math.isfinite(x) and x > 0 else original_text[i]
                for i, x in enumerate(corrected_linear[:, j])
            ]
    frame_log.to_csv(folder / "pg_matrix_corrected_log2.tsv", sep="\t", index=False)
    frame_linear.to_csv(folder / "pg_matrix_corrected_linear.tsv", sep="\t", index=False)
    return frame_linear


def protein_cv_summary(raw, corrected, pg, samples, labels, folder):
    summaries = []
    groups_data = {}
    per_protein = []
    for group in GROUPS:
        cols = [i for i, sid in enumerate(samples) if labels[sid] == group]
        before, n_detected, nmin = cv_by_feature(raw[:, cols])
        after, n_detected2, _ = cv_by_feature(corrected[:, cols])
        if not np.array_equal(n_detected, n_detected2):
            raise AssertionError(f"Detection count changed after correction, {group}")
        eligible = np.isfinite(before) & np.isfinite(after)
        pre, post = before[eligible], after[eligible]
        if not len(pre):
            raise ValueError("No eligible proteins for group " + group)
        groups_data[group] = (pre, post)
        summaries.append({
            "group": group,
            "sample_n": len(cols),
            "min_detected_for_cv": nmin,
            "cv_eligible_proteins": int(len(pre)),
            "median_cv_before_pct": float(np.median(pre)),
            "median_cv_after_pct": float(np.median(post)),
            "median_paired_change_pp": float(np.median(post - pre)),
            "fraction_proteins_cv_lower_after_pct": float(100 * np.mean(post < pre)),
        })
        per_protein.append(pd.DataFrame({
            "protein_row_index": np.arange(1, len(pg)+1),
            "Protein.Group": pg["Protein.Group"].to_numpy(),
            "Protein.Names": pg["Protein.Names"].to_numpy(),
            "group": group,
            "sample_n": len(cols),
            "detected_n": n_detected,
            "cv_before_pct": before,
            "cv_after_pct": after,
        }))
    if not np.array_equal(groups_data["Internal-QC"][0], groups_data["Internal-QC"][1]):
        raise AssertionError("LM CV values are not unchanged")
    if not np.array_equal(groups_data["pooled QC"][0], groups_data["pooled QC"][1]):
        raise AssertionError("Pooled QC CV values are not unchanged")
    pd.concat(per_protein, ignore_index=True).to_csv(folder / "per_protein_cv.csv", index=False)
    pd.DataFrame(summaries).to_csv(folder / "group_cv_summary.csv", index=False)
    return groups_data, summaries


def subset_indices(samples, labels, groups):
    selected = [sid for group in groups for sid in samples if labels[sid] == group]
    return selected, [samples.index(sid) for sid in selected]


def plot_sample_boxes(raw, corrected, samples, labels, groups, path):
    selected, idx = subset_indices(samples, labels, groups)
    arrays = [log2_matrix(raw[:, idx]), log2_matrix(corrected[:, idx])]
    low = min(np.nanpercentile(x, 2) for x in arrays) - 1.
    high = max(np.nanpercentile(x, 98) for x in arrays) + 1.
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(max(14, len(idx)*.25), 9.5),
                             layout="constrained")
    for ax, mat, title in zip(axes, arrays, ["Original", "Corrected"]):
        seq = [mat[:, i][np.isfinite(mat[:, i])] for i in range(mat.shape[1])]
        boxes = ax.boxplot(seq, showfliers=False, patch_artist=True, widths=.55)
        for patch, sid in zip(boxes["boxes"], selected):
            patch.set_facecolor(COLOR[labels[sid]])
            patch.set_alpha(.5)
        offset = 0
        for group in groups[:-1]:
            offset += sum(labels[s] == group for s in selected)
            ax.axvline(offset+.5, ls=":", lw=.8, c="#778391")
        ax.set_ylim(low, high)
        ax.set_ylabel("Protein log2 intensity")
        ax.set_title(title)
    axes[-1].set_xticks(range(1, len(selected)+1), selected,
                        rotation=90, fontsize=6)
    fig.suptitle("Per-sample protein intensity boxplots (same y scale)")
    save_plot(fig, path)


def plot_group_medians(raw, corrected, samples, labels, groups, view, path):
    chosen, idx = subset_indices(samples, labels, groups)
    arrays = [np.nanmedian(log2_matrix(mat[:, idx]), axis=0)
              for mat in (raw, corrected)]
    rows = []
    lo = min(np.min(vals) for vals in arrays)-.3
    hi = max(np.max(vals) for vals in arrays)+.3
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(max(11, len(groups)*3.2), 5.8),
                             layout="constrained")
    for ax, vals, stage in zip(axes, arrays, ("Before", "After")):
        by_group = [[vals[j] for j, sid in enumerate(chosen) if labels[sid] == g]
                    for g in groups]
        b = ax.boxplot(by_group, patch_artist=True, showfliers=False, widths=.45)
        for patch, g in zip(b["boxes"], groups):
            patch.set_facecolor(COLOR[g])
            patch.set_alpha(.4)
        for i, (g, ys) in enumerate(zip(groups, by_group), 1):
            ax.scatter(i + np.linspace(-.15, .15, len(ys)), ys,
                       color=COLOR[g], s=18, alpha=.65, zorder=3)
        ax.set_xticks(range(1, len(groups)+1), groups)
        ax.set_ylim(lo, hi)
        ax.set_ylabel("Sample median log2 protein intensity")
        ax.set_title(stage)
        rows.extend([
            {"view": view, "sample_id": sid, "group": labels[sid],
             "stage": stage, "median_log2_intensity": float(v)}
            for sid, v in zip(chosen, vals)
        ])
    fig.suptitle("Group-wise distribution of per-sample protein median log2 intensities")
    save_plot(fig, path)
    return rows


def plot_cv(cv_data, groups, path):
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(max(12, 3.8*len(groups)), 6.4),
                             layout="constrained")
    for ax, stage, ix in zip(axes, ("Before", "After"), (0, 1)):
        for j, group in enumerate(groups, 1):
            v = cv_data[group][ix]
            capped = np.minimum(v, CV_VISUAL_CAP)
            if len(capped) > 1 and np.ptp(capped) > 0:
                violin = ax.violinplot([capped], positions=[j], widths=.75,
                                       showextrema=False)
                for patch in violin["bodies"]:
                    patch.set_facecolor(COLOR[group])
                    patch.set_alpha(.4)
            box = ax.boxplot([capped], positions=[j], widths=.16,
                             patch_artist=True, showfliers=False, manage_ticks=False)
            for patch in box["boxes"]:
                patch.set_facecolor(COLOR[group])
                patch.set_alpha(.65)
            rng = np.random.default_rng(20261010 + j)
            x = j + rng.uniform(-.22, .22, len(capped))
            ax.scatter(x, capped, c=COLOR[group], s=3, alpha=.06,
                       linewidths=0, rasterized=True, zorder=3)
            ax.text(j, CV_VISUAL_CAP+6, f"Median={np.median(v):.2f}%",
                    fontsize=9, ha="center")
        ax.set_xticks(range(1, len(groups)+1), groups)
        ax.set_xlim(.4, len(groups)+.6)
        ax.set_ylim(0, 222)
        ax.set_ylabel("Per-protein CV on linear intensities (%)")
        ax.set_title(stage)
        ax.grid(axis="y", alpha=.18)
    fig.suptitle("Per-protein CV: violin + box + jitter (raw intensity CV)")
    axes[1].text(.99, .01, "CV >200% clipped for plotting only",
                 transform=axes[1].transAxes, fontsize=8, ha="right")
    save_plot(fig, path)


def pca(raw, corrected, samples, labels, groups, view, path):
    selected, idx = subset_indices(samples, labels, groups)
    pre = raw[:, idx]
    post = corrected[:, idx]
    feature = finite_positive(pre).sum(axis=1) >= math.ceil(PCA_MIN_DETECTION * len(selected))
    if feature.sum() < 100:
        raise ValueError("Too few PCA eligible proteins")
    transformed = []
    for mat in (pre[feature], post[feature]):
        log = log2_matrix(mat)
        per_protein_median = np.nanmedian(log, axis=1)
        transformed.append(np.where(np.isfinite(log), log, per_protein_median[:, None]))
    variable = (np.std(transformed[0], axis=1) > 1e-10) & (
        np.std(transformed[1], axis=1) > 1e-10)
    if variable.sum() < 100:
        raise ValueError("Too few variable PCA proteins")
    rows, overview = [], {"view": view, "n_samples": len(selected),
                          "n_pca_proteins": int(variable.sum())}
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), layout="constrained")
    for ax, mat, stage in zip(axes, transformed, ("Before", "After")):
        fit = PCA(n_components=2, svd_solver="full")
        coords = fit.fit_transform(mat[variable].T)
        ev = fit.explained_variance_ratio_ * 100
        overview[f"PC1_{stage.lower()}_pct"] = float(ev[0])
        overview[f"PC2_{stage.lower()}_pct"] = float(ev[1])
        rows.extend([
            {"view": view, "sample_id": sid, "group": labels[sid],
             "stage": stage, "PC1": float(pt[0]), "PC2": float(pt[1]),
             "PC1_variance_pct": float(ev[0]), "PC2_variance_pct": float(ev[1])}
            for sid, pt in zip(selected, coords)
        ])
        for g in groups:
            points = np.array([p for s, p in zip(selected, coords) if labels[s] == g])
            ax.scatter(points[:, 0], points[:, 1], label=f"{g} (n={len(points)})",
                       s=48, color=COLOR[g], alpha=.85, linewidth=.6,
                       edgecolors="white")
        ax.set_xlabel(f"PC1 ({ev[0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)")
        ax.set_title(stage)
        ax.axhline(0, c="#bbc4ce", lw=.6)
        ax.axvline(0, c="#bbc4ce", lw=.6)
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle(f"PCA of log2 proteins, same {variable.sum()} features before/after")
    save_plot(fig, path)
    return rows, overview


def plot_offsets(factor_df, labels, groups, scenario, path):
    sub = pd.concat([factor_df[factor_df.group == g] for g in groups],
                    ignore_index=True)
    total = sub["applied_total_offset_log2"].to_numpy()
    fig, ax = plt.subplots(figsize=(max(11, len(sub)*.18), 5.3), layout="constrained")
    ax.bar(np.arange(len(sub)), total, width=.8,
           color=[COLOR[g] for g in sub.group])
    ax.axhline(0, color="#313843", ls="--", lw=1)
    ax.set_ylabel("Total applied additive offset in log2 space")
    ax.set_title(f"{scenario}: QC zero, study corrected")
    ax.set_xticks(np.arange(len(sub)), sub.sample_id, fontsize=6, rotation=90)
    n = 0
    for g in groups[:-1]:
        n += sum(sub.group == g)
        ax.axvline(n-.5, ls=":", color="#a3aeb6", lw=.8)
    save_plot(fig, path)


def make_scenario(sc, method, with_help, raw, pg, samples, labels, factors):
    folder = OUT / sc
    folder.mkdir(parents=True, exist_ok=True)
    off = []
    records = []
    for sid in samples:
        group = labels[sid]
        study = group in STUDY_GROUPS
        f1 = factors[sid][f"factor1_{method}"] if study else 0.
        f2 = factors[sid]["factor2_candidate_log2"] if study and with_help else 0.
        combined = f1 + f2
        off.append(combined)
        records.append({
            **factors[sid], "method": method, "stage": "factor1+factor2" if with_help else "factor1",
            "applied_factor1_log2": f1,
            "applied_factor2_log2": f2,
            "applied_total_offset_log2": combined,
            "applied_linear_multiplier": float(np.exp2(combined)),
            "corrected": study,
        })
    df = pd.DataFrame(records)
    df.to_csv(folder / "applied_factors.csv", index=False)
    offset = np.asarray(off)
    before_log = log2_matrix(raw)
    corrected_log = before_log + offset[None, :]
    corrected_linear = raw.copy()
    detected = finite_positive(raw)
    scaled = np.exp2(corrected_log)
    study_columns = np.asarray([labels[sid] in STUDY_GROUPS for sid in samples])
    study_detected = detected & study_columns[None, :]
    # Never roundtrip any QC value through log2/exp2: preserve QC raw numbers exactly.
    corrected_linear[study_detected] = scaled[study_detected]
    if not np.array_equal(finite_positive(corrected_linear), detected):
        raise AssertionError("Detection masks changed")
    for j, sid in enumerate(samples):
        if labels[sid] not in STUDY_GROUPS:
            if not np.allclose(corrected_linear[:, j], raw[:, j], rtol=0, atol=0, equal_nan=True):
                raise AssertionError("QC intensities were modified: " + sid)
    export_protein_matrices(folder, pg, samples, raw, corrected_log,
                            corrected_linear, labels)
    cv_data, cv_summary = protein_cv_summary(raw, corrected_linear, pg, samples,
                                             labels, folder)
    score_rows, pca_info, sample_med_rows = [], [], []
    for view, viewgroups in VIEWS.items():
        figs = folder / "figures" / view
        plot_sample_boxes(raw, corrected_linear, samples, labels, viewgroups,
                          figs / "01_sample_boxplots.png")
        sample_med_rows += plot_group_medians(raw, corrected_linear, samples, labels,
                                              viewgroups, view, figs / "02_group_boxplots.png")
        srows, info = pca(raw, corrected_linear, samples, labels, viewgroups,
                          view, figs / "03_pca.png")
        score_rows.extend(srows)
        pca_info.append(info)
        plot_cv(cv_data, viewgroups, figs / "04_cv_violin_box.png")
        plot_offsets(df, labels, viewgroups, sc, figs / "05_applied_offsets.png")
    write_csv(folder / "pca_scores.csv", score_rows)
    write_csv(folder / "pca_summary.csv", pca_info)
    write_csv(folder / "sample_log2_medians.csv", sample_med_rows)
    return cv_summary, pca_info


def write_report(yea, yeb, ref_a, ref_b, constant, ref_help, results):
    lines = [
        "# PACS Trial 2 - yeast spike-in + log2 HELP offsets", "",
        "This repository experiment is exploratory. It does NOT select a final PACS method.",
        "Only study samples (HC/S; 54) are corrected. All 6 Internal-QC (LM) and",
        "15 pooled QC intensities are unchanged across all four correction scenarios.",
        "Trial 1A / 1B and all original input files remain untouched.", "",
        "## Design and formulas (log2 scale)", "",
        "Let y1_s = log2(AADALLLK intensity), y2_s = log2(VNQIGTLSESIK intensity).",
        "Both peptides are required to be finite and positive in all 54 study samples.",
        "Define mean_ye_s = (y1_s + y2_s)/2.", "",
        "**Method a:** ref_ye_a = median_study(mean_ye_s); factor1_a_s = ref_ye_a - mean_ye_s.",
        "**Method b:** ref_ye_b = (median_study(y1_s)+median_study(y2_s))/2;",
        "factor1_b_s = mean(ref_ye_b-y1_s, ref_ye_b-y2_s).",
        f"ref_ye_a = {ref_a:.12g}, ref_ye_b = {ref_b:.12g},",
        f"and factor1_b - factor1_a = {constant:.12g} for **every** study sample.",
        "Because both yeast peptides are detected for every study sample, the two",
        "factors differ by a common additive constant: their within-group protein",
        "CV and centered PCA are expected to match; the absolute intensity level differs.", "",
        "For every sample, mR_s = median over 97 HELP **log2(L/H)** ratios (positive,",
        "finite L/H only). ref_mR = median(mR_s for 15 pooled QC).",
        f"ref_mR = {ref_help:.12g} in log2 ratio units.",
        "factor2_s = ref_mR - mR_s, applied only to study samples in combined scenarios.",
        "",
        "Study sample protein log2 intensity is corrected as:",
        "- factor1 only: log2(I') = log2(I) + factor1_a or factor1_b.",
        "- factor1+factor2: log2(I') = log2(I) + factor1_a/b + factor2.",
        "- Both QC types: no applied offset (0), so original linear intensities remain exact.",
        "",
        "Each scenario exports corrected protein intensities on log2 and back-",
        "transformed linear scales. The linear scale is only for ordinary CV and",
        "downstream raw-scale comparisons; it is NOT subjected to a second correction.",
        "",
        "## Four correction scenarios x two chart views", "",
        "There are FOUR correction scenarios (a F1, a F1+F2, b F1, b F1+F2).",
        "Each generates an all-groups version (75 samples) and a study-only",
        "version (54 HC/S samples): **8 figure collections, each 5 PNG images**.",
        "No SVG files are generated.",
        "",
    ]
    for sc, method, with_help in SCENARIOS:
        summ, pca_s = results[sc]
        lines += [f"### {sc}", "",
                  "| Group | n | Eligible proteins | Median CV original | Median CV corrected | Median paired Δ (pp) |",
                  "| --- | ---: | ---: | ---: | ---: | ---: |"]
        for item in summ:
            lines.append(
                f'| {item["group"]} | {item["sample_n"]} | '
                f'{item["cv_eligible_proteins"]} | '
                f'{item["median_cv_before_pct"]:.2f}% | '
                f'{item["median_cv_after_pct"]:.2f}% | '
                f'{item["median_paired_change_pp"]:+.2f} |')
        lines.extend(["", "PCA feature counts/variance:", ""])
        for p in pca_s:
            lines.append(
                f'- {p["view"]}: n={p["n_samples"]}, proteins={p["n_pca_proteins"]}; '
                f'before PC1/PC2 {p["PC1_before_pct"]:.1f}%/{p["PC2_before_pct"]:.1f}%; '
                f'after {p["PC1_after_pct"]:.1f}%/{p["PC2_after_pct"]:.1f}%.')
        lines += [""]
        for view in VIEWS:
            lines.append(f'**{view}:**')
            for file, alt in [
                ("01_sample_boxplots.png", "Sample protein distributions"),
                ("02_group_boxplots.png", "Group boxplots"),
                ("03_pca.png", "PCA"),
                ("04_cv_violin_box.png", "CV violin box jitter"),
                ("05_applied_offsets.png", "Applied log2 sample offsets"),
            ]:
                lines += [f"![{alt}]({sc}/figures/{view}/{file})", ""]
    lines += [
        "## Numeric results", "",
        "- yeast_factor_details.csv: both peptide log2 intensities, per-sample means,",
        "  peptide medians, method-a and method-b reference levels/factor1 values.",
        "- help_factor_details.csv: per-sample HELP log2 median, pooled-QC ref_mR,",
        "  candidate factor2 and number of observed HELP ratios.",
        "- factor_summary.csv: per-sample candidate factor1a/b, factor2, and actual",
        "  offset for each of the four scenarios; zero offsets for both QC types.",
        "- comparisons.csv: median protein CV per group per scenario (paired proteins).",
        "- <scenario>/applied_factors.csv: log2 and linear applied multipliers.",
        "- <scenario>/pg_matrix_corrected_log2.tsv: log2 corrected protein matrix.",
        "- <scenario>/pg_matrix_corrected_linear.tsv: back-transformed matrix for CV.",
        "- <scenario>/per_protein_cv.csv: per-protein/group original and corrected CV.",
        "- <scenario>/group_cv_summary.csv: CV median comparison.",
        "- <scenario>/pca_scores.csv and pca_summary.csv: per-view PCA results.",
        "- <scenario>/sample_log2_medians.csv: values underlying group boxplots.",
        "",
        "## Statistical definitions and safeguards", "",
        "- All offsets, yeast factors, HELP reference/factor2, boxplots and PCA",
        "  calculations are on log2 intensities/ratios.",
        "- CV (%) = sample SD(ddof=1) / arithmetic mean of **linear** positive intensities ×100;",
        "  require detected n >= max(3, ceil(70% of group n)) and mean >= 1e-8.",
        "- Missing/nonpositive protein observations are not imputed for CV and",
        "  positive detections are unchanged by the multiplicative correction.",
        "- PCA: choose >=70%-detected proteins within each view; use identical",
        "  protein set before and after in that view, median-impute missing log2",
        "  values for PCA only, center but do not z-scale features.",
        "- all_groups PCA and study_only PCA are **separately fitted** and do not",
        "  have directly comparable PC coordinate axes.",
        "- Both QC groups retain original linear intensity values. LM duplicate",
        "  pairs are intentionally preserved; nominal n=6 is not independent n=6.",
        "- Biological variation influences HC and S CV as well as technical factors.",
        "- This design cannot establish biological validity of either factor without",
        "  further evaluation (QC stability, yeast spike-in, HC/S biology).", "",
    ]
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    labels = metadata()
    pg = pd.read_csv(PG, sep="\t", dtype=str, keep_default_na=False, encoding="utf-8-sig")
    samples = list(pg.columns[5:])
    if len(samples) != 75 or len(set(samples)) != 75 or set(samples) != set(labels):
        raise ValueError("pg_matrix.tsv must contain exactly 75 metadata-matched samples")
    if list(pg.columns[:5]) != ANN:
        raise ValueError("Unexpected protein annotation columns")
    counts = {group: sum(labels[s] == group for s in samples) for group in GROUPS}
    if counts != {"Internal-QC": 6, "pooled QC": 15, "HC": 27, "S": 27}:
        raise ValueError("Unexpected sample count: " + str(counts))
    raw = pg[samples].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    factor_a, factor_b, ref_a, ref_b, constant = compute_yeast_factors(samples, labels)
    factor2, help_medians, help_counts, ref_help = compute_help_factor2(samples, labels)
    factors = {}
    for sid in samples:
        study = labels[sid] in STUDY_GROUPS
        factors[sid] = {
            "sample_id": sid, "group": labels[sid],
            "study_corrected": study,
            "factor1_a": factor_a[sid] if study else 0.,
            "factor1_b": factor_b[sid] if study else 0.,
            "factor2_candidate_log2": factor2[sid],
            "help_log2_median": help_medians[sid],
            "ref_help_log2_median": ref_help,
            "n_valid_help_peptides": help_counts[sid],
        }
    write_csv(OUT / "factor_summary.csv", list(factors.values()))
    write_csv(OUT / "help_factor_details.csv",
              [{"sample_id": sid, "group": labels[sid],
                "n_valid_help_peptides": help_counts[sid],
                "sample_mR_log2": help_medians[sid],
                "ref_mR_log2": ref_help,
                "candidate_factor2_log2": factor2[sid],
                "applied_factor2_in_combined": labels[sid] in STUDY_GROUPS}
               for sid in samples])
    results, comparison = {}, []
    for sc, method, with_help in SCENARIOS:
        cv_summary, pca_summary = make_scenario(
            sc, method, with_help, raw, pg, samples, labels, factors)
        results[sc] = (cv_summary, pca_summary)
        for row in cv_summary:
            comparison.append({"scenario": sc, "method": method,
                               "include_factor2": with_help, **row})
        print("SCENARIO", sc)
        print(pd.DataFrame(cv_summary).to_string(index=False))
        print(pd.DataFrame(pca_summary).to_string(index=False))
    write_csv(OUT / "comparisons.csv", comparison)
    # Crucial property of complete 2-peptide data: a vs b differ only by
    # a constant global shift in log2 Study sample intensities. Consequently
    # CV and centered PCA are invariant up to floating point tolerance.
    for suffix in ["factor1", "factor1_plus_factor2"]:
        a = pd.read_csv(OUT / f"a_{suffix}" / "group_cv_summary.csv")
        b = pd.read_csv(OUT / f"b_{suffix}" / "group_cv_summary.csv")
        if not np.allclose(a["median_cv_after_pct"], b["median_cv_after_pct"],
                           atol=1e-8, rtol=0):
            raise AssertionError("Methods a/b should have matching protein CVs")
    write_report(factor_a, factor_b, ref_a, ref_b, constant, ref_help, results)
    print(f"PASS Trial 2: four scenarios x two views, PNG only, "
          f"QC unchanged; {len(pg)} proteins, 54 study samples; "
          f"ref_a={ref_a:.10f}, ref_b={ref_b:.10f}, "
          f"factor1_b_minus_a={constant:.10f}, HELP_log2_ref={ref_help:.10f}")


if __name__ == "__main__":
    main()
