#!/usr/bin/env python3
"""PACS trial 1: global sample factor from the median of 97 HELP L/H ratios.

factor_s = median_HELP_ratio_s / median_QC1(median_HELP_ratio_q)
corrected_protein_intensity_ps = raw_protein_intensity_ps / factor_s

Dividing is a directional hypothesis, not a validated normalization result.
No biological group labels are used to estimate factors. The raw pg_matrix.tsv,
HELP workbook, and metadata remain unchanged.
"""
import csv
import math
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from filter_samples import ci, readbook, tab
from yeast_protein_qc import COLORS, GROUPS, ROOT, metadata

OUT = ROOT / "results" / "pacs_trial1"
FIG = OUT / "figures"
PG = ROOT / "pg_matrix.tsv"
HELP = ROOT / "LH_ratio_tables_2026.9.22—ALL-HELP.xlsx"
ANN = ["Protein.Group", "Protein.Ids", "Protein.Names",
       "Genes", "First.Protein.Description"]
GROUP_COLOR = dict(zip(GROUPS, COLORS))
PCA_DETECT_FRAC = 0.70
CV_MIN_DETECT_FRAC = 0.70
CV_MIN_N = 3
CV_MIN_MEAN = 1e-8


def save_figure(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg"):
        fig.savefig(FIG / (name + "." + extension), dpi=210,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def load_ratio(samples):
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    with zipfile.ZipFile(HELP) as z:
        book = ET.fromstring(z.read("xl/workbook.xml"))
        names = [node.get("name") for node in book.find(ns + "sheets")]
    if names != ["L", "H", "ratio_LH"]:
        raise ValueError("Expected HELP sheets L, H, ratio_LH; saw " + str(names))
    shared, sheets = readbook(HELP)
    if len(sheets) != 3:
        raise ValueError("HELP sheet count differs from three")
    rows = tab(sheets[2], shared)
    if len(rows) != 98:
        raise ValueError("Expected 97 peptide rows; got " + str(len(rows)-1))
    header = rows[0]
    ids = {name: col for col, name in header.items() if col != "A"}
    if len(ids) != 75 or set(ids) != set(samples):
        raise ValueError("HELP ratio columns differ from metadata IDs")
    seq = [r.get("A", "") for r in rows[1:]]
    if len(set(seq)) != 97 or any(not x for x in seq):
        raise ValueError("HELP peptide sequences not unique or missing")
    medians, n_valid = {}, {}
    for sample, col in ids.items():
        values = []
        for row in rows[1:]:
            try:
                x = float(row.get(col, ""))
                if math.isfinite(x) and x > 0:
                    values.append(x)
            except (TypeError, ValueError):
                pass
        if not values:
            raise ValueError("No valid HELP ratios for " + sample)
        medians[sample] = float(np.median(values))
        n_valid[sample] = len(values)
    return medians, n_valid


def write_factors(samples, medians, n_valid, sample_order):
    qc_medians = [medians[s] for s in sample_order if samples[s] == "pooled QC"]
    if len(qc_medians) != 15:
        raise ValueError("Expected 15 pooled QC samples for reference")
    reference = float(np.median(qc_medians))
    if not np.isfinite(reference) or reference <= 0:
        raise ValueError("Reference median is nonpositive or nonfinite")
    records = []
    for sid in sample_order:
        f = medians[sid] / reference
        if not math.isfinite(f) or f <= 0:
            raise ValueError("Invalid factor for " + sid)
        records.append({
            "sample_id": sid, "group": samples[sid],
            "n_valid_HELP_ratios": n_valid[sid],
            "sample_mR": medians[sid], "ref_mR": reference,
            "factor_sample_over_ref": f,
            "scale_applied_ref_over_sample": 1 / f,
            "log2_factor": math.log2(f),
        })
    pd.DataFrame(records).to_csv(OUT / "normalization_factors.csv", index=False)
    return pd.DataFrame(records), reference


def finite_positive(a):
    return np.isfinite(a) & (a > 0)


def cv_by_feature(matrix):
    """CV for each protein within one sample group, as in supplied R defaults."""
    good = finite_positive(matrix)
    n = good.sum(axis=1)
    valid = np.where(good, matrix, 0.0)
    total = valid.sum(axis=1)
    avg = np.divide(total, n, out=np.full(len(n), np.nan), where=n > 0)
    squares = (valid * valid).sum(axis=1)
    sample_var = np.divide(
        np.maximum(squares - np.divide(total * total, n,
                                      out=np.zeros_like(total), where=n > 0), 0.0),
        n - 1,
        out=np.full(len(n), np.nan), where=n > 1)
    cv = np.divide(100 * np.sqrt(sample_var), avg,
                   out=np.full(len(n), np.nan), where=avg >= CV_MIN_MEAN)
    needed = max(CV_MIN_N, math.ceil(matrix.shape[1] * CV_MIN_DETECT_FRAC))
    cv[(n < needed) | (~np.isfinite(cv))] = np.nan
    return cv, n, needed


def plot_factors(factors, reference):
    arranged = pd.concat([factors[factors.group == g] for g in GROUPS])
    fig, ax = plt.subplots(figsize=(16, 5.5), layout="constrained")
    ax.bar(range(len(arranged)), arranged.factor_sample_over_ref,
           color=[GROUP_COLOR[g] for g in arranged.group], width=0.85)
    ax.axhline(1.0, color="black", ls="--", lw=1, label="pooled-QC reference factor = 1")
    boundaries = np.cumsum([sum(arranged.group == g) for g in GROUPS])[:-1]
    for boundary in boundaries:
        ax.axvline(boundary-0.5, color="#999999", ls=":", lw=0.8)
    ax.set_xticks(range(len(arranged)), arranged.sample_id, rotation=90, fontsize=6)
    ax.set_ylabel("factor = sample_mR / ref_mR")
    ax.set_title("Trial 1: HELP ratio factors by sample group")
    ax.legend(loc="upper right", frameon=False)
    save_figure(fig, "01_sample_factors")


def plot_sample_boxes(before, after, sample_ids, groups):
    """Each sample's distribution of protein log2 intensities, fixed axis across states."""
    import matplotlib.patches as patches
    ids = [sid for g in GROUPS for sid in sample_ids if groups[sid] == g]
    ind = [sample_ids.index(sid) for sid in ids]
    arrays = []
    for mat in (before, after):
        data = np.where(finite_positive(mat), np.log2(np.where(finite_positive(mat), mat, np.nan)), np.nan)
        arrays.append([data[:, j][np.isfinite(data[:, j])] for j in ind])
    ymin = min(np.percentile(a, 5) for batch in arrays for a in batch if len(a))
    ymax = max(np.percentile(a, 95) for batch in arrays for a in batch if len(a))
    fig, axs = plt.subplots(2, 1, figsize=(19, 10), sharex=True, layout="constrained")
    for ax, batch, title in zip(axs, arrays, ["Before correction", "After: divide intensity by factor"]):
        positions = np.arange(1, len(ids)+1)
        bp = ax.boxplot(batch, positions=positions, widths=.6, patch_artist=True,
                        manage_ticks=False, showfliers=False,
                        medianprops={"color": "#172433", "linewidth": .8},
                        whiskerprops={"linewidth": .5}, capprops={"linewidth": .5})
        for patch, sid in zip(bp["boxes"], ids):
            patch.set_facecolor(GROUP_COLOR[groups[sid]])
            patch.set_alpha(.6)
        offset = 0
        for g in GROUPS[:-1]:
            offset += sum(groups[sid] == g for sid in ids)
            ax.axvline(offset+.5, color="#7a8190", ls=":", lw=1)
        ax.set_ylabel("Protein log2 intensity")
        ax.set_title(title)
        ax.set_ylim(ymin-2, ymax+2)
    axs[-1].set_xticks(np.arange(1, len(ids)+1), ids, rotation=90, fontsize=6)
    fig.suptitle("Per-sample protein intensity distributions (same y-axis scale)", fontsize=14)
    save_figure(fig, "02_sample_protein_boxplots")


def plot_group_sample_medians(before, after, sample_ids, groups):
    """Group-level comparison of sample median log2 intensities."""
    median_arrays = []
    rows = []
    for stage, mat in [("Before", before), ("After", after)]:
        log = np.where(finite_positive(mat), np.log2(np.where(finite_positive(mat), mat, np.nan)), np.nan)
        sample_medians = np.nanmedian(log, axis=0)
        median_arrays.append(sample_medians)
        for sid, value in zip(sample_ids, sample_medians):
            rows.append({"sample_id": sid, "group": groups[sid], "stage": stage,
                         "median_log2_protein_intensity": value})
    pd.DataFrame(rows).to_csv(OUT / "sample_protein_medians.csv", index=False)
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True, layout="constrained")
    lo = min(float(np.min(v)) for v in median_arrays)-.3
    hi = max(float(np.max(v)) for v in median_arrays)+.3
    for ax, medians, stage in zip(axs, median_arrays, ["Before", "After"]):
        batches = [[medians[j] for j, sid in enumerate(sample_ids) if groups[sid] == g]
                   for g in GROUPS]
        bp = ax.boxplot(batches, patch_artist=True, widths=.4, showfliers=False)
        for box, color in zip(bp["boxes"], COLORS):
            box.set_facecolor(color)
            box.set_alpha(.4)
        for gi, (vals, col) in enumerate(zip(batches, COLORS), 1):
            jitter = np.linspace(-.17, .17, len(vals))
            ax.scatter(gi + jitter, vals, c=col, s=16, alpha=.65, zorder=3)
        ax.set_xticks(range(1, 5), GROUPS)
        ax.set_title(stage)
        ax.set_ylim(lo, hi)
        ax.set_ylabel("Median log2 protein intensity per sample")
    fig.suptitle("Group comparison of sample median protein abundance", fontsize=14)
    save_figure(fig, "03_group_sample_median_boxes")
    return pd.DataFrame(rows)


def pca_analysis(before, after, sample_ids, groups):
    """Fixed feature set. Per-protein log2 median imputation ONLY for PCA."""
    detection = finite_positive(before)
    feature_mask = detection.sum(axis=1) >= math.ceil(PCA_DETECT_FRAC * len(sample_ids))
    if feature_mask.sum() < 100:
        raise ValueError("Insufficient detected protein features for PCA")
    stage_arrays = []
    explained = []
    rows = []
    for stage, mat in [("Before", before), ("After", after)]:
        picked = mat[feature_mask]
        logged = np.where(finite_positive(picked),
                          np.log2(np.where(finite_positive(picked), picked, np.nan)), np.nan)
        feature_medians = np.nanmedian(logged, axis=1)
        imputed = np.where(np.isfinite(logged), logged, feature_medians[:, None])
        variation = np.std(imputed, axis=1, ddof=0)
        stage_arrays.append(imputed)
    constant = (np.std(stage_arrays[0], axis=1) > 1e-10) & (np.std(stage_arrays[1], axis=1) > 1e-10)
    if constant.sum() < 100:
        raise ValueError("Too few variable proteins for PCA")
    for stage, arr in zip(["Before", "After"], stage_arrays):
        X = arr[constant].T
        pca = PCA(n_components=2, svd_solver="full")
        coordinates = pca.fit_transform(X)
        explained.append(pca.explained_variance_ratio_ * 100)
        for sid, point in zip(sample_ids, coordinates):
            rows.append({"sample_id": sid, "group": groups[sid], "stage": stage,
                         "PC1": point[0], "PC2": point[1],
                         "PC1_explained_pct": explained[-1][0],
                         "PC2_explained_pct": explained[-1][1]})
    pd.DataFrame(rows).to_csv(OUT / "pca_scores.csv", index=False)

    fig, axs = plt.subplots(1, 2, figsize=(14, 6.3), layout="constrained")
    for ax, stage, ev in zip(axs, ["Before", "After"], explained):
        p = [r for r in rows if r["stage"] == stage]
        for g in GROUPS:
            subset = [r for r in p if r["group"] == g]
            ax.scatter([r["PC1"] for r in subset], [r["PC2"] for r in subset],
                       s=45, alpha=.8, color=GROUP_COLOR[g], label=f"{g} (n={len(subset)})",
                       edgecolors="white", linewidths=.5)
        ax.set_xlabel(f"PC1 ({ev[0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)")
        ax.set_title(stage)
        ax.axhline(0, color="#b9c0c8", linewidth=.5)
        ax.axvline(0, color="#b9c0c8", linewidth=.5)
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle(f"Protein abundance PCA, same {int(constant.sum())} proteins before/after", fontsize=14)
    save_figure(fig, "04_pca_before_after")
    return int(constant.sum()), explained


def protein_cv_analysis(before, after, annotations, sample_ids, groups):
    rows = []
    summaries = []
    for g in GROUPS:
        indices = [j for j, sid in enumerate(sample_ids) if groups[sid] == g]
        rawcv, raw_n, need = cv_by_feature(before[:, indices])
        newcv, new_n, _ = cv_by_feature(after[:, indices])
        if not np.array_equal(raw_n, new_n):
            raise AssertionError("Correction changed protein detection counts")
        for i in range(len(annotations)):
            rows.append({
                "row_index_1based": i+1,
                "Protein.Group": annotations.iloc[i]["Protein.Group"],
                "Protein.Names": annotations.iloc[i]["Protein.Names"],
                "group": g, "sample_n": len(indices),
                "cv_min_detected_n": need,
                "detected_n": int(raw_n[i]),
                "cv_before_pct": rawcv[i] if np.isfinite(rawcv[i]) else np.nan,
                "cv_after_pct": newcv[i] if np.isfinite(newcv[i]) else np.nan,
            })
        paired = np.isfinite(rawcv) & np.isfinite(newcv)
        before_values, after_values = rawcv[paired], newcv[paired]
        if not paired.any():
            raise ValueError("No group CVs estimable for " + g)
        summaries.append({
            "group": g, "n_samples": len(indices),
            "cv_min_detected_n": need,
            "eligible_proteins": int(paired.sum()),
            "median_cv_before_pct": float(np.median(before_values)),
            "median_cv_after_pct": float(np.median(after_values)),
            "median_paired_cv_change_pct_points": float(np.median(after_values-before_values)),
            "fraction_proteins_cv_lower_after_pct": float(100 * np.mean(after_values < before_values)),
        })
    pd.DataFrame(rows).to_csv(OUT / "per_protein_cv.csv", index=False)
    sums = pd.DataFrame(summaries)
    sums.to_csv(OUT / "group_cv_summary.csv", index=False)
    return rows, sums


def plot_cv(rows, summary):
    """Violin + box + faint per-protein points, as in reference R plot."""
    fig, axs = plt.subplots(1, 2, figsize=(15, 6.5), sharey=True, layout="constrained")
    for ax, key, title in zip(axs, ["cv_before_pct", "cv_after_pct"],
                              ["Before correction", "After: divide by factor"]):
        for i, (g, color) in enumerate(zip(GROUPS, COLORS), 1):
            v = np.asarray([r[key] for r in rows
                            if r["group"] == g and np.isfinite(r[key])], dtype=float)
            # Display cap does not alter exported numerical CVs or reported medians.
            capped = np.minimum(v, 200)
            if len(capped) >= 2 and np.ptp(capped) > 0:
                parts = ax.violinplot([capped], positions=[i], widths=.75,
                                      showmeans=False, showextrema=False)
                for patch in parts["bodies"]:
                    patch.set_facecolor(color)
                    patch.set_alpha(.4)
            box = ax.boxplot([capped], positions=[i], widths=.18,
                             patch_artist=True, showfliers=False, manage_ticks=False)
            for patch in box["boxes"]:
                patch.set_facecolor(color)
                patch.set_alpha(.65)
            offset = np.linspace(-.22, .22, len(capped))
            rng = np.random.default_rng(20261009 + i)
            rng.shuffle(offset)
            ax.scatter(i + offset, capped, color=color, alpha=.055, s=3,
                       linewidths=0, rasterized=True, zorder=3)
            ax.text(i, 207, f"Median={np.median(v):.1f}%", ha="center", fontsize=9)
        ax.set_xticks(range(1, 5), GROUPS, rotation=10)
        ax.set_xlim(.4, 4.6)
        ax.set_ylim(0, 221)
        ax.set_title(title)
        ax.set_ylabel("Per-protein raw-intensity CV (%)")
        ax.grid(axis="y", alpha=.17)
    fig.suptitle("Protein CV distribution by sample group (violin + box + points)", fontsize=14)
    axs[1].text(.99, .01, "Values >200% are capped at 200% for display only",
                transform=axs[1].transAxes, ha="right", va="bottom", fontsize=8)
    save_figure(fig, "05_protein_cv_violin_box")


def compare_directions(raw, multiplied, corrected, samples, groups, current_summary):
    """Compare direction empirically; keep divide as trial-1 primary result."""
    items = []
    for g in GROUPS:
        ix = [i for i, sid in enumerate(samples) if groups[sid] == g]
        a, na, _ = cv_by_feature(raw[:, ix])
        b, nb, _ = cv_by_feature(corrected[:, ix])
        m, nm, _ = cv_by_feature(multiplied[:, ix])
        if not np.array_equal(na, nb) or not np.array_equal(na, nm):
            raise AssertionError("Direction scaling unexpectedly changed detection")
        matched = np.isfinite(a) & np.isfinite(b) & np.isfinite(m)
        a, b, m = a[matched], b[matched], m[matched]
        items.append({
            "group": g,
            "n_samples": len(ix),
            "eligible_proteins_paired": int(len(a)),
            "median_cv_uncorrected_pct": float(np.median(a)),
            "median_cv_divide_pct": float(np.median(b)),
            "median_cv_multiply_pct": float(np.median(m)),
            "median_divide_minus_original_pp": float(np.median(b - a)),
            "median_multiply_minus_original_pp": float(np.median(m - a)),
            "fraction_cv_lower_divide_pct": float(100 * np.mean(b < a)),
            "fraction_cv_lower_multiply_pct": float(100 * np.mean(m < a)),
        })
    table = pd.DataFrame(items)
    table.to_csv(OUT / "direction_sensitivity_cv.csv", index=False)
    fig, ax = plt.subplots(figsize=(11.5, 6), layout="constrained")
    xs = np.arange(len(GROUPS))
    for shift, name, key, fill in [
        (-0.25, "Original", "median_cv_uncorrected_pct", "#7a8592"),
        (0, "Divide by factor", "median_cv_divide_pct", "#bf5562"),
        (0.25, "Multiply by factor", "median_cv_multiply_pct", "#48a48a"),
    ]:
        values = table[key].to_numpy()
        bars = ax.bar(xs + shift, values, width=.24, color=fill, label=name)
        ax.bar_label(bars, fmt="%.1f%%", padding=3, fontsize=9)
    ax.set_xticks(xs, GROUPS)
    ax.set_ylim(0, max(table[["median_cv_uncorrected_pct",
                             "median_cv_divide_pct",
                             "median_cv_multiply_pct"]].max()) * 1.25)
    ax.set_ylabel("Median across per-protein within-group CVs (%)")
    ax.set_title("Sensitivity check: correction direction matters")
    ax.legend(frameon=False)
    save_figure(fig, "06_cv_direction_sensitivity")
    return table


def report(factors, ref, pg_n, cv_summary, n_pca_features, explained, directions):
    med_qc = factors[factors.group == "pooled QC"]["factor_sample_over_ref"].median()
    lines = [
        "# PACS Trial 1 — 97 HELP median L/H-ratio factor", "",
        "This is an exploratory, **unvalidated** factor correction, not a recommended final PACS method.",
        "Original pg_matrix.tsv and HELP/metadata input files are NOT modified.", "",
        "## Definition", "",
        "For sample s, sample_mR = median of 97 HELP raw L/H ratios (excluding missing/nonpositive).",
        "ref_mR = median of the 15 pooled QC sample_mR values (each QC has one median).",
        "factor_s = sample_mR / ref_mR.",
        "**Applied correction: protein corrected intensity = original intensity / factor_s.**",
        "Thus applied multiplier = ref_mR / sample_mR. This assumes a higher HELP L/H",
        "ratio indicates higher multiplicative protein intensity bias. This direction",
        "cannot be proven from L/H ratios alone (e.g., H-driven drift would break it).",
        "No HC, S or Internal-QC values contribute to ref_mR.",
        "",
        "## Data and assumptions", "",
        "- HELP: 97 sequences, raw ratio_LH sheet; medians taken on raw ratio values, not log2.",
        "- Protein matrix: 5 annotation columns plus exactly 75 metadata-matched samples.",
        f"- Protein Group rows: {pg_n}; all rows (human + yeast) retained.",
        "- Groups: Internal-QC n=6; pooled QC n=15; HC n=27; S n=27.",
        "- Intentional duplicated Internal-QC pairs retained as given; not 6 independent values.",
        "- No factor fitting to protein data, no biological-group-based factor selection.",
        "- Missing protein intensities stay missing; zero stays zero.",
        "- CV uses detected finite positive raw intensities; no log-transform or imputation.",
        "- CV requires detected n >= max(3, ceil(0.70*n_group)) and mean >= 1e-8;",
        "  sample SD (ddof=1) / mean * 100, matching the reference R CV rule.",
        "- CV of HC/S may include biological variation, not solely technical noise.",
        "- PCA uses log2 positive intensities; >=70% overall detection filter.",
        "- For PCA **only**, remaining missing log2 values are imputed with per-protein",
        "  median within each stage; identical selected features before/after; features",
        "  are mean-centered, not z-score standardized; samples are PCA observations.",
        f"- PCA features used: {n_pca_features}. Before PC1/PC2: "
        f"{explained[0][0]:.1f}%/{explained[0][1]:.1f}%; after: "
        f"{explained[1][0]:.1f}%/{explained[1][1]:.1f}%.",
        f"- ref_mR = {ref:.10g}; median pooled QC factor = {med_qc:.6g}.",
        "", "## Group-wise per-protein CV", "",
        "| Group | n | Eligible protein groups | Median CV before | Median CV after | "
        "Median paired CV change (percentage points) | Proteins with lower CV after |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in cv_summary.to_dict("records"):
        lines.append(
            f'| {r["group"]} | {r["n_samples"]} | {r["eligible_proteins"]} | '
            f'{r["median_cv_before_pct"]:.2f}% | {r["median_cv_after_pct"]:.2f}% | '
            f'{r["median_paired_cv_change_pct_points"]:+.2f} | '
            f'{r["fraction_proteins_cv_lower_after_pct"]:.1f}% |'
        )
    lines.extend([
        "", "## Multiplication-direction sensitivity (not the primary corrected matrix)", "",
        "Both candidate directions use exactly the same sample factor:",
        "divide = intensity / factor; multiply = intensity * factor.",
        "This comparison uses the same CV detection threshold and the same",
        "paired eligible proteins in each group; neither direction is selected",
        "as a validated normalization method.",
        "",
        "| Group | Original median CV | Divide median CV | Multiply median CV |",
        "| --- | ---: | ---: | ---: |",
    ])
    for r in directions.to_dict("records"):
        lines.append(
            f'| {r["group"]} | {r["median_cv_uncorrected_pct"]:.2f}% | '
            f'{r["median_cv_divide_pct"]:.2f}% | '
            f'{r["median_cv_multiply_pct"]:.2f}% |'
        )
    lines.extend([
        "", "The multiplier-direction check is a diagnostic, not a third-party",
        "validated quality benchmark. Further assess yeast spike-ins and",
        "biological preservation before committing to a direction.",
    ])
    lines.extend([
        "", "## Figures", "",
        "Group colors throughout: Internal-QC blue, pooled QC orange, HC green, S purple.",
        "![Sample factors](figures/01_sample_factors.png)", "",
        "![Protein distributions per sample](figures/02_sample_protein_boxplots.png)", "",
        "![Group sample median boxplots](figures/03_group_sample_median_boxes.png)", "",
        "![PCA](figures/04_pca_before_after.png)", "",
        "![Violin/box/points per-protein CV](figures/05_protein_cv_violin_box.png)", "",
        "![CV direction comparison](figures/06_cv_direction_sensitivity.png)", "",
        "Figures are saved in PNG and SVG formats.",
        "", "## Data files", "",
        "- normalization_factors.csv: all sample medians, common pooled QC reference,",
        "  factor, actual applied multiplier and log2 factor.",
        "- pg_matrix_corrected.tsv: protein intensities divided by each sample factor.",
        "- sample_protein_medians.csv: before/after per-sample median log2 protein intensity.",
        "- pca_scores.csv: sample PCA coordinates and variance percentages.",
        "- per_protein_cv.csv: each protein-group row × group, before/after raw-intensity CV.",
        "- group_cv_summary.csv: four-group protein CV comparison.",
        "- direction_sensitivity_cv.csv: original vs divide vs multiply per-protein CV.",
        "- pg_matrix_multiply_sensitivity.tsv: alternative direction, not the primary matrix.",
        "- scripts/pacs_trial1.py: reproducible calculation and plots.",
        "- .github/workflows/pacs-trial1.yml: auto-regenerate on input/script changes.",
        "", "## Interpretation warning", "",
        "A shift in median CV or PCA is descriptive; this does not establish whether",
        "HELP ratio tracks protein abundance drift or whether signal should instead",
        "be multiplied, left unmodified, or modeled in other ways. Use pooled-QC",
        "repeatability, yeast spike-in controls and biological-group preservation",
        "to judge any future PACS method.", "",
    ])
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    groups = metadata()
    if len(groups) != 75:
        raise ValueError("Expected exactly 75 metadata sample IDs")
    pg = pd.read_csv(PG, sep="\t", encoding="utf-8-sig", dtype={x: str for x in ANN})
    if list(pg.columns[:5]) != ANN or len(pg.columns) != 80:
        raise ValueError("pg_matrix.tsv annotation/sample column schema unexpected")
    samples = list(pg.columns[5:])
    if len(set(samples)) != 75 or set(samples) != set(groups):
        raise ValueError("pg_matrix.tsv sample IDs do not match metadata exactly")
    raw = pg.loc[:, samples].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(raw).any():
        raise ValueError("No numeric intensities present")
    medians, nvalid = load_ratio(groups)
    factors, reference = write_factors(groups, medians, nvalid, samples)
    multiplier = factors["scale_applied_ref_over_sample"].to_numpy()
    corrected = raw * multiplier[None, :]
    corrected_pg = pg.loc[:, ANN].copy()
    corrected_pg[samples] = corrected
    corrected_pg.to_csv(OUT / "pg_matrix_corrected.tsv", sep="\t", index=False,
                        na_rep="", float_format="%.10g")
    if not np.allclose(corrected[finite_positive(raw)],
                       (raw * multiplier[None, :])[finite_positive(raw)], rtol=1e-10):
        raise AssertionError("Factor arithmetic failed")
    if not np.array_equal(finite_positive(raw), finite_positive(corrected)):
        raise AssertionError("Correction changed detection mask")
    plot_factors(factors, reference)
    plot_sample_boxes(raw, corrected, samples, groups)
    sample_medians = plot_group_sample_medians(raw, corrected, samples, groups)
    nfeatures, explained = pca_analysis(raw, corrected, samples, groups)
    protein_cv, summary = protein_cv_analysis(raw, corrected, pg.loc[:, ANN], samples, groups)
    plot_cv(protein_cv, summary)
    multiplied = raw * factors["factor_sample_over_ref"].to_numpy()[None, :]
    alt_pg = pg.loc[:, ANN].copy()
    alt_pg[samples] = multiplied
    alt_pg.to_csv(OUT / "pg_matrix_multiply_sensitivity.tsv", sep="\t",
                  index=False, na_rep="", float_format="%.10g")
    directions = compare_directions(raw, multiplied, corrected, samples, groups, summary)
    report(factors, reference, len(pg), summary, nfeatures, explained, directions)
    print("Correction direction sensitivity:")
    print(directions.to_string(index=False))
    print(f"PASS trial1: samples={len(samples)}, proteins={len(pg)}, "
          f"HELP_peptides=97, pooledQC=15, reference_mR={reference:.10g}, "
          f"PCA_features={nfeatures}")
    print(summary.to_string(index=False))
    for group in GROUPS:
        vals = sample_medians[sample_medians["group"] == group]
        before = vals[vals["stage"] == "Before"].median_log2_protein_intensity
        after = vals[vals["stage"] == "After"].median_log2_protein_intensity
        print(f"{group}: sample median log2 ranges "
              f"before={before.min():.4f}..{before.max():.4f}, "
              f"after={after.min():.4f}..{after.max():.4f}")


if __name__ == "__main__":
    main()
