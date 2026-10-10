#!/usr/bin/env python3
"""PACS Trial 1B: pooled-QC referenced HELP ratio correction, LM untouched.

S_mR = median(valid positive L/H ratio among 97 HELP peptides, for each sample)
ref_mR = median(15 pooled-QC S_mR)
candidate_factor = S_mR / ref_mR
corrected protein intensity = original / candidate_factor for pooled-QC, HC, S
Internal-QC (LM) intensities remain EXACTLY as in the input pg_matrix.tsv.

Two independent plotting views: all four groups and study-only HC/S.
Only PNG plots are saved. No imputation is used for CV, and PCA-only
imputation is documented. Original files and Trial 1A are left untouched.
"""
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from pacs_trial1 import (
    ANN, COLORS, GROUPS, PG, ROOT, cv_by_feature, finite_positive, load_ratio
)
from yeast_protein_qc import metadata

OUT = ROOT / "results" / "pacs_trial1b"
PALETTE = dict(zip(GROUPS, COLORS))
VIEWS = {
    "all_groups": list(GROUPS),
    "study_only": ["HC", "S"],
}
STAGE_NAMES = ("Before", "After")
MIN_PCA_DETECTION = 0.70
CV_PLOT_CAP = 200.0


def save(fig, view, stem):
    target = OUT / "figures" / view
    target.mkdir(parents=True, exist_ok=True)
    fig.savefig(target / (stem + ".png"), dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def factor_table(labels, samples):
    medians, counts = load_ratio(labels)  # same validated 97-HELP reader as Trial 1A
    qcs = [medians[s] for s in samples if labels[s] == "pooled QC"]
    if len(qcs) != 15:
        raise ValueError("Expected exactly 15 pooled-QC samples")
    ref = float(np.median(qcs))
    if not np.isfinite(ref) or ref <= 0:
        raise ValueError("Invalid pooled QC reference median")
    records = []
    for sid in samples:
        group = labels[sid]
        candidate = medians[sid] / ref
        if not np.isfinite(candidate) or candidate <= 0:
            raise ValueError("Invalid HELP factor for " + sid)
        applied = group != "Internal-QC"
        records.append({
            "sample_id": sid,
            "group": group,
            "n_valid_HELP_ratios": counts[sid],
            "sample_mR": medians[sid],
            "ref_mR": ref,
            "candidate_factor_sample_over_ref": candidate,
            "correction_applied": applied,
            "applied_factor": candidate if applied else 1.0,
            "applied_multiplier": 1.0 / candidate if applied else 1.0,
        })
    table = pd.DataFrame(records)
    table.to_csv(OUT / "normalization_factors.csv", index=False)
    return table, ref


def protein_cv(raw, corrected, annotations, samples, labels):
    """Use reference R defaults for CV, individually for each Protein.Group row."""
    all_rows, group_rows, arrays = [], [], {}
    for group in GROUPS:
        cols = [j for j, s in enumerate(samples) if labels[s] == group]
        pre, npre, min_n = cv_by_feature(raw[:, cols])
        post, npost, _ = cv_by_feature(corrected[:, cols])
        if not np.array_equal(npre, npost):
            raise AssertionError("Detection counts changed unexpectedly")
        if group == "Internal-QC" and not np.allclose(pre, post, equal_nan=True):
            raise AssertionError("Internal-QC CV changed despite leaving LM unchanged")
        paired = np.isfinite(pre) & np.isfinite(post)
        if not paired.any():
            raise ValueError("No eligible CV values: " + group)
        frame = pd.DataFrame({
            "row_index_1based": np.arange(1, len(annotations)+1),
            "Protein.Group": annotations["Protein.Group"].to_numpy(),
            "Protein.Names": annotations["Protein.Names"].to_numpy(),
            "group": group,
            "sample_n": len(cols),
            "cv_min_detected_n": min_n,
            "detected_n": npre,
            "cv_before_pct": pre,
            "cv_after_pct": post,
        })
        all_rows.append(frame)
        baseline, changed = pre[paired], post[paired]
        group_rows.append({
            "group": group, "sample_n": len(cols),
            "cv_min_detected_n": min_n,
            "eligible_proteins_paired": int(np.sum(paired)),
            "median_cv_before_pct": float(np.median(baseline)),
            "median_cv_after_pct": float(np.median(changed)),
            "median_paired_change_pp": float(np.median(changed - baseline)),
            "proteins_cv_lower_after_pct": float(100 * np.mean(changed < baseline)),
        })
        arrays[group] = (baseline, changed)
    pd.concat(all_rows, ignore_index=True).to_csv(OUT / "per_protein_cv.csv", index=False)
    result = pd.DataFrame(group_rows)
    result.to_csv(OUT / "group_cv_summary.csv", index=False)
    return result, arrays


def plot_factors(factors, view, chosen):
    data = pd.concat([factors[factors.group == g] for g in chosen], ignore_index=True)
    applied = data["applied_multiplier"].to_numpy()
    fig, ax = plt.subplots(figsize=(max(11, len(data)*.18), 5.5), layout="constrained")
    ax.bar(np.arange(len(data)), applied,
           color=[PALETTE[g] for g in data.group], width=.8)
    ax.axhline(1, lw=1, ls="--", color="#303744")
    ax.set_xticks(np.arange(len(data)), data.sample_id, rotation=90, fontsize=6)
    ax.set_ylabel("Applied intensity multiplier")
    ax.set_title("Trial 1B - actual multiplier (Internal-QC is exactly 1)")
    offset = 0
    for g in chosen[:-1]:
        offset += sum(data.group == g)
        ax.axvline(offset-.5, lw=.8, ls=":", color="#89939c")
    save(fig, view, "01_applied_multipliers")


def log2_matrix(mat):
    good = finite_positive(mat)
    return np.log2(np.where(good, mat, np.nan))


def plot_per_sample_boxes(raw, corrected, samples, labels, view, chosen):
    selected = [s for g in chosen for s in samples if labels[s] == g]
    ix = [samples.index(s) for s in selected]
    pre = log2_matrix(raw[:, ix])
    post = log2_matrix(corrected[:, ix])
    shared_lim = (
        min(float(np.nanpercentile(pre, 2)), float(np.nanpercentile(post, 2))) - 1,
        max(float(np.nanpercentile(pre, 98)), float(np.nanpercentile(post, 98))) + 1,
    )
    fig, axes = plt.subplots(2, 1, figsize=(max(14, len(ix)*.23), 10),
                             sharex=True, layout="constrained")
    for ax, mat, title in zip(axes, (pre, post), STAGE_NAMES):
        series = [mat[:, k][np.isfinite(mat[:, k])] for k in range(len(ix))]
        boxes = ax.boxplot(series, patch_artist=True, widths=.64, showfliers=False)
        for patch, sid in zip(boxes["boxes"], selected):
            patch.set_facecolor(PALETTE[labels[sid]])
            patch.set_alpha(.55)
        ax.set_ylabel("Protein log2 intensity")
        ax.set_ylim(*shared_lim)
        ax.set_title(title + " - per-sample protein distribution")
        offset = 0
        for g in chosen[:-1]:
            offset += sum(labels[sid] == g for sid in selected)
            ax.axvline(offset+.5, ls=":", color="#87909e")
    axes[-1].set_xticks(range(1, len(selected)+1), selected, rotation=90, fontsize=6)
    fig.suptitle("PACS Trial 1B - per-sample protein boxplots (same y axis)")
    save(fig, view, "02_per_sample_boxplots")


def plot_group_medians(raw, corrected, samples, labels, view, chosen):
    selected = [s for g in chosen for s in samples if labels[s] == g]
    ix = [samples.index(s) for s in selected]
    medians = [np.nanmedian(log2_matrix(m[:, ix]), axis=0) for m in (raw, corrected)]
    lo = min(float(np.min(v)) for v in medians) - .4
    hi = max(float(np.max(v)) for v in medians) + .4
    fig, axes = plt.subplots(1, 2, figsize=(max(12, len(chosen)*3.2), 5.8),
                             sharey=True, layout="constrained")
    rows = []
    for ax, stage, vals in zip(axes, STAGE_NAMES, medians):
        bins = [[vals[k] for k, sid in enumerate(selected) if labels[sid] == g] for g in chosen]
        boxes = ax.boxplot(bins, showfliers=False, patch_artist=True, widths=.45)
        for patch, group in zip(boxes["boxes"], chosen):
            patch.set_facecolor(PALETTE[group])
            patch.set_alpha(.4)
        for j, (g, arr) in enumerate(zip(chosen, bins), 1):
            ax.scatter(j + np.linspace(-.16, .16, len(arr)), arr,
                       c=PALETTE[g], s=17, alpha=.7, zorder=3)
        ax.set_ylim(lo, hi)
        ax.set_xticks(range(1, len(chosen)+1), chosen)
        ax.set_ylabel("Median log2 protein intensity per sample")
        ax.set_title(stage)
        for sid, value in zip(selected, vals):
            rows.append({"view": view, "sample_id": sid, "group": labels[sid],
                         "stage": stage, "median_log2_intensity": float(value)})
    fig.suptitle("PACS Trial 1B - sample median protein intensities by group")
    save(fig, view, "03_group_median_boxplots")
    return rows


def run_pca(raw, corrected, samples, labels, view, chosen):
    """Fit a separate PCA per view, using identical features pre/post within view.

    Features detected in >=70% of view samples are selected from original data.
    Missing log2 intensities are filled with per-feature median for PCA ONLY.
    PCA is feature-mean-centered, no feature variance scaling.
    """
    selected = [s for g in chosen for s in samples if labels[s] == g]
    ix = [samples.index(s) for s in selected]
    pre_raw, post_raw = raw[:, ix], corrected[:, ix]
    chosen_features = finite_positive(pre_raw).sum(axis=1) >= math.ceil(len(ix)*MIN_PCA_DETECTION)
    if chosen_features.sum() < 100:
        raise ValueError(f"Not enough PCA proteins in {view}")
    transformed = []
    for mat in [pre_raw[chosen_features], post_raw[chosen_features]]:
        log = log2_matrix(mat)
        row_med = np.nanmedian(log, axis=1)
        filled = np.where(np.isfinite(log), log, row_med[:, None])
        transformed.append(filled)
    variable = (np.std(transformed[0], axis=1) > 1e-10) & (
        np.std(transformed[1], axis=1) > 1e-10)
    count = int(variable.sum())
    if count < 100:
        raise ValueError("Not enough variable proteins for PCA")
    scores = []
    explained = {}
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.4), layout="constrained")
    for ax, stage, arr in zip(axes, STAGE_NAMES, transformed):
        pca = PCA(n_components=2, svd_solver="full")
        xy = pca.fit_transform(arr[variable].T)
        explained[stage] = pca.explained_variance_ratio_*100
        for sid, point in zip(selected, xy):
            scores.append({
                "view": view, "sample_id": sid, "group": labels[sid], "stage": stage,
                "PC1": float(point[0]), "PC2": float(point[1]),
                "PC1_variance_pct": float(explained[stage][0]),
                "PC2_variance_pct": float(explained[stage][1]),
            })
        for g in chosen:
            p = np.array([pt for sid, pt in zip(selected, xy) if labels[sid] == g])
            ax.scatter(p[:, 0], p[:, 1], s=46, color=PALETTE[g],
                       edgecolors="white", linewidths=.55, alpha=.85,
                       label=f"{g} (n={len(p)})")
        ax.set_title(stage)
        ax.set_xlabel(f"PC1 ({explained[stage][0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({explained[stage][1]:.1f}%)")
        ax.axhline(0, color="#b9c0c8", lw=.6)
        ax.axvline(0, color="#b9c0c8", lw=.6)
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle(f"PACS Trial 1B PCA - {view.replace('_', ' ')} ({count} proteins)")
    save(fig, view, "04_pca")
    return scores, {"view": view, "samples": len(selected), "pca_proteins": count,
                    "pc1_before_pct": explained["Before"][0],
                    "pc2_before_pct": explained["Before"][1],
                    "pc1_after_pct": explained["After"][0],
                    "pc2_after_pct": explained["After"][1]}


def plot_cv(cv_arrays, summary, view, chosen):
    """Reference R style violin + box + individual faint points with median labels."""
    fig, axes = plt.subplots(1, 2, figsize=(max(12, 3.5*len(chosen)), 6.5),
                             sharey=True, layout="constrained")
    for ax, stage, value_ix in zip(axes, STAGE_NAMES, (0, 1)):
        for i, group in enumerate(chosen, 1):
            values = cv_arrays[group][value_ix]
            plotted = np.minimum(values, CV_PLOT_CAP)
            if len(plotted) > 1 and np.ptp(plotted) > 0:
                violin = ax.violinplot([plotted], positions=[i], widths=.75,
                                       showextrema=False)
                for body in violin["bodies"]:
                    body.set_facecolor(PALETTE[group])
                    body.set_alpha(.4)
            box = ax.boxplot([plotted], positions=[i], widths=.16,
                             patch_artist=True, showfliers=False, manage_ticks=False)
            for patch in box["boxes"]:
                patch.set_facecolor(PALETTE[group])
                patch.set_alpha(.6)
            rng = np.random.default_rng(20261010 + i)
            jitter = rng.uniform(-.22, .22, len(plotted))
            ax.scatter(i+jitter, plotted, color=PALETTE[group], s=3,
                       alpha=.055, linewidths=0, rasterized=True, zorder=3)
            ax.text(i, CV_PLOT_CAP+7, f"Median={np.median(values):.2f}%",
                    ha="center", fontsize=9)
        ax.set_xlim(.4, len(chosen)+.6)
        ax.set_ylim(0, 222)
        ax.set_xticks(range(1, len(chosen)+1), chosen)
        ax.set_ylabel("Per-protein CV of raw intensities (%)")
        ax.set_title(stage)
        ax.grid(axis="y", alpha=.2)
    fig.suptitle(f"PACS Trial 1B - {view.replace('_',' ')} protein CV distribution")
    axes[1].text(.99, .01, "CV >200% capped for display only",
                 transform=axes[1].transAxes, ha="right", fontsize=8)
    save(fig, view, "05_cv_violin_box")


def write_report(ref, factors, summary, pca_summary, n_proteins):
    lines = [
        "# PACS Trial 1B - Internal-QC remains unchanged", "",
        "This exploratory test uses the **same division direction** as Trial 1A,",
        "but does NOT correct any Internal-QC (LM) sample. Trial 1A and source files are unchanged.",
        "All generated figures are **PNG only; no SVG**.", "",
        "## Definition", "",
        "- For each of the 75 samples, sample_mR = median of valid positive raw L/H ratios for 97 HELP peptides.",
        "- ref_mR = median of the **15 pooled QC sample_mR values**, not a pooled median of individual peptide ratios.",
        "- candidate_factor = sample_mR / ref_mR.",
        "- For pooled QC, HC and S: corrected intensity = original / candidate_factor.",
        "- For Internal-QC: corrected intensity = **original exactly** (applied multiplier = 1).",
        "- LM diagnostic sample_mR and candidate_factor are included in the factor table but NOT applied.",
        f"- ref_mR = **{ref:.12g}**; corrected matrix includes {n_proteins} protein-group rows × 75 samples.",
        "- n: Internal-QC 6; pooled QC 15; HC 27; S 27.",
        "- The intentional LM_2_1/LM_2_2 and LM_3_2/LM_3_3 duplicates remain as supplied.",
        "- No original pg_matrix.tsv, HELP workbook, metadata or Trial 1A result was edited.",
        "", "## Two visualization versions", "",
        "**All groups**: Internal-QC, pooled QC, HC and S (75 samples).",
        "The LM panel/points are identical before and after correction.",
        "**Study only**: HC and S (54 samples), excluding both Internal-QC and pooled QC.",
        "PCA is re-fitted separately for each view; in each view, identical protein",
        "features are used before and after. PCA axes are not directly comparable across views.",
        "", "### All four groups", "",
    ]
    for name, title in [
        ("01_applied_multipliers", "Applied multipliers"),
        ("02_per_sample_boxplots", "Per-sample protein intensity boxplots"),
        ("03_group_median_boxplots", "Group-level boxes of sample median log2 intensities"),
        ("04_pca", "PCA before and after"),
        ("05_cv_violin_box", "Per-protein CV violin + box + jitter"),
    ]:
        lines += [f"![{title}](figures/all_groups/{name}.png)", ""]
    lines += ["### HC and S only", ""]
    for name, title in [
        ("01_applied_multipliers", "Applied multipliers"),
        ("02_per_sample_boxplots", "Per-sample protein intensity boxplots"),
        ("03_group_median_boxplots", "Group-level boxes of sample median log2 intensities"),
        ("04_pca", "PCA before and after"),
        ("05_cv_violin_box", "Per-protein CV violin + box + jitter"),
    ]:
        lines += [f"![{title}](figures/study_only/{name}.png)", ""]
    lines += [
        "## Per-protein CV summary", "",
        "CV = sample SD (ddof=1) / mean of positive finite raw protein intensities * 100%.",
        "Only proteins with at least max(3, ceil(70% of group n)) detected",
        "and mean >= 1e-8 have a reported CV; all other CVs are missing, not zero.",
        "There is no CV imputation or log transform; HC/S variation includes biological variation.",
        "", "| Group | n | CV-eligible proteins | Median before | Median after | Median paired change (pp) |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for item in summary.to_dict("records"):
        lines.append(f'| {item["group"]} | {item["sample_n"]} | '
                     f'{item["eligible_proteins_paired"]} | '
                     f'{item["median_cv_before_pct"]:.2f}% | '
                     f'{item["median_cv_after_pct"]:.2f}% | '
                     f'{item["median_paired_change_pp"]:+.2f} |')
    lines += [
        "", "## PCA methodology", "",
        "- For each view, select protein groups detected in >=70% of samples in that view.",
        "- Retain an identical selected protein feature set across before/after for that view.",
        "- PCA uses log2 positive protein intensities with per-protein median imputation",
        "  of missing log2 values **only for PCA**, separately for before/after.",
        "- Remove invariant protein features, center columns of sample × protein input,",
        "  and do NOT scale individual proteins to unit variance.",
        "", "| PCA view | Samples | Proteins | PC1 before | PC2 before | PC1 after | PC2 after |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for item in pca_summary:
        lines.append(
            f'| {item["view"]} | {item["samples"]} | {item["pca_proteins"]} | '
            f'{item["pc1_before_pct"]:.2f}% | {item["pc2_before_pct"]:.2f}% | '
            f'{item["pc1_after_pct"]:.2f}% | {item["pc2_after_pct"]:.2f}% |'
        )
    lines += [
        "", "## Outputs", "",
        "- pg_matrix_corrected.tsv: corrected matrix, LM columns byte-text preserved from the source table.",
        "- normalization_factors.csv: sample medians, candidate factors, applied factors and multipliers.",
        "- per_protein_cv.csv: protein × group CV before and after.",
        "- group_cv_summary.csv: CV statistics for four groups (study-only uses same HC/S rows).",
        "- pca_scores_all_groups.csv and pca_scores_study_only.csv: separate PCA scores.",
        "- pca_explained_variance.csv: protein feature and variance information.",
        "- sample_log2_medians.csv: per-sample box-plot aggregation data for both views.",
        "- figures/all_groups/: 5 PNG figures with four groups.",
        "- figures/study_only/: 5 PNG figures with HC/S only.",
        "- scripts/pacs_trial1b.py and .github/workflows/pacs-trial1b.yml: reproducibility.",
        "", "## Important caveat", "",
        "This is a requested exploratory formula, NOT validation of its suitability.",
        "The HELP L/H ratio cannot alone prove that division is the correct correction direction.",
        "No selection or optimization of factor direction has been performed for Trial 1B.", "",
    ]
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    labels = metadata()
    pg = pd.read_csv(PG, sep="\t", encoding="utf-8-sig", dtype=str, keep_default_na=False)
    if list(pg.columns[:5]) != ANN or len(pg.columns) != 80:
        raise ValueError("Unexpected pg_matrix.tsv columns")
    samples = list(pg.columns[5:])
    if set(samples) != set(labels) or len(samples) != 75:
        raise ValueError("Metadata and pg_matrix.tsv sample columns disagree")
    by_group = {g: [s for s in samples if labels[s] == g] for g in GROUPS}
    if {g: len(v) for g, v in by_group.items()} != {
        "Internal-QC": 6, "pooled QC": 15, "HC": 27, "S": 27
    }:
        raise ValueError("Wrong sample group sizes")
    raw = pg[samples].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(raw).any():
        raise ValueError("No numeric protein intensities")
    factors, reference = factor_table(labels, samples)
    multipliers = factors["applied_multiplier"].to_numpy()
    corrected = raw * multipliers[None, :]
    lm_ix = [i for i, s in enumerate(samples) if labels[s] == "Internal-QC"]
    if not np.array_equal(raw[:, lm_ix], corrected[:, lm_ix], equal_nan=True):
        raise AssertionError("LM values were changed")
    if not np.array_equal(finite_positive(raw), finite_positive(corrected)):
        raise AssertionError("Protein detection mask was changed")
    # Preserve original LM text values, including full numeric formatting.
    output = pg.copy()
    for j, sid in enumerate(samples):
        if labels[sid] == "Internal-QC":
            assert output[sid].equals(pg[sid])
            continue
        original = pg[sid].to_numpy()
        output[sid] = [format(v, ".16g") if math.isfinite(v) else original[i]
                       for i, v in enumerate(corrected[:, j])]
    for sid in by_group["Internal-QC"]:
        assert output[sid].equals(pg[sid])
    output.to_csv(OUT / "pg_matrix_corrected.tsv", sep="\t", index=False)
    summary, cvs = protein_cv(raw, corrected, pg[ANN], samples, labels)
    all_pca, pca_overview, sample_medians = {}, [], []
    for view, chosen in VIEWS.items():
        plot_factors(factors, view, chosen)
        plot_per_sample_boxes(raw, corrected, samples, labels, view, chosen)
        sample_medians.extend(plot_group_medians(raw, corrected, samples, labels,
                                                 view, chosen))
        scores, pca_info = run_pca(raw, corrected, samples, labels, view, chosen)
        all_pca[view] = scores
        pca_overview.append(pca_info)
        plot_cv(cvs, summary, view, chosen)
    for view, scores in all_pca.items():
        pd.DataFrame(scores).to_csv(OUT / f"pca_scores_{view}.csv", index=False)
    pd.DataFrame(pca_overview).to_csv(OUT / "pca_explained_variance.csv", index=False)
    pd.DataFrame(sample_medians).to_csv(OUT / "sample_log2_medians.csv", index=False)
    write_report(reference, factors, summary, pca_overview, len(pg))
    print("SUCCESS PACS Trial 1B: LM unchanged, QC/HC/S divided by HELP factor.")
    print(f"Protein groups={len(pg)}, samples=75, LM=6, QC=15, HC=27, S=27, ref_mR={reference:.12g}")
    print(summary.to_string(index=False))
    print(pd.DataFrame(pca_overview).to_string(index=False))


if __name__ == "__main__":
    main()
