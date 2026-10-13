#!/usr/bin/env python3
"""Trial 3: EN01.xlsx yeast protein factor1 replaces Trial 2 yeast peptide factor1.

Two independent strategies:
  single_max_rowsum: select the protein row with the largest SUM of original
      positive linear intensities across all 75 metadata-retained samples.
      factor1_s = median_54Study(log2 selected protein) - log2 selected protein_s.
  sum_two: for each sample SUM the two protein linear intensities, then log2.
      factor1_s = median_54Study(log2(row1+row2)) - log2(row1+row2)_s.

Then optionally add original Trial 2 factor2:
    median_15PooledQC(median_97HELP(log2(L/H))) - median_97HELP(log2(L/H))_s.

Correct HC and S ONLY; leave both QC groups unchanged to exact raw linear
values. All non-CV computations are in log2 space. CV uses original linear
intensities with ddof=1. Reuse Trial 2's diagnostic and 3-stage chart functions
so the ONLY normalization-method change is the definition of factor1.
"""
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pacs_trial2_stepwise as base
from filter_samples import readbook, tab
from pacs_trial1 import ROOT, PG
from yeast_protein_qc import GROUPS

OUT = ROOT / "results"
EN01 = ROOT / "EN01.xlsx"
HELP = ROOT / "LH_ratio_tables_2026.10.13—ALL-HELP.xlsx"
METHODS = ("single_max_rowsum", "sum_two")
STUDY = {"HC", "S"}
EXPECTED_ROWS = ("P00924_ENO1_YEAST", "P00925_ENO1_YEAST")

# Override only this run's variant labels and output readers; old results
# and scripts remain committed as historical references.
base.METHODS = METHODS
base.OUT = OUT
base.old.OUT = OUT
base.old.HELP = HELP


def read_en01(samples):
    strings, sheets = readbook(EN01)
    if len(sheets) != 1:
        raise ValueError(f"Expected one EN01 worksheet, saw {len(sheets)}")
    rows = tab(sheets[0], strings)
    if len(rows) != 3:
        raise ValueError(f"Expected 2 yeast protein rows, saw {len(rows)-1}")
    ids = tuple(row.get("A", "").strip() for row in rows[1:])
    if ids != EXPECTED_ROWS:
        raise ValueError(f"EN01 protein identifiers changed: {ids}")
    columns = {name: col for col, name in rows[0].items() if col != "A"}
    if len(columns) != 75 or set(columns) != set(samples):
        raise ValueError("EN01 sample columns must match exactly the 75 metadata samples")
    arr = np.empty((2, len(samples)), dtype=float)
    for i, row in enumerate(rows[1:]):
        for j, sample in enumerate(samples):
            value = row.get(columns[sample], "")
            try:
                arr[i, j] = float(value)
            except (ValueError, TypeError):
                arr[i, j] = np.nan
    if not np.all(np.isfinite(arr) & (arr > 0)):
        invalid = [(ids[i], samples[j], arr[i, j])
                   for i, j in zip(*np.where(~(np.isfinite(arr) & (arr > 0))))]
        raise ValueError(f"EN01 requires two positive measured intensities in every sample: {invalid[:30]}")
    return ids, arr


def compute_protein_factors(samples, labels):
    if not EN01.is_file() or not HELP.is_file():
        raise FileNotFoundError("EN01.xlsx or corrected 2026.10.13 HELP workbook missing")
    ids, arr = read_en01(samples)
    study = np.array([labels[s] in STUDY for s in samples], dtype=bool)
    if study.sum() != 54:
        raise ValueError(f"Expected 54 Study samples, got {study.sum()}")
    row_sums = arr.sum(axis=1)  # unlogged, across all 75 included samples
    if np.isclose(row_sums[0], row_sums[1], rtol=1e-12):
        raise ValueError("Cannot choose a strictly larger EN01 row sum")
    win = int(np.argmax(row_sums))
    single_linear = arr[win]
    combined_linear = arr[0] + arr[1]   # sum BEFORE log2; NOT mean of log2s
    signals = {
        "single_max_rowsum": np.log2(single_linear),
        "sum_two": np.log2(combined_linear),
    }
    ref = {method: float(np.median(logged[study]))
           for method, logged in signals.items()}
    raw_factors = {method: ref[method] - logged
                   for method, logged in signals.items()}
    f2, help_medians, counts, help_ref = base.old.compute_help_factor2(samples, labels)

    rows = []
    for j, sid in enumerate(samples):
        rows.append({
            "sample_id": sid, "group": labels[sid],
            "P00924_linear_intensity": arr[0, j],
            "P00925_linear_intensity": arr[1, j],
            "two_protein_linear_sum": combined_linear[j],
            "selected_protein": ids[win],
            "single_selected_linear_intensity": single_linear[j],
            "single_selected_log2": signals["single_max_rowsum"][j],
            "sum_two_log2": signals["sum_two"][j],
            "candidate_factor1_single": raw_factors["single_max_rowsum"][j],
            "candidate_factor1_sum_two": raw_factors["sum_two"][j],
            "applied_factor1_single": raw_factors["single_max_rowsum"][j] if study[j] else 0.0,
            "applied_factor1_sum_two": raw_factors["sum_two"][j] if study[j] else 0.0,
            "candidate_factor2": f2[sid],
            "applied_factor2": f2[sid] if study[j] else 0.0,
            "sample_HELP_mR_log2": help_medians[sid],
            "reference_HELP_mR_log2": help_ref,
            "n_valid_HELP_ratios": counts[sid],
        })
    base.save_csv(OUT / "sample_factors.csv", rows)
    base.save_csv(OUT / "eno1_protein_details.csv", rows)
    base.save_csv(OUT / "eno1_rowsum_ranking.csv", [{
        "protein_id": ids[i],
        "row_sum_linear_all75": float(row_sums[i]),
        "row_sum_linear_study54": float(arr[i, study].sum()),
        "nonmissing_positive_all75": int((np.isfinite(arr[i]) & (arr[i] > 0)).sum()),
        "used_for_single_protein_correction": bool(i == win)
    } for i in range(2)])

    factors = {
        sid: {
            method: float(raw_factors[method][j]) if study[j] else 0.0
            for method in METHODS
        } | {"factor2": float(f2[sid]) if study[j] else 0.0}
        for j, sid in enumerate(samples)
    }
    return factors, {
        "selected_protein": ids[win],
        "row_ids": ids,
        "row_sums": row_sums.tolist(),
        "references": ref,
        "ref_HELP_log2": help_ref,
    }


def offset_plot(samples, labels, groups, factors, method, folder):
    selected, _ = base.select(samples, labels, groups)
    x = np.arange(len(selected))
    f1 = np.array([factors[s][method] for s in selected])
    f2 = np.array([factors[s]["factor2"] for s in selected])
    fig, axes = plt.subplots(2, 1,
                            figsize=(max(12, len(selected) * .23), 8.2),
                            sharex=True, layout="constrained")
    axes[0].bar(x-.2, f1, width=.38, label="ENO1 protein factor1")
    axes[0].bar(x+.2, f2, width=.38, label="97-HELP factor2")
    axes[0].legend(frameon=False)
    axes[0].axhline(0, c="black", lw=.8)
    axes[0].set_ylabel("Individual log2 offsets")
    axes[1].bar(x, f1 + f2, color=[base.PALETTE[labels[s]] for s in selected])
    axes[1].axhline(0, c="black", lw=.8)
    axes[1].set_ylabel("Combined log2 offset")
    axes[1].set_xticks(x, selected, rotation=90, fontsize=6)
    fig.suptitle(f"PACS yeast protein replacement trial: {method} (both QC groups uncorrected)")
    base.save_fig(fig, folder / "06_factor_components.png")


def create_report(audit, refs, cvrows, pcarows):
    lines = [
        "# PACS new attempt: EN01 yeast protein replaces yeast peptide factor1", "",
        f"Archived previous results/ as **{os.environ.get('PACS_RESULTS_ARCHIVE', 'results_vN')}/**.",
        "This new results/ contains two yeast-protein-based factor1 strategies;",
        "the former Trial 2 factor2 and other calculation/plotting rules are unchanged.", "",
        "## Source data", "",
        f"- Protein matrix: pg_matrix.tsv, SHA256 {hashlib.sha256(PG.read_bytes()).hexdigest()}.",
        f"- Yeast protein input: EN01.xlsx, SHA256 {hashlib.sha256(EN01.read_bytes()).hexdigest()}.",
        f"- HELP workbook: {HELP.name}, SHA256 {hashlib.sha256(HELP.read_bytes()).hexdigest()}.",
        "- The previous yeast **peptide** workbook is NOT used to calculate either factor1.",
        f"- Protein rows: {audit['input_protein_rows']} (yeast P00924 and P00925 removed from main matrix).",
        f"- Samples: {audit['retained_sample_columns']} (LM 6, pooled QC 15, HC 27, S 27).",
        "- No input or older result files are modified by calculation.", "",
        "## Rank two yeast proteins by row sum", "",
        "Ranking uses the ordinary sum of original **linear** intensities for",
        "**all 75 retained samples**, not the sum of log2 intensities. Both rows",
        "are present and strictly positive in all 75 samples.", "",
        "| Source row identifier (verbatim) | Linear row sum (75 samples) | Chosen for single method |",
        "| --- | ---: | --- |",
    ]
    for ident, sm in zip(refs["row_ids"], refs["row_sums"]):
        lines.append(
            f"| {ident} | {sm:.8g} | {'Yes' if ident == refs['selected_protein'] else 'No'} |"
        )
    lines += [
        "", f"Selected larger-row-sum row: **{refs['selected_protein']}**.",
        "The two source rows are P00924_ENO1_YEAST and P00925_ENO1_YEAST;",
        "the source names are kept as uploaded, including the P00925 naming discrepancy",
        "documented in ENO1_DATA_CLEANUP_REPORT.md.", "",
        "## Corrections in log2", "",
        "Let x1_s and x2_s be the two original positive **linear** EN01 protein intensities.",
        "All yeast references are the median across **54 Study samples** (HC and S).", "",
        "**Variant 1: single_max_rowsum**",
        f"- Select the row with larger 75-sample row sum: {refs['selected_protein']}.",
        "- y_s = log2(x_selected,s).",
        "- factor1_s = median_Study(y) - y_s.", "",
        "**Variant 2: sum_two**",
        "- y_s = log2(x1_s + x2_s); sum linear intensities BEFORE log2.",
        "- factor1_s = median_Study(y) - y_s.", "",
        "**Unchanged Trial 2 HELP correction**",
        "- mR_s = median_97HELP(log2(L/H) for valid positive ratios).",
        "- ref_mR = median_15_pooled_QC(mR_s).",
        "- factor2_s = ref_mR - mR_s.", "",
        "Stages within each method:",
        "- Raw: original protein intensity.",
        "- Factor1: log2 corrected protein = log2(raw) + factor1.",
        "- Factor1+Factor2: log2 corrected protein = log2(raw) + factor1 + factor2.",
        "- **Only HC and S are corrected.** Internal-QC and pooled QC raw linear data are exactly unchanged.",
        "- Factor1 and factor2 are added in log2; CV is sample SD(ddof=1)/arithmetic mean",
        "  of original-scale positive intensities, times 100, with >=70% group detection, >=3 and mean >=1e-8.",
        "- PCA and boxplots are log2; PCA imputes only for projection, uses >=70% observed",
        "  common proteins per view; PCA refitted at each stage, so PC coordinates are not directly aligned.", "",
        "## Computed references", "",
        f"- Larger-row-sum protein Study log2 median: {refs['references']['single_max_rowsum']:.12g}.",
        f"- Two-protein linear sum Study log2 median: {refs['references']['sum_two']:.12g}.",
        f"- Pooled QC HELP ref_mR log2: {refs['ref_HELP_log2']:.12g}.", "",
        "## Median per-protein CV (%)", "",
        "| Method | Stage | Internal-QC | pooled QC | HC | S |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    cv = pd.DataFrame(cvrows)
    titles = dict(zip(base.STAGES, base.TITLES))
    for method in METHODS:
        for stage in base.STAGES:
            values = [float(cv[(cv.method==method)&(cv.stage==stage)&(cv.group==g)]
                           .median_cv_pct.iloc[0]) for g in GROUPS]
            lines.append("| " + method + " | " + titles[stage] + " | " +
                         " | ".join(f"{v:.2f}%" for v in values) + " |")
    lines += [
        "", "## Visualizations: all three stages together", "",
        "2 methods × 2 views (all_groups / study_only) × 6 PNG = **24 PNG**; no SVG.",
        "Each plot places Raw / Factor1 / Factor1+Factor2 next to each other.", "",
    ]
    chartfiles = [
        ("01_sample_boxplots_3stages.png", "Sample protein intensity boxplots"),
        ("02_group_boxes_3stages.png", "Per-sample median intensity by group"),
        ("03_pca_3stages.png", "PCA across 3 stages"),
        ("04_cv_3stages.png", "Protein-level CV distribution"),
        ("05_cv_median_comparison.png", "Median protein CV by stage"),
        ("06_factor_components.png", "Factor1 and factor2 component offsets"),
    ]
    for method in METHODS:
        for view in base.VIEWS:
            lines += [f"### {method} / {view}", ""]
            for name, desc in chartfiles:
                lines += [f"![{desc}](method_{method}/figures/{view}/{name})", ""]
    lines += [
        "## PCA explained variance", "",
        "| Method | View | Stage | Proteins | PC1 | PC2 |",
        "| --- | --- | --- | ---: | ---: | ---: |",
    ]
    for row in pcarows:
        lines.append(f"| {row['method']} | {row['view']} | {row['stage']} | "
                     f"{row['n_pca_proteins']} | {row['PC1_pct']:.2f}% | {row['PC2_pct']:.2f}% |")
    lines += [
        "", "## Data files", "",
        "- eno1_rowsum_ranking.csv: both original linear protein row sums and selected row.",
        "- eno1_protein_details.csv: both raw protein intensities, two-protein sum, references and factors by sample.",
        "- sample_factors.csv: both factor1 candidates and applied factors plus unchanged HELP factor2.",
        "- pg_matrix_filtered.tsv, screening_report.json: unchanged source protein screening.",
        "- all_stage_cv_summary.csv, all_stage_pca_summary.csv: comparable method/stage summaries.",
        "- method_single_max_rowsum/ and method_sum_two/: corrected linear/log2 matrices",
        "  at Factor1 and Factor1+Factor2 stages, per-protein CV, PCA data, offsets",
        "  and six joint-stage PNGs for all_groups and study_only.", "",
        "## Interpretation caveats", "",
        "An improvement in HC/S CV alone is insufficient to prove batch-effect removal:",
        "biological variability contributes to study CV. Neither QC group is corrected",
        "in the main outputs, so their unchanged CV cannot serve as independent validation.",
        "This is a method comparison, not a claim of successful batch normalization.", "",
    ]
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")


def run_comparison():
    """Same stage-by-stage pipeline as Trial 2, without the obsolete assumption
    that the two yeast strategies differ only by a global constant."""
    OUT.mkdir(parents=True, exist_ok=True)
    df, labels, samples, raw, audit = base.prepare_input()
    if len(df) != 6675:
        raise ValueError(f"Expected 6675 endogenous protein groups, got {len(df)}")
    factors, refs = compute_protein_factors(samples, labels)
    cv_rows, pca_rows = [], []
    for method in METHODS:
        folder = OUT / ("method_" + method)
        folder.mkdir(parents=True, exist_ok=True)
        mats, offsets = base.correction_stages(raw, samples, labels, factors, method)
        for stage in base.STAGES[1:]:
            base.write_matrix(folder / ("pg_matrix_" + stage + "_linear.tsv"),
                              df, samples, labels, mats[stage])
            base.write_matrix(folder / ("pg_matrix_" + stage + "_log2.tsv"),
                              df, samples, labels, mats[stage], log2_scale=True)
        base.save_csv(folder / "stage_offsets.csv", [
            {"sample_id": sid, "group": labels[sid], "stage": stage,
             "applied_log2_offset": float(offset),
             "applied_linear_multiplier": float(np.exp2(offset))}
            for stage in base.STAGES for sid, offset in zip(samples, offsets[stage])
        ])
        summaries, per_protein = base.calculate_cv(mats, df, samples, labels, method, folder)
        cv_rows.extend(summaries)
        pca_scores, pca_summaries, medians = [], [], []
        for view, groups in base.VIEWS.items():
            figures = folder / "figures" / view
            base.plot_sample_boxes(mats, samples, labels, groups, figures)
            medians.extend(base.plot_group_medians(mats, samples, labels, groups, view, figures))
            scores, pca_info = base.pca_three(mats, samples, labels, groups, view, figures)
            pca_scores.extend(scores)
            pca_summaries.extend(pca_info)
            pca_rows.extend([{"method":method, **record} for record in pca_info])
            base.plot_cv(per_protein, groups, figures)
            base.plot_cv_medians(summaries, groups, figures)
            offset_plot(samples, labels, groups, factors, method, figures)
        base.save_csv(folder / "pca_scores.csv", pca_scores)
        base.save_csv(folder / "pca_summary.csv", pca_summaries)
        base.save_csv(folder / "sample_log2_medians.csv", medians)
    base.save_csv(OUT / "all_stage_cv_summary.csv", cv_rows)
    base.save_csv(OUT / "all_stage_pca_summary.csv", pca_rows)
    create_report(audit, refs, cv_rows, pca_rows)


def main():
    run_comparison()
    print("TRIAL3 SUCCESS: EN01 protein single-max-row-sum and linear-two-row-sum, Trial2 factor2 unchanged")


if __name__ == "__main__":
    main()
