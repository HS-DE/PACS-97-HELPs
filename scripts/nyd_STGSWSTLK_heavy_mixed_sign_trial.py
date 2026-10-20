#!/usr/bin/env python3
"""NYD 39 samples: STGSWSTLK Heavy mixed-sign calibration.

For both global and batch-local Heavy reference, f_i = median_ref - log2(H_i).
Apply -f_i to raw log2 PG in T001/T002/T003 and +f_i in T004.
Keep the original input matrices, missingness, sample IDs and all old trials.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import nyd_STGSWSTLK_heavy_exclude_T003_outlier as prev
import nyd_STGSWSTLK_heavy_batch_internal_trial as batch
import nyd_STGSWSTLK_heavy_reverse_factor_trial as reverse
import nyd_client_boxplots_correlation_heatmaps as client
import nyd_STGSWSTLK_subtract_factor_CV_violin as violins
from nyd_pg_batch_boxplots import BATCHES

OUT = prev.BASE / "STGSWSTLK_Heavy_mixed_sign_39samples"
VARIANTS = ("overall_reference", "within_batch_reference")
SIGN = {"T001": -1., "T002": -1., "T003": -1., "T004": +1.}


def csv_write(path, records):
    if not records:
        raise ValueError(f"Empty CSV table: {path}")
    with path.open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def draw_shift(samples, records, shift, folder, title):
    fig, ax = plt.subplots(figsize=(17.5, 5.9), layout="constrained")
    color = {"T001": "#6B9CBF", "T002": "#67A68D",
             "T003": "#B38CC3", "T004": "#D2A36A"}
    xx = np.arange(1, len(samples)+1)
    for b in BATCHES:
        idx = np.array([row["batch"] == b for row in records])
        ax.bar(xx[idx], shift[idx], color=color[b], width=.75,
               alpha=.8, label=f"{b} ({'-fi' if SIGN[b] < 0 else '+fi'})")
    for cut in (16.5, 24.5, 31.5):
        ax.axvline(cut, color="#71879A", linewidth=.9)
    ax.axhline(0, color="#384E62", linewidth=1)
    ax.set_xticks(xx, [s[-6:] for s in samples], rotation=90, fontsize=7.3)
    ax.set_ylabel("Applied log2 shift to every valid PG")
    ax.set_title(f"{title}: T001-T003 -fi, T004 +fi", fontweight="bold")
    ax.legend(frameon=False, ncol=4, fontsize=8, loc="upper right")
    ax.grid(axis="y", alpha=.14)
    reverse.savefig(fig, folder/"04_sample_applied_log2_shifts.png")


def cv_report_plot(compare):
    x = np.arange(4)
    raw = [next(row["raw_median_cv_pct"] for row in compare if row["batch"] == b)
           for b in BATCHES]
    corrected = [next(row["overall_mixed_median_cv_pct"] for row in compare if row["batch"] == b)
                 for b in BATCHES]
    fig, ax = plt.subplots(figsize=(10.5, 5.5), layout="constrained")
    ax.bar(x-.2, raw, .38, color="#8798A5", label="Raw")
    ax.bar(x+.2, corrected, .38, color="#347CA5", label="Mixed-sign corrected")
    for pos, cv in zip(list(x-.2)+list(x+.2), raw+corrected):
        ax.text(pos, cv+.32, f"{cv:.2f}%", ha="center", fontsize=9)
    ax.set_ylim(0, max(raw+corrected)*1.24)
    ax.set_xticks(x, [f"{b}\n{'-fi' if SIGN[b]<0 else '+fi'}" for b in BATCHES])
    ax.set_ylabel("Within-batch median per-protein CV (%)")
    ax.set_title("NYD STGSWSTLK Heavy | mixed directions vs Raw", fontweight="bold")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=.14)
    ax.set_axisbelow(True)
    reverse.savefig(fig, OUT/"00_mixed_sign_CV_vs_Raw.png")


def main():
    if OUT.exists():
        raise FileExistsError("Independent output already exists: " + str(OUT))

    records, names, annotations, raw, targets, heavy, _ = prev.inputs()
    if len(names) != 39 or prev.EXCLUDED_ID in names:
        raise AssertionError("Expected exactly 39 retained samples")
    if {b: sum(r["batch"] == b for r in records) for b in BATCHES} != prev.BATCH_COUNTS:
        raise AssertionError("Batch sizes changed")

    h = np.log2(heavy)
    ref_global = float(np.median(h))
    h_again, idx, refs_batch, f_batch, _ = batch.factors(records, heavy, raw)
    assert np.allclose(h, h_again, atol=0, rtol=0)
    f_global = ref_global - h

    # All factors must match the previously validated +fi trials exactly.
    plus_global = pd.read_csv(
        reverse.PLUS_GLOBAL/"sample_STGSWSTLK_Heavy_factors.csv"
    ).set_index("sample_id").loc[names]
    plus_batch = pd.read_csv(
        reverse.PLUS_BATCH/"sample_batch_internal_factors.csv"
    ).set_index("sample_id").loc[names]
    assert np.allclose(f_global, plus_global.fi_log2.to_numpy(float), atol=1e-11, rtol=0)
    assert np.allclose(f_batch, plus_batch.fi_log2.to_numpy(float), atol=1e-11, rtol=0)

    if not np.isclose(ref_global, 21.085857242488732, atol=1e-8, rtol=0):
        raise AssertionError("Global Heavy reference differs from previous trial")
    known = {"T001":21.122744,"T002":21.163281,"T003":20.918815,"T004":21.230995}
    for b in BATCHES:
        if not np.isclose(refs_batch[b], known[b], atol=1e-6, rtol=0):
            raise AssertionError(f"Batch reference changed: {b}: {refs_batch[b]}")

    sign = np.array([SIGN[r["batch"]] for r in records])
    factors = {"overall_reference": f_global,
               "within_batch_reference": f_batch}
    matrices = {variant: raw+(sign*f)[None, :]
                for variant, f in factors.items()}

    # Verify the per-sample algebra and the batch-wise constant-offset invariant.
    for variant, m in matrices.items():
        if not np.array_equal(np.isfinite(m), np.isfinite(raw)):
            raise AssertionError(f"Missingness changed in {variant}")
        if not np.allclose(m-(sign*factors[variant])[None,:], raw,
                           atol=1e-11, rtol=0, equal_nan=True):
            raise AssertionError(f"Incorrect per-sample correction in {variant}")
    for b, sel in idx.items():
        observed=(sign*(f_global-f_batch))[sel]
        expected=SIGN[b]*(ref_global-refs_batch[b])
        if not np.allclose(observed, expected, atol=1e-12, rtol=0):
            raise AssertionError(f"Within-batch constant offset failed for {b}")

    # Same Raw PCA fit and plot bounds for both versions.
    pca_results, pca_bounds = reverse.common_raw_fit_pca(raw, matrices, names, records)
    sample_short = client.sample_order(records, names)
    box_bounds = client.box_extent([raw, *matrices.values()])
    OUT.mkdir(parents=True, exist_ok=False)
    reverse.pg_matrix(OUT/"00_Raw_NYD_pg_matrix2_log2_39samples.tsv",
                      annotations, names, raw)
    cvs = {}
    for variant in VARIANTS:
        folder = OUT / variant
        folder.mkdir(exist_ok=False)
        applied = sign*factors[variant]
        corrected = matrices[variant]
        refs = np.full(39, ref_global) if variant == "overall_reference" else np.array(
            [refs_batch[r["batch"]] for r in records])

        audit_rows = []
        for j, sid in enumerate(names):
            audit_rows.append({
                "sample_id": sid, "batch": records[j]["batch"],
                "correction_direction": "-fi" if sign[j] < 0 else "+fi",
                "STGSWSTLK_Heavy_linear": float(heavy[j]),
                "STGSWSTLK_Heavy_log2": float(h[j]),
                "reference_median_log2": float(refs[j]),
                "fi_reference_minus_Heavy_log2": float(factors[variant][j]),
                "actual_PG_applied_shift_log2": float(applied[j]),
                "actual_PG_multiplier_linear": float(np.exp2(applied[j])),
                "Raw_PG_sample_median_log2": float(np.nanmedian(raw[:,j])),
                "Corrected_PG_sample_median_log2": float(np.nanmedian(corrected[:,j])),
            })
        csv_write(folder/"sample_factors_internal.csv", audit_rows)
        reverse.pg_matrix(folder/"corrected_NYD_pg_matrix2_log2.tsv",
                          annotations, names, corrected)
        summary, long_cv, _ = prev.stats_cv(raw, corrected, records, annotations)
        csv_write(folder/"batch_protein_CV_summary.csv", summary)
        csv_write(folder/"per_protein_batch_CV_long.csv", long_cv)
        cvs[variant] = summary
        csv_write(folder/"PCA_Raw_fitted_sample_coordinates.csv", pca_results[variant][0])

        client.box_single(raw, records, names, sample_short, box_bounds,
                          "Before normalization", folder/"01_protein_boxplots_before.png")
        client.box_single(corrected, records, names, sample_short, box_bounds,
                          "After mixed-sign normalization", folder/"02_protein_boxplots_after.png")
        client.box_combined(raw, corrected, records, names, sample_short, box_bounds,
                            folder/"03_protein_boxplots_before_after.png")
        draw_shift(names, records, applied, folder, variant)
        reverse.pca_two_separate(*pca_results[variant], pca_bounds, folder, variant)
        arrays, labels = violins.data(folder)
        violins.plot(folder, f"{variant} | T001-T003 -fi; T004 +fi", arrays, labels)
        (folder/"reference_and_processing.json").write_text(json.dumps({
            "variant":variant, "reference_log2":ref_global if variant=="overall_reference" else refs_batch,
            "correction_direction":SIGN,
            "PG_log2_correction":"Raw + sign[batch] * (reference - log2(STGSWSTLK Heavy))",
            "factors_equal_earlier_positive_direction_trials":True,
            "n_protein_groups":4264, "n_samples":39,
            "no_missing_imputation_except_PCA":True,
            "PCA_raw_fit":pca_results[variant][1],
            "PCA_common_axis_bounds":pca_bounds
        }, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")

    baseline = {"T001":(22.55927460503227,26.431399617918338),
                "T002":(16.57163077669857,18.676613123603047),
                "T003":(16.236967208760674,17.782982420792013),
                "T004":(17.776581368637007,21.95522315225137)}
    comparison = []
    for b in BATCHES:
        row_a = next(r for r in cvs["overall_reference"] if r["stage"] == "Corrected" and r["batch"] == b)
        row_b = next(r for r in cvs["within_batch_reference"] if r["stage"] == "Corrected" and r["batch"] == b)
        raw_row = next(r for r in cvs["overall_reference"] if r["stage"] == "Raw" and r["batch"] == b)
        raw_cv = float(raw_row["median_protein_cv_pct"])
        a = float(row_a["median_protein_cv_pct"])
        bb = float(row_b["median_protein_cv_pct"])
        if not np.isclose(a, bb, atol=1e-10, rtol=0):
            raise AssertionError(f"Global and batch CV differ in {b}")
        if not np.isclose(raw_cv, baseline[b][0], atol=1e-9, rtol=0):
            raise AssertionError(f"Raw CV changed from previous trial in {b}")
        if not np.isclose(a, baseline[b][1], atol=1e-9, rtol=0):
            raise AssertionError(f"Mixed CV differs from original matching-sign trial in {b}")
        comparison.append({
            "batch": b, "n_samples":prev.BATCH_COUNTS[b],
            "correction_direction": "-fi" if SIGN[b] < 0 else "+fi",
            "raw_median_cv_pct":raw_cv,
            "overall_mixed_median_cv_pct":a,
            "within_batch_mixed_median_cv_pct":bb,
            "delta_mixed_minus_raw_percentage_points":a-raw_cv,
            "n_proteins_CV":row_a["n_proteins_with_valid_cv"]
        })
    csv_write(OUT/"00_mixed_sign_CV_comparison.csv", comparison)
    cv_report_plot(comparison)

    manifest = {
        "n_original":40,"n_retained":39,"n_protein_groups":4264,
        "excluded_sample_id":prev.EXCLUDED_ID,
        "batch_sizes":prev.BATCH_COUNTS,
        "stgswstlk_only":True,
        "directions":{b:"-fi" if SIGN[b]<0 else "+fi" for b in BATCHES},
        "fi_definition":"reference_median_log2 - log2(original STGSWSTLK Heavy)",
        "reference_log2":{"overall_reference":ref_global,
                          "within_batch_reference":refs_batch},
        "correction_rule":"PG_corrected_log2 = raw_PG_log2 + sign_by_batch * fi",
        "original_inputs_sha256":{p.name:prev.sha256(p) for p in (prev.PG,prev.HELP,prev.META)},
        "raw_pg_valid_values":int(np.isfinite(raw).sum()),
        "no_changes_to_previous_experiments":True,
        "within_batch_cv_identical_between_two_reference_methods":True,
        "same_Raw_PCA_basis_and_axes":True,
        "no_imputation_for_CV_or_matrices":True,
    }
    (OUT/"run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+"\n",encoding="utf-8")

    lines = [
        "# NYD 39 repeated measurements | STGSWSTLK Heavy mixed-sign normalization",
        "",
        "T001/T002/T003: subtract the unchanged original factor fi; T004: add fi.",
        "The 39 samples are repeat measurements of the same specimen.",
        "QC-excluded sample MAP26092326E01 is omitted from the analysis but retained in source.",
        "No original matrix or historical experiment is changed.",
        "",
        "## Definition",
        "",
        "fi = median_ref - log2(original STGSWSTLK Heavy).",
        "PG_corrected_log2 = PG_raw_log2 + sign(batch)*fi.",
        "sign(T001/T002/T003)=-1; sign(T004)=+1.",
        "The overall_reference version uses one 39-sample Heavy reference.",
        "The within_batch_reference version uses each batch's own Heavy median.",
        "All computed references use full floating-point precision.",
        "",
        f"- Global Heavy log2 reference: {ref_global:.15f}",
        *[f"- {b} Heavy log2 reference: {refs_batch[b]:.15f}" for b in BATCHES],
        "",
        "## Within-batch median per-protein CV (%)",
        "",
        "| Batch | n | Sign | Raw | Overall mixed | Within-batch mixed | Difference vs Raw (pp) |",
        "| --- | ---: | --- | ---: | ---: | ---: | ---: |",
        *[f"| {r['batch']} | {r['n_samples']} | {r['correction_direction']} | {r['raw_median_cv_pct']:.2f} | "
          f"{r['overall_mixed_median_cv_pct']:.2f} | {r['within_batch_mixed_median_cv_pct']:.2f} | "
          f"{r['delta_mixed_minus_raw_percentage_points']:+.2f} |" for r in comparison],
        "",
        "CV is sample SD(ddof=1)/arithmetic mean x100 on back-transformed linear intensities.",
        "For each protein within a batch require >=max(3,ceil(.7*n)) valid positive intensities.",
        "Both methods have the SAME within-batch protein CV by construction:",
        "the difference between global and batch medians is a constant within each batch.",
        "This does not make absolute intensities or between-batch positions identical.",
        "Mixed-sign calibration does not beat Raw in these four batch CV medians.",
        "",
        "![Raw vs mixed CV](00_mixed_sign_CV_vs_Raw.png)",
        "",
        "## Outputs",
        "",
        "- 00_Raw_NYD_pg_matrix2_log2_39samples.tsv: 4,264 PG x 39 sample Raw log2 matrix",
        "- 00_mixed_sign_CV_comparison.csv: Raw/overall/batch median CV comparison.",
        "- overall_reference/ and within_batch_reference/: each has",
        "  - corrected_NYD_pg_matrix2_log2.tsv: full 4,264 PG x 39 sample corrected log2 matrix",
        "  - sample_factors_internal.csv: Heavy, original fi, actual direction and applied shift.",
        "  - batch_protein_CV_summary.csv and per_protein_batch_CV_long.csv",
        "  - PCA_Raw_fitted_sample_coordinates.csv and reference_and_processing.json",
        "  - 01/02/03: separate before/after and paired protein intensity boxplots (no targeted highlights)",
        "  - 04: applied factor chart, for internal audit only",
        "  - 04a_PCA_Raw.png and 04b_PCA_Corrected.png (both on one Raw-fitted basis)",
        "  - 08_per_protein_CV_violin_Raw_vs_Corrected_median_labels.png and CSV",
        "",
        "## Before/after plots",
        "",
        "### Overall reference",
        "![Overall boxplots](overall_reference/03_protein_boxplots_before_after.png)",
        "![Overall CV violin](overall_reference/08_per_protein_CV_violin_Raw_vs_Corrected_median_labels.png)",
        "",
        "### Within-batch reference",
        "![Within-batch boxplots](within_batch_reference/03_protein_boxplots_before_after.png)",
        "![Within-batch CV violin](within_batch_reference/08_per_protein_CV_violin_Raw_vs_Corrected_median_labels.png)",
        "",
        "## Interpretation",
        "",
        "This is an exploratory mixture of two directions, motivated by earlier data.",
        "Do not interpret selecting the preferred signs on this SAME 39-sample dataset",
        "as out-of-sample proof of improved technical precision or physical validity.",
        "Within-sample cross-protein Pearson heatmaps would be unchanged by these shifts.",
        "No missing-value imputation is used outside the standard Raw-fitted PCA procedure.",
    ]
    (OUT/"README.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    print("REFERENCES", json.dumps(manifest["reference_log2"]))
    print("CV",json.dumps(comparison))
    print("PASS mixed signs: 2 corrected 4264x39 matrices, fixed Raw PCA, CV, clean PG boxes, per-protein CV violins", flush=True)


if __name__ == "__main__":
    main()
