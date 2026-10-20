#!/usr/bin/env python3
"""Compare the five user-chosen PACS HELP Heavy intensities with NYD PG medians.

39 NYD repeated measurements; exclude MAP26092326E01. Two independent PG
summaries: sample-wise all detected PGs vs the fixed intersection across all 39.
Only read original Heavy.xlsx / NYD_pg_matrix2.tsv / NYD_sampleinfo2.xlsx.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

import nyd_STGSWSTLK_heavy_exclude_T003_outlier as base
from nyd_HELP5_Heavy_CV_MAD_normalization import heavy_97

PEPTIDES = ("FSVVYAK", "STGSWSTLK", "VTGWGNLK", "VVEESELAR", "VYFAGFPR")
BATCHES = ("T001", "T002", "T003", "T004")
OUT = base.BASE / "user_selected5_HELP_Heavy_PG_Pearson_39samples"
COLORS = dict(zip(BATCHES, ("#337DB2", "#E39C43", "#329A7D", "#B06691")))


def save_csv(path: Path, records: list[dict]) -> None:
    if not records:
        raise ValueError("No records: " + path.name)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def pearson_valid(x, y):
    keep = np.isfinite(x) & np.isfinite(y)
    a, b = x[keep], y[keep]
    if a.size < 3 or np.std(a) <= 1e-12 or np.std(b) <= 1e-12:
        return int(a.size), None, None
    result = pearsonr(a, b)
    return int(a.size), float(result.statistic), float(result.pvalue)


def draw_scatter(names, batches, heavy_log, medians, correlations):
    fig, axes = plt.subplots(5, 2, figsize=(13.5, 19), layout="constrained")
    for k, sequence in enumerate(PEPTIDES):
        for col, (definition, yy) in enumerate(medians.items()):
            ax = axes[k, col]
            for batch in BATCHES:
                mask = np.array([b == batch for b in batches])
                a, b = heavy_log[k, mask], yy[mask]
                good = np.isfinite(a) & np.isfinite(b)
                ax.scatter(a[good], b[good], s=36, color=COLORS[batch],
                           alpha=0.8, label=batch, edgecolors="white", linewidths=0.35)
            summary = next(r for r in correlations
                           if r["scope"] == "ALL_39"
                           and r["HELP"] == sequence
                           and r["PG_basis"] == definition)
            rv = summary["pearson_r"]
            annotation = "undefined" if rv is None else f"{rv:+.3f}"
            ax.set_title(f"{sequence} | {definition} | overall r = {annotation}", fontsize=10.5)
            ax.set_xlabel("HELP Heavy intensity (log2)")
            ax.set_ylabel("PG median intensity (log2)")
            ax.grid(alpha=0.17)
            if k == 0 and col == 1:
                ax.legend(frameon=False, fontsize=8, ncol=4)
    fig.suptitle("NYD 39 repeated measurements | five user-selected HELP Heavy vs PG median",
                 fontsize=14, fontweight="bold")
    fig.savefig(OUT / "01_scatter_5HELP_two_PG_definitions.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


def draw_heatmap(correlations):
    scopes = ("ALL_39", *BATCHES)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), layout="constrained")
    for ax, definition in zip(axes, ("all_detected_PG", "common_39_PG")):
        value = np.full((len(scopes), len(PEPTIDES)), np.nan)
        for i, scope in enumerate(scopes):
            for j, sequence in enumerate(PEPTIDES):
                match = next(r for r in correlations
                             if r["scope"] == scope
                             and r["HELP"] == sequence
                             and r["PG_basis"] == definition)
                if match["pearson_r"] is not None:
                    value[i, j] = match["pearson_r"]
        im = ax.imshow(value, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        ax.set_xticks(range(len(PEPTIDES)), PEPTIDES, rotation=40, ha="right")
        ax.set_yticks(range(len(scopes)), scopes)
        ax.set_title("PG median: " + definition)
        for i in range(len(scopes)):
            for j in range(len(PEPTIDES)):
                if np.isfinite(value[i, j]):
                    ax.text(j, i, f"{value[i,j]:+.2f}", va="center", ha="center",
                            color="white" if abs(value[i,j]) > 0.65 else "#202C35",
                            fontsize=10, fontweight="bold")
    fig.colorbar(im, ax=list(axes), fraction=0.024, label="Pearson r")
    fig.savefig(OUT / "02_Pearson_r_by_batch_heatmaps.png", dpi=190, bbox_inches="tight")
    plt.close(fig)


def main():
    if OUT.exists():
        raise FileExistsError("Results directory already exists: " + str(OUT))
    records, names, annotations, raw_pg, _, _, _ = base.inputs()
    assert raw_pg.shape == (4264, 39) and len(set(names)) == 39
    batches = [row["batch"] for row in records]
    assert {b: batches.count(b) for b in BATCHES} == base.BATCH_COUNTS
    assert base.EXCLUDED_ID not in names

    sequences, all_heavy = heavy_97(names)
    missing = set(PEPTIDES) - set(sequences)
    if missing:
        raise ValueError("Missing requested HELP sequences in NYD Heavy.xlsx: " + repr(missing))
    selected = all_heavy[[sequences.index(p) for p in PEPTIDES], :]
    good = np.isfinite(selected) & (selected > 0)
    heavy_log = np.full(selected.shape, np.nan)
    heavy_log[good] = np.log2(selected[good])
    n_valid_pg = np.isfinite(raw_pg).sum(axis=0)
    shared_mask = np.isfinite(raw_pg).all(axis=1)
    n_shared = int(shared_mask.sum())
    if n_shared < 100:
        raise AssertionError(f"Unexpectedly few shared PG: {n_shared}")
    all_median = np.nanmedian(raw_pg, axis=0)
    shared_median = np.median(raw_pg[shared_mask], axis=0)
    if not (np.isfinite(all_median).all() and np.isfinite(shared_median).all()):
        raise AssertionError("Sample PG medians have non-finite values")

    # The five-peptide composite is extra context, NOT a new calibration factor.
    # It requires all five Heavy observations; per-peptide centering removes only
    # peptide-specific constants and preserves each peptide's Pearson r.
    per_peptide_ref = np.nanmedian(heavy_log, axis=1)
    centered = heavy_log - per_peptide_ref[:, None]
    composite = np.full(39, np.nan)
    all_five_good = np.isfinite(centered).all(axis=0)
    composite[all_five_good] = centered[:, all_five_good].mean(axis=0)
    features = {p: heavy_log[k] for k, p in enumerate(PEPTIDES)}
    features["5HELP_mean_centered_log2"] = composite
    medians = {"all_detected_PG": all_median, "common_39_PG": shared_median}

    details = []
    for j, sid in enumerate(names):
        row = {"sample_id": sid, "batch": batches[j],
               "n_detected_PG": int(n_valid_pg[j]),
               "n_shared_PG": n_shared,
               "all_detected_PG_log2_median": float(all_median[j]),
               "common_39_PG_log2_median": float(shared_median[j])}
        for k, seq in enumerate(PEPTIDES):
            row[seq + "_Heavy_linear"] = float(selected[k, j]) if good[k, j] else None
            row[seq + "_Heavy_log2"] = float(heavy_log[k, j]) if good[k, j] else None
        row["n_valid_selected_HELP"] = int(np.isfinite(heavy_log[:, j]).sum())
        row["5HELP_mean_centered_log2"] = float(composite[j]) if all_five_good[j] else None
        details.append(row)

    correlations = []
    for scope in ("ALL_39", *BATCHES):
        idx = np.array([scope == "ALL_39" or b == scope for b in batches])
        for label, x in features.items():
            for pg_basis, y in medians.items():
                n, r, p = pearson_valid(x[idx], y[idx])
                correlations.append({
                    "scope": scope, "n_samples_scope": int(idx.sum()), "HELP": label,
                    "PG_basis": pg_basis, "n_pairs": n,
                    "pearson_r": r, "p_value_two_sided": p
                })

    OUT.mkdir(parents=True)
    save_csv(OUT / "sample_level_log2_intensities.csv", details)
    save_csv(OUT / "pearson_summary_all_and_by_batch.csv", correlations)
    shared = []
    for i, ann in enumerate(annotations):
        if shared_mask[i]:
            shared.append({"source_PG_row_number": i + 2,
                           **dict(zip(("Protein.Group","Protein.Ids","Protein.Names",
                                       "Genes","First.Protein.Description"), ann))})
    save_csv(OUT / "common_39_PG_protein_list.csv", shared)
    draw_scatter(names, batches, heavy_log, medians, correlations)
    draw_heatmap(correlations)

    source_paths = (base.PG, base.HELP, base.META)
    manifest = {
        "n_samples": len(names), "excluded_sample": base.EXCLUDED_ID,
        "batch_sizes": base.BATCH_COUNTS, "n_PG_total": len(annotations),
        "n_PG_common_39": n_shared,
        "HELP_selected_manually_from_PACS": list(PEPTIDES),
        "HELP_missing_by_peptide": {
            p: int((~np.isfinite(heavy_log[k])).sum()) for k, p in enumerate(PEPTIDES)
        },
        "data": "NYD original Heavy.xlsx H intensity and raw PG linear intensities",
        "PG_aggregation": "median(log2(valid positive PG intensity))",
        "PG_definitions": {
            "all_detected_PG": "Each sample's OWN valid positive PG set",
            "common_39_PG": "The SAME PG rows finite positive in all 39"
        },
        "correlation": "Pearson across sample observations, pairwise complete",
        "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    }
    (OUT / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = [
        "# NYD 39 sample | five user-selected HELP Heavy vs raw PG medians",
        "",
        "Uses the five MANUALLY specified PACS HELP peptides: " + ", ".join(PEPTIDES) + ".",
        "Do not confuse them with the separate NYD CV/MAD-selected five.",
        "",
        "All values come from original Heavy.xlsx and NYD_pg_matrix2.tsv.",
        "The T003 sample MAP26092326E01 is excluded without changing input files.",
        "The 39 samples are repeated measurements of the same specimen.",
        "",
        "## The two PG median definitions",
        "",
        "1. all_detected_PG: for each sample independently, take median of log2 positive",
        "   protein-group intensities among the PGs detected in THAT sample.",
        "   Protein composition can differ by sample due to missingness.",
        "2. common_39_PG: select once the PG rows with positive finite intensity",
        "   in every one of the 39 samples, then take per-sample log2 medians",
        "   using precisely the SAME protein rows for all samples.",
        f"   Common set size: {n_shared} PG out of 4264.",
        "",
        "For each HELP, Pearson correlates its log2 raw Heavy intensity across samples",
        "with each of the two per-sample PG log2 medians, overall and within batch.",
        "n_pairs counts the samples usable for both values. No imputation.",
        "Two-sided p-values are exploratory and unadjusted for multiple comparisons.",
        "A five-HELP mean of peptide-wise centered log2 Heavy values is also shown",
        "as a secondary composite; it is NOT used to normalize PGs.",
        "",
        "## Outputs",
        "",
        "- pearson_summary_all_and_by_batch.csv: 6 HELP representations x 2 PG definitions x 5 scopes.",
        "- sample_level_log2_intensities.csv: 39 auditable sample-level records.",
        "- common_39_PG_protein_list.csv: annotation list of fixed shared PG set.",
        "- 01_scatter_5HELP_two_PG_definitions.png: two scatter comparisons per peptide.",
        "- 02_Pearson_r_by_batch_heatmaps.png: each selected peptide's Pearson r.",
        "- run_manifest.json: provenance and exact processing definitions.",
        "",
        "![Selected HELP vs two PG medians](01_scatter_5HELP_two_PG_definitions.png)",
        "",
        "![Per-batch Pearson r](02_Pearson_r_by_batch_heatmaps.png)",
        "",
        "## Interpretation",
        "",
        "These are correlations across repeated sample measurements, NOT sample-to-sample",
        "PG correlation heatmaps. Negative r does not by itself validate a -fi correction.",
        "A change in sign between the two PG definitions suggests composition/missingness",
        "effects may affect sample-wise medians. Batch sample sizes are small.",
        "No historical files or source matrices were modified."
    ]
    (OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({"n_PG_common_39": n_shared, "HELP_missing": manifest["HELP_missing_by_peptide"],
                      "overall": [r for r in correlations if r["scope"] == "ALL_39"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
