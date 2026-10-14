#!/usr/bin/env python3
"""Repeat existing HELP-dependent PACS experiments with exactly five peptides.

Original 97-HELP workbooks and historic scripts are left untouched. This
orchestrator patches ONLY the HELP reader for the previously implemented
methods; each method retains its exact original formula and QC treatment.

Reruns: Trial 1A, Trial 1B, Trial 2 stepwise yeast-peptide a/b, and
the two yeast-protein EN01 strategies. Input pg_matrix.tsv is already
metadata-filtered with both yeast protein rows removed (6675 proteins).
"""
import hashlib
import json
import math
import os
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

import pacs_trial1 as t1
import pacs_trial1b as t1b
import pacs_trial2 as t2
import pacs_trial2_stepwise as t2step
from filter_samples import readbook, tab
from yeast_protein_qc import GROUPS, metadata

ROOT = t1.ROOT
OUT = ROOT / "results"
HELP = ROOT / "LH_ratio_tables_2026.10.13—ALL-HELP.xlsx"
CHOSEN = (
    "FSVVYAK",
    "STGSWSTLK",
    "VTGWGNLK",
    "VVEESELAR",
    "VYFAGFPR",
)
METHODS = ("trial1a", "trial1b", "trial2_peptide", "trial3_eno1")
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
help_cache = None


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def prepare_five_help():
    """Validate the full source workbook, select exactly the requested 5 rows
    from ratio_LH, and preserve positive-ratio/missingness semantics."""
    labels = metadata()
    if len(labels) != 75:
        raise ValueError("Expected exactly 75 metadata samples")
    with zipfile.ZipFile(HELP) as z:
        book = ET.fromstring(z.read("xl/workbook.xml"))
        sheet_names = [node.get("name") for node in book.find(NS + "sheets")]
    if sheet_names != ["L", "H", "ratio_LH"]:
        raise ValueError(f"Unexpected HELP sheet layout: {sheet_names}")
    ss, sheets = readbook(HELP)
    if len(sheets) != 3:
        raise ValueError(f"Unexpected HELP sheet count: {len(sheets)}")
    entries = [tab(sheet, ss) for sheet in sheets]
    for sheet_name, rows in zip(sheet_names, entries):
        if len(rows) != 98 or len(set(r.get("A", "") for r in rows[1:])) != 97:
            raise ValueError(f"Expected 97 unique sequences in sheet {sheet_name}")
        if {r.get("A", "") for r in rows[1:]} != {r.get("A", "") for r in entries[2][1:]}:
            raise ValueError("HELP sheet peptides differ")
    ratio = entries[2]
    colmap = {sid: col for col, sid in ratio[0].items() if col != "A"}
    if len(colmap) != 75 or set(colmap) != set(labels):
        raise ValueError("HELP sample IDs differ from 75 metadata sample IDs")
    by_name = {row["A"]: row for row in ratio[1:]}
    if len(CHOSEN) != len(set(CHOSEN)) or not set(CHOSEN).issubset(by_name):
        raise ValueError(f"Missing/duplicated requested HELP sequence: {set(CHOSEN)-set(by_name)}")
    chosen_data = [by_name[seq] for seq in CHOSEN]
    sample_values = {}
    rows = []
    for sid in colmap:
        values = []
        for seq, peptide in zip(CHOSEN, chosen_data):
            cell = peptide.get(colmap[sid], "")
            try:
                v = float(cell)
            except (ValueError, TypeError):
                v = float("nan")
            if math.isfinite(v) and v > 0:
                values.append(v)
        if not values:
            raise ValueError(f"All five HELP ratios missing/nonpositive for sample {sid}; refusing to continue")
        sample_values[sid] = values
        rows.append({
            "sample_id": sid,
            "group": labels[sid],
            "n_requested_peptides": len(CHOSEN),
            "n_valid_selected_HELP_ratios": len(values),
            "raw_ratio_median": float(np.median(values)),
            "log2_ratio_median": float(np.median(np.log2(values))),
            **{
                "ratio_" + seq: next(
                    (float(peptide[colmap[sid]])
                     for peptide in [by_name[seq]]
                     if peptide.get(colmap[sid], "") not in ("", None)
                     and _positive_numeric(peptide.get(colmap[sid], ""))),
                    np.nan,
                )
                for seq in CHOSEN
            },
        })
    frame = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT / "five_help_sample_diagnostics.csv", index=False)
    return {
        "labels": labels,
        "samples": sample_values,
        "raw_medians": {row["sample_id"]: row["raw_ratio_median"] for row in rows},
        "log2_medians": {row["sample_id"]: row["log2_ratio_median"] for row in rows},
        "counts": {row["sample_id"]: row["n_valid_selected_HELP_ratios"] for row in rows},
        "raw_ref": float(frame.loc[frame.group == "pooled QC", "raw_ratio_median"].median()),
        "log2_ref": float(frame.loc[frame.group == "pooled QC", "log2_ratio_median"].median()),
        "n_missing_ratios": int(5 * 75 - frame.n_valid_selected_HELP_ratios.sum()),
        "qc_reference_sample_n": int(sum(frame.group == "pooled QC")),
    }


def _positive_numeric(value):
    try:
        n = float(value)
        return math.isfinite(n) and n > 0
    except (TypeError, ValueError):
        return False


def five_raw_median(requested_labels):
    global help_cache
    if set(requested_labels) != set(help_cache["labels"]):
        raise ValueError("HELP/raw ratio sample IDs not matched to metadata")
    return help_cache["raw_medians"].copy(), help_cache["counts"].copy()


def five_log2_factor(samples, labels):
    global help_cache
    if len(samples) != 75 or set(samples) != set(help_cache["labels"]):
        raise ValueError("HELP/log2 ratio sample IDs not matched to metadata")
    ref = help_cache["log2_ref"]
    if not math.isfinite(ref):
        raise ValueError("HELP pooled QC reference invalid")
    medians = help_cache["log2_medians"]
    return (
        {sid: float(ref - medians[sid]) for sid in samples},
        medians.copy(),
        help_cache["counts"].copy(),
        ref,
    )


def update_legacy_readme(path, label):
    """Keep old script prose but correct its literal hardcoded 97 count."""
    original = path.read_text(encoding="utf-8")
    corrected = re.sub(r"\b97\b", "5", original)
    notice = (
        f"> **Five-HELP rerun — {label}.** The historical procedure has been "
        "re-executed with exactly the following five HELP sequences: "
        + ", ".join(CHOSEN) + ". The original 97-row workbook was not modified. "
        "Only the HELP ratio reader was replaced; all other algorithm choices "
        "and QC handling are inherited from the historical script. "
        "The protein source is the current 6675-row matrix.\n\n"
    )
    heading, body = corrected.split("\n", 1)
    path.write_text(heading + "\n\n" + notice + body, encoding="utf-8")


def ensure_expected_preset():
    # Both source row and sample count are checked before any run.
    with (ROOT / "pg_matrix.tsv").open("r", encoding="utf-8-sig") as f:
        header = f.readline().rstrip("\r\n").split("\t")
        if len(header) != 80 or header[:5] != t1.ANN:
            raise ValueError("pg_matrix.tsv must have 5 annotations and 75 samples")
        if set(header[5:]) != set(help_cache["labels"]):
            raise ValueError("Protein matrix sample columns != metadata")
        count = 0
        for line in f:
            if line.split("\t")[2] in ("ENO1_YEAST", "ENO2_YEAST"):
                raise ValueError("Yeast protein row leaked into endogenous matrix")
            count += 1
    if count != 6675:
        raise ValueError(f"Expected 6675 endogenous protein rows, got {count}")


def run_1a():
    t1.OUT = OUT / "trial1a"
    t1.FIG = t1.OUT / "figures"
    t1.HELP = HELP
    t1.load_ratio = five_raw_median
    print("=== RUN 1A: raw L/H ratio median, correct all groups ===", flush=True)
    t1.main()
    update_legacy_readme(t1.OUT / "README.md", "Trial 1A")


def run_1b():
    t1b.OUT = OUT / "trial1b"
    t1b.load_ratio = five_raw_median
    print("=== RUN 1B: raw L/H ratio median, LM unchanged ===", flush=True)
    t1b.main()
    update_legacy_readme(t1b.OUT / "README.md", "Trial 1B")


def run_2():
    t2step.METHODS = ("a", "b")
    t2step.OUT = OUT / "trial2_peptide"
    t2step.old.OUT = t2step.OUT
    t2step.old.HELP = HELP
    t2step.old.compute_help_factor2 = five_log2_factor
    print("=== RUN 2: yeast peptide a,b + five-HELP log2 Factor2 ===", flush=True)
    t2step.main()
    update_legacy_readme(t2step.OUT / "README.md", "Trial 2 stepwise, yeast peptide a/b")


def run_3():
    import pacs_eno1_protein_trial as t3
    # Importing t3 changes base.METHODS, so import AFTER Trial 2 is finished.
    # base is t2step and t3's two factors come from unchanged EN01.xlsx.
    t3.OUT = OUT / "trial3_eno1"
    t3.base.OUT = t3.OUT
    t3.base.old.OUT = t3.OUT
    t3.base.old.HELP = HELP
    t3.base.old.compute_help_factor2 = five_log2_factor
    print("=== RUN 3: EN01 row-sum/sum-two + five-HELP log2 Factor2 ===", flush=True)
    t3.main()
    update_legacy_readme(t3.OUT / "README.md", "Trial 3 yeast protein replacement")


def read_historical_comparison():
    """Prior runs are displayed for context, never overwritten. Their matrix
    predates yeast row removal, so cross-version CV is descriptive only."""
    map_archives = {
        "trial2_peptide": ROOT / "results_v3" / "all_stage_cv_summary.csv",
        "trial3_eno1": ROOT / "results_v4" / "all_stage_cv_summary.csv",
    }
    records = []
    for name, old in map_archives.items():
        new = OUT / name / "all_stage_cv_summary.csv"
        if not new.exists():
            raise FileNotFoundError(new)
        current = pd.read_csv(new)
        if not old.is_file():
            raise FileNotFoundError(f"Expected historical 97-HELP CV reference {old}")
        previous = pd.read_csv(old)
        for _, rec in current.iterrows():
            key = ((previous.method == rec["method"]) &
                   (previous.stage == rec["stage"]) &
                   (previous.group == rec["group"]))
            found = previous.loc[key]
            if len(found) != 1:
                raise ValueError("Historical comparison method/stage/group missing")
            prior = found.iloc[0]
            records.append({
                "trial": name, "method": rec["method"],
                "stage": rec["stage"], "group": rec["group"],
                "median_cv_5HELP_pct": float(rec["median_cv_pct"]),
                "median_cv_97HELP_pct_historical": float(prior["median_cv_pct"]),
                "delta_5minus97_percentage_points":
                    float(rec["median_cv_pct"] - prior["median_cv_pct"]),
                "new_protein_rows": 6675,
                "historical_protein_rows": 6677,
                "comparison_warning": "Historic 97-HELP matrix included 2 yeast proteins; descriptive only",
            })
    pd.DataFrame(records).to_csv(OUT / "five_vs_97_historical_cv.csv", index=False)
    return records


def validate_main_outputs():
    factor = pd.read_csv(OUT / "trial2_peptide" / "sample_factors.csv")
    factor3 = pd.read_csv(OUT / "trial3_eno1" / "sample_factors.csv")
    if factor.n_valid_HELP_ratios.max() > 5 or factor.n_valid_HELP_ratios.min() < 1:
        raise AssertionError("Trial2 HELP factor used something other than requested five")
    if factor3.n_valid_HELP_ratios.max() > 5 or factor3.n_valid_HELP_ratios.min() < 1:
        raise AssertionError("Trial3 HELP factor used something other than requested five")
    for f in [factor, factor3]:
        for _, r in f.iterrows():
            expected = (help_cache["log2_ref"] -
                        help_cache["log2_medians"][r["sample_id"]])
            if abs(float(r["candidate_factor2"]) - expected) > 1e-10:
                raise AssertionError("HELP Factor2 mismatch from exactly-five computation")
    for name in ("trial1a", "trial1b"):
        fs = pd.read_csv(OUT / name / "normalization_factors.csv")
        if fs.n_valid_HELP_ratios.min() < 1 or fs.n_valid_HELP_ratios.max() > 5:
            raise AssertionError(f"Incorrect five HELP count in {name}")
        for _, r in fs.iterrows():
            if abs(float(r["sample_mR"]) - help_cache["raw_medians"][r["sample_id"]]) > 1e-9:
                raise AssertionError(f"Incorrect five HELP raw median in {name}")


def summary_report(comparison):
    samples = pd.read_csv(OUT / "five_help_sample_diagnostics.csv")
    cv2 = pd.read_csv(OUT / "trial2_peptide" / "all_stage_cv_summary.csv")
    cv3 = pd.read_csv(OUT / "trial3_eno1" / "all_stage_cv_summary.csv")
    archive = os.environ.get("PACS_RESULTS_ARCHIVE", "results_vN")
    lines = [
        "# PACS five-HELP reruns — same existing methods, only HELP set reduced", "",
        f"Previous results/ archived intact as **{archive}/**; all previous results_vN/ retained.",
        "The five peptides were selected explicitly by the user, not algorithmically optimized on these results.",
        "This is a controlled re-run of the existing four types of experiments,",
        "NOT a new batch-correction model or batch-effect audit.", "",
        "## Exactly five HELP sequences, all from ratio_LH sheet", "",
        *[f"{i}. {peptide}" for i, peptide in enumerate(CHOSEN, 1)],
        "",
        "The original 97-row 3-sheet workbook is unchanged. Exactly these five sequence",
        "rows are selected from the existing 97, for every 75-sample HELP calculation.",
        "Only valid finite positive ratios are included in each sample median.",
        f"The 75 × 5 ratio grid contains {help_cache['n_missing_ratios']} missing/nonpositive values.",
        f"Valid HELP counts per sample range {samples.n_valid_selected_HELP_ratios.min()}",
        f"to {samples.n_valid_selected_HELP_ratios.max()}; none has zero usable HELP ratios.",
        f"All five-HELP raw-ratio pooled-QC reference: {help_cache['raw_ref']:.12g}.",
        f"All five-HELP log2-ratio pooled-QC reference: {help_cache['log2_ref']:.12g}.", "",
        "## Input provenance", "",
        f"- pg_matrix.tsv 6675 protein groups: sha256 {sha256(ROOT / 'pg_matrix.tsv')}.",
        f"- Full HELP workbook used only as row source: {HELP.name}: sha256 {sha256(HELP)}.",
        f"- Yeast peptide source, unchanged: sha256 {sha256(ROOT / 'sample-QC-InternalQC-HELP_yeast.xlsx')}.",
        f"- Yeast protein source, unchanged: sha256 {sha256(ROOT / 'EN01.xlsx')}.",
        "- Metadata: 75 samples, HC=27, S=27, Internal-QC=6, pooled QC=15.", "",
        "## Four historical algorithms re-executed", "",
        "| Folder | Exact inherited procedure | QC treatment |",
        "| --- | --- | --- |",
        "| [trial1a](trial1a/README.md) | median(raw L/H), pooled-QC ref; divide protein by sample/ref factor; includes multiplied-direction sensitivity | All groups corrected, including QC |",
        "| [trial1b](trial1b/README.md) | Same as 1A, two PCA views | Internal-QC unchanged, pooled QC corrected |",
        "| [trial2_peptide](trial2_peptide/README.md) | Yeast peptide Factor1 methods a/b (unchanged) + median(log2 HELP L/H) Factor2; Raw/F1/F1+F2 | Both QC groups unchanged |",
        "| [trial3_eno1](trial3_eno1/README.md) | EN01 protein Factor1: higher row sum vs linear sum of both (unchanged) + median(log2 HELP L/H) Factor2 | Both QC groups unchanged |",
        "",
        "Do not equate raw-ratio medians used in Trial 1A/1B with median(log2 ratio)",
        "used in Trial2/3; the old procedures explicitly differed in this respect.", "",
        "## Trial 2 and EN01 median protein CV (%)", "",
        "CV on linear protein intensities, HC/S include biological variation;",
        "unchanged QC CV in Trial2/3 cannot independently validate normalization.", "",
        "| Trial | Method | Stage | HC | S |",
        "| --- | --- | --- | ---: | ---: |",
    ]
    for name, cv in [("trial2_peptide", cv2), ("trial3_eno1", cv3)]:
        for method in cv.method.drop_duplicates():
            for stage in ("raw", "factor1", "factor1_factor2"):
                v = {g: cv.loc[(cv.method == method) & (cv.stage == stage) &
                               (cv.group == g), "median_cv_pct"].iloc[0]
                     for g in ["HC", "S"]}
                lines.append(f"| {name} | {method} | {stage} | {v['HC']:.2f}% | {v['S']:.2f}% |")
    lines += [
        "", "## 5-vs-97 historical context — interpretation limit", "",
        "See [five_vs_97_historical_cv.csv](five_vs_97_historical_cv.csv).",
        "Historical Trial2 97-HELP data are read only from results_v3/;",
        "historical EN01 97-HELP data are read only from results_v4/.",
        "**Important:** Those older outputs contain 6677 proteins; the current matrix",
        "contains 6675 after removing two yeast proteins. Therefore the historical",
        "CV differences are not a perfectly isolated 5-vs-97 HELP comparison.",
        "Factor1-only results are not expected to change because HELP is used at Factor2,",
        "apart from any changes in the protein set. Do not claim causal batch reduction.",
        "Raw, Factor1, and Factor1+Factor2 are shown together for both Trial2/3.", "",
        "## Files and reproducibility", "",
        "- [five_help_sample_diagnostics.csv](five_help_sample_diagnostics.csv): 5 requested",
        "  ratios, valid-count, raw/log2 medians, and QC/reference per 75 samples.",
        "- [five_vs_97_historical_cv.csv](five_vs_97_historical_cv.csv): prior-run comparison",
        "  with clear protein-row caveat.",
        "- [run_manifest.json](run_manifest.json): sources, selected sequences, methods.",
        "- Each trial directory has its complete matrices, factors, CVs, PCA and original plots.",
        "- Replay only via scripts/pacs_help5_reruns.py; the historical scripts remain",
        "  configured for all 97 as a preserved historical baseline.", "",
    ]
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")
    manifest = {
        "archive_previous_results": archive,
        "selected_HELP_sequences": list(CHOSEN),
        "selected_count": len(CHOSEN),
        "original_HELP_count": 97,
        "sample_count": 75,
        "protein_rows": 6675,
        "valid_peptide_count_distribution": {
            str(k): int(v) for k,v in samples.n_valid_selected_HELP_ratios.value_counts().sort_index().items()
        },
        "input_sha256": {
            "pg_matrix.tsv": sha256(ROOT / "pg_matrix.tsv"),
            "HELP": sha256(HELP),
            "yeast_peptides": sha256(ROOT / "sample-QC-InternalQC-HELP_yeast.xlsx"),
            "EN01.xlsx": sha256(ROOT / "EN01.xlsx"),
        },
        "methods": list(METHODS),
        "other_parameters_unchanged": True,
        "description": "Trial1A/1B: raw L/H median. Trial2/3: log2 L/H median.",
        "comparison_has_protein_set_change": True,
    }
    (OUT / "run_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")


def main():
    global help_cache
    if OUT.exists() and any(OUT.iterdir()):
        raise ValueError("Target results/ must be new and empty, to protect archived results")
    if not HELP.exists():
        raise FileNotFoundError(HELP)
    help_cache = prepare_five_help()
    ensure_expected_preset()
    print("FIVE PEPTIDES:", CHOSEN, flush=True)
    print("Valid selected HELP counts:", help_cache["counts"], flush=True)
    print("HELP pooled-QC reference raw/log2:",
          help_cache["raw_ref"], help_cache["log2_ref"], flush=True)
    run_1a()
    run_1b()
    run_2()
    run_3()
    validate_main_outputs()
    comparison = read_historical_comparison()
    summary_report(comparison)
    print("PASS: complete five-HELP reruns, 1A 1B Trial2 peptide and Trial3 ENO1",
          "archive =", os.environ.get("PACS_RESULTS_ARCHIVE", "unknown"), flush=True)


if __name__ == "__main__":
    main()
