#!/usr/bin/env python3
"""PACS yeast spike-in: metadata-group detection and raw-intensity CV.

Select rows where Protein.Names contains YEAST (case-insensitive).
Detected = finite positive intensity. CV = sample SD / mean * 100,
provided >= max(3, ceil(0.70 * group size)) detections and mean >= 1e-8.
Missing/zero values do not contribute to CV. No imputation or normalization.
"""
import csv
import math
import re
import statistics as st
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

GROUPS = ["Internal-QC", "pooled QC", "HC", "S"]
COLORS = ["#4777B5", "#E69F45", "#48A48A", "#B26ABB"]
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "yeast_qc"
ANN = ["Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description"]


def excel_rows(path):
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            ss = ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared = ["".join(t.text or "" for t in si.iter(NS + "t"))
                      for si in ss.findall(NS + "si")]
        doc = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        rows = []
        for row in doc.find(NS + "sheetData").findall(NS + "row"):
            cells = {}
            for c in row.findall(NS + "c"):
                letters = re.match(r"[A-Z]+", c.get("r")).group()
                i = 0
                for letter in letters:
                    i = i * 26 + ord(letter) - 64
                v = c.find(NS + "v")
                if c.get("t") == "inlineStr":
                    val = "".join(t.text or "" for t in c.iter(NS + "t"))
                elif v is None:
                    val = ""
                elif c.get("t") == "s":
                    val = shared[int(v.text)]
                else:
                    val = v.text or ""
                cells[i - 1] = str(val).strip()
            if cells:
                rows.append([cells.get(i, "") for i in range(max(cells) + 1)])
        return [dict(zip(rows[0], row)) for row in rows[1:]]


def metadata():
    result = {}
    for file, qc in (("all_sample_metadata.xlsx", False),
                     ("all_QC_metadata.xlsx", True)):
        for row in excel_rows(ROOT / file):
            labels = {row.get("group1", "").upper(), row.get("group2", "").upper()}
            if "BLANK" in labels:
                continue
            possible = ({"INTERNAL-QC": "Internal-QC", "QC1": "pooled QC"}
                        if qc else {"HC": "HC", "S": "S"})
            groups = {possible[x] for x in labels if x in possible}
            sid = row["sample_id"]
            if len(groups) != 1 or sid in result or not sid:
                raise ValueError("Ambiguous or duplicate metadata sample: " + str(row))
            result[sid] = groups.pop()
    return result


def intensity(value):
    try:
        value = float(value)
        return value if math.isfinite(value) and value > 0 else None
    except (ValueError, TypeError):
        return None


def cv(values, n):
    if len(values) < max(3, math.ceil(n * 0.70)):
        return None
    avg = st.mean(values)
    return 100 * st.stdev(values) / avg if avg >= 1e-8 else None


def save_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def figures(proteins, summary, sample_stats):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    dest = OUT / "figures"
    dest.mkdir(parents=True, exist_ok=True)

    def save(fig, name):
        for extension in ("svg", "png"):
            fig.savefig(dest / (name + "." + extension),
                        bbox_inches="tight", dpi=220, facecolor="white")
        plt.close(fig)

    def empty(ax):
        ax.text(0.5, 0.5, "No YEAST entries in Protein.Names",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=13, color="#64748b")
        ax.set_xticks([])
        ax.set_yticks([])

    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    counts = [s["detected_protein_groups"] for s in summary]
    bars = ax.bar(range(4), counts, color=COLORS, width=0.6)
    ax.bar_label(bars, padding=4)
    ax.set_xticks(range(4), [f'{g} (n={summary[i]["sample_n"]})'
                             for i, g in enumerate(GROUPS)])
    ax.set_ylim(0, max([1] + counts) * 1.2)
    ax.set_ylabel("Protein groups detected in at least one sample")
    ax.set_title("Yeast proteins detected by group")
    save(fig, "01_group_detected_counts")

    # Protein-level rates and CVs; paginate so all identities remain legible.
    ordered = sorted(proteins, key=lambda r: (-r["total_detected"], r["Protein.Names"]))
    for kind, prefix, title, xlabel in (
        ("detection_rate_pct", "02_per_protein_detection_rates",
         "Yeast protein detection rate by group", "Detection rate (%)"),
        ("CV_pct", "04_per_protein_cv",
         "Yeast protein intensity CV by group", "CV (%)"),
    ):
        batches = [ordered[i:i + 30] for i in range(0, len(ordered), 30)] or [[]]
        for page, batch in enumerate(batches):
            fig, ax = plt.subplots(figsize=(12, max(5, len(batch) * 0.45 + 2)),
                                   layout="constrained")
            ax.set_title(title)
            if batch:
                ys = list(range(len(batch)))
                for j, (g, color) in enumerate(zip(GROUPS, COLORS)):
                    vals = [r[g + "_" + kind] for r in batch]
                    yy = [v + (j - 1.5) * 0.19 for v in ys]
                    if kind == "detection_rate_pct":
                        ax.barh(yy, [x if x is not None else 0 for x in vals],
                                height=0.18, color=color, label=g)
                    else:
                        ax.scatter([v for v in vals if v is not None],
                                   [y for y, v in zip(yy, vals) if v is not None],
                                   s=34, color=color, label=g)
                labels = [r["Protein.Names"].split(";")[0] + " [" + r["Protein.Group"] + "]"
                          for r in batch]
                ax.set_yticks(ys, labels, fontsize=8)
                ax.invert_yaxis()
                ax.legend(loc="upper left", bbox_to_anchor=(1, 1), frameon=False)
                ax.set_xlabel(xlabel)
                if kind == "detection_rate_pct":
                    ax.set_xlim(0, 105)
                ax.grid(axis="x", alpha=0.2)
            else:
                empty(ax)
            extra = f"_{page+1:02d}" if len(batches) > 1 else ""
            save(fig, prefix + extra)

    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    ax.set_title("Yeast protein intensity CV distribution")
    vals = [[r[g + "_CV_pct"] for r in proteins if r[g + "_CV_pct"] is not None]
            for g in GROUPS]
    if any(vals):
        for i, (v, color) in enumerate(zip(vals, COLORS), 1):
            if len(v) >= 2 and max(v) > min(v):
                parts = ax.violinplot([v], positions=[i], widths=0.75,
                                      showmeans=False, showextrema=False)
                for body in parts["bodies"]:
                    body.set_facecolor(color)
                    body.set_alpha(0.4)
            if v:
                bx = ax.boxplot([v], positions=[i], widths=0.18,
                                showfliers=False, patch_artist=True,
                                manage_ticks=False)
                for b in bx["boxes"]:
                    b.set_facecolor(color)
                    b.set_alpha(0.6)
                jitter = [((j * 37) % 23 - 11) / 85 for j in range(len(v))]
                ax.scatter([i + a for a in jitter], v, color=color,
                           alpha=0.65, s=15, zorder=3)
                ax.annotate(f"Median={st.median(v):.2f}%",
                            (i, max(v)), xytext=(0, 12),
                            textcoords="offset points", ha="center", fontsize=9)
        ax.set_xlim(0.4, 4.6)
        ax.set_xticks(range(1, 5),
                      [f"{g}\n(n={summary[i]['sample_n']}, CV={len(vals[i])})"
                       for i, g in enumerate(GROUPS)])
    else:
        empty(ax)
        ax.text(0.5, 0.35, "CV undefined; no imputed zero CV",
                transform=ax.transAxes, ha="center", color="#64748b")
    ax.set_ylabel("CV (%) = sample SD / mean × 100")
    save(fig, "03_group_cv_violin_box")


    # Supplementary identification plots (bars); CV keeps the violin and scatter views.
    fig, ax = plt.subplots(figsize=(9.5, 5.2), layout="constrained")
    means = [r["mean_detected_per_sample"] for r in summary]
    bars = ax.bar(range(4), means, color=COLORS, width=0.6)
    ax.bar_label(bars, labels=[f"{v:.1f}" for v in means], padding=5)
    ax.set_xticks(range(4), [f'{g} (n={summary[i]["sample_n"]})'
                             for i, g in enumerate(GROUPS)])
    ax.set_ylim(0, max(means) * 1.2 if any(means) else 1)
    ax.set_ylabel("Mean detected yeast protein groups / sample")
    ax.set_title("Mean yeast protein identification per sample")
    save(fig, "05_group_mean_identification")

    ordered_samples = [r for g in GROUPS for r in sample_stats if r["group"] == g]
    fig, ax = plt.subplots(figsize=(16, 5.5), layout="constrained")
    indexes = list(range(len(ordered_samples)))
    color_values = [COLORS[GROUPS.index(r["group"])] for r in ordered_samples]
    ax.bar(indexes, [r["detected_yeast_protein_groups"] for r in ordered_samples],
           color=color_values, width=0.85)
    ax.set_xticks(indexes, [r["sample_id"] for r in ordered_samples],
                  rotation=85, fontsize=6)
    ax.set_ylabel("Detected yeast protein groups")
    ax.set_title("Yeast protein identification by sample")
    ax.legend(handles=[Patch(color=c, label=f'{g} (n={summary[i]["sample_n"]})')
                       for i, (g, c) in enumerate(zip(GROUPS, COLORS))],
              loc="upper right", frameon=False)
    offset = 0
    for g in GROUPS[:-1]:
        offset += sum(r["group"] == g for r in ordered_samples)
        ax.axvline(offset - 0.5, color="#94a3b8", lw=0.9, ls="--")
    save(fig, "06_per_sample_identification")

    fig, ax = plt.subplots(figsize=(9.5, 5.2), layout="constrained")
    counts = [r["CV_eligible_protein_groups"] for r in summary]
    bars = ax.bar(range(4), counts, color=COLORS, width=0.6)
    ax.bar_label(bars, padding=5)
    ax.set_xticks(range(4), [f'{g} (n={summary[i]["sample_n"]})'
                             for i, g in enumerate(GROUPS)])
    ax.set_ylim(0, max(counts) * 1.2 if any(counts) else 1)
    ax.set_ylabel("Protein groups with estimable CV")
    ax.set_title("Yeast proteins eligible for CV (70% detected)")
    save(fig, "07_cv_eligible_counts")


def main():
    samples = metadata()
    with (ROOT / "pg_matrix.tsv").open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        cols = reader.fieldnames
        if not cols or cols[:5] != ANN or len(cols[5:]) != len(set(cols[5:])):
            raise ValueError("Unexpected pg_matrix.tsv header")
        if set(cols[5:]) != set(samples):
            raise ValueError(f"pg_matrix sample IDs differ from metadata: "
                             f"extra={set(cols[5:]) - set(samples)}, "
                             f"missing={set(samples) - set(cols[5:])}")
        cols_by_group = {g: [s for s in cols[5:] if samples[s] == g] for g in GROUPS}
        if any(not x for x in cols_by_group.values()):
            raise ValueError("Missing one of the four analysis groups")
        yeast = [row for row in reader if "YEAST" in row["Protein.Names"].upper()]
    if not yeast:
        raise ValueError("No YEAST protein names found; check the new pg_matrix.tsv")
    if len({r["Protein.Group"] for r in yeast}) != len(yeast):
        raise ValueError("Duplicate Protein.Group among yeast rows")

    protein_fields = ANN + ["total_detected"]
    for g in GROUPS:
        protein_fields += [g + "_detected_n", g + "_detection_rate_pct", g + "_CV_pct"]
    sample_stats = [{"sample_id": s, "group": samples[s], "detected_yeast_protein_groups": 0}
                    for s in cols[5:]]
    stats_index = {r["sample_id"]: r for r in sample_stats}
    proteins = []
    for row in yeast:
        out = {a: row[a] for a in ANN}
        out["total_detected"] = 0
        for g, ids in cols_by_group.items():
            observed = [(s, intensity(row[s])) for s in ids]
            values = [v for _, v in observed if v is not None]
            for s, v in observed:
                if v is not None:
                    stats_index[s]["detected_yeast_protein_groups"] += 1
            out[g + "_detected_n"] = len(values)
            out[g + "_detection_rate_pct"] = 100 * len(values) / len(ids)
            out[g + "_CV_pct"] = cv(values, len(ids))
            out["total_detected"] += len(values)
        proteins.append(out)

    summary = []
    for g, ids in cols_by_group.items():
        matched = [r for r in proteins if r[g + "_detected_n"]]
        cv_vals = [r[g + "_CV_pct"] for r in proteins if r[g + "_CV_pct"] is not None]
        per_sample = [stats_index[s]["detected_yeast_protein_groups"] for s in ids]
        summary.append({
            "group": g, "sample_n": len(ids),
            "detected_protein_groups": len(matched),
            "distinct_Protein.Names": len({r["Protein.Names"] for r in matched}),
            "mean_detected_per_sample": st.mean(per_sample),
            "median_detected_per_sample": st.median(per_sample),
            "CV_eligible_protein_groups": len(cv_vals),
            "median_protein_CV_pct": st.median(cv_vals) if cv_vals else None,
            "identification_count_CV_pct": (
                100 * st.stdev(per_sample) / st.mean(per_sample)
                if len(per_sample) >= 2 and st.mean(per_sample) != 0 else None
            ),
        })
    save_csv(OUT / "group_summary.csv", list(summary[0]), summary)
    save_csv(OUT / "per_protein.csv", protein_fields, proteins)
    save_csv(OUT / "per_sample.csv", list(sample_stats[0]), sample_stats)
    figures(proteins, summary, sample_stats)
    detected = sum(p["total_detected"] > 0 for p in proteins)
    print(f"YEAST annotated protein groups={len(proteins)}, detected in >=1 sample={detected}, unique names={len({r['Protein.Names'] for r in proteins})}")
    for r in summary:
        print(f"{r['group']}: n={r['sample_n']}, detected={r['detected_protein_groups']}, CV eligible={r['CV_eligible_protein_groups']}")
    if not proteins:
        print("SOURCE DATA ISSUE: NO YEAST IN Protein.Names; CV is UNDEFINED.")


if __name__ == "__main__":
    main()
