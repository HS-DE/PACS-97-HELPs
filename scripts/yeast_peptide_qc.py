#!/usr/bin/env python3
"""Yeast spike-in peptides: metadata sample filtering and four-group QC."""
import statistics as stats
from filter_samples import T, ci, cn, readbook, rewrite, tab
from yeast_protein_qc import COLORS, GROUPS, ROOT, cv, intensity, metadata, save_csv

SOURCE = ROOT / "sample-QC-InternalQC-HELP_yeast.xlsx"
OUT = ROOT / "results" / "yeast_peptide_qc"


def load_and_filter(accepted):
    ss, sheets = readbook(SOURCE)
    assert len(sheets) == 1
    root = sheets[0]
    raw = tab(root, ss)
    assert len(raw) == 3, "Expected two peptide rows"
    peptides = [r["A"] for r in raw[1:]]
    assert len(set(peptides)) == 2
    columns = sorted([(ci(c), sid) for c, sid in raw[0].items() if c != "A"])
    source_ids = [sid for _, sid in columns]
    assert len(set(source_ids)) == len(source_ids), "Duplicate source sample ID"
    assert set(accepted).issubset(source_ids), "Metadata samples missing"
    keep, next_col = {1: 1}, 2
    for old_col, sid in columns:
        if sid in accepted:
            keep[old_col] = next_col
            next_col += 1
    assert next_col == 77, "Expected exactly 75 retained columns"
    removed = [sid for sid in source_ids if sid not in accepted]
    if removed:
        for row in root.find(T("sheetData")).findall(T("row")):
            rn = row.get("r")
            row.set("spans", "1:76")
            for cell in list(row.findall(T("c"))):
                old = ci("".join(x for x in cell.get("r") if x.isalpha()))
                if old not in keep:
                    row.remove(cell)
                else:
                    cell.set("r", cn(keep[old]) + rn)
        dim = root.find(T("dimension"))
        if dim is not None:
            dim.set("ref", "A1:BX3")
        rewrite(SOURCE, {"xl/worksheets/sheet1.xml": root})
    ss, sheets = readbook(SOURCE)
    clean = tab(sheets[0], ss)
    samples = [sid for c, sid in sorted(clean[0].items(), key=lambda x: ci(x[0])) if c != "A"]
    assert len(samples) == 75 and set(samples) == set(accepted)
    assert [r["A"] for r in clean[1:]] == peptides
    intensities = {r["A"]: {sid: intensity(r.get(cn(j), ""))
                             for j, sid in enumerate(samples, 2)} for r in clean[1:]}
    return peptides, samples, intensities, source_ids, removed


def figure(stats_rows, summaries, peptides):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dest = OUT / "figures"
    dest.mkdir(parents=True, exist_ok=True)

    def save(fig, name):
        for ext in ("png", "svg"):
            fig.savefig(dest / (name + "." + ext), dpi=240, bbox_inches="tight", facecolor="white")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    bars = ax.bar(GROUPS, [s["detected_peptide_species"] for s in summaries], color=COLORS)
    ax.bar_label(bars, padding=4)
    ax.set_ylim(0, 2.7)
    ax.set_ylabel("Peptide species detected (of 2)")
    ax.set_title("Yeast spike-in peptide counts by group")
    save(fig, "01_group_detected_counts")

    fig, ax = plt.subplots(figsize=(10, 5.5), layout="constrained")
    for j, (g, color) in enumerate(zip(GROUPS, COLORS)):
        values = [next(r["detection_rate_pct"] for r in stats_rows
                       if r["peptide"] == p and r["group"] == g) for p in peptides]
        xs = [i + (j - 1.5) * .19 for i in range(2)]
        b = ax.bar(xs, values, .18, color=color, label=g)
        ax.bar_label(b, fmt="%.1f%%", padding=2, fontsize=8, rotation=90)
    ax.set_xticks([0, 1], peptides)
    ax.set_ylim(0, 120)
    ax.set_ylabel("Detection rate (%)")
    ax.set_title("Per-peptide detection rates")
    ax.legend(ncol=4, bbox_to_anchor=(.5, -.10), loc="upper center", frameon=False)
    save(fig, "02_per_peptide_detection_rates")

    fig, ax = plt.subplots(figsize=(10, 5.5), layout="constrained")
    for j, (g, color) in enumerate(zip(GROUPS, COLORS)):
        for i, p in enumerate(peptides):
            r = next(r for r in stats_rows if r["peptide"] == p and r["group"] == g)
            v = r["cv_pct"]
            if v is None:
                continue
            x = i + (j - 1.5) * .17
            ax.scatter(x, v, s=95, color=color, edgecolor="white", zorder=3,
                       label=g if i == 0 else None)
            ax.annotate(f"{v:.2f}%", (x, v), xytext=(0, 9), textcoords="offset points",
                        ha="center", fontsize=8)
    ax.set_xticks([0, 1], peptides)
    ax.set_xlim(-.5, 1.5)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("Raw-intensity CV (%)")
    ax.set_title("Per-peptide CV (sample SD / mean)")
    ax.grid(axis="y", alpha=.2)
    ax.legend(ncol=4, bbox_to_anchor=(.5, -.10), loc="upper center", frameon=False)
    save(fig, "03_per_peptide_cv")

    fig, ax = plt.subplots(figsize=(9, 5.5), layout="constrained")
    for i, (g, color) in enumerate(zip(GROUPS, COLORS), 1):
        values = [r["cv_pct"] for r in stats_rows if r["group"] == g and r["cv_pct"] is not None]
        if not values:
            continue
        b = ax.boxplot([values], positions=[i], widths=.18, patch_artist=True,
                       showfliers=False, manage_ticks=False)
        for patch in b["boxes"]:
            patch.set_facecolor(color)
            patch.set_alpha(.4)
        for j, v in enumerate(values):
            ax.scatter(i + (j-(len(values)-1)/2)*.12, v, color=color, s=70,
                       edgecolor="white", zorder=3)
        ax.annotate(f"Median={stats.median(values):.2f}%", (i, max(values)),
                    xytext=(0, 12), textcoords="offset points", ha="center", fontsize=9)
    ax.set_xticks(range(1, 5), [f"{g}\n(n={summaries[i]['sample_n']})"
                                for i, g in enumerate(GROUPS)])
    ax.set_xlim(.4, 4.6)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("Raw-intensity CV (%)")
    ax.set_title("CV summary (2 peptides/group; box and points)")
    ax.grid(axis="y", alpha=.2)
    save(fig, "04_group_cv_box_points")


def main():
    labels = metadata()
    peptides, samples, intensities, original, excluded = load_and_filter(labels)
    group_samples = {g: [s for s in samples if labels[s] == g] for g in GROUPS}
    assert {g: len(v) for g, v in group_samples.items()} == {
        "Internal-QC": 6, "pooled QC": 15, "HC": 27, "S": 27}
    per_peptide, per_sample, group_summary = [], [], []
    for p in peptides:
        for g in GROUPS:
            vs = [intensities[p][s] for s in group_samples[g] if intensities[p][s] is not None]
            per_peptide.append({
                "peptide": p, "group": g, "sample_n": len(group_samples[g]),
                "detected_n": len(vs), "detection_rate_pct": 100 * len(vs) / len(group_samples[g]),
                "cv_pct": cv(vs, len(group_samples[g])),
                "mean_detected_intensity": stats.mean(vs) if vs else None,
            })
        for s in samples:
            v = intensities[p][s]
            per_sample.append({"peptide": p, "sample_id": s, "group": labels[s],
                               "raw_intensity": v, "detected": int(v is not None)})
    for g in GROUPS:
        rows = [r for r in per_peptide if r["group"] == g]
        cvs = [r["cv_pct"] for r in rows if r["cv_pct"] is not None]
        per_sample_n = [sum(intensities[p][s] is not None for p in peptides)
                        for s in group_samples[g]]
        group_summary.append({
            "group": g, "sample_n": len(group_samples[g]),
            "detected_peptide_species": sum(r["detected_n"] > 0 for r in rows),
            "mean_detected_per_sample": stats.mean(per_sample_n),
            "median_cv_pct": stats.median(cvs) if cvs else None,
            "cv_eligible_peptides": len(cvs),
        })
    OUT.mkdir(parents=True, exist_ok=True)
    save_csv(OUT / "filtered_intensities.csv", ["peptide", *samples],
             [{"peptide": p, **intensities[p]} for p in peptides])
    save_csv(OUT / "per_peptide.csv", list(per_peptide[0]), per_peptide)
    save_csv(OUT / "per_sample.csv", list(per_sample[0]), per_sample)
    save_csv(OUT / "group_summary.csv", list(group_summary[0]), group_summary)
    save_csv(OUT / "excluded_samples.csv", ["sample_id"],
             [{"sample_id": s} for s in excluded])
    figure(per_peptide, group_summary, peptides)

    lines = [
        "# Yeast spike-in peptide QC", "",
        "The raw Excel was filtered in place to metadata-listed sample IDs only.",
        f"Original sample columns: {len(original)}; retained: 75; removed: {len(excluded)}.",
        f"Two peptide sequences: {', '.join(peptides)}.",
        "No intensity normalization, log transform or imputation was performed.", "",
        "## Detection rates and CV by group", "",
        "| Peptide | Group | n | Detected | Rate | CV (%) |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for r in per_peptide:
        formatted = "NA" if r["cv_pct"] is None else f'{r["cv_pct"]:.2f}'
        lines.append(f'| {r["peptide"]} | {r["group"]} | {r["sample_n"]} | '
                     f'{r["detected_n"]} | {r["detection_rate_pct"]:.2f}% | {formatted} |')
    lines.extend([
        "", "## Method", "",
        "Detected = finite, positive raw intensity. Group rate = detected_n / all group samples * 100.",
        "CV = sample standard deviation (ddof=1) / arithmetic mean * 100 of the positive raw intensities.",
        "As in the attached R script, CV requires at least max(3, ceil(70% * group n)) observations",
        "and mean >= 1e-8. Minimum CV counts: Internal-QC 5/6, pooled QC 11/15, HC 19/27, S 19/27.",
        "Two pairs of Internal-QC columns are intentionally duplicated; nominal n=6 is not independent n=6.",
        "With only two peptides per group, boxplots are descriptive, not stable distribution estimates.",
        "", "## Figures", "",
        "Group colors: Internal-QC blue, pooled QC orange, HC green, S purple.", "",
        "![Group detection count](figures/01_group_detected_counts.png)", "",
        "![Detection rates](figures/02_per_peptide_detection_rates.png)", "",
        "![Per-peptide CV](figures/03_per_peptide_cv.png)", "",
        "![Group CV](figures/04_group_cv_box_points.png)", "",
        "Each figure is available as PNG and SVG.",
        "", "## CSV outputs", "",
        "- filtered_intensities.csv: two peptides by 75 samples",
        "- per_peptide.csv: per-peptide, four-group rates/CVs",
        "- per_sample.csv: raw peptide intensities and detection flags",
        "- group_summary.csv: four groups",
        "- excluded_samples.csv: removed sample IDs",
        "",
    ])
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"PASS: {len(peptides)} peptides, {len(samples)} samples, {len(excluded)} samples removed")
    for r in per_peptide:
        print(f'{r["peptide"]} {r["group"]}: {r["detected_n"]}/{r["sample_n"]}, '
              f'detection={r["detection_rate_pct"]:.2f}%, CV={r["cv_pct"]}')


if __name__ == "__main__":
    main()
