#!/usr/bin/env python3
"""Plot NYD 40-sample raw PG matrix by batch, no normalization or imputation."""
import csv
import json
import hashlib
import math
from pathlib import Path
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from filter_samples import readbook,tab

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"NYD_data"/"40samples"
PG=BASE/"NYD_pg_matrix2.tsv"
META=BASE/"NYD_sampleinfo2.xlsx"
OUT=BASE/"PG_batch_boxplot_results"
GROUPS={"T001":16,"T002":8,"T003":8,"T004":8}
BATCHES=list(GROUPS)
COLORS={"T001":"#4E79A7","T002":"#F28E2B","T003":"#59A14F","T004":"#E15759"}
ANNOT=("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description")

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def writecsv(name,rows):
    with (OUT/name).open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def meta():
    shared,sheets=readbook(META)
    if len(sheets)!=1:raise ValueError("Expected one metadata sheet")
    r=tab(sheets[0],shared)
    if r[0].get("A")!="样品编号" or r[0].get("C")!="batch":
        raise ValueError("Unexpected metadata headers")
    info=[]
    for row in r[1:]:
        sid=row.get("A","").strip()
        if not sid:continue
        batch=row.get("C","").strip()
        if batch not in GROUPS or sid!=row.get("B","").strip():
            raise ValueError("Metadata sample/batch mismatch "+sid)
        if row.get("D")!=batch or row.get("E")!=batch:
            raise ValueError("Metadata group columns differ from batch "+sid)
        info.append({"sample_id":sid,"batch":batch})
    if len(info)!=40 or len({x["sample_id"] for x in info})!=40:
        raise ValueError("Expected 40 unique metadata samples")
    if Counter(x["batch"] for x in info)!=GROUPS:
        raise ValueError("Unexpected batch counts")
    return info

def numeric(raw):
    try:
        v=float(raw)
        return v if math.isfinite(v) and v>0 else np.nan
    except (ValueError,TypeError):
        return np.nan

def pg_matrix(info):
    ids=[r["sample_id"] for r in info]
    with PG.open("r",encoding="utf-8-sig",newline="") as f:
        rd=csv.reader(f,delimiter="\t")
        header=next(rd)
        if len(set(header))!=len(header):raise ValueError("Duplicate PG column")
        if not set(ANNOT).issubset(header):
            raise ValueError("Required PG protein annotation fields missing "+str(header[:10]))
        if set(header)-set(ANNOT)!=set(ids):
            raise ValueError("PG intensity column names do not exactly match sampleinfo IDs; "
                             +repr(sorted(set(ids)-set(header))))
        inds=[header.index(s) for s in ids]
        raw=[]
        for lineno,row in enumerate(rd,start=2):
            if not row or all(c=="" for c in row):continue
            if len(row)!=len(header):
                raise ValueError(f"Malformed PG line {lineno}: {len(row)} cells vs {len(header)}")
            raw.append([numeric(row[j]) for j in inds])
    a=np.asarray(raw,dtype=float)
    if a.ndim!=2 or a.shape[1]!=40 or a.shape[0]<100:
        raise ValueError("Unexpected PG matrix shape "+str(a.shape))
    good=np.isfinite(a)&(a>0)
    if np.min(good.sum(axis=0))<20:
        raise ValueError("Some sample has unexpectedly few quantified protein groups")
    logged=np.full_like(a,np.nan)
    logged[good]=np.log2(a[good])
    return logged,int(a.shape[0]),int(good.sum())

def clean(x):
    arr=np.asarray(x)
    return arr[np.isfinite(arr)]

def box(ax,values,positions,colors,width=.6):
    raw=[clean(x) for x in values]
    bp=ax.boxplot(raw,positions=positions,widths=width,showfliers=False,
                  patch_artist=True,manage_ticks=False,
                  medianprops={"color":"#203345","linewidth":1.8},
                  boxprops={"color":"#758696","linewidth":.9},
                  whiskerprops={"color":"#758696","linewidth":.9},
                  capprops={"color":"#758696","linewidth":.9})
    for artist,color in zip(bp["boxes"],colors):
        artist.set_facecolor(color);artist.set_alpha(.75)

def style(ax):
    ax.grid(axis="y",color="#7F8A94",alpha=.15)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_axisbelow(True)

def save(fig,name):
    fig.savefig(OUT/name,dpi=220,bbox_inches="tight",facecolor="white")
    plt.close(fig)

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError("Refusing to overwrite an existing results directory")
    info=meta()
    log,proteins,npositive=pg_matrix(info)
    by={g:[i for i,x in enumerate(info) if x["batch"]==g] for g in BATCHES}
    rows=[]
    for j,r in enumerate(info):
        v=clean(log[:,j])
        rows.append(dict(sample_id=r["sample_id"],batch=r["batch"],
           n_protein_groups=len(log),n_positive=len(v),n_not_positive_or_missing=len(log)-len(v),
           positive_fraction=len(v)/len(log),
           median_log2_raw=float(np.median(v)),q1_log2_raw=float(np.percentile(v,25)),
           q3_log2_raw=float(np.percentile(v,75)),
           mean_log2_raw=float(np.mean(v)),sd_log2_raw=float(np.std(v,ddof=1))))
    batches=[]
    for b,ind in by.items():
        values=clean(log[:,ind].ravel())
        sm=np.array([rows[i]["median_log2_raw"] for i in ind])
        batches.append(dict(batch=b,n_samples=len(ind),n_protein_groups=proteins,
             possible_pg_intensity_values=proteins*len(ind),
             valid_positive_pg_values=len(values),
             positive_fraction=len(values)/(proteins*len(ind)),
             median_log2_all_pg_values=float(np.median(values)),
             q1_log2_all_pg_values=float(np.percentile(values,25)),
             q3_log2_all_pg_values=float(np.percentile(values,75)),
             median_of_sample_medians_log2=float(np.median(sm)),
             min_sample_median_log2=float(np.min(sm)),
             max_sample_median_log2=float(np.max(sm))))
    OUT.mkdir(parents=True,exist_ok=True)
    writecsv("sample_pg_intensity_summary.csv",rows)
    writecsv("batch_pg_intensity_summary.csv",batches)

    fig,ax=plt.subplots(figsize=(10.8,7.5),layout="constrained")
    box(ax,[log[:,by[b]].ravel() for b in BATCHES],range(1,5),
        [COLORS[b] for b in BATCHES],width=.64)
    for idx,b in enumerate(BATCHES,1):
        sm=[rows[j]["median_log2_raw"] for j in by[b]]
        ax.scatter(idx+np.linspace(-.19,.19,len(sm)),sm,
                   s=28,color="#17374A",edgecolor="white",linewidth=.4,zorder=5)
    ax.set_xticks(range(1,5),[f"{b}\n(n={GROUPS[b]} samples)" for b in BATCHES])
    ax.set_ylabel("log2 raw protein-group intensity")
    ax.set_title("NYD PG | raw intensities by batch",fontsize=16,
                 loc="left",fontweight="bold",pad=12)
    ax.text(.02,.02,"Each box pools all positive protein x sample values in one batch.\n"
            "Black dots represent the median PG intensity of individual samples.",
            fontsize=9,color="#44596A",va="bottom",transform=ax.transAxes)
    style(ax);save(fig,"01_batch_all_pg_boxplots.png")

    order=[i for b in BATCHES for i in by[b]]
    fig,ax=plt.subplots(figsize=(19.5,8.5),layout="constrained")
    points=np.arange(1,41)
    box(ax,[log[:,j] for j in order],points,
        [COLORS[info[j]["batch"]] for j in order],width=.66)
    left=1
    for b in BATCHES:
        n=GROUPS[b]
        if left>1:ax.axvline(left-.5,color="#8B9CAB",linewidth=1.3)
        ax.text(left+(n-1)/2,1.025,f"{b} (n={n})",ha="center",va="bottom",
                fontweight="bold",color=COLORS[b],fontsize=12,
                transform=ax.get_xaxis_transform())
        left+=n
    ax.set_xticks(points,[info[j]["sample_id"] for j in order],
                  rotation=90,fontsize=7)
    ax.set_xlim(.35,40.65)
    ax.set_ylabel("log2 raw protein-group intensity")
    ax.set_title("NYD PG | individual samples ordered by batch",
                 fontsize=16,loc="left",fontweight="bold",pad=28)
    style(ax);save(fig,"02_40_sample_pg_boxplots.png")

    fig,ax=plt.subplots(figsize=(10.8,7.5),layout="constrained")
    data=[[rows[j]["median_log2_raw"] for j in by[b]] for b in BATCHES]
    box(ax,data,range(1,5),[COLORS[b] for b in BATCHES],width=.57)
    for i,(b,ys) in enumerate(zip(BATCHES,data),1):
        ax.scatter(i+np.linspace(-.18,.18,len(ys)),ys,
                   s=42,color=COLORS[b],edgecolor="white",linewidth=.6,zorder=6)
    ax.set_xticks(range(1,5),[f"{b}\n(n={GROUPS[b]})" for b in BATCHES])
    ax.set_ylabel("Sample median log2 raw PG intensity")
    ax.set_title("NYD PG | sample medians by batch",fontsize=16,
                 loc="left",fontweight="bold",pad=12)
    style(ax);save(fig,"03_batch_sample_median_boxplots.png")

    manifest={"source_pg":"NYD_data/40samples/NYD_pg_matrix2.tsv",
              "source_meta":"NYD_data/40samples/NYD_sampleinfo2.xlsx",
              "pg_sha256":sha(PG),"metadata_sha256":sha(META),
              "protein_groups":proteins,"n_samples":40,"batch_sizes":GROUPS,
              "positive_intensity_values":npositive,
              "invalid_or_nonpositive_values":proteins*40-npositive,
              "scale":"log2(raw linear PG intensity), positive finite entries only",
              "no_imputation":True,"no_normalization":True,
              "figures":["01_batch_all_pg_boxplots.png",
                         "02_40_sample_pg_boxplots.png",
                         "03_batch_sample_median_boxplots.png"]}
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    readme=[
       "# Raw NYD PG protein intensity boxplots grouped by batch","",
       "The source NYD_pg_matrix2.tsv contains "+str(proteins)+" protein groups x 40 samples.",
       "Sample assignment is taken from NYD_sampleinfo2.xlsx (sample ID joins, not column position).",
       "T001 n=16; T002, T003 and T004 n=8 each.",
       "No imputation, calibration, normalization or batch correction applied.",
       "Finite positive PG intensities only are log2 transformed for display.","",
       "## 1. Protein intensity distribution by batch","",
       "Each colored box pools all valid protein x sample intensities within one batch.",
       "Dark points show individual sample medians, one per sample.",
       "Unequal batch sizes affect the number of observations in the pooled background.","",
       "![All PG by batch](01_batch_all_pg_boxplots.png)","",
       "## 2. Forty per-sample PG boxplots, ordered by batch","",
       "A single sample box summarizes its measured proteins; samples are grouped by batch.","",
       "![All 40 samples](02_40_sample_pg_boxplots.png)","",
       "## 3. Sample median distribution by batch","",
       "Each observation in this panel is a sample-level median PG intensity.",
       "Use it to compare batch location/dispersion without mistaking proteins",
       "within a sample for independent samples.","",
       "![Per-sample median by batch](03_batch_sample_median_boxplots.png)","",
       "## Descriptive summary","",
       "| Batch | Samples | Positive/possible PG values | Detection rate | All-intensity median log2 | Sample median of medians log2 |",
       "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for item in batches:
        readme.append(f'| {item["batch"]} | {item["n_samples"]} | '
             f'{item["valid_positive_pg_values"]}/{item["possible_pg_intensity_values"]} | '
             f'{item["positive_fraction"]*100:.2f}% | '
             f'{item["median_log2_all_pg_values"]:.3f} | '
             f'{item["median_of_sample_medians_log2"]:.3f} |')
    readme += ["","## Exported audit tables","",
       "- sample_pg_intensity_summary.csv: 40 sample rows, protein detection, medians and quartiles.",
       "- batch_pg_intensity_summary.csv: pooled protein values and batch medians.",
       "- run_manifest.json: original input SHA-256, matrix size, positive counts and plotting rules.",
       "- scripts/nyd_pg_batch_boxplots.py: reproducible source code.","",
       "These figures describe raw signal distributions. They do not by themselves",
       "prove the presence, absence or removal of an experimental batch effect.",""]
    (OUT/"README.md").write_text("\n".join(readme),encoding="utf-8")
    print("PG_BATCH_AUDIT",json.dumps({
        "proteins":proteins,"samples":40,"positive_values":npositive,
        "batch_sample_medians":{x["batch"]:round(x["median_of_sample_medians_log2"],4) for x in batches}
        }),flush=True)
    print("SUCCESS: 3 PNG, 2 CSV and reproducible source metadata",flush=True)

if __name__=="__main__":
    main()
