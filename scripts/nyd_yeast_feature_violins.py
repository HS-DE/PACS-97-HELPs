#!/usr/bin/env python3
"""NYD YEAST: correctly draw one intensity violin per peptide/protein.

For each of the 3 distinct analytes, each batch violin consists of that
analyte's measured log2 sample intensities (n=16/8/8/8). Batch CV labels
use previously verified unlogged-linear sample SD/mean x100. Never build
a violin out of the 3 unrelated CV estimates within a batch.

Style follows the user-supplied R plot function: violin as background,
boxplot centered inside, sample dots jittered, explicit batch n and top
CV text. The plots are intended to show sample distributions, not CV KDE.
"""
from pathlib import Path
import csv
import json
import hashlib
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/"NYD_data"/"40samples"/"YEAST_batch_CV_results"
FEATURES=("AADALLLK","VNQIGTLSESIK","P00924_protein")
GROUPS=("T001","T002","T003","T004")
EXPECTED={"T001":16,"T002":8,"T003":8,"T004":8}
BATCH_COLORS={"T001":"#4E79A7","T002":"#59A14F","T003":"#F28E2B","T004":"#B07AA1"}
LABELS={"AADALLLK":"AADALLLK (z1 + z2)",
        "VNQIGTLSESIK":"VNQIGTLSESIK (z2 + z3)",
        "P00924_protein":"P00924 protein"}
FIELD={"AADALLLK":"AADALLLK_sum_linear",
       "VNQIGTLSESIK":"VNQIGTLSESIK_sum_linear",
       "P00924_protein":"P00924_protein_linear"}
OUTFILES={f:"violin_"+f+".png" for f in FEATURES}
COMBINED="02_feature_specific_violin_panels.png"
OLD="02_batch_cv_violin_box_scatter.png"

def read_csv(filename,output=OUTPUT):
    with (output/filename).open("r",encoding="utf-8-sig",newline="") as fp:
        return list(csv.DictReader(fp))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def records_and_cv(samples,stats):
    if len(samples)!=40 or len(stats)!=12:
        raise ValueError("Expected 40 samples and 12 analyte-batch CV records")
    if len({s["sample_id"] for s in samples})!=40:
        raise ValueError("Duplicate sample IDs")
    for g in GROUPS:
        subset=[r for r in samples if r["batch"]==g]
        if len(subset)!=EXPECTED[g]:
            raise ValueError(f"Unexpected sample count in {g}")
    cv={}
    values={}
    for feat in FEATURES:
        for group in GROUPS:
            measured=np.asarray([float(row[FIELD[feat]])
                     for row in samples if row["batch"]==group],dtype=float)
            if len(measured)!=EXPECTED[group] or not np.all(np.isfinite(measured)&(measured>0)):
                raise ValueError(f"Nonpositive/missing {feat}/{group} values")
            match=[r for r in stats if r["batch"]==group and r["analyte"]==feat]
            if len(match)!=1:
                raise ValueError(f"Ambiguous CV row for {feat}/{group}")
            v=float(match[0]["CV_pct"])
            independently=100*np.std(measured,ddof=1)/np.mean(measured)
            if not np.isclose(v,independently,rtol=1e-10,atol=1e-10):
                raise AssertionError(f"CV from report disagrees with raw sample intensities for {feat}/{group}")
            cv[group,feat]=v
            values[group,feat]=np.log2(measured)   # only for chart distribution
    return values,cv

def style(ax,feature,values,cv,compact):
    logged=[values[g,feature] for g in GROUPS]
    all_y=np.concatenate(logged)
    delta=np.ptp(all_y)
    pad=max(delta*.09,.035)
    ymin=float(np.min(all_y)-pad)
    ymax=float(np.max(all_y)+max(delta*.37,.14))
    for group in GROUPS:
        pos=GROUPS.index(group)+1
        measured=values[group,feature]
        color=BATCH_COLORS[group]
        # violin summarizes samples of ONE analyte in ONE batch
        parts=ax.violinplot([measured],positions=[pos],widths=.74,
                            showmeans=False,showmedians=False,showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor(color);body.set_alpha(.57)
            body.set_edgecolor(color);body.set_linewidth(1.0)
        box=ax.boxplot([measured],positions=[pos],widths=.15,patch_artist=True,
                       showfliers=False,manage_ticks=False,
                       medianprops=dict(color="#0A273D",linewidth=2.0),
                       whiskerprops=dict(color="#344A5A",linewidth=1.1),
                       capprops=dict(color="#344A5A",linewidth=1.1))
        for patch in box["boxes"]:
            patch.set_facecolor("white");patch.set_alpha(.8)
            patch.set_edgecolor("#364B60")
        # deterministically spaced x offsets are actual samples, not random
        dx=np.linspace(-.16,.16,len(measured))
        ax.scatter(np.full(len(measured),pos)+dx,measured,s=20 if compact else 34,
                   facecolor=color,edgecolor="white",linewidth=.45,alpha=.94,zorder=5)
        ax.text(pos,float(np.max(all_y)+max(delta*.15,.075)),
                f"CV={cv[group,feature]:.2f}%",ha="center",va="center",
                fontsize=9 if compact else 11,fontweight="bold",color="#21364A")
    ax.set_ylim(ymin,ymax)
    ax.set_xlim(.45,4.55)
    ax.set_xticks(range(1,5),[f"{g}\n(n={EXPECTED[g]})" for g in GROUPS],
                  fontsize=9 if compact else 11)
    ax.set_ylabel("log2 summed peptide intensity" if feature!="P00924_protein"
                  else "log2 protein intensity",fontsize=10 if compact else 12)
    ax.set_title(LABELS[feature],fontsize=12 if compact else 15,
                 fontweight="bold",loc="left",pad=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y",alpha=.18,color="#536775")
    ax.set_axisbelow(True)

def generate(samples,stats,out=OUTPUT):
    out=Path(out)
    out.mkdir(parents=True,exist_ok=True)
    values,cv=records_and_cv(samples,stats)
    for feature in FEATURES:
        fig,ax=plt.subplots(figsize=(9.8,7),layout="constrained")
        style(ax,feature,values,cv,compact=False)
        fig.suptitle("NYD YEAST · raw intensity distribution by batch",
                     fontsize=16,fontweight="bold")
        fig.savefig(out/OUTFILES[feature],dpi=230,bbox_inches="tight",facecolor="white")
        plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(20.6,7.2),layout="constrained")
    for ax,f in zip(axes,FEATURES):
        style(ax,f,values,cv,compact=True)
    fig.suptitle("NYD YEAST: two unique peptides and one protein · independent batch violins",
                 fontsize=17,fontweight="bold")
    fig.savefig(out/COMBINED,dpi=220,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    return [COMBINED,*OUTFILES.values()]

def main():
    if not (OUTPUT/"yeast_per_sample_charge_sums.csv").exists():
        raise FileNotFoundError("Existing CV report not found")
    files_before={name:sha(OUTPUT/name) for name in
                  ("yeast_per_sample_charge_sums.csv",
                   "yeast_batch_cv_long.csv","yeast_batch_cv_wide.csv")}
    rows=read_csv("yeast_per_sample_charge_sums.csv")
    cvs=read_csv("yeast_batch_cv_long.csv")
    outputs=generate(rows,cvs)
    wrong=OUTPUT/OLD
    if not wrong.exists():
        raise FileNotFoundError("Expected former incorrect CV violin to replace")
    wrong.unlink()
    report=OUTPUT/"README.md"
    text=report.read_text(encoding="utf-8")
    original="""### 2. CV violin / boxplot / scatter

![Batch CV violin](02_batch_cv_violin_box_scatter.png)

Note: Each violin summarizes just three CV values per batch. Its smoothed
shape is not a statistically robust distribution estimate."""
    updated="""### 2. Corrected: one intensity violin per unique peptide / protein

Each violin is based on **sample-level log2 intensity for a single analyte in a single batch**;
the color indicates batch. Each chart overlays the narrow boxplot and every measured sample,
as in the supplied R visualization. The labels **CV=...%** are the independently verified
**linear-intensity** CV estimates; each x-axis label gives the batch sample size.
These are NOT violins of the four batch CV values or a pooled mixture of three analytes.

![Three analytes separately, each across four batches](02_feature_specific_violin_panels.png)

Individual PNG figures:

- [AADALLLK charge-summed peptide](violin_AADALLLK.png)
- [VNQIGTLSESIK charge-summed peptide](violin_VNQIGTLSESIK.png)
- [P00924 protein](violin_P00924_protein.png)"""
    if original not in text:
        raise ValueError("Original figure section not found; refused silent README rewrite")
    report.write_text(text.replace(original,updated),encoding="utf-8")
    manifest_path=OUTPUT/"run_manifest.json"
    manifest=json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["violin_definition"]="One violin per batch per analyte, from sample-level log2 intensity; annotated with CV calculated on linear sample intensities"
    manifest["figures"]=["01_batch_cv_bars_heatmap.png",COMBINED,
                         "03_batch_log2_intensity_boxplots.png",*OUTFILES.values()]
    manifest_path.write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    for name,old_digest in files_before.items():
        if sha(OUTPUT/name)!=old_digest:raise AssertionError(f"Input result changed: {name}")
    if len(list(OUTPUT.glob("*.png")))!=6:
        raise AssertionError("Exactly 6 figures expected after replacement")
    print("SUCCESS: replaced invalid 3-CV violin with 3 distinct analyte intensity violin plots")
    print("FIGURES",outputs)

if __name__=="__main__":
    main()
