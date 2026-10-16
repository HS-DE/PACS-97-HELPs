#!/usr/bin/env python3
"""Violin+box CV distributions with median labels, NYD 39-sample subtract-factor trials."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE=Path(__file__).resolve().parents[1]/"NYD_data"/"40samples"/"STGSWSTLK_Heavy_subtract_factor_39samples"
BRANCHES={"global_subtract_fi":"Global reference -fi",
          "batch_subtract_fi":"Batch-internal reference -fi"}
BATCHES=("T001","T002","T003","T004")
N={"T001":16,"T002":8,"T003":7,"T004":8}
COLORS={"T001":"#669AC4","T002":"#63A68F","T003":"#C08BC5","T004":"#D2A36A"}
RAW="#8997A3"
PNG="08_per_protein_CV_violin_Raw_vs_Corrected_median_labels.png"
CSV="08_violin_median_CV_values.csv"

def data(folder):
    a=pd.read_csv(folder/"per_protein_batch_CV_long.csv",encoding="utf-8-sig")
    s=pd.read_csv(folder/"batch_protein_CV_summary.csv",encoding="utf-8-sig")
    if len(a)!=2*4*4264 or len(s)!=8:
        raise AssertionError("Unexpected existing CV export dimensions")
    if set(a.stage)!=set(("Raw","Corrected")) or set(a.batch)!=set(BATCHES):
        raise AssertionError("Unexpected CV stage or batch")
    arrays={}
    labels=[]
    for b in BATCHES:
        for stage in ("Raw","Corrected"):
            sub=a[(a.batch==b)&(a.stage==stage)].sort_values("protein_row_1based")
            if len(sub)!=4264 or not np.array_equal(sub.protein_row_1based.to_numpy(),np.arange(1,4265)):
                raise AssertionError("Protein row identity mismatch")
            allcv=sub.cv_pct.to_numpy(float)
            good=allcv[np.isfinite(allcv)]
            row=s[(s.batch==b)&(s.stage==stage)].iloc[0]
            median=float(np.median(good))
            if len(good)!=int(row.n_proteins_with_valid_cv) or not np.isclose(median,row.median_protein_cv_pct,atol=1e-10):
                raise AssertionError("CV median mismatches source")
            for q,field in ((25,"q1_protein_cv_pct"),(75,"q3_protein_cv_pct")):
                if not np.isclose(np.percentile(good,q),row[field],atol=1e-10):
                    raise AssertionError("CV quartile mismatch")
            arrays[(stage,b)]=allcv
            labels.append({"batch":b,"stage":stage,"n_samples":N[b],
                           "n_proteins_with_valid_CV":len(good),
                           "median_CV_pct":median,
                           "q1_CV_pct":float(np.percentile(good,25)),
                           "q3_CV_pct":float(np.percentile(good,75))})
    return arrays,labels

def plot(folder,method,arr,labels):
    groups=[v[np.isfinite(v)] for v in arr.values()]
    limit=float(np.percentile(np.concatenate(groups),99))*1.10
    fig,ax=plt.subplots(figsize=(14.8,7.8))
    fig.subplots_adjust(left=.09,right=.98,bottom=.16,top=.77)
    lookup={(r["stage"],r["batch"]):r for r in labels}
    xpos=np.arange(1,5,dtype=float)
    for stage,offset in (("Raw",-.19),("Corrected",.19)):
        xx=xpos+offset
        data=[arr[(stage,b)][np.isfinite(arr[(stage,b)])] for b in BATCHES]
        violin=ax.violinplot(data,positions=xx,widths=.34,
                             showmeans=False,showmedians=False,
                             showextrema=False,points=140)
        for patch,b in zip(violin["bodies"],BATCHES):
            patch.set_facecolor(RAW if stage=="Raw" else COLORS[b])
            patch.set_edgecolor("#5C6D7A")
            patch.set_linewidth(.85)
            patch.set_alpha(.32 if stage=="Raw" else .52)
        boxes=ax.boxplot(data,positions=xx,widths=.13,
                         patch_artist=True,showfliers=False,manage_ticks=False,
                         medianprops={"color":"#1A3548","linewidth":1.7},
                         boxprops={"color":"#4F6170","linewidth":.8},
                         whiskerprops={"color":"#647785","linewidth":.8},
                         capprops={"color":"#647785","linewidth":.8})
        for patch,b in zip(boxes["boxes"],BATCHES):
            patch.set_facecolor(RAW if stage=="Raw" else COLORS[b])
            patch.set_alpha(.65)
        for x,b in zip(xx,BATCHES):
            median=lookup[(stage,b)]["median_CV_pct"]
            ax.text(x,.955,f"{median:.2f}%",
                    transform=ax.get_xaxis_transform(),ha="center",va="top",
                    fontsize=10.2,fontweight="bold",
                    color="#3E4D59" if stage=="Raw" else "#173D57",
                    bbox={"facecolor":"white","edgecolor":"none","alpha":.83,"pad":1.8})
    ax.set_ylim(0,limit)
    ax.set_xlim(.5,4.5)
    ax.set_xticks(xpos,[f"{b}\nn={N[b]}" for b in BATCHES],fontsize=11)
    ax.set_ylabel("Per-protein CV (%) on linear intensities",fontsize=11)
    ax.set_title("Raw (left) vs Corrected (right) · violin and embedded boxplot",
                 loc="left",fontweight="bold",fontsize=13.5,pad=9)
    ax.grid(axis="y",alpha=.14)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(.01,.995,"Median CV (%) above each violin",transform=ax.transAxes,
            fontsize=10,color="#4C6170",va="bottom")
    fig.suptitle(f"NYD STGSWSTLK Heavy · {method}\n39 retained samples; T003 outlier excluded",
                 fontsize=15.5,fontweight="bold",y=.992)
    fig.legend(handles=[
        Patch(facecolor=RAW,edgecolor="#5C6D7A",alpha=.6,label="Raw (left)"),
        Patch(facecolor="#76A8B7",edgecolor="#5C6D7A",alpha=.65,
              label="Corrected (right, colored by batch)")],
        loc="upper center",bbox_to_anchor=(.5,.855),ncol=2,frameon=False)
    fig.text(.5,.055,
             "Display capped at pooled 99th percentile × 1.10; all valid CVs retained for calculations.",
             ha="center",fontsize=9,color="#5E707B")
    file=folder/PNG
    fig.savefig(file,dpi=210,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    if file.stat().st_size<30000:raise AssertionError("Violin image unexpectedly small")
    for r in labels:r["plot_y_max_CV_pct"]=limit
    with (folder/CSV).open("w",encoding="utf-8-sig",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(labels[0]))
        writer.writeheader()
        writer.writerows(labels)
    return limit

def main():
    for branch in BRANCHES:
        folder=BASE/branch
        if not folder.is_dir():raise FileNotFoundError(folder)
        if (folder/PNG).exists() or (folder/CSV).exists():
            raise FileExistsError("Refuse to replace existing files in "+str(folder))
    x={b:data(BASE/b) for b in BRANCHES}
    for batch in BATCHES:
        for stage in ("Raw","Corrected"):
            global_=x["global_subtract_fi"][0][(stage,batch)]
            perbatch=x["batch_subtract_fi"][0][(stage,batch)]
            if not np.allclose(global_,perbatch,atol=1e-10,rtol=0,equal_nan=True):
                raise AssertionError("Global/batch subtract factor CV should match per-protein")
    for branch,method in BRANCHES.items():
        arr,labels=x[branch]
        limit=plot(BASE/branch,method,arr,labels)
        print(branch,"y-limit",limit,
              "medians",[(r["batch"],r["stage"],round(r["median_CV_pct"],2)) for r in labels],
              flush=True)
    readme=BASE/"README.md"
    source=readme.read_text(encoding="utf-8")
    section="## Added violin plots with median CV labels"
    if section in source:raise FileExistsError("README section already exists")
    lines=[section,"",
         "Following the original 39-sample single-Heavy CV figure RIGHT-hand violin panel:",
         "paired Raw (left) and Corrected (right) violins with embedded boxplots,",
         "each with its median per-protein CV (%) labeled at the very top.",
         "Axes display up to the pooled 99th percentile x 1.10, matching",
         "the prior style. All finite protein CV values are retained for",
         "distributions, medians, and quartiles.","",
         "### Global -fi violin plot","",
         f"![Global minus factor violin](global_subtract_fi/{PNG})","",
         "### Batch-internal -fi violin plot","",
         f"![Batch-local minus factor violin](batch_subtract_fi/{PNG})","",
         "Both methods yield the same within-batch CV distribution, because",
         "their factors differ only by a constant within each batch.",
         "Corresponding 08_violin_median_CV_values.csv files record the plotted medians.",
         ""]
    readme.write_text(source.rstrip()+"\n\n"+"\n".join(lines),encoding="utf-8")
    print("PASS: two independent CV violin PNG, eight median annotations each, README updated",flush=True)

if __name__=="__main__":
    main()
