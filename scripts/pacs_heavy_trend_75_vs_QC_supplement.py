#!/usr/bin/env python3
"""Add direct full-75 PG-trend Pearson rankings and QC contrast plots.

This supplement retains the QC-first ranking produced by the primary study.
All files are written ONLY under new distinct names in that completed folder;
existing tables/figures are read-only and never overwritten.
"""
from __future__ import annotations

import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pacs_heavy_97_trend_screening as previous

ROOT=previous.OUT
PLATES=previous.PLATES
COLORS=previous.GROUP_COLOR


def main():
    original=ROOT/"01_full_detect_HELP_QC_trend_ranking.csv"
    if not original.is_file():raise FileNotFoundError(original)
    rank=pd.read_csv(original)
    if len(rank)<5 or not rank.HELP_sequence.is_unique:raise AssertionError("Bad original ranking")
    allpg=rank.sort_values("overall75_allPG_r",ascending=False,kind="stable").copy()
    fixed=rank.sort_values("overall75_fixedPG_r",ascending=False,kind="stable").copy()
    qc=rank.sort_values("QC_within_plate_fixedPG_r",ascending=False,kind="stable").copy()
    for df, col, fn in (
        (allpg,"overall75_Pearson_rank","06_all75_ownPG_Pearson_rank.csv"),
        (fixed,"all75_fixedPG_Pearson_rank","07_all75_fixedPG_Pearson_rank.csv"),
        (qc,"QC_within_plate_Pearson_rank","08_QC_within_plate_Pearson_rank.csv")
    ):
        if (ROOT/fn).exists():raise FileExistsError(ROOT/fn)
        df.insert(0,col,range(1,len(df)+1))
        df.to_csv(ROOT/fn,index=False)

    pg,names,info,raw,seqs,H,valid,group,plates=previous.raw_inputs()
    if len(names)!=75 or int(valid.sum())!=len(rank):
        raise ValueError("Original source differed from ranking")
    pgmed, pgmask=previous.median_reference_pg(raw,group,plates)
    order=previous.sample_order(names,group,plates)
    qidx=np.where(group=="pooled QC")[0]
    q_order=np.array([j for j in order if group[j]=="pooled QC"],int)
    rorder={j:k for k,j in enumerate(qidx)}
    q_position=np.array([rorder[j] for j in q_order],int)
    pg_own=pgmed["all_detected_PG"]
    pg_q=previous.qc_plate_residual(pgmed["QC_common_15_PG"],plates,group=="pooled QC")
    sequences=allpg.HELP_sequence.head(5).tolist()

    fig,axs=plt.subplots(5,2,figsize=(23,18.2))
    fig.subplots_adjust(left=.065,right=.985,bottom=.11,top=.90,hspace=.43,wspace=.17)
    x=np.arange(1,76)
    xq=np.arange(1,16)
    deltaPG=pg_own-np.median(pg_own)
    for j,seq in enumerate(sequences):
        a,b=axs[j]
        hv=H[seqs.index(seq)]
        r=float(rank.set_index("HELP_sequence").loc[seq,"overall75_allPG_r"])
        rqc=float(rank.set_index("HELP_sequence").loc[seq,"QC_within_plate_fixedPG_r"])
        a.plot(x,deltaPG[order],color="#283E51",lw=1.4,marker=".",ms=3,label="Original PG median (own detected)")
        a.plot(x,(hv-np.median(hv))[order],color="#D08048",lw=1.25,marker=".",ms=3,
               label="Original HELP Heavy")
        previous.axis_groups(a,order,group,plates)
        a.axhline(0,color="#AAB7C2",lw=.85)
        a.set_title(f"{seq}: all75 Pearson r={r:+.3f}",fontsize=11)
        a.set_ylabel("Overall-median centered log2")
        a.set_xlim(.4,75.6)
        b.plot(xq,pg_q[q_position],color="#283E51",lw=1.5,marker="o",ms=3.8,
               label="QC shared-PG median")
        hvq=previous.qc_plate_residual(hv,plates,group=="pooled QC")
        b.plot(xq,hvq[q_position],color="#D08048",lw=1.35,marker="o",ms=3.8,
               label="QC HELP Heavy")
        b.axhline(0,color="#AAB7C2",lw=.85)
        for cut in (5.5,10.5):b.axvline(cut,color="#8799A8",lw=.85)
        b.set_title(f"Same {seq}: within-plate pooled-QC r={rqc:+.3f}",fontsize=11)
        b.set_ylabel("Plate-QC-median centered log2")
        b.set_xlim(.3,15.7);b.grid(axis="y",alpha=.12)
    axs[-1,0].set_xticks(x,[str(i) for i in x],rotation=90,fontsize=5.4)
    axs[-1,0].set_xlabel("75 sample index (see 03_sample_order_and_medians.csv)")
    axs[-1,1].set_xticks(xq,[f"QC{i}" for i in xq],fontsize=8)
    axs[-1,1].set_xlabel("15 QC ordered by plate (5 per plate)")
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,["PG median","HELP Heavy"],loc="lower center",
               bbox_to_anchor=(.5,.02),ncol=2,frameon=False)
    fig.suptitle("Top five Heavy HELP matching 75-sample Raw PG median | does the match hold in pooled QC?",
                 fontsize=15,fontweight="bold")
    path=ROOT/"07_top5_overall75_rawPG_trend_vs_QC_check.png"
    if path.exists():raise FileExistsError(path)
    previous.save(fig,path,dpi=180)

    meta={
       "criterion":"Five highest Pearson correlations to per-sample all detected PG log2 median across all 75",
       "all75_top5":sequences,
       "major_caveat":"Heterogeneous biology and PG missingness mean an overall trend match does not validate a technical calibration",
       "QC_check":"Fixed 15-QC-complete PG median and Heavy log2, both median-centered within each plate",
       "outputs":["06_all75_ownPG_Pearson_rank.csv","07_all75_fixedPG_Pearson_rank.csv",
                  "08_QC_within_plate_Pearson_rank.csv",path.name],
    }
    target=ROOT/"07_overall75_ranking_supplement.json"
    if target.exists():raise FileExistsError(target)
    target.write_text(json.dumps(meta,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("TOP_ALL75",json.dumps({
        "all75_top":allpg[["HELP_sequence","overall75_allPG_r",
                          "QC_within_plate_fixedPG_r",
                          "QC_within_plate_slope_PG_per_HELP"]].head(10).to_dict("records"),
        "QC_top_by_r":qc[["HELP_sequence","QC_within_plate_fixedPG_r",
                          "QC_within_plate_slope_PG_per_HELP"]].head(10).to_dict("records")
        },ensure_ascii=False))
    print("PASS overall75 rank tables and paired all-sample/QC candidate trace chart")


if __name__=="__main__":
    main()
