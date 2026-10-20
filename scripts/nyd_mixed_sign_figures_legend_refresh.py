#!/usr/bin/env python3
"""Refresh legend placement for recent NYD mixed-sign figures ONLY.

Read previously completed NYD run's audit tables and PCA coordinates;
never recompute/replace numerical results or original matrices.
Replaces overlapping legends in factor-shift, sample-median and PCA PNGs
by one figure-level legend BELOW the axes.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]/"NYD_data"/"40samples"/"STGSWSTLK_Heavy_mixed_sign_39samples"
BRANCHES=("overall_reference","within_batch_reference")
BATCHES=("T001","T002","T003","T004")
COLORS={"T001":"#6B9CBF","T002":"#67A68D","T003":"#B38CC3","T004":"#D2A36A"}


def save(fig,path):
    fig.savefig(path,dpi=195,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    if path.stat().st_size<25000:raise AssertionError(f"Small figure {path.name}")


def figure_legend(fig,ax,ncol=4):
    handles,labels=ax.get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.5,-.035),
               ncol=ncol,fontsize=9,frameon=False)


def plot_shifts(folder,records,branch):
    fig,ax=plt.subplots(figsize=(18,6),layout="constrained")
    x=np.arange(1,40)
    for b in BATCHES:
        index=np.array([s==b for s in records.batch])
        ax.bar(x[index],records.loc[index,"actual_PG_applied_shift_log2"],
               color=COLORS[b],alpha=.85,width=.78,
               label=f"{b} ({'-fi' if b!='T004' else '+fi'})")
    for cut in (16.5,24.5,31.5):ax.axvline(cut,color="#94A4B2",lw=.8)
    ax.axhline(0,color="#405564",lw=.9)
    ax.set_xlim(.1,39.9)
    ax.set_xticks(x,[str(s)[-6:] for s in records.sample_id],rotation=90,fontsize=7.1)
    ax.set_ylabel("Applied per-sample log2 shift")
    ax.set_title(f"{branch}: actual PG shifts (mixed Heavy factor signs)",pad=13)
    ax.grid(axis="y",alpha=.15)
    figure_legend(fig,ax)
    save(fig,folder/"04_sample_applied_log2_shifts.png")


def plot_medians(folder,records,branch):
    fig,ax=plt.subplots(figsize=(18,6),layout="constrained")
    x=np.arange(1,40)
    raw=records.Raw_PG_sample_median_log2.to_numpy(float)
    corrected=records.Corrected_PG_sample_median_log2.to_numpy(float)
    ax.plot(x,raw,c="#697E91",lw=1.5,marker=".",markersize=5,label="Raw PG median")
    ax.plot(x,corrected,c="#2F8F91",lw=1.6,marker=".",markersize=5,label="Corrected PG median")
    for b in BATCHES:
        idx=np.array(records.batch==b)
        ax.scatter(x[idx],corrected[idx],c=COLORS[b],s=15,alpha=.8)
    for cut in (16.5,24.5,31.5):ax.axvline(cut,color="#94A4B2",lw=.8)
    ax.set_xlim(.1,39.9)
    ax.set_xticks(x,[str(s)[-6:] for s in records.sample_id],rotation=90,fontsize=7.1)
    ax.set_ylabel("Per-sample median detected PG intensity (log2)")
    ax.set_title(f"{branch}: median PG distribution before / after calibration",pad=13)
    ax.grid(axis="y",alpha=.14)
    figure_legend(fig,ax,ncol=2)
    save(fig,folder/"06_sample_PG_medians_Raw_vs_corrected.png")


def plot_pca(folder,points,bounds,raw_meta,branch):
    stages=(
        ("Raw_PC1","Raw_PC2","Raw PCA · fitted only on original data","04a_PCA_Raw.png"),
        ("Corrected_projected_PC1","Corrected_projected_PC2",
         "Corrected PCA · projected into Raw PCs","04b_PCA_Corrected.png")
    )
    for xkey,ykey,title,filename in stages:
        fig,ax=plt.subplots(figsize=(9.4,7.3),layout="constrained")
        for b in BATCHES:
            sub=points[points.batch==b]
            ax.scatter(sub[xkey],sub[ykey],s=50,alpha=.85,
                       edgecolor="white",lw=.6,color=COLORS[b],
                       label=f"{b} (n={len(sub)})")
        ax.set_xlim(bounds[0],bounds[1]);ax.set_ylim(bounds[2],bounds[3])
        ax.set_xlabel(f"Raw PC1 ({raw_meta['PC1_raw_explained_variance_pct']:.2f}%)")
        ax.set_ylabel(f"Raw PC2 ({raw_meta['PC2_raw_explained_variance_pct']:.2f}%)")
        ax.set_title(f"{title} | {branch}",fontsize=12,pad=12)
        ax.axhline(0,c="#B2BFC9",lw=.8);ax.axvline(0,c="#B2BFC9",lw=.8)
        ax.grid(alpha=.12)
        figure_legend(fig,ax)
        save(fig,folder/filename)


def main():
    manifest=json.loads((ROOT/"run_manifest.json").read_text())
    assert manifest["n_retained"]==39 and manifest["n_protein_groups"]==4264
    bounds=None
    for branch in BRANCHES:
        folder=ROOT/branch
        factors=pd.read_csv(folder/"sample_factors_internal.csv")
        points=pd.read_csv(folder/"PCA_Raw_fitted_sample_coordinates.csv")
        meta=json.loads((folder/"reference_and_processing.json").read_text())
        if len(factors)!=39 or len(points)!=39:
            raise ValueError(f"Incomplete 39-sample tables: {branch}")
        if list(points.sample_id)!=list(factors.sample_id):
            raise ValueError(f"Sample IDs not aligned: {branch}")
        branch_bounds=np.array(meta["PCA_common_axis_bounds"],float)
        if bounds is None:bounds=branch_bounds
        else:
            if not np.allclose(bounds,branch_bounds,rtol=0,atol=1e-12):
                raise AssertionError("PCA axes must match across variants")
        plot_shifts(folder,factors,branch)
        plot_medians(folder,factors,branch)
        plot_pca(folder,points,bounds,meta["PCA_raw_fit"],branch)
    print("PASS: 8 NYD mixed-sign PNGs refreshed with figure legends outside axes; no numerical files changed.")

if __name__=="__main__":
    main()
