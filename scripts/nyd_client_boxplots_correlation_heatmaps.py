#!/usr/bin/env python3
"""Client-ready figures for two existing NYD 39-sample normalization results.

Input: unchanged original NYD PG matrix, original sample metadata, and the
two independently saved corrected log2 PG matrices.

Create for EACH variant:
  - 39-sample raw / normalized protein intensity boxplots (no markers)
  - 39x39 raw / normalized sample Pearson correlations (no factor values)
  - paired before/after overview images for boxplots and heatmaps.

Sample ordering follows original metadata T001, T002, T003, T004, with
previously QC-excluded T003 sample omitted. No factor calculation, no
selected peptide/protein annotation and no factor metadata in images.
For a sample-wise constant log2 shift, Pearson correlation between samples
is exactly invariant; verify ALL 39x39 pairwise coefficients accordingly.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mplcolors
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

import nyd_STGSWSTLK_heavy_exclude_T003_outlier as prev
from nyd_pg_batch_boxplots import BATCHES

BASE=prev.BASE
ROOT=BASE/"client_visualizations_39samples"
SOURCES={
    "overall_reference":{
        "source":BASE/"STGSWSTLK_Heavy_subtract_factor_39samples"/"global_subtract_fi"/"corrected_NYD_pg_matrix2_log2.tsv",
        "name":"Overall-reference normalization"
    },
    "within_batch_reference":{
        "source":BASE/"STGSWSTLK_Heavy_subtract_factor_39samples"/"batch_subtract_fi"/"corrected_NYD_pg_matrix2_log2.tsv",
        "name":"Within-batch normalization"
    }
}
GROUP_COLORS={"T001":"#6495B7","T002":"#60A58C","T003":"#AA82BE","T004":"#D6A15D"}
FILE_NAMES=(
 "01_protein_boxplots_before.png",
 "02_protein_boxplots_after.png",
 "03_protein_boxplots_before_after.png",
 "04_sample_correlation_before.png",
 "05_sample_correlation_after.png",
 "06_sample_correlation_before_after.png",
)
FIELDS=("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description")

def assert_file_fresh(path):
    if path.exists():
        raise FileExistsError("Refusing to overwrite client figures: "+str(path))

def save(fig,path,dpi=190):
    assert_file_fresh(path)
    fig.savefig(path,dpi=dpi,facecolor="white",bbox_inches="tight")
    plt.close(fig)
    if path.stat().st_size<25000:
        raise AssertionError("PNG unusually small: "+str(path))

def load_corrected(path,names,anno,raw):
    frame=pd.read_csv(path,sep="\t",encoding="utf-8-sig",
                      dtype=str,keep_default_na=False,low_memory=False)
    if len(frame)!=raw.shape[0] or list(frame.columns[:5])!=list(FIELDS):
        raise ValueError("Unexpected corrected protein source shape/annotations")
    if list(frame.columns[5:])!=names:
        raise ValueError("Corrected sample columns differ from common metadata order")
    for i,key in enumerate(FIELDS):
        if not np.array_equal(frame[key].astype(str).to_numpy(),
                              np.array([str(v[i]) for v in anno])):
            raise AssertionError("Corrected protein annotations have changed")
    matrix=frame[names].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    if matrix.shape!=raw.shape:raise AssertionError("Corrected PG source dimension")
    if not np.array_equal(np.isfinite(matrix),np.isfinite(raw)):
        raise AssertionError("Sample/protein positive-detection pattern changed")
    if not np.all(np.isfinite(matrix[np.isfinite(matrix)])):
        raise AssertionError("Invalid corrected protein log2 value")
    return matrix

def sample_order(records,names):
    if len(records)!=39 or len(names)!=39 or prev.EXCLUDED_ID in names:
        raise AssertionError("Must use 39 QC-screened samples")
    if {b:sum(r["batch"]==b for r in records) for b in BATCHES}!=prev.BATCH_COUNTS:
        raise ValueError("Wrong batch populations")
    short=[sid[-6:] for sid in names]
    if len(set(short))!=39:
        short=names
    return short

def compute_correlations(matrix,names):
    # A column is one sample and a row is one protein.
    # Pairwise complete cases (no sample-wide missing-value imputation).
    data=pd.DataFrame(matrix,columns=names)
    corr=data.corr(method="pearson",min_periods=500).to_numpy(float)
    good=np.isfinite(matrix)
    overlap=good.astype(np.int64).T@good.astype(np.int64)
    if corr.shape!=(39,39) or not np.isfinite(corr).all():
        raise AssertionError("Expected complete 39 x 39 correlations")
    if np.min(overlap)<500:
        raise ValueError(f"Too few pairwise quantified proteins: {overlap.min()}")
    if not np.allclose(corr,corr.T,atol=2e-12):
        raise AssertionError("Pearson matrix must be symmetric")
    if not np.allclose(np.diag(corr),1,atol=1e-12):
        raise AssertionError("Pearson self-correlation is not 1")
    return corr,overlap

def box_extent(mats):
    lows=[];highs=[]
    for mat in mats:
        for j in range(mat.shape[1]):
            x=mat[:,j][np.isfinite(mat[:,j])]
            q1,q3=np.percentile(x,[25,75])
            iqr=q3-q1
            left=x[x>=q1-1.5*iqr]
            right=x[x<=q3+1.5*iqr]
            lows.append(float(left.min() if len(left) else q1))
            highs.append(float(right.max() if len(right) else q3))
    low=float(min(lows));hi=float(max(highs))
    pad=max(1.0,(hi-low)*.06)
    return low-pad,hi+pad

def frame(ax):
    ax.grid(axis="y",alpha=.14,color="#90A1AC")
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def setup_box(ax,mat,records,names,short,y_range,title,show_names=True):
    position=np.arange(1,len(names)+1)
    contents=[v[np.isfinite(v)] for v in mat.T]
    boxes=ax.boxplot(contents,positions=position,widths=.66,
        patch_artist=True,showfliers=False,manage_ticks=False,
        boxprops={"edgecolor":"#607687","linewidth":.85},
        medianprops={"color":"#233D4E","linewidth":1.55},
        whiskerprops={"color":"#647989","linewidth":.78},
        capprops={"color":"#647989","linewidth":.78})
    for item,record in zip(boxes["boxes"],records):
        item.set_facecolor(GROUP_COLORS[record["batch"]])
        item.set_alpha(.58)
    edge=0
    for b in BATCHES:
        n=prev.BATCH_COUNTS[b]
        if edge:ax.axvline(edge+.5,color="#8898A5",linewidth=.9)
        ax.text(edge+(n+1)/2,1.026,f"{b}  (n={n})",
                color=GROUP_COLORS[b],fontweight="bold",
                transform=ax.get_xaxis_transform(),ha="center",
                fontsize=10,clip_on=False)
        edge+=n
    ax.set_xlim(.35,len(names)+.65)
    ax.set_ylim(*y_range)
    ax.set_xticks(position,short if show_names else [""]*len(names),
                  rotation=90,fontsize=7)
    ax.set_ylabel("Protein-group intensity (log2)")
    ax.set_title(title,loc="left",fontsize=13.5,fontweight="bold",pad=29)
    frame(ax)

def box_single(matrix,records,names,short,limit,title,path):
    fig,ax=plt.subplots(figsize=(19.5,7.3),layout="constrained")
    setup_box(ax,matrix,records,names,short,limit,title)
    save(fig,path,dpi=220)

def box_combined(raw,corrected,records,names,short,limit,path):
    fig,ax=plt.subplots(2,1,figsize=(19.5,12.8),
                        sharex=True,sharey=True,layout="constrained")
    setup_box(ax[0],raw,records,names,short,limit,
              "Before normalization",show_names=False)
    setup_box(ax[1],corrected,records,names,short,limit,
              "After normalization",show_names=True)
    ax[0].tick_params(axis="x",labelbottom=False)
    fig.suptitle("Protein intensity distribution · 39 samples",
                 fontsize=16,fontweight="bold")
    save(fig,path,dpi=190)

def heatmap_single(matrix,short,records,limits,title,path):
    fig,ax=plt.subplots(figsize=(14.5,12.8))
    fig.subplots_adjust(left=.15,right=.9,bottom=.15,top=.89)
    draw_heatmap(ax,matrix,short,records,limits,title)
    save(fig,path,dpi=190)

def draw_heatmap(ax,matrix,short,records,limits,title):
    cmap=plt.get_cmap("viridis")
    im=ax.imshow(matrix,cmap=cmap,interpolation="nearest",
                 vmin=limits[0],vmax=limits[1],aspect="equal")
    ax.set_xticks(np.arange(len(short)),short,rotation=90,fontsize=6.5)
    ax.set_yticks(np.arange(len(short)),short,fontsize=6.5)
    ax.xaxis.tick_bottom()
    ax.tick_params(length=0,pad=3)
    ax.set_xlim(-1.15,len(short)-.5)
    ax.set_ylim(len(short)-.5,-1.15)
    # Batch segmentation and small colored group strips.
    start=0
    for batch in BATCHES:
        n=prev.BATCH_COUNTS[batch]
        center=start+(n-1)/2
        ax.add_patch(Rectangle((start-.5,-1.0),n,.34,
                    fc=GROUP_COLORS[batch],ec="none",clip_on=False))
        ax.add_patch(Rectangle((-1.0,start-.5),.34,n,
                    fc=GROUP_COLORS[batch],ec="none",clip_on=False))
        ax.text(center,-1.56,batch,fontsize=9,fontweight="bold",
                color=GROUP_COLORS[batch],ha="center",va="bottom",clip_on=False)
        start+=n
        if start<len(short):
            ax.axhline(start-.5,color="#F9FBFC",lw=1.9)
            ax.axvline(start-.5,color="#F9FBFC",lw=1.9)
    ax.set_title(title,loc="left",fontsize=13.5,fontweight="bold",pad=26)
    ax.set_xlabel("Samples (batch order)")
    ax.set_ylabel("Samples (batch order)")
    colorbar=ax.figure.colorbar(im,ax=ax,fraction=.032,pad=.03,shrink=.87)
    colorbar.set_label("Pearson correlation (r)",fontsize=10)
    colorbar.ax.tick_params(labelsize=9)
    for spine in ax.spines.values():spine.set_visible(False)

def heatmap_combined(rawcorr,corr,short,records,limits,path):
    # Two panels are intended for client-friendly side-by-side comparison.
    fig,axs=plt.subplots(1,2,figsize=(25,11.5))
    fig.subplots_adjust(left=.055,right=.98,top=.88,bottom=.14,wspace=.28)
    draw_heatmap(axs[0],rawcorr,short,records,limits,"Before normalization")
    draw_heatmap(axs[1],corr,short,records,limits,"After normalization")
    fig.suptitle("Sample-to-sample protein correlations · 39 samples",
                 fontsize=16,fontweight="bold")
    save(fig,path,dpi=180)

def main():
    # Keep originals and algorithm files untouched; output only an isolated package.
    if ROOT.exists():
        raise FileExistsError("Client presentation directory exists; won't overwrite: "+str(ROOT))
    records,names,annotations,raw,targets,heavy,older_ref=prev.inputs()
    short=sample_order(records,names)
    mats={v:load_corrected(spec["source"],names,annotations,raw)
          for v,spec in SOURCES.items()}
    rawcorr,overlap=compute_correlations(raw,names)
    correlations={}
    max_diffs={}
    for variant,matrix in mats.items():
        c,o=compute_correlations(matrix,names)
        if not np.array_equal(o,overlap):
            raise AssertionError("Correlations use different pairwise protein counts")
        error=float(np.max(np.abs(c-rawcorr)))
        if error>1e-10:
            raise AssertionError(f"Pearson invariance violated for {variant}: {error}")
        correlations[variant]=c
        max_diffs[variant]=error
    boxlimits=box_extent([raw,*mats.values()])
    offdiag=rawcorr[~np.eye(len(names),dtype=bool)]
    minimum=float(offdiag.min())
    corrmin=max(-1.,np.floor(minimum*20)/20) # common 0.05-bin minimum
    if corrmin>=1.:corrmin=.95
    corrlimits=(corrmin,1.)
    ROOT.mkdir(parents=True,exist_ok=False)
    # No graphs of factors or special protein/peptide markers.
    csv_map=[{"batch":r["batch"],"sample_label":small,"sample_id":name}
             for r,small,name in zip(records,short,names)]
    pd.DataFrame(csv_map).to_csv(ROOT/"sample_label_key.csv",index=False)
    for variant,spec in SOURCES.items():
        out=ROOT/variant
        out.mkdir(parents=True,exist_ok=False)
        before=raw
        after=mats[variant]
        box_single(before,records,names,short,boxlimits,
                   "Before normalization · protein intensity distribution",
                   out/FILE_NAMES[0])
        box_single(after,records,names,short,boxlimits,
                   "After normalization · protein intensity distribution",
                   out/FILE_NAMES[1])
        box_combined(before,after,records,names,short,boxlimits,out/FILE_NAMES[2])
        heatmap_single(rawcorr,short,records,corrlimits,
                       "Before normalization · Pearson correlations",
                       out/FILE_NAMES[3])
        heatmap_single(correlations[variant],short,records,corrlimits,
                       "After normalization · Pearson correlations",
                       out/FILE_NAMES[4])
        heatmap_combined(rawcorr,correlations[variant],
                         short,records,corrlimits,out/FILE_NAMES[5])
        pd.DataFrame(rawcorr,index=names,columns=names).to_csv(
            out/"Pearson_correlation_before.csv",encoding="utf-8-sig")
        pd.DataFrame(correlations[variant],index=names,columns=names).to_csv(
            out/"Pearson_correlation_after.csv",encoding="utf-8-sig")
        if len(list(out.glob("*.png")))!=6:
            raise AssertionError("Expected six client-ready figures per variant")
    report=[
      "# NYD protein-intensity and correlation figures","",
      "## Scope","",
      "39 samples following sample-level quality control.",
      "Four acquisition batches (T001: 16, T002: 8, T003: 7, T004: 8).",
      "All available positive protein-group measurements are visualized",
      "on the log2 intensity scale. Samples remain in the same order in",
      "all panels; boxplots share one intensity scale.","",
      "## Processing variants","",
      "### Overall-reference normalization","",
      "![Overall reference protein distributions](overall_reference/03_protein_boxplots_before_after.png)","",
      "![Overall reference sample correlations](overall_reference/06_sample_correlation_before_after.png)","",
      "### Within-batch normalization","",
      "![Within-batch protein distributions](within_batch_reference/03_protein_boxplots_before_after.png)","",
      "![Within-batch sample correlations](within_batch_reference/06_sample_correlation_before_after.png)","",
      "## Notes on interpretation","",
      "Each heatmap is a 39 x 39 sample-to-sample Pearson correlation",
      "matrix across quantified log2 protein intensities. A pair of samples",
      "uses proteins detected in both (pairwise complete observations).",
      "All heatmaps share the same correlation color scale and sample ordering.",
      "**Pearson correlations are unchanged by a constant log2 shift of",
      "every protein measurement in a sample.** Therefore the before and",
      "after heatmaps are identical by design for these normalization variants.",
      "The heatmap similarity is not evidence that nothing was corrected;",
      "protein-intensity boxplots provide the appropriate shift comparison.",
      "These are presentation figures, not an independent accuracy assessment.",
      "",
      "## Figures supplied for each variant","",
      "- Protein intensity distributions: individual before and after PNGs.",
      "- Protein intensity distributions: before-and-after paired PNG.",
      "- Sample Pearson correlations: individual before and after PNGs.",
      "- Sample Pearson correlations: before-and-after paired PNG.",
      "- Supporting sample-to-sample correlation CSVs and sample-label key.",
      "",
      "All figures are free of reference-target highlights and numeric",
      "correction-factor annotations.",
      ""]
    (ROOT/"README.md").write_text("\n".join(report),encoding="utf-8")
    manifest={
      "sample_count":len(names),
      "n_protein_groups":raw.shape[0],
      "n_valid_protein_intensities":int(np.isfinite(raw).sum()),
      "batches":prev.BATCH_COUNTS,
      "boxplot_ylim":[float(v) for v in boxlimits],
      "pearson_color_limits":[float(v) for v in corrlimits],
      "min_common_protein_count":int(np.min(overlap)),
      "max_common_protein_count":int(np.max(overlap)),
      "max_abs_pearson_change":max_diffs,
      "raw_pearson_offdiagonal_min":minimum,
      "note":"Pearson sample correlation is mathematically invariant to per-sample constant log2 shifts.",
      "figures_by_variant":list(FILE_NAMES),
      "source_PG_sha256":prev.sha256(prev.PG),
      "original_processing_results_retained":True,
    }
    (ROOT/"verification.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding="utf-8")
    print("BOXPLOT_Y_LIMITS",boxlimits,flush=True)
    print("PEARSON_COLOR_LIMITS",corrlimits,flush=True)
    print("CORRELATION_MAX_ABS_CHANGE",max_diffs,flush=True)
    print("PAIRWISE_PROTEIN_MIN",int(overlap.min()),flush=True)
    print("PASS: 2 variants x (2 raw/corrected boxplot PNG + paired + 2 Pearson heatmap PNG + paired) = 12 customer PNGs",flush=True)

if __name__=="__main__":
    main()
