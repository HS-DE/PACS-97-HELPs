#!/usr/bin/env python3
"""PACS 75 samples: user-selected five HELP Heavy vs per-sample PG median.

Use raw original PG and the 'H' sheet of original 2026-10-13 HELP workbook.
Pearson across 75 sample observations and within each of three plates.
Compare per-sample ALL detected PG versus the same 75-sample-common PG set.
Stratify by sample class as a supplementary sensitivity check; no normalization.
"""
from __future__ import annotations
import csv, hashlib, json
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pacs_two_reference_plate_trial as old
import pacs_help5_heavy_only_trial as help5
from yeast_protein_qc import GROUPS

OUT=old.ROOT/"results"/"pacs_selected5_HELP_Heavy_PG_Pearson_75samples"
PEPTIDES=("FSVVYAK","STGSWSTLK","VTGWGNLK","VVEESELAR","VYFAGFPR")
PLATES=("1","2","3")
PLOT_COLORS={"1":"#2B70B1","2":"#1D947A","3":"#CC8050"}
SCOPES=("ALL_75","Plate_1","Plate_2","Plate_3","HC","S","pooled_QC","Internal_QC")


def csvwrite(path, rows):
    if not rows:raise ValueError(f"Empty file {path}")
    with path.open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def corr(x,y):
    ok=np.isfinite(x)&np.isfinite(y)
    a,b=x[ok],y[ok]
    if len(a)<3 or np.std(a)<=1e-12 or np.std(b)<=1e-12:
        return int(len(a)),None,None
    result=pearsonr(a,b)
    return int(len(a)),float(result.statistic),float(result.pvalue)


def plot_scatter(heavy, medians, info, samples, results):
    fig,axs=plt.subplots(5,2,figsize=(13.4,19.4),layout="constrained")
    for j,seq in enumerate(PEPTIDES):
        for k,(basis,y) in enumerate(medians.items()):
            ax=axs[j,k]
            for p in PLATES:
                take=np.array([info[s]["plate"]==p for s in samples])
                x=heavy[j,take]; yy=y[take]
                finite=np.isfinite(x)&np.isfinite(yy)
                ax.scatter(x[finite],yy[finite],s=29,alpha=.69,
                           color=PLOT_COLORS[p],edgecolor="white",linewidths=.4,
                           label=f"Plate {p}")
            r=next(v["pearson_r"] for v in results if v["scope"]=="ALL_75"
                   and v["HELP"]==seq and v["PG_basis"]==basis)
            ax.set_title(f"{seq} | {basis} | all-75 r={r:+.3f}",fontsize=10.1,pad=8)
            ax.set_xlabel("Original HELP Heavy log2 intensity")
            ax.set_ylabel("Sample PG median log2")
            ax.grid(alpha=.13)
            ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.5,-.005),
               ncol=3,frameon=False,fontsize=10)
    fig.suptitle("PACS 75 samples | manually selected five HELP Heavy vs raw PG median",
                 fontsize=14,fontweight="bold")
    fig.savefig(OUT/"01_scatter_selected5_Heavy_vs_PG_medians.png",dpi=170,
                bbox_inches="tight",facecolor="white")
    plt.close(fig)


def plot_heatmaps(results):
    scopes=SCOPES
    fig,axs=plt.subplots(1,2,figsize=(14.5,7.7),layout="constrained")
    im=None
    for ax,basis in zip(axs,("all_detected_PG","common_75_PG")):
        matrix=np.full((len(scopes),len(PEPTIDES)),np.nan)
        for row,scope in enumerate(scopes):
            for col,seq in enumerate(PEPTIDES):
                val=next(r["pearson_r"] for r in results if r["scope"]==scope
                         and r["HELP"]==seq and r["PG_basis"]==basis)
                if val is not None:matrix[row,col]=val
        im=ax.imshow(matrix,vmin=-1,vmax=1,cmap="RdBu_r",aspect="auto")
        ax.set_yticks(range(len(scopes)),scopes)
        ax.set_xticks(range(len(PEPTIDES)),PEPTIDES,rotation=38,ha="right")
        ax.set_title("Sample PG median: "+basis,fontsize=11,pad=13)
        for i in range(len(scopes)):
            for j in range(len(PEPTIDES)):
                if np.isfinite(matrix[i,j]):
                    ax.text(j,i,f"{matrix[i,j]:+.2f}",ha="center",va="center",
                            color="white" if abs(matrix[i,j])>.67 else "#213342",
                            fontsize=10,fontweight="bold")
    fig.colorbar(im,ax=list(axs),location="right",fraction=.025,pad=.03,label="Pearson r")
    fig.suptitle("PACS five HELP Heavy correlations | overall, by plate, by group",
                 fontsize=14,fontweight="bold")
    fig.savefig(OUT/"02_Pearson_r_overall_plate_group_heatmaps.png",
                dpi=190,bbox_inches="tight",facecolor="white")
    plt.close(fig)


def main():
    if OUT.exists():raise FileExistsError(f"Result directory already exists: {OUT}")
    pg,samples,info,raw,_,_=old.load()
    if list(pg.columns[:5])!=old.ANN or raw.shape!=(6675,75):
        raise AssertionError("Original PACS protein matrix does not match handoff")
    if tuple(help5.SEQS)!=PEPTIDES:
        raise AssertionError("Manual HELP list differs from previous selected-five trial")
    sheets,_=help5.sheet_data(samples)
    heavy=sheets["H"]
    if heavy.shape!=(5,75) or not np.isfinite(heavy).all():
        raise AssertionError("Five manually chosen HELP Heavy must exist in 75 samples")
    detection=np.isfinite(raw)
    shared=detection.all(axis=1)
    n_shared=int(np.count_nonzero(shared))
    if n_shared<100:raise AssertionError(f"Too few common PG: {n_shared}")
    means={
      "all_detected_PG":np.nanmedian(raw,axis=0),
      "common_75_PG":np.median(raw[shared],axis=0)
    }
    if not all(np.isfinite(x).all() for x in means.values()):
        raise AssertionError("Invalid PG medians")
    ref=np.median(heavy,axis=1)
    composite=(heavy-ref[:,None]).mean(axis=0)
    feature={p:heavy[j] for j,p in enumerate(PEPTIDES)}
    feature["5HELP_mean_centered_log2"]=composite
    sample_details=[]
    for j,s in enumerate(samples):
        row={"sample_id":s,"group":info[s]["group"],"plate":info[s]["plate"],
             "run":info[s]["run"],
             "n_detected_PG":int(detection[:,j].sum()),
             "n_common_75_PG":n_shared,
             "all_detected_PG_log2_median":float(means["all_detected_PG"][j]),
             "common_75_PG_log2_median":float(means["common_75_PG"][j])}
        for k,p in enumerate(PEPTIDES):
            row[p+"_Heavy_log2"]=float(heavy[k,j])
            row[p+"_Heavy_linear"]=float(np.exp2(heavy[k,j]))
        row["5HELP_mean_centered_log2"]=float(composite[j])
        sample_details.append(row)
    correlations=[]
    def choose(scope):
        if scope=="ALL_75":return np.ones(75,dtype=bool)
        if scope.startswith("Plate_"):
            return np.array([info[s]["plate"]==scope.split("_")[1] for s in samples])
        g={"pooled_QC":"pooled QC","Internal_QC":"Internal-QC"}.get(scope,scope)
        return np.array([info[s]["group"]==g for s in samples])
    for scope in SCOPES:
        mask=choose(scope)
        assert mask.sum()>=6
        for pep,x in feature.items():
            for basis,y in means.items():
                n,r,p=corr(x[mask],y[mask])
                correlations.append({"scope":scope,
                   "scope_type":"overall" if scope=="ALL_75" else
                                "plate" if scope.startswith("Plate_") else "sample_group",
                   "n_samples_scope":int(mask.sum()),"HELP":pep,"PG_basis":basis,
                   "n_pairs":n,"pearson_r":r,"p_value_two_sided":p})
    # Correlation between detected count and PG median as a missingness diagnostic
    qc_composition=[]
    for scope in SCOPES:
        mask=choose(scope)
        for basis,y in means.items():
            n,r,p=corr(detection.sum(axis=0)[mask].astype(float),y[mask])
            qc_composition.append({"scope":scope,"PG_basis":basis,"n_samples_scope":int(mask.sum()),
                                   "n_pairs":n,"pearson_r_detected_count_vs_PG_median":r,
                                   "p_value_two_sided":p})
    OUT.mkdir(parents=True)
    csvwrite(OUT/"pearson_summary_overall_plate_group.csv",correlations)
    csvwrite(OUT/"sample_75_Heavy_and_PG_log2_medians.csv",sample_details)
    csvwrite(OUT/"PG_detected_count_vs_medians.csv",qc_composition)
    annotations=[{"source_PG_row_number":i+2, **dict(zip(old.ANN,pg.iloc[i][old.ANN].tolist()))}
                 for i in np.where(shared)[0]]
    csvwrite(OUT/"common_75_PG_protein_list.csv",annotations)
    plot_scatter(heavy,means,info,samples,correlations)
    plot_heatmaps(correlations)
    manifest={
        "dataset":"PACS main, not NYD",
        "n_samples":75,"n_raw_protein_groups":6675,
        "n_common_75_PG":n_shared,
        "n_valid_PG_by_sample_range":[int(detection.sum(axis=0).min()),
                                      int(detection.sum(axis=0).max())],
        "manual_5_HELP":list(PEPTIDES),
        "source_Heavy":"Original H sheet only, not L/H ratio",
        "raw_log2_PG":"Log2 of each positive finite original PG intensity",
        "PG_summary":"Samplewise median of log2 PG intensity",
        "all_detected_PG":"Own sample positive finite PGs, different sets possible",
        "common_75_PG":"Fixed 75-sample intersection, positive finite in all 75 samples",
        "scope_counts":{s:int(choose(s).sum()) for s in SCOPES},
        "extra_group_scopes":"HC/S/pooled QC/Internal-QC separate diagnostic, no cross-group correlation interpreted as technical precision",
        "statistics":"Pearson across samples, no imputation; two-sided unadjusted p-values",
        "composite_HELP":"Mean of five peptide-wise QC-independent centered log2 Heavy signals",
        "source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                         (old.PG,help5.HELP,old.META_QC,old.META_STUDY)}
    }
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    readme=[
        "# PACS 75 samples | user-selected five HELP Heavy vs sample raw PG medians",
        "",
        "Uses 75 original samples and these five MANUALLY chosen HELP peptide sequences:",
        ", ".join(PEPTIDES),". Heavy intensity is from H sheet, NOT ratio_LH.",
        "No normalization/factor is applied here. Study, pooled-QC and Internal-QC",
        "samples are all included. This is not NYD; biological sample heterogeneity matters.",
        "",
        "## Matched protein sets",
        "",
        "1. all_detected_PG: each sample independently uses all its own valid PGs.",
        "2. common_75_PG: the exact same rows valid across all 75 samples.",
        f"Common set size is {n_shared} out of 6675 original PGs.",
        "",
        "Both sample summaries take median(log2(original linear PG intensity)).",
        "Pearson is across sample observations, overall and within each plate.",
        "Group-stratified results are supplementary; unequal group composition",
        "across plates is a potential confounder of overall/plate correlation.",
        "The 5HELP composite is the mean of the five peptide-wise centered Heavy log2",
        "intensities; it is not a correction factor and is included only as an extra diagnostic.",
        "No missing-value imputation; p-values unadjusted.",
        "",
        "## Overall Pearson r across 75 samples",
        "",
        "| HELP | All detected PG median | Fixed 75-common PG median |",
        "| --- | ---: | ---: |"
    ]
    for pep in list(PEPTIDES)+["5HELP_mean_centered_log2"]:
        vals={x["PG_basis"]:x["pearson_r"] for x in correlations
              if x["scope"]=="ALL_75" and x["HELP"]==pep}
        readme.append(f"| {pep} | {vals['all_detected_PG']:+.3f} | {vals['common_75_PG']:+.3f} |")
    readme.extend([
       "", "## By-plate Pearson r",
       "", "Each cell displays all_detected_PG / common_75_PG.",
       "", "| HELP | Plate 1 (n=27) | Plate 2 (n=27) | Plate 3 (n=21) |",
       "| --- | ---: | ---: | ---: |"
    ])
    for pep in list(PEPTIDES)+["5HELP_mean_centered_log2"]:
        vals=[]
        for plate in PLATES:
            d={x["PG_basis"]:x["pearson_r"] for x in correlations
               if x["scope"]=="Plate_"+plate and x["HELP"]==pep}
            vals.append(f"{d['all_detected_PG']:+.3f} / {d['common_75_PG']:+.3f}")
        readme.append("| "+pep+" | "+" | ".join(vals)+" |")
    readme.extend([
       "",
       "## Figures",
       "![Five HELP scatters](01_scatter_selected5_Heavy_vs_PG_medians.png)",
       "",
       "![Correlation heatmaps](02_Pearson_r_overall_plate_group_heatmaps.png)",
       "",
       "## Files",
       "- pearson_summary_overall_plate_group.csv: overall, 3 plates, 4 sample groups.",
       "- sample_75_Heavy_and_PG_log2_medians.csv: all 75 records with source group/plate.",
       "- PG_detected_count_vs_medians.csv: diagnostic of missingness/composition impact.",
       "- common_75_PG_protein_list.csv: exact common protein set, with annotations.",
       "- run_manifest.json: precise source hashes and methodology.",
       "",
       "## Interpretation",
       "HC and S are biologically different samples. Their combined correlations",
       "are descriptive only; they are not technical repeatability metrics.",
       "All charts place their legend OUTSIDE the plotting axes.",
       "Correlation does not demonstrate that either +fi or -fi improves data.",
    ])
    (OUT/"README.md").write_text("\n".join(readme)+"\n",encoding="utf-8")
    print("SUMMARY",json.dumps({"common":n_shared,"counts":manifest["scope_counts"],
                 "overall":[x for x in correlations if x["scope"]=="ALL_75"],
                 "plate":[x for x in correlations if x["scope"].startswith("Plate_")]},ensure_ascii=False))
    print("PASS PACS five-HELP Pearson 75 samples, 3 plates, 4 study/QC groups, two PG sets.")


if __name__=="__main__":
    main()
