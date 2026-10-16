#!/usr/bin/env python3
"""NYD 39-sample reverse-direction STGSWSTLK Heavy tests.

Keep the exact 39-sample global and per-batch factors f_i = reference - log2(H_i)
used by the prior +f experiments; now apply PG_corrected = log2(PG_raw) - f_i.
Two independent variants:
  global_subtract_fi : one median across all 39 retained Heavy observations;
  batch_subtract_fi  : own-batch median across 16/8/7/8 Heavy observations.
Exclude previously confirmed bad T003 sample before medians, CV, PCA, exports.
DO NOT re-normalize Heavy.xlsx, cap factors, or modify previous results.
PCA fitted on Raw ONCE and both corrected variants projected to the same
Raw-only feature selection, imputation medians, centering and PC loadings.
Raw and Corrected PCA are separate PNG files WITHOUT arrows.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import nyd_STGSWSTLK_heavy_exclude_T003_outlier as prev
import nyd_STGSWSTLK_heavy_batch_internal_trial as batch
from nyd_pg_batch_boxplots import BATCHES, COLORS, ANNOT

ROOT=prev.BASE
OUT=ROOT/"STGSWSTLK_Heavy_subtract_factor_39samples"
PLUS_GLOBAL=ROOT/"STGSWSTLK_Heavy_single_reference_normalization_without_T003_outlier"
PLUS_BATCH=ROOT/"STGSWSTLK_Heavy_batch_internal_normalization_without_T003_outlier"
VARIANTS={"global_subtract_fi":"Global 39-sample median",
          "batch_subtract_fi":"Batch-internal median"}
PNGS=[
    "01_PG_Raw_vs_subtract_factor_blue_CFB_yellow_PROS1.png",
    "02_Heavy_reference_and_reverse_sample_factors.png",
    "03_batch_protein_CV_Raw_vs_subtract_factor.png",
    "04a_PCA_Raw.png",
    "04b_PCA_Corrected.png",
    "05_CFB_PROS1_by_batch_Raw_vs_corrected.png",
    "06_sample_PG_medians_Raw_vs_corrected.png"
]

def store_csv(path,rows):
    if not rows:raise ValueError("Empty table "+str(path))
    with path.open("w",encoding="utf-8-sig",newline="") as handle:
        w=csv.DictWriter(handle,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def pg_matrix(path,ann,names,log):
    with path.open("w",encoding="utf-8",newline="") as handle:
        w=csv.writer(handle,delimiter="\t")
        w.writerow(list(ANNOT)+names)
        for i,row in enumerate(ann):
            w.writerow(list(row)+[float(v) if np.isfinite(v) else "" for v in log[i]])

def savefig(fig,path):
    fig.savefig(path,dpi=190,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    if path.stat().st_size<15000:raise ValueError("Figure may be empty: "+str(path))

def pg_boxes(raw,corrected,records,names,targets,dir_,variant):
    fig,axs=plt.subplots(2,1,figsize=(19.5,14),sharex=True,sharey=True,layout="constrained")
    prev.boxplot(axs[0],raw,records,names,targets,"Raw · 4,264 protein groups",full_x=False)
    prev.boxplot(axs[1],corrected,records,names,targets,
        "Corrected = original protein log2 intensity MINUS fi · "+variant)
    fig.suptitle("NYD 39 samples · blue CFB/P00751 · yellow PROS1/P07225",
                 fontsize=15,fontweight="bold")
    savefig(fig,dir_/PNGS[0])

def factor_plot(h,ref_per_sample,f,names,records,dir_,variant,ref_by_batch=None):
    pos=np.arange(1,len(names)+1)
    fig,axs=plt.subplots(2,1,figsize=(18,9.4),sharex=True,layout="constrained")
    for b in BATCHES:
        idx=np.array([j for j,r in enumerate(records) if r["batch"]==b])
        axs[0].scatter(pos[idx],h[idx],s=48,c=COLORS[b],
                       edgecolors="white",lw=.5,label=b)
        axs[1].bar(pos[idx],-f[idx],color=COLORS[b],width=.73,alpha=.83)
        if ref_by_batch:
            axs[0].plot([pos[idx].min()-.45,pos[idx].max()+.45],
                        [ref_by_batch[b]]*2,c=COLORS[b],ls="--",lw=1.5)
    if not ref_by_batch:
        axs[0].axhline(float(ref_per_sample[0]),c="#224D86",lw=1.5,ls="--",
                       label="Global median ref")
    # Draw applied sign ONLY: raw H - f; when f=ref-H this moves AWAY.
    axs[0].scatter(pos,h-f,c="#323C43",marker="x",s=26,alpha=.8,
                   label="Theoretical Heavy minus fi (not changed in source)")
    axs[0].set_ylabel("log2 STGSWSTLK Heavy")
    axs[0].set_title("Original Heavy reference signal / original reference medians",
                     loc="left",fontweight="bold")
    axs[0].legend(frameon=False,ncol=4,fontsize=8)
    axs[1].set_ylabel("Actual protein additive shift = -fi (log2)")
    axs[1].set_title("Reverse-direction offset applied to EVERY protein in sample",
                     loc="left",fontweight="bold")
    axs[1].axhline(0,color="#32485A",lw=.7)
    for ax in axs:
        for cut in (16.5,24.5,31.5):ax.axvline(cut,c="#8498A6",lw=.8)
        ax.grid(axis="y",alpha=.12)
        ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    axs[-1].set_xticks(pos,names,rotation=90,fontsize=6.4)
    fig.suptitle("Reverse factor sign | "+variant,fontweight="bold",fontsize=15)
    savefig(fig,dir_/PNGS[1])

def cv_chart(stats,dir_):
    x=np.arange(4)
    fig,axs=plt.subplots(1,2,figsize=(15,6.1),layout="constrained")
    raw=[next(r["median_protein_cv_pct"] for r in stats
              if r["stage"]=="Raw" and r["batch"]==b) for b in BATCHES]
    corr=[next(r["median_protein_cv_pct"] for r in stats
               if r["stage"]=="Corrected" and r["batch"]==b) for b in BATCHES]
    a=axs[0]
    a.bar(x-.19,raw,width=.36,color="#8193A1",label="Raw")
    a.bar(x+.19,corr,width=.36,color="#3379AD",label="Raw - fi")
    for xx,v in zip(list(x-.19)+list(x+.19),raw+corr):
        a.text(xx,v+.35,f"{v:.2f}",ha="center",fontsize=9)
    a.set_xticks(x,[f"{b}\nn={prev.BATCH_COUNTS[b]}" for b in BATCHES])
    a.set_ylim(0,max(raw+corr)*1.26)
    a.set_ylabel("Within-batch per-protein median CV (%)")
    a.set_title("Batch median CV (linear-scale intensities)",loc="left",fontweight="bold")
    a.legend(frameon=False)
    a=axs[1]
    diff=np.array(corr)-np.array(raw)
    a.bar(x,diff,color=["#238674" if y<0 else "#B76542" for y in diff],width=.59)
    for xx,d in zip(x,diff):
        a.text(xx,d+(.20 if d>=0 else -.30),f"{d:+.2f} pp",
               ha="center",va="bottom" if d>=0 else "top",fontsize=10)
    a.axhline(0,c="#455B6E",lw=.8)
    a.set_xticks(x,list(BATCHES))
    a.set_ylabel("Corrected CV − Raw CV (percentage points)")
    a.set_title("Change after reversing the original factor",loc="left",fontweight="bold")
    for ax in axs:
        ax.grid(axis="y",alpha=.15)
        ax.set_axisbelow(True)
        ax.spines["right"].set_visible(False);ax.spines["top"].set_visible(False)
    fig.suptitle("STGSWSTLK Heavy sign-reversal trial",fontweight="bold",fontsize=14)
    savefig(fig,dir_/PNGS[2])

def common_raw_fit_pca(raw,corrected_by_variant,names,records):
    # The same existing Raw-only eigenvector routine is reused for both.
    result={}
    for label,corrected in corrected_by_variant.items():
        rows,audit=prev.raw_fit_pca(raw,corrected,names,records)
        result[label]=(rows,audit)
    a=list(result.values())[0]
    other=list(result.values())[1]
    cols=("Raw_PC1","Raw_PC2")
    for col in cols:
        if not np.allclose([r[col] for r in a[0]],
                           [r[col] for r in other[0]],rtol=0,atol=1e-11):
            raise AssertionError("Raw PCA must be identical between both corrected variants")
    if a[1]!=other[1]:raise AssertionError("Different Raw PCA fitting metadata")
    arrays=[]
    for rows,_ in result.values():
        arrays.extend([
            np.array([[r["Raw_PC1"],r["Raw_PC2"]] for r in rows]),
            np.array([[r["Corrected_projected_PC1"],r["Corrected_projected_PC2"]] for r in rows])
        ])
    full=np.vstack(arrays)
    lower=full.min(axis=0)
    upper=full.max(axis=0)
    pad=np.maximum((upper-lower)*.12,1.)
    bounds=[float(lower[0]-pad[0]),float(upper[0]+pad[0]),
            float(lower[1]-pad[1]),float(upper[1]+pad[1])]
    return result,bounds

def pca_two_separate(rows,audit,bounds,dir_,variant):
    raw=np.array([[r["Raw_PC1"],r["Raw_PC2"]] for r in rows])
    corr=np.array([[r["Corrected_projected_PC1"],r["Corrected_projected_PC2"]]
                   for r in rows])
    for val,title,file in (
        (raw,"Raw PCA · fitted only on 39 original samples",PNGS[3]),
        (corr,"Corrected PCA · projected on identical Raw axes",PNGS[4])
    ):
        fig,ax=plt.subplots(figsize=(8.7,7.25),layout="constrained")
        for batch in BATCHES:
            mask=np.array([r["batch"]==batch for r in rows])
            ax.scatter(val[mask,0],val[mask,1],color=COLORS[batch],
                       s=60,lw=.6,edgecolor="white",
                       label=f"{batch} (n={prev.BATCH_COUNTS[batch]})")
        ax.set_xlim(bounds[0],bounds[1]);ax.set_ylim(bounds[2],bounds[3])
        ax.set_aspect("equal",adjustable="box")
        ax.set_xlabel(f"Raw PC1 ({audit['PC1_raw_explained_variance_pct']:.2f}%)")
        ax.set_ylabel(f"Raw PC2 ({audit['PC2_raw_explained_variance_pct']:.2f}%)")
        ax.set_title(title+"\n"+variant,loc="left",fontsize=13,fontweight="bold")
        ax.legend(frameon=False,ncol=2,fontsize=9)
        ax.axhline(0,c="#AEB9C2",lw=.75);ax.axvline(0,c="#AEB9C2",lw=.75)
        ax.grid(alpha=.12)
        ax.spines["right"].set_visible(False);ax.spines["top"].set_visible(False)
        savefig(fig,dir_/file)

def summary_plot(rows,dir_):
    fig,ax=plt.subplots(figsize=(14.5,6.9),layout="constrained")
    x=np.arange(4)
    variants=[
        ("Raw","#8799A5",-.28),
        ("Global + fi","#D9A85B",-.14),
        ("Global - fi","#2A88B8",0.),
        ("Batch + fi","#D7975A",.14),
        ("Batch - fi","#397F98",.28),
    ]
    for method,col,shift in variants:
        nums=[next(r[method] for r in rows if r["batch"]==b) for b in BATCHES]
        bars=ax.bar(x+shift,nums,width=.125,color=col,alpha=.9,label=method)
    ax.set_xticks(x,[f"{b}\nn={prev.BATCH_COUNTS[b]}" for b in BATCHES])
    ax.set_ylim(0,max(r[v[0]] for r in rows for v in variants)*1.24)
    ax.set_ylabel("Within-batch median protein CV (%)")
    ax.set_title("Same 39 samples · compare original and reversed factor direction",
                 loc="left",fontweight="bold",fontsize=14)
    ax.legend(frameon=False,ncol=5,fontsize=8)
    ax.grid(axis="y",alpha=.14);ax.set_axisbelow(True)
    ax.spines["right"].set_visible(False);ax.spines["top"].set_visible(False)
    savefig(fig,dir_/"00_Raw_plus_factor_minus_factor_CV_comparison.png")

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError("Independent reverse-factor experiment directory already exists")
    for folder in (PLUS_GLOBAL,PLUS_BATCH):
        if not (folder/"batch_protein_CV_summary.csv").is_file():
            raise FileNotFoundError(f"Comparison baseline missing: {folder}")
    records,names,annotations,raw,targets,heavy,ref40=prev.inputs()
    if len(names)!=39 or prev.EXCLUDED_ID in names:
        raise AssertionError("The confirmed T003 outlier must be excluded")
    h=np.log2(heavy)
    global_ref=float(np.median(h))
    f_global=global_ref-h
    hb,batch_idx,batch_refs,f_batch,positive_batch=batch.factors(records,heavy,raw)
    if not np.allclose(h,hb,atol=1e-12):raise AssertionError("Batch Heavy mismatch")
    prev_global_f=pd.read_csv(PLUS_GLOBAL/"sample_STGSWSTLK_Heavy_factors.csv").set_index("sample_id").loc[names]
    prev_batch_f=pd.read_csv(PLUS_BATCH/"sample_batch_internal_factors.csv").set_index("sample_id").loc[names]
    if not np.allclose(f_global,prev_global_f.fi_log2,atol=1e-12):
        raise AssertionError("Reverse global factor must EXACTLY reuse original global f")
    if not np.allclose(f_batch,prev_batch_f.fi_log2,atol=1e-12):
        raise AssertionError("Reverse batch factor must EXACTLY reuse original batch f")
    reverse_factors={"global_subtract_fi":f_global,"batch_subtract_fi":f_batch}
    matrices={label:raw-f[None,:] for label,f in reverse_factors.items()}
    for label,after in matrices.items():
        if not np.array_equal(np.isfinite(after),np.isfinite(raw)):
            raise AssertionError("Reverse calibration cannot change missingness")
        if not np.allclose(after+reverse_factors[label][None,:],raw,atol=1e-11,equal_nan=True):
            raise AssertionError("Reverse application formula incorrect")
    pca_results,pca_bounds=common_raw_fit_pca(raw,matrices,names,records)
    OUT.mkdir(parents=True,exist_ok=False)
    cv_compare={}
    for label,title in VARIANTS.items():
        path=OUT/label
        path.mkdir(parents=True,exist_ok=False)
        prev.OUT=path
        f=reverse_factors[label]
        after=matrices[label]
        refs=np.array([global_ref]*len(names)) if label=="global_subtract_fi" else np.array([batch_refs[r["batch"]] for r in records])
        tables=[]
        for j,sid in enumerate(names):
            b=records[j]["batch"]
            item={"sample_id":sid,"batch":b,
                  "raw_STGSWSTLK_Heavy_intensity_linear":float(heavy[j]),
                  "median_hi_log2_single_peptide":float(h[j]),
                  "median_ref_log2":float(refs[j]),
                  "fi_original_log2":float(f[j]),
                  "actual_subtracted_shift_log2":float(-f[j]),
                  "applied_multiplier_linear":float(np.exp2(-f[j])),
                  "theoretical_Heavy_after_subtract_log2":float(h[j]-f[j]),
                  "Raw_PG_median_log2":float(np.nanmedian(raw[:,j])),
                  "Corrected_PG_median_log2":float(np.nanmedian(after[:,j]))}
            for acc,p in prev.PROTEINS.items():
                ix=targets[acc]
                item[p["gene"]+"_Raw_log2"]=float(raw[ix,j])
                item[p["gene"]+"_Corrected_log2"]=float(after[ix,j])
            tables.append(item)
        store_csv(path/"sample_factors_subtract_fi.csv",tables)
        pg_matrix(path/"corrected_NYD_pg_matrix2_log2.tsv",annotations,names,after)
        cv_rows,detail,vec=prev.stats_cv(raw,after,records,annotations)
        store_csv(path/"batch_protein_CV_summary.csv",cv_rows)
        store_csv(path/"per_protein_batch_CV_long.csv",detail)
        cv_compare[label]=cv_rows
        store_csv(path/"PCA_Raw_fitted_sample_coordinates.csv",pca_results[label][0])
        pg_boxes(raw,after,records,names,targets,path,title)
        factor_plot(h,refs,f,names,records,path,title,
                    batch_refs if label=="batch_subtract_fi" else None)
        cv_chart(cv_rows,path)
        pca_two_separate(*pca_results[label],pca_bounds,path,title)
        prev.plot_target_proteins(raw,after,records,targets)
        prev.plot_sample_medians(raw,after,names,records)
        meta={
            "variant":label,"reference":title,"sample_n":39,
            "excluded_sample":prev.EXCLUDED_ID,
            "fi":"median reference (global or batch-local) - raw log2 Heavy",
            "application":"corrected_PG_log2 = raw_PG_log2 - fi",
            "unmodified_reference_factor_source":str(PLUS_GLOBAL if label=="global_subtract_fi" else PLUS_BATCH),
            "reference_log2":float(global_ref) if label=="global_subtract_fi" else batch_refs,
            "n_protein_groups":4264,"n_valid_pg_values":int(np.isfinite(after).sum()),
            "PCA_Raw_fit_only":pca_results[label][1],
            "PCA_axis_bounds_identical_across_BOTH_variants":pca_bounds,
            "separate_PCA_figures":True,
            "no_imputation_in_CV_or_protein_matrix":True,
            "figures":PNGS
        }
        (path/"factor_reference.json").write_text(json.dumps(meta,indent=2,ensure_ascii=False),encoding="utf-8")
    prev_plus={}
    for label,source in (("Global + fi",PLUS_GLOBAL),("Batch + fi",PLUS_BATCH)):
        prior=pd.read_csv(source/"batch_protein_CV_summary.csv")
        prev_plus[label]={(r.stage,r.batch):float(r.median_protein_cv_pct)
                          for r in prior.itertuples()}
    comparison=[]
    for b in BATCHES:
        rawcv=next(r["median_protein_cv_pct"] for r in cv_compare["global_subtract_fi"]
                   if r["stage"]=="Raw" and r["batch"]==b)
        globminus=next(r["median_protein_cv_pct"] for r in cv_compare["global_subtract_fi"]
                       if r["stage"]=="Corrected" and r["batch"]==b)
        batchminus=next(r["median_protein_cv_pct"] for r in cv_compare["batch_subtract_fi"]
                        if r["stage"]=="Corrected" and r["batch"]==b)
        if not np.isclose(globminus,batchminus,rtol=1e-12,atol=1e-10):
            raise AssertionError("Within-batch reverse-factor CV must be invariant to constant batch shifts")
        if not np.isclose(rawcv,prev_plus["Global + fi"][("Raw",b)],atol=1e-11):
            raise AssertionError("Raw group CV changed from 39-sample reference")
        comparison.append({"batch":b,"n_samples":prev.BATCH_COUNTS[b],
            "Raw":rawcv,
            "Global + fi":prev_plus["Global + fi"][("Corrected",b)],
            "Global - fi":globminus,
            "Batch + fi":prev_plus["Batch + fi"][("Corrected",b)],
            "Batch - fi":batchminus,
            "Reverse_minus_Raw_percentage_points":globminus-rawcv})
    store_csv(OUT/"00_global_batch_plus_minus_CV_comparison.csv",comparison)
    summary_plot(comparison,OUT)
    audit={
      "inputs_sha256":{file.name:prev.sha256(file) for file in (prev.PG,prev.HELP,prev.META)},
      "excluded_sample_id":prev.EXCLUDED_ID,
      "removed_before_all_references_and_analysis":True,
      "n_original":40,"n_retained":39,"batch_sizes":prev.BATCH_COUNTS,
      "references_log2":{"global":global_ref,"batch":batch_refs},
      "factors_exactly_match_original_plus_f_trials":True,
      "reverse_factor_application":"log2 PG corrected = raw log2 PG - (reference - log2 H)",
      "corrected_does_not_align_STGSWSTLK_Heavy_to_median":True,
      "previous_results_and_origins_untouched":True,
      "PCA_method":"Raw PCA fitted once and reused for both corrected branches",
      "same_PCA_limits_all_four_images":pca_bounds,
      "Raw_PCA_details":pca_results["global_subtract_fi"][1],
      "CV_evaluated_in_original_linear_intensity":True,
      "per_batch_reverse_CV_global_equals_batch_internal":True,
      "summary_figures":["00_Raw_plus_factor_minus_factor_CV_comparison.png"]}
    (OUT/"run_manifest.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding="utf-8")
    doc=[
      "# NYD STGSWSTLK Heavy · reverse direction (-fi) trials","",
      "39 retained NYD clinical samples (T001=16, T002=8, T003=7, T004=8);",
      "confirmed bad T003 sample MAP26092326E01 is excluded from BOTH trials.",
      "Raw input files and both earlier positive-direction corrections untouched.","",
      "## Only change is the sign on the protein correction","",
      "Original reference definition is unchanged:",
      "fi = median_ref − log2(H_STGSWSTLK,sample_i).",
      "**Old:** corrected PG log2 = Raw PG log2 + fi.",
      "**NEW:** corrected PG log2 = Raw PG log2 − fi.",
      "Global version uses same reference median of all 39 samples.",
      "Batch-internal version independently uses T001/T002/T003/T004 reference medians.",
      "No other HELP peptide, weighting, limiting, outlier redefinition,",
      "or other normalization steps. Heavy.xlsx is not modified.","",
      "These reverse shifts mathematically move STGSWSTLK Heavy *away*",
      "from its reference if applied to Heavy. They are exploratory,",
      "not a physically justified calibration.","",
      "## Median within-batch protein CV (%)","",
      "| Batch | n | Raw | Global + fi | Global − fi | Batch + fi | Batch − fi |",
      "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for r in comparison:
        doc.append(f'| {r["batch"]} | {r["n_samples"]} | {r["Raw"]:.2f} | {r["Global + fi"]:.2f} | {r["Global - fi"]:.2f} | {r["Batch + fi"]:.2f} | {r["Batch - fi"]:.2f} |')
    doc += ["",
      "![All CV branches](00_Raw_plus_factor_minus_factor_CV_comparison.png)","",
      "**Within-batch CV is the same for global -fi and batch -fi**",
      "because the two reversed corrections differ by only one constant",
      "within each batch. The corresponding +fi approaches likewise",
      "have identical within-batch CV. This does not imply the two",
      "methods yield the same across-batch intensity or PCA.",
      "",
      "## Variant directories","",
      "### global_subtract_fi","",
      "![Raw versus reverse-corrected PG](global_subtract_fi/01_PG_Raw_vs_subtract_factor_blue_CFB_yellow_PROS1.png)","",
      "![Raw PCA](global_subtract_fi/04a_PCA_Raw.png)","",
      "![Corrected PCA](global_subtract_fi/04b_PCA_Corrected.png)","",
      "### batch_subtract_fi","",
      "![Raw versus reverse-corrected PG](batch_subtract_fi/01_PG_Raw_vs_subtract_factor_blue_CFB_yellow_PROS1.png)","",
      "![Raw PCA](batch_subtract_fi/04a_PCA_Raw.png)","",
      "![Corrected PCA](batch_subtract_fi/04b_PCA_Corrected.png)","",
      "## Files per variant","",
      "Seven PNGs including separate Raw and Corrected PCA images;",
      "full 4264 x 39 corrected log2 PG matrix, per-sample original fi,",
      "actual applied (-fi), per-batch protein CV, per-protein CV,",
      "Raw-fitted PCA scores, and factor reference metadata.",
      "All four PCA images from both variants use the exact same Raw",
      "PC components, raw-based center, missingness treatment and axis limits.",
      "Blue CFB/P00751 and yellow PROS1/P07225 are marked on PG boxes.","",
      "Because samples within each batch are different clinical specimens,",
      "within-batch protein CV includes biological variation. A lower CV",
      "is not by itself evidence of more accurate instrument calibration."
    ]
    (OUT/"README.md").write_text("\n".join(doc)+"\n",encoding="utf-8")
    print("REFERENCES:",json.dumps({"global":global_ref,"batch":batch_refs}),flush=True)
    print("CV:",json.dumps(comparison),flush=True)
    print("PASS: 2 reversed factors, 39 samples, 2 PG matrices, 14 method PNGs + 1 joint CV, common Raw PCA",flush=True)

if __name__=="__main__":
    main()
