#!/usr/bin/env python3
"""NYD 39 retained samples: STGSWSTLK Heavy batch-specific median correction.

Remove confirmed bad T003 sample MAP26092326E01 before analysis.
For each batch b: median_ref_b = median(log2 STGSWSTLK Heavy of that batch).
For each sample i in batch b: f_i = median_ref_b - log2(H_i);
corrected log2 PG_i = raw log2 PG_i + f_i.

Separate Raw and Corrected PCA figures, both in one Raw-fitted space.
Compare to previous 39-sample global-median result without overwriting it.
"""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import nyd_STGSWSTLK_heavy_exclude_T003_outlier as prev
from nyd_pg_batch_boxplots import BATCHES, COLORS, ANNOT

OUT=prev.BASE/"STGSWSTLK_Heavy_batch_internal_normalization_without_T003_outlier"
OLD=prev.BASE/"STGSWSTLK_Heavy_single_reference_normalization_without_T003_outlier"

def factors(records,heavy,raw):
    h=np.log2(heavy)
    idx={b:np.array([i for i,r in enumerate(records) if r["batch"]==b],dtype=int)
         for b in BATCHES}
    assert {b:len(j) for b,j in idx.items()}==prev.BATCH_COUNTS
    med={b:float(np.median(h[j])) for b,j in idx.items()}
    f=np.empty(len(h))
    for b,j in idx.items():
        f[j]=med[b]-h[j]
        assert np.allclose(h[j]+f[j],med[b],atol=1e-12,rtol=0)
        assert abs(np.median(f[j]))<1e-12
    corrected=raw+f[None,:]
    assert np.array_equal(np.isfinite(corrected),np.isfinite(raw))
    return h,idx,med,f,corrected

def plot_heavy(h,idx,med,f,names):
    fig,axs=plt.subplots(2,1,figsize=(18.5,9.2),sharex=True,layout="constrained")
    x=np.arange(1,len(names)+1)
    for b,sel in idx.items():
        axs[0].scatter(x[sel],h[sel],color=COLORS[b],s=46,
                       edgecolors="white",linewidths=.6,label=b)
        axs[0].plot([x[sel].min()-.4,x[sel].max()+.4],[med[b]]*2,
                    color=COLORS[b],lw=1.5,ls="--",label=b+" median_ref")
        axs[0].scatter(x[sel],h[sel]+f[sel],marker="_",s=65,c="#293A46")
        axs[1].bar(x[sel],f[sel],color=COLORS[b],alpha=.78,width=.75)
    axs[0].set_title("Original STGSWSTLK Heavy vs batch median_ref (dashed)",
                     loc="left",fontweight="bold")
    axs[0].set_ylabel("log2 Heavy intensity")
    axs[0].legend(frameon=False,ncol=4,fontsize=8)
    axs[1].axhline(0,color="#43596C",lw=.8)
    axs[1].set_title("Own-batch correction factors",loc="left",fontweight="bold")
    axs[1].set_ylabel("fi = batch median_ref - log2(H_i)")
    for ax in axs:
        for cut in (16.5,24.5,31.5):ax.axvline(cut,c="#9BA8B3",lw=.9)
        ax.grid(axis="y",alpha=.13)
        ax.spines["right"].set_visible(False);ax.spines["top"].set_visible(False)
    axs[1].set_xticks(x,names,rotation=90,fontsize=6.5)
    fig.suptitle("39 NYD samples, four independent Heavy reference medians",
                 fontweight="bold",fontsize=15)
    prev.save(fig,"02_STGSWSTLK_Heavy_batch_references_and_factors.png")

def separate_common_pca(scores,audit):
    raw=np.array([[r["Raw_PC1"],r["Raw_PC2"]] for r in scores])
    corr=np.array([[r["Corrected_projected_PC1"],r["Corrected_projected_PC2"]] for r in scores])
    allx=np.vstack([raw,corr])
    mi=allx.min(axis=0);ma=allx.max(axis=0)
    padding=np.maximum((ma-mi)*.12,1.)
    bounds=[float(mi[0]-padding[0]),float(ma[0]+padding[0]),
            float(mi[1]-padding[1]),float(ma[1]+padding[1])]
    for arr,title,filename in (
      (raw,"Raw PCA (fitted on Raw)","04a_PCA_Raw.png"),
      (corr,"Corrected PCA (projected onto Raw axes)","04b_PCA_Corrected.png")):
        fig,ax=plt.subplots(figsize=(9,7.5),layout="constrained")
        for batch in BATCHES:
            mask=np.array([x["batch"]==batch for x in scores])
            ax.scatter(arr[mask,0],arr[mask,1],s=54,color=COLORS[batch],
                       edgecolors="white",linewidths=.6,
                       label=f"{batch} (n={prev.BATCH_COUNTS[batch]})")
        ax.set_xlim(bounds[0],bounds[1])
        ax.set_ylim(bounds[2],bounds[3])
        ax.set_aspect("equal",adjustable="box")
        ax.set_xlabel(f"Raw PC1 ({audit['PC1_raw_explained_variance_pct']:.2f}%)")
        ax.set_ylabel(f"Raw PC2 ({audit['PC2_raw_explained_variance_pct']:.2f}%)")
        ax.set_title(title+"\nSame 39 samples, same axis limits",loc="left",fontweight="bold")
        ax.legend(frameon=False,ncol=2,fontsize=9)
        ax.axhline(0,color="#B6C3CB",lw=.7)
        ax.axvline(0,color="#B6C3CB",lw=.7)
        ax.grid(alpha=.14)
        ax.spines["right"].set_visible(False);ax.spines["top"].set_visible(False)
        prev.save(fig,filename)
    return {"x_min":bounds[0],"x_max":bounds[1],"y_min":bounds[2],"y_max":bounds[3]}

def plot_comparison(comp):
    x=np.arange(4)
    fig,ax=plt.subplots(figsize=(11,6.1),layout="constrained")
    stage=[
        ("Raw_median_cv_pct","Raw","#8695A3",-.26),
        ("Previous_global_median_cv_pct","Previous global median","#D69D55",0),
        ("Batch_internal_median_cv_pct","New batch-internal median","#2591A1",.26),
    ]
    for key,label,col,shift in stage:
        vals=[r[key] for r in comp]
        bars=ax.bar(x+shift,vals,width=.245,color=col,label=label)
        for bar,v in zip(bars,vals):
            ax.text(bar.get_x()+bar.get_width()/2,v+.25,f"{v:.2f}",
                    ha="center",fontsize=9)
    ax.set_xticks(x,[f"{b}\nn={prev.BATCH_COUNTS[b]}" for b in BATCHES])
    ax.set_ylim(0,max(r["Batch_internal_median_cv_pct"] for r in comp)*1.28)
    ax.set_ylabel("Within-batch median protein CV (%)")
    ax.set_title("Same 39 samples: Raw / global / batch-local reference",
                 loc="left",fontweight="bold",fontsize=14)
    ax.legend(frameon=False,ncol=3,fontsize=9)
    ax.grid(axis="y",alpha=.15)
    ax.spines["right"].set_visible(False);ax.spines["top"].set_visible(False)
    prev.save(fig,"07_Raw_global_batch_internal_CV.png")

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError("Existing batch-internal output directory protected")
    if not (OLD/"batch_protein_CV_summary.csv").is_file():
        raise FileNotFoundError("Prior 39-sample global trial unavailable")
    records,names,ann,raw,targets,heavy,ref40=prev.inputs()
    assert len(names)==39 and prev.EXCLUDED_ID not in names
    h,idx,refs,f,corrected=factors(records,heavy,raw)
    if np.max(np.abs(f))>=1.:raise ValueError("Extreme unexpected factor remains")
    OUT.mkdir(parents=True,exist_ok=False)
    prev.OUT=OUT
    reference=[{"batch":b,"n_samples":len(sel),
       "median_ref_log2_STGSWSTLK_Heavy":refs[b],
       "median_ref_linear_STGSWSTLK_Heavy":float(np.exp2(refs[b])),
       "fi_min_log2":float(f[sel].min()),
       "fi_max_log2":float(f[sel].max()),
       "fi_median_log2":float(np.median(f[sel]))} for b,sel in idx.items()]
    prev.writecsv("batch_STGSWSTLK_Heavy_median_references.csv",reference)
    factor_rows=[]
    for j,s in enumerate(names):
        batch=records[j]["batch"]
        row={"sample_id":s,"batch":batch,
            "Heavy_STGSWSTLK_linear":float(heavy[j]),
            "median_hi_log2":float(h[j]),
            "median_ref_batch_log2":refs[batch],
            "fi_log2":float(f[j]),
            "multiplier_linear":float(np.exp2(f[j])),
            "Heavy_after_factor_diagnostic_log2":float(h[j]+f[j]),
            "PG_Raw_sample_median_log2":float(np.nanmedian(raw[:,j])),
            "PG_Corrected_sample_median_log2":float(np.nanmedian(corrected[:,j]))}
        for accession,params in prev.PROTEINS.items():
            k=targets[accession]
            row[params["gene"]+"_Raw_log2"]=float(raw[k,j])
            row[params["gene"]+"_Corrected_log2"]=float(corrected[k,j])
        factor_rows.append(row)
    prev.writecsv("sample_batch_internal_factors.csv",factor_rows)
    with (OUT/"corrected_NYD_pg_matrix2_log2.tsv").open("w",encoding="utf-8",newline="") as fh:
        w=csv.writer(fh,delimiter="\t")
        w.writerow(list(ANNOT)+names)
        for i,anno in enumerate(ann):
            w.writerow(list(anno)+[float(v) if np.isfinite(v) else "" for v in corrected[i]])
    cvsummary,cvdetail,arrays=prev.stats_cv(raw,corrected,records,ann)
    prev.writecsv("batch_protein_CV_summary.csv",cvsummary)
    prev.writecsv("per_protein_batch_CV_long.csv",cvdetail)
    scores,pca_audit=prev.raw_fit_pca(raw,corrected,names,records)
    prev.writecsv("PCA_Raw_fitted_sample_coordinates.csv",scores)
    # Original charting routine: raw/corrected log2 PG boxplots and two target proteins.
    prev.plot_pgs(raw,corrected,records,names,targets)
    plot_heavy(h,idx,refs,f,names)
    prev.plot_cv(cvsummary,arrays)
    limits=separate_common_pca(scores,pca_audit)
    prev.plot_target_proteins(raw,corrected,records,targets)
    prev.plot_sample_medians(raw,corrected,names,records)
    prior=pd.read_csv(OLD/"batch_protein_CV_summary.csv").set_index(["stage","batch"])
    comparison=[]
    for b in BATCHES:
        rawcv=next(r["median_protein_cv_pct"] for r in cvsummary
                   if r["stage"]=="Raw" and r["batch"]==b)
        newcv=next(r["median_protein_cv_pct"] for r in cvsummary
                   if r["stage"]=="Corrected" and r["batch"]==b)
        globalcv=float(prior.loc[("Corrected",b),"median_protein_cv_pct"])
        earlier_raw=float(prior.loc[("Raw",b),"median_protein_cv_pct"])
        if not np.isclose(rawcv,earlier_raw,atol=1e-10,rtol=0):
            raise AssertionError("Raw baseline changed relative to 39-sample global trial")
        if not np.isclose(newcv,globalcv,atol=1e-10,rtol=1e-12):
            raise AssertionError("Within-batch CV must equal prior global-median CV")
        comparison.append({"batch":b,"n_samples":len(idx[b]),
          "Raw_median_cv_pct":float(rawcv),
          "Previous_global_median_cv_pct":globalcv,
          "Batch_internal_median_cv_pct":float(newcv),
          "Batch_minus_previous_global_pp":float(newcv-globalcv),
          "Batch_minus_Raw_pp":float(newcv-rawcv)})
    prev.writecsv("Raw_global_batch_internal_CV_comparison.csv",comparison)
    plot_comparison(comparison)
    images=[
        "01_PG_39samples_Raw_vs_corrected_highlight_CFB_PROS1.png",
        "02_STGSWSTLK_Heavy_batch_references_and_factors.png",
        "03_batch_protein_CV_raw_vs_corrected.png",
        "04a_PCA_Raw.png","04b_PCA_Corrected.png",
        "05_CFB_PROS1_by_batch_Raw_vs_corrected.png",
        "06_sample_PG_medians_Raw_vs_corrected.png",
        "07_Raw_global_batch_internal_CV.png"]
    manifest={
        "n_original_samples":40,"n_retained_samples":39,
        "excluded_sample_id":prev.EXCLUDED_ID,"excluded_before_medians":True,
        "source_sha256":{p.name:prev.sha256(p) for p in (prev.PG,prev.HELP,prev.META)},
        "previous_global_CV_sha256":prev.sha256(OLD/"batch_protein_CV_summary.csv"),
        "n_protein_groups":len(ann),"batch_counts":{b:len(sel) for b,sel in idx.items()},
        "selected_HELP":"STGSWSTLK Heavy only",
        "median_ref":"median(log2 STGSWSTLK Heavy), independently calculated within each batch",
        "fi":"median_ref_batch - log2 STGSWSTLK Heavy per sample",
        "corrected":"raw PG log2 intensity + fi for 39 retained samples",
        "batch_references_log2":refs,
        "reference_not_global":True,
        "median_of_each_batch_factors_is_zero":True,
        "no_other_HELP_or_weight":True,"no_imputation_for_PG_or_CV":True,
        "original_input_files_and_old_results_unchanged":True,
        "CV_measured_on_linear_intensity":True,
        "within_batch_CV_equal_to_prior_global_result":True,
        "PCA_Raw_fit":pca_audit,"common_PCA_axis_limits":limits,
        "PCA_display":"Two separate PNG images; no combined raw/corrected plot and no arrows",
        "figures":images,
        "valid_PG_cells":int(np.isfinite(corrected).sum())}
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    doc=[
      "# NYD batch-internal STGSWSTLK Heavy normalization · 39 samples","",
      "Exclude MAP26092326E01 (T003) from analysis before recomputing factors.",
      "Retained T001=16, T002=8, T003=7, T004=8. Original files unchanged.","",
      "## Method","",
      "One Heavy STGSWSTLK intensity per sample. Transform to log2 first.",
      "For each batch, reference = median of only that batch's log2 Heavy values.",
      "Each sample's fi = median_ref_own_batch - log2 Heavy_i.",
      "Add fi to every observed log2 PG protein intensity in the same sample.",
      "No sample gets a reference from another batch. Missing PGs remain missing.",
      "This method aligns Heavy measurements WITHIN batches; it does not",
      "force the four batch median references to be equal.","",
      "## Batch-specific reference medians","",
      "| Batch | n | median_ref log2 | median_ref linear | fi min | fi max |",
      "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for r in reference:
        doc.append(f"| {r['batch']} | {r['n_samples']} | {r['median_ref_log2_STGSWSTLK_Heavy']:.6f} | {r['median_ref_linear_STGSWSTLK_Heavy']:.1f} | {r['fi_min_log2']:.4f} | {r['fi_max_log2']:.4f} |")
    doc += ["","## Batch protein median CV (%)","",
       "| Batch | Raw | Global 39-sample median correction | Batch-internal correction | Batch - Raw (pp) |",
       "| --- | ---: | ---: | ---: | ---: |"]
    for r in comparison:
        doc.append(f"| {r['batch']} | {r['Raw_median_cv_pct']:.2f} | {r['Previous_global_median_cv_pct']:.2f} | {r['Batch_internal_median_cv_pct']:.2f} | {r['Batch_minus_Raw_pp']:+.2f} |")
    doc += ["",
       "Because the difference between global and batch reference is a",
       "constant offset within each batch, within-batch protein CV is",
       "mathematically identical for both approaches (up to precision).",
       "They differ in the BETWEEN-batch intensity alignment and PCA positions.","",
       "## PCA: two separate figures with the same Raw-fitted coordinates","",
       "Raw fit performed exactly once. Corrected projected to the Raw PCA",
       "loadings, raw-selected proteins, Raw missingness imputation medians",
       "and original Raw centering. Both figures have identical axis ranges.",
       "No connecting arrows or before/after overlay.","",
       "### Raw PCA","",
       "![Raw PCA](04a_PCA_Raw.png)","",
       "### Corrected PCA","",
       "![Corrected PCA](04b_PCA_Corrected.png)","",
       "## Other diagnostic visualizations","",
       "![Raw vs Corrected PG](01_PG_39samples_Raw_vs_corrected_highlight_CFB_PROS1.png)","",
       "![Heavy references and sample factors](02_STGSWSTLK_Heavy_batch_references_and_factors.png)","",
       "![CV](03_batch_protein_CV_raw_vs_corrected.png)","",
       "![Blue CFB and yellow PROS1](05_CFB_PROS1_by_batch_Raw_vs_corrected.png)","",
       "![PG sample medians](06_sample_PG_medians_Raw_vs_corrected.png)","",
       "![Global vs batch-local reference CV](07_Raw_global_batch_internal_CV.png)","",
       "## Audit tables","",
       "batch_STGSWSTLK_Heavy_median_references.csv: four batch references.",
       "sample_batch_internal_factors.csv: 39 exact Heavy factors and two source proteins.",
       "corrected_NYD_pg_matrix2_log2.tsv: 4,264 PG proteins x 39 samples.",
       "batch_protein_CV_summary.csv, per_protein_batch_CV_long.csv,",
       "Raw_global_batch_internal_CV_comparison.csv and",
       "PCA_Raw_fitted_sample_coordinates.csv include detailed checks.",
       "These batch groups may include true biological differences; batch CV",
       "alone is not independent technical-quality validation.",""]
    (OUT/"README.md").write_text("\n".join(doc),encoding="utf-8")
    print("BATCH_REFS",json.dumps(refs),flush=True)
    print("CV",json.dumps(comparison),flush=True)
    print("PCA",json.dumps({**pca_audit,**limits}),flush=True)
    print("PASS: 39 samples, four batch medians, 8 PNG including TWO separate shared-coordinate PCA",flush=True)

if __name__=="__main__":
    main()
