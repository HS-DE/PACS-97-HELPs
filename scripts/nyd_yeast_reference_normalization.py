#!/usr/bin/env python3
"""Experimental yeast-reference correction. Not MacCoss median normalization.

PG factor: 40-sample median YEAST P00924 protein / sample P00924 protein.
PR factor: geometric mean of reference/sample ratios for TWO DIFFERENT
peptides, AADALLLK and VNQIGTLSESIK. Only charge states of EACH peptide
are added in linear intensity. Distinct peptides are NOT summed together.
Yeast-derived post-CV is circular; endogenous PG/PR CV is independently
described but biological variance may also be present.
"""
import csv,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nyd_sample_median_normalization as core
from nyd_pg_batch_boxplots import meta,GROUPS,BATCHES,COLORS
from nyd_yeast_batch_cv import load as load_yeast

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"NYD_data"/"40samples"
OUT=BASE/"yeast_reference_normalization_40samples"
ANALYTES=("AADALLLK","VNQIGTLSESIK","P00924_protein")
FIELDS=("AADALLLK_sum_linear","VNQIGTLSESIK_sum_linear","P00924_protein_linear")
PGFILE="NYD_pg_matrix2_yeast_protein_normalized.tsv"
PRFILE="NYD_pr_matrix2_yeast_two_peptide_normalized.tsv"

def writecsv(name,records):
    with (OUT/name).open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]))
        w.writeheader();w.writerows(records)

def save(fig,name):
    fig.savefig(OUT/name,dpi=200,bbox_inches="tight",facecolor="white")
    plt.close(fig)

def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y",alpha=.15)
    ax.set_axisbelow(True)

def yeast_factors(info):
    r,missing=load_yeast()
    if missing:raise ValueError("Yeast spike-in intensities contain missing/nonpositive cells")
    idmap={x["sample_id"]:x for x in r}
    if set(idmap)!={x["sample_id"] for x in info}:
        raise ValueError("Yeast/sample metadata ID mismatch")
    arr=np.array([[float(idmap[x["sample_id"]][k]) for k in FIELDS] for x in info])
    if arr.shape!=(40,3) or np.any(~np.isfinite(arr)) or np.any(arr<=0):
        raise ValueError("Expected 40x3 finite-positive yeast values")
    log=np.log2(arr)
    refs=np.median(log,axis=0)
    by_feat=np.exp2(refs[None,:]-log)
    # Two distinct peptides remain separate references; aggregate their
    # scale ratios on log scale (equivalent to their geometric mean).
    pg=by_feat[:,2]
    pr=np.exp2(np.mean(refs[None,:2]-log[:,:2],axis=1))
    if not np.allclose(pr,np.sqrt(by_feat[:,0]*by_feat[:,1]),rtol=1e-12):
        raise AssertionError("Two-reference geometric mean mismatch")
    factors=[]
    for i,x in enumerate(info):
        factors.append(dict(sample_id=x["sample_id"],batch=x["batch"],
          AADALLLK_raw_linear=float(arr[i,0]),
          VNQIGTLSESIK_raw_linear=float(arr[i,1]),
          P00924_protein_raw_linear=float(arr[i,2]),
          AADALLLK_reference_ratio=float(by_feat[i,0]),
          VNQIGTLSESIK_reference_ratio=float(by_feat[i,1]),
          PG_yeast_protein_factor=float(pg[i]),
          PR_two_distinct_peptides_factor=float(pr[i])))
    return factors,pg,pr,arr,refs

def plot_boxes(raw,after,info,level):
    fig,axs=plt.subplots(2,1,figsize=(20,11.6),sharex=True,
                         gridspec_kw={"hspace":.25})
    for ax,mat,label in zip(axs,(raw,after),("Raw","Yeast reference corrected")):
        data=[mat[:,i][np.isfinite(mat[:,i])] for i in range(40)]
        bp=ax.boxplot(data,positions=np.arange(1,41),widths=.62,
                      patch_artist=True,showfliers=False,manage_ticks=False,
                      medianprops={"color":"#203A4D","linewidth":1.6})
        for p,x in zip(bp["boxes"],info):
            p.set_facecolor(COLORS[x["batch"]]);p.set_alpha(.7)
        for boundary in (16.5,24.5,32.5):
            ax.axvline(boundary,color="#8293A4",lw=1,alpha=.65)
        ax.set_xlim(.5,40.5)
        ax.set_ylabel("log2 original / scaled endogenous intensity")
        ax.set_title(label,loc="left",fontweight="bold",fontsize=14)
        style(ax)
    axs[1].set_xticks(np.arange(1,41),[x["sample_id"] for x in info],
                       rotation=90,fontsize=7)
    fig.suptitle(level+" | 40 original vs yeast-reference-normalized samples",
                 fontsize=16,fontweight="bold",y=.98)
    save(fig,"01_pg_before_after.png" if level=="PG" else "02_pr_before_after.png")

def plot_factor(record):
    fig,axs=plt.subplots(1,2,figsize=(15.3,6),layout="constrained")
    for ax,field,label in zip(axs,
         ("PG_yeast_protein_factor","PR_two_distinct_peptides_factor"),
         ("PG: one yeast protein","PR: geometric mean of two yeast-peptide ratios")):
        ys=[[float(x[field]) for x in record if x["batch"]==b] for b in BATCHES]
        bp=ax.boxplot(ys,positions=range(1,5),patch_artist=True,widths=.52)
        for p,b in zip(bp["boxes"],BATCHES):
            p.set_facecolor(COLORS[b]);p.set_alpha(.52)
        for i,(b,v) in enumerate(zip(BATCHES,ys),1):
            ax.scatter(i+np.linspace(-.14,.14,len(v)),v,s=34,
                       color=COLORS[b],edgecolor="white",lw=.5,zorder=4)
        ax.axhline(1,ls="--",lw=.9,color="#596E7C")
        ax.set_xticks(range(1,5),[f"{b} (n={GROUPS[b]})" for b in BATCHES])
        ax.set_ylabel("Multiplicative factor")
        ax.set_title(label,loc="left",fontweight="bold")
        style(ax)
    save(fig,"03_yeast_scaling_factors_by_batch.png")

def yeast_check(info,arr,pg,pr):
    per=[];summ=[]
    for i,x in enumerate(info):
        for k,name in enumerate(ANALYTES):
            f=float(pg[i] if k==2 else pr[i])
            per.append(dict(sample_id=x["sample_id"],batch=x["batch"],
               analyte=name,source="PG" if k==2 else "PR",
               raw_linear=float(arr[i,k]),factor=f,
               corrected_linear=float(arr[i,k]*f)))
    for group in ("All_40",*BATCHES):
        for name in ANALYTES:
            ss=[r for r in per if r["analyte"]==name and
                (group=="All_40" or r["batch"]==group)]
            raw=np.array([r["raw_linear"] for r in ss])
            new=np.array([r["corrected_linear"] for r in ss])
            summ.append(dict(batch=group,analyte=name,n_samples=len(ss),
               raw_yeast_CV_pct=float(100*np.std(raw,ddof=1)/np.mean(raw)),
               corrected_yeast_CV_pct=float(100*np.std(new,ddof=1)/np.mean(new)),
               caution="CIRCULAR: analyzed yeast is input to the scaling factor"))
    return per,summ

def yeast_figure(rows):
    fig,axs=plt.subplots(1,3,figsize=(16.5,6),layout="constrained")
    for ax,name in zip(axs,ANALYTES):
        a=[next(r for r in rows if r["batch"]==b and r["analyte"]==name)
           for b in BATCHES]
        x=np.arange(4)
        ax.bar(x-.18,[r["raw_yeast_CV_pct"] for r in a],.34,
               color="#C1CCD6",label="Raw")
        ax.bar(x+.18,[r["corrected_yeast_CV_pct"] for r in a],.34,
               color="#397EA6",label="Yeast scaled")
        ax.set_title(name,fontweight="bold")
        ax.set_xticks(x,BATCHES)
        ax.set_ylabel("CV (%) · linear intensity")
        ax.legend(frameon=False,fontsize=9);style(ax)
    fig.suptitle("Yeast CV: circular diagnostic, NOT independent correction validation",
                 fontsize=14,fontweight="bold")
    save(fig,"04_yeast_cv_CIRCULAR.png")

def endogenous_cv_figure(records):
    fig,axs=plt.subplots(1,2,figsize=(15.5,6),layout="constrained")
    for ax,name in zip(axs,("pg","pr")):
        a=[next(r for r in records if r["matrix"]==name and
                r["batch"]==b and r["stage"]=="raw") for b in BATCHES]
        b=[next(r for r in records if r["matrix"]==name and
                r["batch"]==bb and r["stage"]=="median_normalized") for bb in BATCHES]
        x=np.arange(4)
        ax.bar(x-.17,[r["median_feature_cv_pct"] for r in a],.32,
               color="#C1CCD6",label="Raw")
        ax.bar(x+.17,[r["median_feature_cv_pct"] for r in b],.32,
               color="#397EA6",label="Yeast factor")
        ax.set_xticks(x,BATCHES)
        ax.set_title(name.upper()+" endogenous feature CVs",fontweight="bold")
        ax.set_ylabel("Within-batch median feature CV (%)")
        ax.legend(frameon=False);style(ax)
    fig.suptitle("Endogenous PG/PR biological-sample variability before vs after yeast scaling",
                 fontsize=14,fontweight="bold")
    save(fig,"05_endogenous_feature_CV_before_after.png")

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError("Existing yeast reference calibration results will not be overwritten")
    info=meta()
    factors,pgf,prf,yeast,refs=yeast_factors(info)
    OUT.mkdir(parents=True,exist_ok=True)
    writecsv("yeast_reference_factors.csv",factors)
    plot_factor(factors)
    core.OUT=OUT
    ids=[x["sample_id"] for x in info]
    sample_stats=[];batch_stats=[];cvstats=[];matrix_meta={}
    for label,src,mults,filename in (
        ("pg",BASE/"NYD_pg_matrix2.tsv",pgf,PGFILE),
        ("pr",BASE/"NYD_pr_matrix2.tsv",prf,PRFILE)):
        h,rows,indices,logged,exclude=core.load_matrix(src,ids)
        result=core.write_normalized(h,rows,indices,mults,filename)
        vh,vr,vi,check,_=core.load_matrix(result,ids)
        corrected=logged+np.log2(mults)[None,:]
        if vh!=h or len(vr)!=len(rows) or not np.array_equal(
              np.isfinite(check),np.isfinite(corrected)):
            raise AssertionError("Corrected matrix schema or missingness changed")
        valid=np.isfinite(check)
        if not np.allclose(check[valid],corrected[valid],atol=2e-10,rtol=0):
            raise AssertionError("Re-imported corrected values disagree with factors")
        # core.matrix_stats explicitly excludes yeast accession rows
        # when calculating endogenous feature CVs.
        _,ss,bb,cv=core.matrix_stats(logged,mults,label,exclude,info)
        sample_stats+=ss;batch_stats+=bb;cvstats+=cv
        plot_boxes(logged[~exclude],corrected[~exclude],info,label.upper())
        matrix_meta[label]={"source_rows":len(rows),"nonyeast_rows":int((~exclude).sum()),
             "excluded_yeast_rows":int(exclude.sum()),"factor_min":float(min(mults)),
             "factor_max":float(max(mults))}
        print("CORRECTED",label,json.dumps(matrix_meta[label]),flush=True)
    writecsv("endogenous_sample_stats.csv",sample_stats)
    writecsv("endogenous_batch_stats.csv",batch_stats)
    writecsv("endogenous_batch_CV.csv",cvstats)
    yrows,ycv=yeast_check(info,yeast,pgf,prf)
    writecsv("yeast_reference_response.csv",yrows)
    writecsv("yeast_CV_CIRCULAR.csv",ycv)
    yeast_figure(ycv)
    endogenous_cv_figure(cvstats)
    manifest=dict(
       method="Experimental yeast spike-in reference normalization, NOT the MacCoss 2026 method",
       source_sha256={x:core.sha(BASE/f) for x,f in (
           ("PG","NYD_pg_matrix2.tsv"),("PR","NYD_pr_matrix2.tsv"),
           ("YEAST","YEAST.xlsx"),("META","NYD_sampleinfo2.xlsx"))},
       sample_n=40,batch_sizes=GROUPS,yeast_used_for_factor=True,
       yeast_CV_after_is_independent=False,peptides_summed_together=False,
       charge_states_summed_within_peptide=True,
       protein_factor="40-sample median yeast P00924 / sample P00924",
       peptide_factor="geometric mean of (40-sample median peptide / sample peptide) for each of two unique yeast peptides",
       ref_yeast_log2=dict(zip(ANALYTES,refs.tolist())),
       no_imputation=True,batch_reference="All 40 samples; no per-batch factor",
       corrected_matrices={"PG":PGFILE,"PR":PRFILE},
       matrix_summary=matrix_meta,
       figures=["01_pg_before_after.png","02_pr_before_after.png",
         "03_yeast_scaling_factors_by_batch.png","04_yeast_cv_CIRCULAR.png",
         "05_endogenous_feature_CV_before_after.png"])
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n")
    summary=[
      "# Experimental yeast internal-standard normalization: 40 NYD samples","",
      "**Not the median normalization method described in MacCoss et al. 2026.**",
      "The paper used proteome/peptide median-normalized data and checked",
      "spiked-in yeast ENO1 as internal QC. The current test instead uses yeast",
      "intensities themselves to derive reference-based scaling factors.","",
      "## Method","",
      "PG factor: median yeast P00924 protein signal across all 40 samples",
      "divided by the sample's yeast P00924 signal.",
      "PR factor: first add different precursor charge states of EACH peptide",
      "on linear scale; then compute the reference/sample ratio independently",
      "for AADALLLK and VNQIGTLSESIK, and take their geometric mean.",
      "**The two DIFFERENT yeast peptides are never summed together.**",
      "All sample intensities in the corresponding matrix are multiplied",
      "by that sample factor. Raw source matrices remain intact.",
      "No missing-value imputation, no batch-specific factor.","",
      "## Files","",
      "- "+PGFILE+": yeast-P00924-normalized PG matrix.",
      "- "+PRFILE+": yeast-two-peptide-normalized PR matrix.",
      "- yeast_reference_factors.csv: two peptide ratios and PG/PR factor per sample.",
      "- endogenous_batch_CV.csv: raw and yeast-corrected endogenous PG/PR CV distributions.",
      "- yeast_CV_CIRCULAR.csv: raw and post yeast CV; cannot independently validate correction.",
      "- endogenous_sample_stats.csv and endogenous_batch_stats.csv: descriptive sample summaries.",
      "- run_manifest.json: source hashes and exact definition.","",
      "## PG: 40 sample boxplots before/after yeast protein reference correction","",
      "![PG before after](01_pg_before_after.png)","",
      "## PR: 40 sample boxplots before/after two yeast peptide reference correction","",
      "![PR before after](02_pr_before_after.png)","",
      "## Factor distributions by batch","",
      "![Factors](03_yeast_scaling_factors_by_batch.png)","",
      "## Yeast CV (CIRCULAR; not independent evidence)","",
      "By construction, the P00924 reference has a post-normalization CV",
      "near zero because its own signal determined the correction factor.",
      "This does NOT prove successful biological or technical correction.",
      "The same caution applies to the two peptides used to calculate the PR factor.","",
      "![Yeast CV circular](04_yeast_cv_CIRCULAR.png)","",
      "## Endogenous PG/PR feature CV by batch","",
      "![Endogenous CV](05_endogenous_feature_CV_before_after.png)","",
      "| Matrix | Batch | Raw median CV (%) | Yeast corrected median CV (%) | n features |",
      "| --- | --- | ---: | ---: | ---: |"]
    for label in ("pg","pr"):
        for batch in BATCHES:
            raw=next(r for r in cvstats if r["matrix"]==label and
                   r["batch"]==batch and r["stage"]=="raw")
            adj=next(r for r in cvstats if r["matrix"]==label and
                   r["batch"]==batch and r["stage"]=="median_normalized")
            summary.append(f'| {label.upper()} | {batch} | {raw["median_feature_cv_pct"]:.2f} | {adj["median_feature_cv_pct"]:.2f} | {raw["n_features_cv"]} |')
    summary += ["","## Limitations","",
       "The yeast spike-in amount must truly be equal across all samples.",
       "The two yeast peptides previously showed limited correlation, so a",
       "two-peptide PR reference can be noisy. Known T003 low-detection sample",
       "MAP26092326E01 is retained. Endogenous CV also contains biological",
       "differences; without repeated pooled QCs, apparent changes are not",
       "proof of batch-effect removal.",""]
    (OUT/"README.md").write_text("\n".join(summary),encoding="utf-8")
    print("CV_SUMMARY",json.dumps(cvstats,ensure_ascii=False),flush=True)
    print("YEAST_CV",json.dumps(ycv,ensure_ascii=False),flush=True)
    print("SUCCESS: two matrices / factors / 5 figures / endogenous assessments",flush=True)

if __name__=="__main__":
    main()
