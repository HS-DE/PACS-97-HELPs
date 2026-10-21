#!/usr/bin/env python3
"""NYD v2: one Pearson per Heavy HELP across samples against PG sample medians.
If Pearson >0, all proteins in each sample receive +fi; if Pearson <0,
all proteins receive -fi. NOT per-protein Pearson as in v1.
Main 40 originals (74 full-detection HELP expected), historical 39 sensitivity.
"""
import json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import nyd_all_HELP_Pearson_sign_directional_normalization_v1 as helper
from nyd_pg_batch_boxplots import PG,META,BATCHES,COLORS

ROOT=Path(__file__).resolve().parents[1]/"NYD_data"/"40samples"
OUT=ROOT/"HELP97_PGsampleMedian_globalPearson_v2"
HEAVY=ROOT/"Heavy.xlsx"

def pearson(h,med):
    if not (np.isfinite(h).all() and np.isfinite(med).all()):
        raise ValueError("Nonfinite Heavy or sample PG medians")
    if np.std(h)<1e-12 or np.std(med)<1e-12:return np.nan,0
    r=float(np.corrcoef(h,med)[0,1])
    return r,int(np.sign(r))

def trendchart(seq,h,med,applied,rec,path):
    x=np.arange(1,len(h)+1)
    fig,axs=plt.subplots(2,1,figsize=(19,9.8),sharex=True)
    fig.subplots_adjust(left=.07,right=.975,top=.85,bottom=.16,hspace=.34)
    for b in BATCHES:
        ix=np.array([i for i,a in enumerate(rec) if a["batch"]==b])
        axs[0].plot(x[ix],(h-np.median(h))[ix],color=COLORS[b],
                    marker="D",ls="--",lw=1.5,ms=4)
        axs[0].plot(x[ix],(med-np.median(med))[ix],color=COLORS[b],
                    marker="o",ls="-",lw=1.5,ms=4.5)
        axs[1].bar(x[ix],applied[ix],color=COLORS[b],alpha=.8)
    for ax in axs:
        helper.frame(ax);helper.add_batch_labels(ax,rec)
        ax.axhline(0,color="#8CA0AF",lw=.8)
    axs[0].set_title("Original PG median (solid circle) vs Heavy (dashed diamond)",
                     loc="left",pad=27,fontsize=11)
    axs[1].set_title("Applied offset = sign(global Pearson) * fi for all PG",
                     loc="left",pad=27,fontsize=11)
    axs[0].set_ylabel("Log2 deviation from own median")
    axs[1].set_ylabel("Applied log2 correction")
    axs[1].set_xticks(x,[r["sample_id"][-6:] for r in rec],rotation=90,fontsize=7)
    r,s=pearson(h,med)
    fig.suptitle("NYD "+seq+" | global Pearson r = {:+.4f} | sign = {:+d}".format(r,s),
                 fontsize=14,weight="bold")
    fig.savefig(path,dpi=120,facecolor="white",bbox_inches="tight")
    plt.close(fig)

def summaryplot(metrics,cohort):
    df=pd.DataFrame(metrics)
    fig,axes=plt.subplots(1,2,figsize=(14.5,6.2))
    fig.subplots_adjust(left=.075,right=.98,top=.85,bottom=.16,wspace=.3)
    axes[0].hist(df["overall_delta_medianCV_pctpt"],bins=20,
                 color="#538CA4",edgecolor="white")
    axes[0].axvline(0,c="#D26B50",ls="--",lw=1.4)
    axes[0].set_title("Distribution of median PG CV changes")
    axes[0].set_xlabel("Corrected - Raw PG CV (percentage points)")
    axes[0].set_ylabel("Number of HELP")
    for direction,color in [(1,"#248C77"),(-1,"#CC8740"),(0,"#8C919E")]:
        mask=df.pearson_direction==direction
        if mask.any():
            axes[1].scatter(df.loc[mask,"Pearson_Heavy_vs_samplePGmedian"],
                            df.loc[mask,"overall_delta_medianCV_pctpt"],
                            color=color,s=43,label="sign {:+d}".format(direction))
    axes[1].axhline(0,c="#D26B50",ls="--",lw=1.4)
    axes[1].set_xlabel("Pearson Heavy vs sample PG median")
    axes[1].set_ylabel("Corrected - Raw PG CV (pp)")
    axes[1].legend(frameon=False,fontsize=9)
    for ax in axes:helper.frame(ax)
    fig.suptitle("NYD "+cohort+" | one Pearson per HELP across clinical samples",
                 weight="bold")
    fig.savefig(OUT/cohort/"00_summary_all_HELP.png",dpi=165,bbox_inches="tight")
    plt.close(fig)

def run(cohort,recs,names,raw,anno,seqs,H,full_outputs):
    d=OUT/cohort
    d.mkdir(parents=True,exist_ok=False)
    n=len(names)
    originalPGmed=np.nanmedian(raw,axis=0)
    rawCV=helper.summaryCV(raw,recs)
    basis=helper.raw_pca_basis(raw) if full_outputs else None
    complete=np.flatnonzero((np.isfinite(H)&(H>0)).all(axis=1))
    results=[];factors=[]
    for num,k in enumerate(complete,1):
        seq=seqs[k];h=np.log2(H[k]);ref=float(np.median(h))
        fi=ref-h
        r,sign=pearson(h,originalPGmed)
        applied=fi*sign
        corrected=raw+applied[None,:]
        assert np.array_equal(np.isfinite(corrected),np.isfinite(raw))
        assert np.allclose(corrected-applied[None,:],raw,equal_nan=True,atol=1e-11)
        newPGmed=np.nanmedian(corrected,axis=0)
        assert np.allclose(newPGmed,originalPGmed+applied,atol=1e-11)
        corrCV=helper.summaryCV(corrected,recs)
        m={"HELP_sequence":seq,"n_samples":n,"global_median_Heavy_log2":ref,
           "Pearson_Heavy_vs_samplePGmedian":r,"pearson_direction":sign,
           "fi_min_log2":float(fi.min()),"fi_max_log2":float(fi.max()),
           "applied_min_log2":float(applied.min()),
           "applied_max_log2":float(applied.max()),
           "Raw_PGmedian_SD":float(np.std(originalPGmed,ddof=1)),
           "Corrected_PGmedian_SD":float(np.std(newPGmed,ddof=1))}
        for b in ["ALL"]+list(BATCHES):
            m[b+"_Raw_medianCV_pct"]=rawCV[b]["median_CV_pct"]
            m[b+"_Corrected_medianCV_pct"]=corrCV[b]["median_CV_pct"]
            m[b+"_delta_medianCV_pctpt"]=corrCV[b]["median_CV_pct"]-rawCV[b]["median_CV_pct"]
            m[b+"_n_PG_for_CV"]=corrCV[b]["n_proteins"]
        m["overall_delta_medianCV_pctpt"]=m["ALL_delta_medianCV_pctpt"]
        results.append(m)
        for i,sid in enumerate(names):
            factors.append({
                "HELP_sequence":seq,"sample_id":sid,"batch":recs[i]["batch"],
                "Heavy_log2":float(h[i]),"median_hi_log2":float(h[i]),
                "median_h_log2":ref,"fi_log2":float(fi[i]),
                "Pearson_global_across_samples":r,"Pearson_direction":sign,
                "applied_offset_log2":float(applied[i]),
                "Raw_sample_PGmedian_log2":float(originalPGmed[i]),
                "Corrected_sample_PGmedian_log2":float(newPGmed[i])})
        if full_outputs:
            target=d/"per_HELP"/seq
            target.mkdir(parents=True,exist_ok=False)
            table=anno.copy()
            for i,sid in enumerate(names):table[sid]=corrected[:,i]
            table.to_csv(target/"corrected_PG_log2.tsv.gz",sep="\t",index=False,
                         na_rep="",compression={"method":"gzip","compresslevel":4})
            helper.chart_box(seq,raw,corrected,recs,
                             target/"01_Raw_vs_Corrected_boxplots.png")
            trendchart(seq,h,originalPGmed,applied,recs,
                       target/"02_Pearson_Heavy_vs_PGmedian_and_factors.png")
            helper.chart_pca(seq,basis,helper.project_pca(corrected,basis),
                             recs,target/"03_Raw_vs_Corrected_common_PCA.png")
        if num%10==0 or num==len(complete):
            print(cohort,"COMPLETED",num,"/",len(complete),flush=True)
    df=pd.DataFrame(results).sort_values("overall_delta_medianCV_pctpt")
    df.to_csv(d/"01_global_Pearson_HELP_ranking_and_CV.csv",index=False)
    pd.DataFrame(factors).to_csv(d/"02_every_HELP_sample_factors.csv.gz",
          index=False,compression={"method":"gzip","compresslevel":4})
    summaryplot(results,cohort)
    return {"n_samples":n,"n_full_HELP":len(complete),
       "n_positive_r":int((df.pearson_direction==1).sum()),
       "n_negative_r":int((df.pearson_direction==-1).sum()),
       "n_other_r":int((df.pearson_direction==0).sum()),
       "n_improved_allPG_medianCV":int((df.overall_delta_medianCV_pctpt<0).sum()),
       "batches":{b:sum(q["batch"]==b for q in recs) for b in BATCHES},
       "raw_CV":rawCV,
       "best8":df.head(8).to_dict(orient="records"),
       "all_sequences":[seqs[k] for k in complete]}

def main():
    if OUT.exists():raise FileExistsError(str(OUT))
    rec,names,raw,anno,seqs,H=helper.load()
    assert raw.shape==(4264,40) and len(seqs)==97
    OUT.mkdir(parents=True,exist_ok=False)
    a=run("full40",rec,names,raw,anno,seqs,H,True)
    ix=[i for i,s in enumerate(names) if s!=helper.EXCLUDE]
    assert len(ix)==39
    b=run("exclude_T003_bad_sample_39",
          [rec[i] for i in ix],[names[i] for i in ix],raw[:,ix],
          anno,seqs,H[:,ix],False)
    summary={
        "method":"One across-sample Pearson PER HELP vs sample PG median, uniform direction per HELP",
        "full40":a,"sensitivity39":b,"n_protein_groups":4264,"input_HELP":97,
        "formula_fi":"median(log2 Heavy_j over cohort) - log2 Heavy_ji",
        "formula_sign":"Pearson(log2 Heavy_j across samples, PG sample median log2 across same samples)",
        "formula_corrected":"PG_log2_gi + sign(r_j)*fi for every protein g in sample i",
        "zero_nan_r":"no PG correction",
        "clinical_biology_caveat":"Pearson against same Study PG medians can encode biology and creates no independent validation",
        "raw_PG_missingness":"unchanged; positive original linear intensities only",
        "sensitivity_note":"40 main; 39 refit excluding historical T003 poor-quality sample",
        "cv_note":"Median per-PG CV on original linear scale, minimum 70% detection in each cohort/batch",
        "PCA_note":"One raw-fit PCA basis; same projection before/after, per sample common shift",
        "source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (PG,HEAVY,META)},
        "v1_per_protein_Pearson":"kept unchanged separately"
    }
    (OUT/"run_manifest.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    readme=["# NYD: one global Pearson per Heavy HELP (v2)","",
      "Each HELP provides a 40-observation Heavy log2 vector. Each sample provides",
      "one original log2 PG median, computed across its own detected PGs.",
      "Pearson across these TWO 40-value vectors yields ONE coefficient per HELP.",
      "It is NOT per-protein Pearson and NOT per-sample Pearson.","",
      "- fi_i = median_all_samples(log2 Heavy_j) - log2 Heavy_ji.",
      "- r_j > 0: ALL proteins in sample i receive +fi_i.",
      "- r_j < 0: ALL proteins in sample i receive -fi_i.",
      "- Undefined/zero r: no correction.","",
      "## Main full40 results","",
      "- 40/40 detected Heavy HELP: **"+str(a["n_full_HELP"])+" / 97**.",
      "- Positive Pearson HELP: "+str(a["n_positive_r"])+
        "; negative: "+str(a["n_negative_r"])+
        "; zero or undefined: "+str(a["n_other_r"])+".",
      "- Better overall clinical median PG CV: "+
        str(a["n_improved_allPG_medianCV"])+" / "+str(a["n_full_HELP"])+".",
      "","| HELP | Pearson r | sign | Raw median CV | Corrected median CV | change pp |",
      "|---|---:|---:|---:|---:|---:|"]
    for t in a["best8"]:
        readme.append("| {} | {:+.4f} | {:+d} | {:.2f}% | {:.2f}% | {:+.2f} |".format(
            t["HELP_sequence"],t["Pearson_Heavy_vs_samplePGmedian"],t["pearson_direction"],
            t["ALL_Raw_medianCV_pct"],t["ALL_Corrected_medianCV_pct"],
            t["overall_delta_medianCV_pctpt"]))
    readme += ["","## Historical exclusion sensitivity (39 samples)","",
      "- Full detection 39/39: "+str(b["n_full_HELP"])+" / 97.",
      "- Better overall median PG CV: "+str(b["n_improved_allPG_medianCV"])+" / "+str(b["n_full_HELP"])+".",
      "","## File layout","",
      "- full40/01_global_Pearson_HELP_ranking_and_CV.csv",
      "- full40/02_every_HELP_sample_factors.csv.gz",
      "- full40/00_summary_all_HELP.png",
      "- full40/per_HELP/SEQUENCE/corrected_PG_log2.tsv.gz",
      "- full40/per_HELP/SEQUENCE/01_Raw_vs_Corrected_boxplots.png",
      "- full40/per_HELP/SEQUENCE/02_Pearson_Heavy_vs_PGmedian_and_factors.png",
      "- full40/per_HELP/SEQUENCE/03_Raw_vs_Corrected_common_PCA.png",
      "- exclude_T003_bad_sample_39/01_global_Pearson_HELP_ranking_and_CV.csv",
      "- exclude_T003_bad_sample_39/02_every_HELP_sample_factors.csv.gz",
      "","## Caveats","",
      "No identical-composition pooled QC in NYD. CV over clinical samples includes biology.",
      "Pearson derived using clinical PG sample medians is not an independent technical reference.",
      "A single per-sample additive log2 shift leaves sample-to-sample protein-profile",
      "Pearson correlations invariant; PCA, sample median and CV are the useful diagnostics.",
      "Previously completed per-protein Pearson v1 remains unchanged."]
    (OUT/"README.md").write_text("\n".join(readme)+"\n")
    print("RESULT",json.dumps({
       "40":{k:a[k] for k in ("n_full_HELP","n_positive_r","n_negative_r","n_improved_allPG_medianCV")},
       "39":{k:b[k] for k in ("n_full_HELP","n_positive_r","n_negative_r","n_improved_allPG_medianCV")},
       "best8":a["best8"]},ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
