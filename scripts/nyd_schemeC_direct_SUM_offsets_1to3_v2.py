#!/usr/bin/env python3
"""NYD Scheme C SUM-ONLY v2.

Goal: match sample PG log2 median deviation by a direct unweighted SUM of
each included HELP peptide's log2 Heavy deviation from its cohort median.

A peptide is either included (+1) or not included (0). No average, no learned
weight, no regression scaling, no per-peptide sign flip, and no recentering
of the final sum. Evaluate 1, 2, and 3 distinct peptides exhaustively.
Separate n=40 and outlier-excluded n=39, both independently recomputed.
No protein correction is applied by this program.
"""
from __future__ import annotations

import gzip
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import nyd_all_HELP_Pearson_sign_directional_normalization_v1 as source
from nyd_pg_batch_boxplots import BATCHES, COLORS

OUT=Path(__file__).resolve().parents[1]/"NYD_data"/"40samples"/"SchemeC_v2_directSUM_Heavy_offsets_1to3"
CHUNK=3000
KMAX=3

def stream_indices(n,k):
    it=itertools.combinations(range(n),k)
    while True:
        b=list(itertools.islice(it,CHUNK))
        if not b:break
        yield np.array(b,dtype=np.int16)

def score(pred,target):
    # pred: (#comb, #samples), sum of centerings, NO final centering.
    d=pred-target[None,:]
    rmse=np.sqrt(np.mean(d*d,axis=1))
    x=pred-pred.mean(axis=1,keepdims=True)
    y=target-target.mean()
    denom=np.sqrt(np.sum(x*x,axis=1)*np.sum(y*y))
    r=np.divide(x@y,denom,out=np.full(len(pred),np.nan),where=denom>1e-10)
    sd=np.std(pred,axis=1,ddof=1)/np.std(target,ddof=1)
    inverted=np.sqrt(np.mean((pred+target[None,:])**2,axis=1))
    med_shift=np.median(pred,axis=1)
    return rmse,r,sd,inverted,med_shift

def search_one(recs,names,raw,seqs,H,cohort):
    n=len(names)
    if n not in (39,40):raise AssertionError("Cohort population not expected")
    sub=OUT/cohort
    sub.mkdir(parents=True,exist_ok=False)
    idx=np.where(np.isfinite(H).all(axis=1)&(H>0).all(axis=1))[0]
    peptides=[seqs[i] for i in idx]
    h=np.log2(H[idx])
    offset=h-np.median(h,axis=1,keepdims=True)
    protein_medians=np.nanmedian(raw,axis=0)
    target=protein_medians-np.median(protein_medians)
    if not (np.isfinite(offset).all() and np.isfinite(target).all()):
        raise AssertionError("Invalid inputs")
    baseline=float(np.sqrt(np.mean(target*target)))
    sample_index=np.arange(n)
    pd.DataFrame({"order":sample_index+1,"sample_id":names,
          "batch":[r["batch"] for r in recs],
          "rawPG_sampleMedian_log2":protein_medians,
          "target_PG_medianDeviation_log2":target}).to_csv(
             sub/"00_sample_PG_median_deviations.csv",index=False)

    # Using compact gzip avoids huge files in GitHub without omitting combos.
    path=sub/"01_all_combinations_1to3_unweighted_SUM.csv.gz"
    saved=[]
    stats=[]
    first=True
    with gzip.open(path,"wt",encoding="utf-8",newline="",compresslevel=5) as handle:
        for k in (1,2,3):
            n_this=0
            for indices in stream_indices(len(peptides),k):
                # Each Heavy offset gets coefficient EXACTLY 1.
                predicted=offset[indices].sum(axis=1)
                rmse,r,sd,inv,med=score(predicted,target)
                out=pd.DataFrame({
                    "n_HELP":np.full(len(indices),k),
                    "HELP_1":[peptides[t[0]] for t in indices],
                    "HELP_2":[peptides[t[1]] if k>=2 else "" for t in indices],
                    "HELP_3":[peptides[t[2]] if k>=3 else "" for t in indices],
                    "Pearson_r":r,
                    "RMSE_log2":rmse,
                    "RMSE_minus_zero_predictor":rmse-baseline,
                    "sum_SD_divided_target_SD":sd,
                    "median_of_combination_sum_log2":med,
                    "inverted_whole_SUM_RMSE_diagnostic":inv
                })
                if not np.isfinite(out[["Pearson_r","RMSE_log2",
                                      "sum_SD_divided_target_SD"]].to_numpy()).all():
                    raise AssertionError("Invalid combination scores")
                out.to_csv(handle,index=False,header=first)
                first=False
                n_this+=len(out)
                # Retain the strongest candidates from every chunk under
                # multiple criteria without selecting weights or directions.
                p=out[out.Pearson_r>0]
                if len(p):
                    saved.extend(p.nsmallest(min(12,len(p)),"RMSE_log2").to_dict("records"))
                    saved.extend(p.nlargest(min(12,len(p)),"Pearson_r").to_dict("records"))
                saved.extend(out.nsmallest(min(6,len(out)),"RMSE_log2").to_dict("records"))
                invgood=out[out.Pearson_r<0]
                if len(invgood):
                    saved.extend(invgood.nsmallest(
                        min(6,len(invgood)),"inverted_whole_SUM_RMSE_diagnostic"
                    ).to_dict("records"))
                # Keep only nondominated candidates at chunk boundaries.
                if len(saved)>5000:
                    cand=pd.DataFrame(saved).drop_duplicates(subset=[
                       "n_HELP","HELP_1","HELP_2","HELP_3"])
                    selected=[]
                    for kk in (1,2,3):
                        part=cand[cand.n_HELP==kk]
                        pos=part[part.Pearson_r>0]
                        neg=part[part.Pearson_r<0]
                        selected.append(part.nsmallest(min(60,len(part)),"RMSE_log2"))
                        selected.append(pos.nsmallest(min(60,len(pos)),"RMSE_log2"))
                        selected.append(pos.nlargest(min(60,len(pos)),"Pearson_r"))
                        selected.append(neg.nsmallest(min(25,len(neg)),
                                                   "inverted_whole_SUM_RMSE_diagnostic"))
                    saved=pd.concat(selected,ignore_index=True).drop_duplicates(
                           subset=["n_HELP","HELP_1","HELP_2","HELP_3"]).to_dict("records")
            expected=math.comb(len(peptides),k)
            assert n_this==expected,(n_this,expected)
            stats.append({"HELP_count":k,"n_combinations":n_this})
            print("SEARCH_FINISHED",cohort,k,n_this,flush=True)
    # ALL combos are saved. Read back to guarantee genuine global rankings
    # rather than using chunk-candidate approximations.
    full=pd.read_csv(path,keep_default_na=False)
    expected_total=sum(math.comb(len(peptides),k) for k in (1,2,3))
    assert len(full)==expected_total
    assert full[["HELP_1","HELP_2","HELP_3"]].apply(
           lambda x:len(set([z for z in x if z]))==len([z for z in x if z]),
           axis=1).all()
    global_positive=full[full.Pearson_r>0].copy()
    global_negative=full[full.Pearson_r<0].copy()
    ranked=[]
    by_k={}
    for k in (1,2,3):
        d=full[full.n_HELP==k]
        p=d[d.Pearson_r>0].copy()
        q=d[d.Pearson_r<0].copy()
        by_k[str(k)]={"tested":len(d),"positive":len(p),"negative":len(q),
                     "RMSE_below_zero_baseline":int((d.RMSE_log2<baseline).sum())}
        table_items={
           "all_signs_min_RMSE":d.sort_values(
                         ["RMSE_log2","Pearson_r"],ascending=[True,False]).head(30),
           "positive_min_RMSE":p.sort_values(
                         ["RMSE_log2","Pearson_r"],ascending=[True,False]).head(30),
           "positive_max_Pearson":p.sort_values(
                         ["Pearson_r","RMSE_log2"],ascending=[False,True]).head(30),
           "negative_whole_sum_inverted_diagnostic":q.sort_values(
                        ["inverted_whole_SUM_RMSE_diagnostic","Pearson_r"],
                        ascending=[True,True]).head(30)
        }
        for name,top in table_items.items():
            top.to_csv(sub/("02_"+str(k)+"HELP_"+name+"_top30.csv"),index=False)
            if len(top):
                r=top.iloc[0].to_dict()
                ranked.append({"criterion":name,**r})
    winners=pd.DataFrame(ranked)
    winners.to_csv(sub/"03_winners_all_sizes_criteria.csv",index=False)

    # Compute a non-optimized control: zero offset (the raw sample medians).
    # Any RMSE >= this baseline has not predicted sample deviations better
    # than doing nothing in a squared-error sense.
    comparison={
       "baseline_zero_prediction_RMSE_log2":baseline,
       "n_eligible":len(peptides),
       "n_samples":n,
       "n_combinations":len(full),
       "by_size":by_k,
       "n_total_positive":len(global_positive),
       "n_total_negative":len(global_negative),
       "best_positive_RMSE":(
           global_positive.nsmallest(10,"RMSE_log2").to_dict("records") if len(global_positive) else []),
       "best_positive_Pearson":(
           global_positive.nlargest(10,"Pearson_r").to_dict("records") if len(global_positive) else []),
       "best_unconstrained_RMSE":full.nsmallest(10,"RMSE_log2").to_dict("records"),
       "best_negative_global_inverted_diagnostic":(
           global_negative.nsmallest(10,"inverted_whole_SUM_RMSE_diagnostic").to_dict("records")
           if len(global_negative) else []),
       "size_winners":winners.to_dict("records"),
       "full_detected_HELP":peptides,
    }
    # Make top result plots per size; separate Pearson-maximum and
    # RMSE-minimum and do not mislabel opposite-trend candidates.
    make_plots(sub,comparison,offset,peptides,target,recs,names)
    (sub/"summary.json").write_text(json.dumps(comparison,ensure_ascii=False,indent=2)+"\n")
    return comparison

def make_plots(out,report,offset,peptides,target,recs,names):
    x=np.arange(1,len(names)+1)
    rank=pd.DataFrame(report["size_winners"])
    fig,axes=plt.subplots(3,1,figsize=(20,14.8),sharex=True)
    fig.subplots_adjust(left=.08,right=.985,top=.85,bottom=.19,hspace=.58)
    for ax,k in zip(axes,(1,2,3)):
        win=rank[(rank.n_HELP==k)&(rank.criterion=="positive_min_RMSE")]
        if len(win)==0:
            ax.text(.5,.5,f"No positive Pearson for k={k}",ha="center",
                    transform=ax.transAxes)
            continue
        row=win.iloc[0]
        ids=[peptides.index(row[col]) for col in ("HELP_1","HELP_2","HELP_3")
             if row[col]!=""]
        z=offset[ids].sum(axis=0)
        # In particular NO subtracting median(z) and NO dividing by k.
        assert np.isclose(np.sqrt(np.mean((z-target)**2)),row.RMSE_log2,
                          atol=1e-9)
        for b in BATCHES:
            ix=np.array([j for j,r in enumerate(recs) if r["batch"]==b])
            ax.plot(x[ix],target[ix],marker="o",ms=4.5,ls="-",lw=1.45,color=COLORS[b])
            ax.plot(x[ix],z[ix],marker="D",ms=3.7,ls="--",lw=1.5,color=COLORS[b])
        for a in range(1,len(recs)):
            if recs[a]["batch"]!=recs[a-1]["batch"]:
                ax.axvline(a+.5,color="#8B9CA9",ls=":",lw=.85)
        ax.axhline(0,color="#A1B2BD",lw=.8)
        ax.set_ylabel("Signed log2 deviation")
        ax.set_title(f"k={k}: "+" + ".join([peptides[j] for j in ids])+
            f"\nPearson r={row.Pearson_r:+.4f}  |  SUM RMSE={row.RMSE_log2:.4f}"+
            f"  |  Zero-offset baseline={report['baseline_zero_prediction_RMSE_log2']:.4f}",
            loc="left",pad=8,fontsize=10.7)
        ax.grid(axis="y",alpha=.12)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
    axes[-1].set_xticks(x,[s[-6:] for s in names],rotation=90,fontsize=7)
    axes[-1].set_xlabel("Samples ordered by T001 T002 T003 T004, then metadata order")
    lg=[Line2D([0],[0],color=COLORS[b],lw=3,label=b) for b in BATCHES]
    lg.extend([Line2D([0],[0],color="#556978",marker="o",ls="-",lw=1.5,
                      label="Sample PG median deviation"),
               Line2D([0],[0],color="#556978",marker="D",ls="--",lw=1.5,
                      label="Direct SUM HELP deviation")])
    fig.legend(handles=lg,loc="lower center",bbox_to_anchor=(.5,.015),
               frameon=False,ncol=6,fontsize=9)
    fig.suptitle(f"NYD Scheme C DIRECT SUM | best positive Pearson combinations by RMSE | n={len(names)}",
                 fontsize=14,weight="bold")
    fig.savefig(out/"04_best_same_direction_1_2_3HELP_direct_SUM.png",
                dpi=175,bbox_inches="tight",facecolor="white")
    plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(15.4,6.2))
    fig.subplots_adjust(left=.075,right=.985,top=.85,bottom=.16,wspace=.27)
    for k in (1,2,3):
        win=rank[(rank.n_HELP==k)&(rank.criterion=="positive_min_RMSE")]
        if len(win):
            row=win.iloc[0]
            axs[0].scatter(row.Pearson_r,row.RMSE_log2,
                s=140,label=f"k={k} selected",edgecolors="white",linewidths=.7)
            axs[0].annotate(str(k),xy=(row.Pearson_r,row.RMSE_log2),
                            xytext=(7,5),textcoords="offset points")
    axs[0].axhline(report["baseline_zero_prediction_RMSE_log2"],
        lw=1.35,color="#B66749",ls="--",label="Zero HELP offset")
    axs[0].set_ylabel("RMSE signed log2")
    axs[0].set_xlabel("Original Pearson r (positive only)")
    axs[0].set_title("Best same-direction combination per size",loc="left")
    axs[0].legend(fontsize=9,frameon=False)
    axs[1].bar([1,2,3],[report["by_size"][str(k)]["positive"] for k in (1,2,3)],
                 color=["#5684A4","#26947E","#D39B5C"],width=.62)
    axs[1].set_xticks([1,2,3])
    axs[1].set_xlabel("Number of directly summed HELP peptides")
    axs[1].set_ylabel("Count")
    axs[1].set_title("Combinations with positive Pearson",loc="left")
    for ax in axs:
        ax.grid(axis="y",alpha=.15)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("NYD Scheme C direct-SUM exhaustive search: shape and amplitude",weight="bold")
    fig.savefig(out/"05_SUM_search_overview.png",dpi=180,bbox_inches="tight",facecolor="white")
    plt.close(fig)


def main():
    if OUT.exists():raise FileExistsError("Do not overwrite an existing run: "+str(OUT))
    rec,names,pg,annotation,seqs,H=source.load()
    assert len(names)==40 and len(seqs)==97 and pg.shape==(4264,40)
    OUT.mkdir(parents=True,exist_ok=False)
    a=search_one(rec,names,pg,seqs,H,"full40")
    keep=[i for i,s in enumerate(names) if s!=source.EXCLUDE]
    assert len(keep)==39
    b=search_one([rec[i] for i in keep],[names[i] for i in keep],pg[:,keep],
                 seqs,H[:,keep],"exclude_T003_outlier_39")
    summary={"method":"direct SUM of individual signed log2 Heavy deviations",
             "formula_for_HELP_j":"delta_ji=log2 Heavy_ji - median_over_samples(log2 Heavy_j)",
             "formula_combination":"C_i = SUM_j(delta_ji) with coefficient +1 for all selected j, no division",
             "formula_target":"Y_i = median_across_detectedPG(log2 PG_i) - median_over_samples(samplePGmedians)",
             "formula_RMSE":"sqrt(mean_i((C_i-Y_i)^2)), no center of C afterwards",
             "pearson":"Pearson_i(C_i,Y_i), sign preserved as originally measured",
             "no_weights":True,"no_averaging":True,"no_amplitude_fitting":True,
             "no_sign_flip_per_HELP":True,"no_combination_recenter":True,
             "no_calibration_applied":True,"max_n_HELP":KMAX,
             "main40":a,"sensitivity39":b,
             "known_T003_outlier":"Main 40 includes original MAP26092326E01; 39 independently excludes it",
             "warning":"Fitting to clinical PG medians is not independent demonstration of technical effect"}
    (OUT/"run_manifest.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    doc=["# NYD Scheme C v2: direct SUM of Heavy deviations (no weights / no mean)",
         "","For each Heavy peptide, subtract its own across-sample log2 median.",
         "For the chosen combination, add peptide deviations directly with coefficient 1.",
         "Neither divide by k nor recenter the final combination, nor fit scale/weights/sign.",
         "Compare to each sample's own detected-PG median deviation (40 or 39 samples).",
         "","Only 1, 2, 3 HELP are exhaustively searched. This is a screening experiment",
         "not proof of independent technical correction.",
         ""]
    for title,res in (("full40",a),("exclude_T003_outlier_39",b)):
        doc.extend(["## "+title,
             "","Full detected: "+str(res["n_eligible"])+"/97; combinations: "+
                str(res["n_combinations"])+".",
             "","Zero-offset baseline RMSE: {:.5f} log2.".format(
                        res["baseline_zero_prediction_RMSE_log2"]),
             "","| # HELP | # tested | same direction # | best positive r | best positive RMSE log2 | winner |",
             "|---:|---:|---:|---:|---:|---|"])
        for k in (1,2,3):
            best=[w for w in res["size_winners"]
                     if w["n_HELP"]==k and w["criterion"]=="positive_min_RMSE"]
            if best:
                w=best[0];names2=" + ".join([w[c] for c in (
                       "HELP_1","HELP_2","HELP_3") if w[c]])
                doc.append("| {} | {} | {} | {:+.4f} | {:.5f} | {} |".format(
                    k,res["by_size"][str(k)]["tested"],
                    res["by_size"][str(k)]["positive"],
                    w["Pearson_r"],w["RMSE_log2"],names2))
            else:
                doc.append("| {} | {} | {} | - | - | (none) |".format(
                  k,res["by_size"][str(k)]["tested"],
                  res["by_size"][str(k)]["positive"]))
        doc.extend(["",f"![Top 1-3 trends]({title}/04_best_same_direction_1_2_3HELP_direct_SUM.png)",
                    f"![Overview]({title}/05_SUM_search_overview.png)",""])
    doc.extend(["## Result files","","- run_manifest.json: formulas and full winners",
         "- COHORT/00_sample_PG_median_deviations.csv",
         "- COHORT/01_all_combinations_1to3_unweighted_SUM.csv.gz",
         "- COHORT/02_*HELP_*_top30.csv (RMSE, Pearson, reverse diagnostic)",
         "- COHORT/03_winners_all_sizes_criteria.csv",
         "- COHORT/04_best_same_direction_1_2_3HELP_direct_SUM.png",
         "- COHORT/05_SUM_search_overview.png",
         "","## Limits","",
         "NYD uses different clinical samples, not pooled technical replicates.",
         "Pearson is shape similarity; RMSE penalizes amplitude and offset.",
         "The score of a zero-combination predictor is reported as an essential sanity check.",
         "Negative-correlation combinations cannot be called same-direction;",
         "the optional opposite-signed diagnostic does not influence primary search.",
         "40 and 39 samples have different eligibility and must not be conflated.",
         "Older Scheme C averaged HELP values. This NEW experiment uses SUM only."])
    (OUT/"README.md").write_text("\n".join(doc)+"\n")
    payload=lambda z:{"sample_count":z["n_samples"],"eligible":z["n_eligible"],
             "combinations":z["n_combinations"],
             "baseline":z["baseline_zero_prediction_RMSE_log2"],
             "by_size":z["by_size"],"winners":[w for w in z["size_winners"]
                 if w["criterion"]=="positive_min_RMSE"]}
    print("RESULT_NYD_SCHEME_C_SUM",json.dumps({
        "full40":payload(a),"outlier_excluded39":payload(b)},ensure_ascii=False),flush=True)

if __name__=="__main__":main()
