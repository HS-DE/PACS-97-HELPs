#!/usr/bin/env python3
"""NYD: screen fully-detected Heavy HELP peptide trajectories against each
sample's own detected raw log2 PG median trajectory.

Pearson and median-centered RMSE over identical samples, not per-PG or
per-sample Pearson. Rank positive (same-direction) HELP separately from
negative (opposite-direction) HELP; reversed candidate curves are diagnostic
only and explicitly labeled. Main original 40 samples and independent
exclude-T003-outlier 39-sample sensitivity. No calibration on this run.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import nyd_all_HELP_Pearson_sign_directional_normalization_v1 as data
from nyd_pg_batch_boxplots import BATCHES, COLORS

OUT=Path(__file__).resolve().parents[1]/"NYD_data"/"40samples"/"HELP97_Heavy_vs_samplePGmedian_Pearson_RMSE_screen"
TARGET_NAME="Each sample median of its own detected original log2 PG"


def metrics(h,y):
    hc=h-np.median(h)
    yc=y-np.median(y)
    if not (np.isfinite(h).all() and np.isfinite(y).all()):
        raise AssertionError("Missing target or HELP")
    rr=float(np.corrcoef(h,y)[0,1])
    rmse=float(np.sqrt(np.mean((hc-yc)**2)))
    inv=float(np.sqrt(np.mean((-hc-yc)**2)))
    return {
        "Pearson_r":rr,
        "median_centered_RMSE_log2":rmse,
        "SD_Heavy_over_PGmedian":float(np.std(h,ddof=1)/np.std(y,ddof=1)),
        "inverted_Heavy_RMSE_log2_diagnostic":inv,
        "Pearson_abs":abs(rr),
        "trend_same_direction":bool(rr>0),
        "median_abs_deviation_between_curves":float(np.median(abs(hc-yc))),
    }


def plots_scatter(df,cohort):
    fig,axs=plt.subplots(1,2,figsize=(15.5,6.2))
    fig.subplots_adjust(left=.08,right=.985,top=.85,bottom=.17,wspace=.29)
    neg=df[df.Pearson_r<0]
    pos=df[df.Pearson_r>0]
    axs[0].scatter(neg.Pearson_r,neg.median_centered_RMSE_log2,
        s=53,alpha=.78,c="#C88548",label="Opposite trend (Pearson <0)")
    axs[0].scatter(pos.Pearson_r,pos.median_centered_RMSE_log2,
        s=57,alpha=.85,c="#268F83",label="Same-direction (Pearson >0)")
    for k in ("Pearson_positive_rank_1","RMSE_positive_rank_1","combined_positive_rank_1"):
        lab=df[df[k]].head(1)
        for _,r in lab.iterrows():
            axs[0].annotate(r.HELP_sequence,(r.Pearson_r,
                   r.median_centered_RMSE_log2),fontsize=8,
                   xytext=(5,5),textcoords="offset points")
    axs[0].axvline(0,c="#8493A2",ls="--",lw=.85)
    axs[0].set_xlabel("Pearson r (Heavy vs sample PG median)")
    axs[0].set_ylabel("Median-centered RMSE (log2)")
    axs[0].set_title("Original Heavy trend (no sign flip)",loc="left")
    axs[0].set_yscale("log")
    axs[0].legend(frameon=False,fontsize=8)
    axs[1].scatter(df.Pearson_abs,
          np.where(df.Pearson_r>0,df.median_centered_RMSE_log2,
                   df.inverted_Heavy_RMSE_log2_diagnostic),
          c=np.where(df.Pearson_r>0,"#268F83","#C88548"),s=52,
          edgecolors="white",linewidths=.5,alpha=.84)
    axs[1].set_xlabel("Absolute Pearson |r| (direction ignored)")
    axs[1].set_ylabel("Aligned (sign-inverted if needed) RMSE log2")
    axs[1].set_title("Diagnostic only: negative Heavy mirrored",loc="left")
    axs[1].set_yscale("log")
    for ax in axs:
        ax.grid(alpha=.15)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(f"NYD {cohort}: screening all fully detected HELP by trend",weight="bold",fontsize=14)
    name="01_scatter_all_HELP_Pearson_RMSE.png"
    fig.savefig(OUT/cohort/name,dpi=180,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    return name


def plot_top_curves(hall,targets,records,names,seqs,chosen,cohort,reversed_mode=False):
    subset=chosen.head(min(5,len(chosen)))
    if subset.empty:return None
    n=len(records);x=np.arange(1,n+1)
    y=targets-np.median(targets)
    fig,axes=plt.subplots(len(subset),1,figsize=(20,3.3*len(subset)+2),sharex=True)
    fig.subplots_adjust(top=.87,bottom=.17,left=.08,right=.985,hspace=.54)
    if len(subset)==1:axes=[axes]
    for ax,(_,r) in zip(axes,subset.iterrows()):
        h=hall[seqs.index(r.HELP_sequence)]
        z=h-np.median(h)
        if reversed_mode:z=-z
        for batch in BATCHES:
            idx=np.array([j for j,rec in enumerate(records) if rec["batch"]==batch])
            ax.plot(x[idx],y[idx],color=COLORS[batch],ls="-",marker="o",ms=4.5,lw=1.45)
            ax.plot(x[idx],z[idx],color=COLORS[batch],ls="--",marker="D",ms=3.9,lw=1.35)
        for p in (16,24,32):
            if p<n:ax.axvline(p+.5,c="#8295A5",ls=":",lw=.9)
        ax.axhline(0,c="#9EAEB9",lw=.7)
        rmse=(r.inverted_Heavy_RMSE_log2_diagnostic if reversed_mode else r.median_centered_RMSE_log2)
        ax.set_title(f"{r.HELP_sequence} | original r={r.Pearson_r:+.4f} | "
                    f"{'inverted' if reversed_mode else 'original'} RMSE={rmse:.4f} log2"
                    f" | Heavy/PG SD={r.SD_Heavy_over_PGmedian:.2f}",
                    fontsize=10.1,loc="left",pad=8)
        ax.set_ylabel("Centered log2")
        ax.grid(axis="y",alpha=.12)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
    axes[-1].set_xticks(x,[rec["sample_id"][-6:] for rec in records],
                       rotation=90,fontsize=7)
    axes[-1].set_xlabel("Samples ordered T001-T004 (metadata order within batch)")
    leg=[Line2D([0],[0],color=COLORS[b],lw=3,label=b) for b in BATCHES]
    leg.extend([Line2D([0],[0],c="#566F80",ls="-",marker="o",lw=1.3,label="Raw PG median"),
                Line2D([0],[0],c="#566F80",ls="--",marker="D",lw=1.3,
                    label="Heavy log2"+(" (inverted)" if reversed_mode else ""))])
    fig.legend(handles=leg,loc="lower center",bbox_to_anchor=(.5,.012),frameon=False,ncol=6,fontsize=9)
    fig.suptitle(f"NYD {cohort}: {'inverted negative-correlation' if reversed_mode else 'positive same-direction'} Heavy vs sample PG median",fontsize=14,weight="bold")
    png="02_top5_opposite_inverted_diagnostic.png" if reversed_mode else "02_top5_positive_same_direction.png"
    fig.savefig(OUT/cohort/png,dpi=175,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    return png


def ranks(df):
    a=df.copy()
    pos=a[a.Pearson_r>0].copy()
    neg=a[a.Pearson_r<0].copy()
    if not pos.empty:
        pos["rank_r"]=pos.Pearson_r.rank(ascending=False,method="min")
        pos["rank_rmse"]=pos.median_centered_RMSE_log2.rank(ascending=True,method="min")
        pos["rank_sum"]=pos.rank_r+pos.rank_rmse
        pos=pos.sort_values(["rank_sum","median_centered_RMSE_log2","Pearson_r"],
                            ascending=[True,True,False],kind="stable")
        map_rank=dict(zip(pos.HELP_sequence,range(1,len(pos)+1)))
        a["positive_joint_rank"]=a.HELP_sequence.map(map_rank)
        by_r=pos.sort_values(["Pearson_r","median_centered_RMSE_log2"],
                             ascending=[False,True],kind="stable")
        by_rmse=pos.sort_values(["median_centered_RMSE_log2","Pearson_r"],
                                ascending=[True,False],kind="stable")
        a["Pearson_positive_rank_1"]=a.HELP_sequence.eq(by_r.iloc[0].HELP_sequence)
        a["RMSE_positive_rank_1"]=a.HELP_sequence.eq(by_rmse.iloc[0].HELP_sequence)
        a["combined_positive_rank_1"]=a.HELP_sequence.eq(pos.iloc[0].HELP_sequence)
    else:
        for name in ("positive_joint_rank","Pearson_positive_rank_1",
                     "RMSE_positive_rank_1","combined_positive_rank_1"):
            a[name]=np.nan if "joint" in name else False
        by_r=by_rmse=pos
    if not neg.empty:
        neg["rank_absr"]=neg.Pearson_abs.rank(ascending=False,method="min")
        neg["rank_inverted_rmse"]=neg.inverted_Heavy_RMSE_log2_diagnostic.rank(
                             ascending=True,method="min")
        neg["rank_sum"]=neg.rank_absr+neg.rank_inverted_rmse
        neg=neg.sort_values(["rank_sum","inverted_Heavy_RMSE_log2_diagnostic",
                             "Pearson_abs"],ascending=[True,True,False],kind="stable")
        by_abs=neg.sort_values("Pearson_abs",ascending=False,kind="stable")
        by_fliprmse=neg.sort_values("inverted_Heavy_RMSE_log2_diagnostic",
                                    ascending=True,kind="stable")
    else:by_abs=by_fliprmse=neg
    return a,pos,by_r,by_rmse,neg,by_abs,by_fliprmse


def one_cohort(cohort,recs,names,raw,seqs,allH):
    assert len(recs)==len(names)==raw.shape[1]==allH.shape[1]
    eligible=np.where(np.isfinite(allH).all(axis=1)&(allH>0).all(axis=1))[0]
    target=np.nanmedian(raw,axis=0)
    assert np.isfinite(target).all()
    d=OUT/cohort;d.mkdir(parents=True,exist_ok=False)
    stats=[];details=[]
    for k in eligible:
        seq=seqs[k]
        h=np.log2(allH[k])
        m=metrics(h,target)
        stats.append({"HELP_sequence":seq,"n_samples":len(names),
                      "full_detected":True,**m})
        for i,rec in enumerate(recs):
            details.append({"HELP_sequence":seq,"sample_id":names[i],
                "batch":rec["batch"],"sample_number":i+1,
                "original_Heavy_log2":h[i],
                "original_samplePGmedian_log2":target[i],
                "Heavy_centered_log2":h[i]-np.median(h),
                "PGmedian_centered_log2":target[i]-np.median(target)})
    df,pos,byr,byrm,neg,byabs,byflip=ranks(pd.DataFrame(stats))
    df.sort_values("median_centered_RMSE_log2").to_csv(
        d/"01_all_full_detection_HELP_Pearson_centeredRMSE.csv",index=False)
    pd.DataFrame(details).to_csv(d/"02_all_HELP_curves_by_sample.csv.gz",
                        index=False,compression={"method":"gzip","compresslevel":4})
    pos.to_csv(d/"03_positive_same_direction_combined_rank.csv",index=False)
    byr.to_csv(d/"04_positive_rank_by_Pearson.csv",index=False)
    byrm.to_csv(d/"05_positive_rank_by_centered_RMSE.csv",index=False)
    neg.to_csv(d/"06_negative_opposite_direction_inverted_diagnostic_rank.csv",index=False)
    byabs.to_csv(d/"07_negative_rank_by_absPearson.csv",index=False)
    byflip.to_csv(d/"08_negative_rank_by_inverted_RMSE.csv",index=False)
    figs=[plots_scatter(df,cohort)]
    p=plot_top_curves(np.log2(allH),target,recs,names,seqs,pos,cohort)
    if p:figs.append(p)
    q=plot_top_curves(np.log2(allH),target,recs,names,seqs,neg,cohort,reversed_mode=True)
    if q:figs.append(q)
    return {"n_samples":len(names),"full_detected":int(len(eligible)),
            "n_positive":int(len(pos)),"n_negative":int(len(neg)),
            "n_zero_or_undefined":int(len(df)-len(pos)-len(neg)),
            "best_joint_positive":pos.head(5).to_dict("records"),
            "best_pearson_positive":byr.head(5).to_dict("records"),
            "best_rmse_positive":byrm.head(5).to_dict("records"),
            "best_joint_negative_if_inverted":neg.head(5).to_dict("records"),
            "plots":figs,
            "median_samplePG_reference":"Original PG matrix, own detected PG in each sample",
            "not_a_correction":True}


def main():
    if OUT.exists():raise FileExistsError("Protect prior results: "+str(OUT))
    recs,names,raw,annotations,seqs,H=data.load()
    assert len(names)==40 and len(seqs)==97 and raw.shape==(4264,40)
    OUT.mkdir(parents=True)
    full=one_cohort("full40",recs,names,raw,seqs,H)
    keep=np.array([sid!=data.EXCLUDE for sid in names],dtype=bool)
    assert int(keep.sum())==39
    sens=one_cohort("exclude_T003_bad_sample_39",
        [r for r,k in zip(recs,keep) if k],
        [s for s,k in zip(names,keep) if k],
        raw[:,keep],seqs,H[:,keep])
    manifest={"method":"Raw each-sample PG median vs individual original Heavy log2 peptide",
        "target":"Per-sample median of all positive detected PG log2 values",
        "pearson":"One correlation across sample dimension for each HELP",
        "rmse":"sqrt(mean((HeavyLog2 - median(HeavyLog2) - PGsampleMedianLog2 + median(PGsampleMedianLog2))**2))",
        "main40":full,"sensitivity39":sens,
        "joint_rank":"Among Pearson-positive HELP: rank descending Pearson + rank ascending centered RMSE; tied by RMSE",
        "negative_analysis":"Separate; flipped Heavy centered trajectory ONLY for comparison, never counted same-direction",
        "no_normalization":True,
        "outlier":"Full 40 includes historical T003 outlier; sensitivity recomputed on remaining 39",
        "original_source_files_untouched":True}
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n")
    doc=["# NYD: Pearson and centered RMSE screen (Heavy vs PG sample medians)",
         "","One Heavy HELP (log2) across N samples vs N per-sample PG log2 medians.",
         "Excludes missing PG intensity values from each median; no calibration.",
         "Pearson computed on raw log2 values; centered RMSE has no free fit or rescaling.",
         "Pearson-positive HELP are the only ones that naturally share the PG trend direction.",
         "Negative Pearson HELP are presented separately. A sign-flipped visualization",
         "is not evidence that their original curves match.",""]
    for n,res in (("full40",full),("exclude_T003_bad_sample_39",sens)):
        doc.extend(["## "+n,"",
         f"- HELP completely detected: {res['full_detected']}/97.",
         f"- Same direction r>0: {res['n_positive']}; opposite r<0: {res['n_negative']}.",
         "","| Combined positive rank | HELP | Pearson r | Centered RMSE log2 | Heavy/PG SD |",
         "|---:|---|---:|---:|---:|"])
        for i,row in enumerate(res["best_joint_positive"],1):
            doc.append(f"| {i} | {row['HELP_sequence']} | {row['Pearson_r']:+.4f} | "
                       f"{row['median_centered_RMSE_log2']:.4f} | {row['SD_Heavy_over_PGmedian']:.3f} |")
        doc += ["",f"![Pearson and RMSE overview]({n}/01_scatter_all_HELP_Pearson_RMSE.png)",
                f"![Same-direction top]({n}/02_top5_positive_same_direction.png)",
                f"![Opposite-direction inverted diagnostic]({n}/02_top5_opposite_inverted_diagnostic.png)",""]
    doc+=["## Tables",
      "- 01_all_full_detection_HELP_Pearson_centeredRMSE.csv: every fully detected HELP",
      "- 03_positive_same_direction_combined_rank.csv",
      "- 04_positive_rank_by_Pearson.csv",
      "- 05_positive_rank_by_centered_RMSE.csv",
      "- 06_negative_opposite_direction_inverted_diagnostic_rank.csv",
      "- 02_all_HELP_curves_by_sample.csv.gz for reconstructing charts",
      "","## Cautions",
      "- Full40 outlier may dominate correlation and RMSE; compare independent 39 set.",
      "- Equal centered-RMSE is sensitive to Heavy/PG amplitudes; Pearson only shape.",
      "- Good fit to clinical PG medians does not demonstrate better technical calibration."]
    (OUT/"README.md").write_text("\n".join(doc)+"\n",encoding="utf-8")
    print("FINAL",json.dumps({
       "full40":{key:full[key] for key in ("full_detected","n_positive","n_negative","best_joint_positive","best_pearson_positive","best_rmse_positive")},
       "exclude39":{key:sens[key] for key in ("full_detected","n_positive","n_negative","best_joint_positive","best_pearson_positive","best_rmse_positive")}},
       ensure_ascii=False),flush=True)


if __name__=="__main__":main()
