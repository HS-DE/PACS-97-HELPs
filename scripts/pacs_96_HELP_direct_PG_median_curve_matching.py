#!/usr/bin/env python3
"""Directly compare every fully detected PACS HELP Heavy curve with the 75-sample
raw PG per-sample median curve. Distinguish level Pearson, amplitude-matched
median-centered deviations, and successive plotted-sample changes.

No normalization. No target-selection by pooled QC. This precisely addresses
the requested all-75 median-curve trend matching, with QC as extra diagnostics.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pacs_heavy_97_trend_screening as core

OUT=core.old.ROOT/"results"/"pacs_96_HELP_direct_raw_PG_median_curve_matching"
PLATES=("1","2","3")
PLOT_GROUPS=core.GROUPS
COLORS={"PG":"#283E54","Heavy":"#C47744"}


def stats(x,y):
    if not(np.isfinite(x).all() and np.isfinite(y).all()):raise ValueError("Nonfinite full-75 curve")
    xp=x-np.median(x); yp=y-np.median(y)
    diff=xp-yp
    r=float(pearsonr(x,y).statistic)
    spearman=float(spearmanr(x,y).statistic)
    cov=float(np.dot(x-x.mean(),y-y.mean()))
    beta=cov/float(np.dot(x-x.mean(),x-x.mean()))
    sd_x=float(np.std(x,ddof=1))
    sd_y=float(np.std(y,ddof=1))
    return {
        "pearson_r":r,"spearman_rho":spearman,
        "median_centered_RMSE_log2":float(np.sqrt(np.mean(diff**2))),
        "median_centered_MAE_log2":float(np.mean(abs(diff))),
        "HELP_SD_log2":sd_x,
        "PG_SD_log2":sd_y,
        "HELP_SD_over_PG_SD":sd_x/sd_y,
        "slope_PG_per_HELP":beta,
        "mean_absolute_PG_deviation":float(np.mean(abs(yp))),
        "mean_absolute_HELP_deviation":float(np.mean(abs(xp))),
    }


def adjacent_stats(x,y,order,group,plate):
    """Only neighboring points within the SAME sample type and plate.
    Sample ID sorting is not necessarily chronological: this is graphical
    curve-shape diagnostics, not validated acquisition-time drift."""
    xx,yy=[],[]
    for k in range(len(order)-1):
        i,j=int(order[k]),int(order[k+1])
        if group[i]==group[j] and plate[i]==plate[j]:
            xx.append(float(x[j]-x[i]))
            yy.append(float(y[j]-y[i]))
    xx=np.array(xx);yy=np.array(yy)
    if len(xx)<3:raise ValueError("Too few comparable within-type-plate neighbor steps")
    return {"adjacent_pairs_same_group_plate":len(xx),
            "neighbor_step_pearson_r":float(pearsonr(xx,yy).statistic),
            "neighbor_step_direction_agreement_pct":float(np.mean(np.sign(xx)==np.sign(yy))*100),
            "neighbor_step_RMSE_log2":float(np.sqrt(np.mean((xx-yy)**2)))}


def annotate_axes(ax,order,group,plates):
    for k in range(1,len(order)):
        i,j=int(order[k-1]),int(order[k])
        if group[i]!=group[j]:ax.axvline(k+.5,c="#788B9E",lw=1.15)
        elif plates[i]!=plates[j]:ax.axvline(k+.5,c="#BBC5CE",lw=.8,ls=":")
    ax.grid(axis="y",alpha=.15)
    ax.axhline(0,c="#B1BEC9",lw=.7)
    ax.set_xlim(.5,75.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def overlay(seqs,H,PG,order,groups,plates,ranking,title,filename):
    top=list(ranking.HELP_sequence.head(5))
    fig,axes=plt.subplots(5,1,figsize=(22,17.5),sharex=True)
    fig.subplots_adjust(left=.06,right=.98,top=.92,bottom=.075,hspace=.36)
    x=np.arange(1,76)
    pgc=PG-np.median(PG)
    for ax,seq in zip(axes,top):
        h=H[seqs.index(seq)]
        hcent=h-np.median(h)
        metrics=ranking.set_index("HELP_sequence").loc[seq]
        ax.plot(x,pgc[order],c=COLORS["PG"],lw=1.55,marker=".",ms=5,
                label="Raw PG sample median deviation")
        ax.plot(x,hcent[order],c=COLORS["Heavy"],lw=1.35,marker=".",ms=4,
                label=seq+" Heavy deviation")
        annotate_axes(ax,order,groups,plates)
        ax.set_ylabel("Log2 median deviation")
        ax.set_title(f"{seq}: r = {metrics['ownPG_Pearson_r']:+.3f}; "
                     f"RMSE = {metrics['ownPG_centered_RMSE_log2']:.3f}; "
                     f"SD(HELP)/SD(PG) = {metrics['ownPG_SD_ratio']:.2f}",
                     loc="left",fontsize=11,pad=10,fontweight="bold")
    axes[-1].set_xticks(x,[str(i) for i in x],rotation=90,fontsize=7)
    axes[-1].set_xlabel("Sample display index, group -> plate -> sample ID")
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,["Raw PG log2 median deviation","Candidate HELP Heavy log2 deviation"],
               ncol=2,loc="lower center",bbox_to_anchor=(.5,.012),
               frameon=False,fontsize=10)
    fig.suptitle(title,fontsize=15,fontweight="bold")
    fig.savefig(OUT/filename,dpi=175,bbox_inches="tight",facecolor="white")
    plt.close(fig)


def diag_scatter(ranked):
    fig,ax=plt.subplots(figsize=(12,7.8))
    fig.subplots_adjust(left=.09,right=.97,top=.89,bottom=.20)
    xx=ranked.ownPG_Pearson_r
    yy=ranked.ownPG_centered_RMSE_log2
    sizes=np.clip(ranked.ownPG_SD_ratio,0,8)
    points=ax.scatter(xx,yy,c=sizes,s=55,cmap="viridis",alpha=.75,edgecolors="white",linewidths=.6)
    highlight=pd.concat([ranked.nlargest(4,"ownPG_Pearson_r"),
                         ranked.nsmallest(4,"ownPG_centered_RMSE_log2")])
    for item in highlight.drop_duplicates("HELP_sequence").itertuples():
        ax.annotate(item.HELP_sequence,(item.ownPG_Pearson_r,item.ownPG_centered_RMSE_log2),
                    xytext=(5,7),textcoords="offset points",fontsize=9)
    ax.axvline(.9,c="#9B3346",lw=1,ls="--",alpha=.8,label="Pearson r = 0.90")
    ax.set_xlabel("Pearson r: original Heavy vs per-sample Raw PG median (75 samples)")
    ax.set_ylabel("RMSE after median-centering each log2 curve (lower = closer)")
    ax.set_title("No standardized curve fitting: matching direction is not equal to matching amplitude",
                 fontsize=12,fontweight="bold")
    fig.colorbar(points,ax=ax,label="Heavy/PG variation SD ratio",fraction=.03,pad=.02)
    h,l=ax.get_legend_handles_labels()
    fig.legend(h,l,loc="lower center",bbox_to_anchor=(.5,.08),frameon=False)
    ax.grid(alpha=.12)
    fig.savefig(OUT/"03_all96_Pearson_vs_centered_RMSE_scatter.png",dpi=180,bbox_inches="tight")
    plt.close(fig)


def group_overlay(seqs,H,PG,order,groups,plates,sequences):
    x=np.arange(1,76)
    fig,axes=plt.subplots(2,1,figsize=(22,9.8),sharex=True)
    fig.subplots_adjust(left=.06,right=.985,bottom=.13,top=.89,hspace=.32)
    for ax,seq in zip(axes,sequences):
        h=H[seqs.index(seq)]
        xh=h-np.median(h)
        xp=PG-np.median(PG)
        ax.plot(x,xp[order],c=COLORS["PG"],marker=".",lw=1.6,label="Raw PG median")
        ax.plot(x,xh[order],c=COLORS["Heavy"],marker=".",lw=1.4,label=seq+" original Heavy")
        annotate_axes(ax,order,groups,plates)
        ax.set_title(f"Previously highlighted reference: {seq}",loc="left",fontsize=12)
        ax.set_ylabel("Median-centered log2 deviation")
    axes[-1].set_xticks(x,[str(i) for i in x],rotation=90,fontsize=6)
    h,l=axes[0].get_legend_handles_labels()
    fig.legend(h,l,loc="lower center",bbox_to_anchor=(.5,.02),ncol=2,frameon=False)
    fig.suptitle("Previously highlighted Heavy candidates vs exact Raw PG median trajectory",
                 fontsize=15,fontweight="bold")
    fig.savefig(OUT/"04_LAVYQAGAR_AVLTIDEK_vs_Raw_PG_median.png",dpi=185,bbox_inches="tight")
    plt.close(fig)


def main():
    if OUT.exists():raise FileExistsError(f"Refuse to overwrite existing analysis {OUT}")
    pg,names,info,logPG,seqs,H,complete,groups,plates=core.raw_inputs()
    medians,masks=core.median_reference_pg(logPG,groups,plates)
    order=core.sample_order(names,groups,plates)
    if sum(complete)!=96:raise AssertionError("Expected 96 of 97 fully detected Heavy")
    y=medians["all_detected_PG"]
    yf=medians["common_75_PG"]
    rows=[]
    for k,seq in enumerate(seqs):
        if not complete[k]:continue
        x=H[k]
        a=stats(x,y)
        f=stats(x,yf)
        d=adjacent_stats(x,y,order,groups,plates)
        row={"HELP_sequence":seq,"n_samples":75,
             "ownPG_Pearson_r":a["pearson_r"],
             "ownPG_Spearman_rho":a["spearman_rho"],
             "ownPG_centered_RMSE_log2":a["median_centered_RMSE_log2"],
             "ownPG_centered_MAE_log2":a["median_centered_MAE_log2"],
             "ownPG_SD_ratio":a["HELP_SD_over_PG_SD"],
             "ownPG_slope_PG_per_HELP":a["slope_PG_per_HELP"],
             "fixed75PG_Pearson_r":f["pearson_r"],
             "fixed75PG_centered_RMSE_log2":f["median_centered_RMSE_log2"],
             "fixed75PG_SD_ratio":f["HELP_SD_over_PG_SD"],
             **d}
        qcidx=(groups=="pooled QC")
        qcH=core.qc_plate_residual(x,plates,qcidx)
        qcPG=core.qc_plate_residual(medians["QC_common_15_PG"],plates,qcidx)
        qcR=stats(qcH,qcPG)
        row["within_plate_QC15_Pearson_r"]=qcR["pearson_r"]
        row["within_plate_QC15_centered_RMSE_log2"]=qcR["median_centered_RMSE_log2"]
        rows.append(row)
    rank=pd.DataFrame(rows)
    rank["rank_by_overall_Pearson"]=rank.ownPG_Pearson_r.rank(method="first",ascending=False).astype(int)
    rank["rank_by_centered_RMSE"]=rank.ownPG_centered_RMSE_log2.rank(method="first").astype(int)
    rank["rank_by_neighbor_changes"]=rank.neighbor_step_pearson_r.rank(method="first",ascending=False).astype(int)
    by_r=rank.sort_values("ownPG_Pearson_r",ascending=False)
    by_rmse=rank.sort_values("ownPG_centered_RMSE_log2")
    by_step=rank.sort_values("neighbor_step_pearson_r",ascending=False)
    if max(rank.ownPG_Pearson_r)>.9:
        print("NOTE strong matching sequence detected, check sample groups")
    if not np.isclose(by_r.iloc[0].ownPG_Pearson_r, .7178607757641912, atol=1e-10):
        raise AssertionError("Data differs from prior PACS Pearson run")
    OUT.mkdir(parents=True)
    by_r.to_csv(OUT/"00_all96_ranking_by_overall_Pearson.csv",index=False)
    by_rmse.to_csv(OUT/"01_all96_ranking_by_centered_RMSE.csv",index=False)
    by_step.to_csv(OUT/"02_all96_ranking_by_neighbor_curve_changes.csv",index=False)
    sample_details=[]
    for k,j in enumerate(order,1):
        item={"plot_index":k,"sample_id":names[j],
              "group":groups[j],"plate":plates[j],
              "pg_log2_median_all_detected":float(y[j]),
              "pg_log2_median_common_75":float(yf[j]),
              "PG_all_detected_centered":float(y[j]-np.median(y))}
        for seq in sorted(set(by_r.HELP_sequence.head(5))|set(by_rmse.HELP_sequence.head(5))|
                          {"LAVYQAGAR","AVLTIDEK"}):
            Hval=float(H[seqs.index(seq),j])
            item[seq+"_Heavy_log2"]=Hval
            item[seq+"_centered_log2"]=Hval-float(np.median(H[seqs.index(seq)]))
        sample_details.append(item)
    pd.DataFrame(sample_details).to_csv(OUT/"03_sample_curve_values.csv",index=False)
    overlay(seqs,H,y,order,groups,plates,by_r,
        "Top 5 individual Heavy by Pearson | original 75-sample PG medians",
        "01_top5_by_Pearson_75sample_median_curve.png")
    overlay(seqs,H,y,order,groups,plates,by_rmse,
        "Top 5 individual Heavy by unscaled centered-log2 RMSE | original PG medians",
        "02_top5_by_real_deviation_RMSE_75sample_median_curve.png")
    diag_scatter(rank)
    group_overlay(seqs,H,y,order,groups,plates,("LAVYQAGAR","AVLTIDEK"))
    statlist={
        "full_75_HELP":96,"reference":"median(log2 positive PG intensity) per sample",
        "reference_median_across_75":float(np.median(y)),
        "reference_PG_sd":float(np.std(y,ddof=1)),
        "no_candidate_with_r_ge_09":int(np.sum(rank.ownPG_Pearson_r>=.9))==0,
        "count_r_ge_09":int(np.sum(rank.ownPG_Pearson_r>=.9)),
        "count_r_ge_08":int(np.sum(rank.ownPG_Pearson_r>=.8)),
        "best_by_r":by_r.HELP_sequence.head(5).tolist(),
        "best_by_rmse":by_rmse.HELP_sequence.head(5).tolist(),
        "best_by_neighbor_steps":by_step.HELP_sequence.head(5).tolist(),
        "n_fixed_75PG":int(masks["common_75"].sum()),
        "analysis_nature":"No calibration. Raw original values; comparing 75-sample curves",
        "caveat":"Adjacent differences depend on chosen group/plate/sample ID ordering; not chronological or causal",
        "note":"Pearson is invariant to offsets and positive scaling, while centered RMSE retains amplitude mismatch",
    }
    (OUT/"summary.json").write_text(json.dumps(statlist,ensure_ascii=False,indent=2)+"\n")
    lines=["# Which Heavy HELP really follows the PG median across all 75 samples?","",
           "One fixed approach: m_i = median(log2(valid positive original PG intensities in sample i)).",
           "For each complete peptide j, H_ji = log2(original Heavy intensity).",
           "To compare ACTUAL amplitude in plot units, subtract median_i(m_i) and",
           "median_i(H_ji) from the two 75-point curves, then compute RMSE.",
           "Pearson only checks linear covariation; it does NOT check amplitude.",
           "",
           "## Screening results",
           "",
           f"- Full-75 Heavy HELP sequences: {len(rank)}.",
           f"- Best Pearson: {by_r.iloc[0].HELP_sequence}, r={by_r.iloc[0].ownPG_Pearson_r:.4f}.",
           f"- No. with r >= 0.90: {statlist['count_r_ge_09']}.",
           f"- No. with r >= 0.80: {statlist['count_r_ge_08']}.",
           f"- Raw sample PG median SD on log2: {statlist['reference_PG_sd']:.4f}.",
           "",
           "### Top five by Pearson (shape, scale ignored)",
           "",
           "| HELP | r | median-centered RMSE log2 | SD(Heavy)/SD(PG) | Neighbor step r |",
           "|---|---:|---:|---:|---:|",
           *[f"| {r.HELP_sequence} | {r.ownPG_Pearson_r:+.3f} | "
             f"{r.ownPG_centered_RMSE_log2:.3f} | {r.ownPG_SD_ratio:.2f} | "
             f"{r.neighbor_step_pearson_r:+.3f} |" for r in by_r.head(5).itertuples()],
           "",
           "### Top five by median-centered RMSE (actual shape AND amplitude)",
           "",
           "| HELP | RMSE log2 | r | SD(Heavy)/SD(PG) | Neighbor step r |",
           "|---|---:|---:|---:|---:|",
           *[f"| {r.HELP_sequence} | {r.ownPG_centered_RMSE_log2:.3f} | "
             f"{r.ownPG_Pearson_r:+.3f} | {r.ownPG_SD_ratio:.2f} | "
             f"{r.neighbor_step_pearson_r:+.3f} |" for r in by_rmse.head(5).itertuples()],
           "",
           "## Charts",
           "",
           "![Best Pearson](01_top5_by_Pearson_75sample_median_curve.png)",
           "",
           "![Best actual RMSE](02_top5_by_real_deviation_RMSE_75sample_median_curve.png)",
           "",
           "![Scatter](03_all96_Pearson_vs_centered_RMSE_scatter.png)",
           "",
           "![Earlier two candidates](04_LAVYQAGAR_AVLTIDEK_vs_Raw_PG_median.png)",
           "",
           "## Caveats",
           "",
           "- Samples are heterogeneous HC/S/QC/Internal-QC, not the same biological sample.",
           "- Median of own detected PG is associated with missingness/composition.",
           "- Pearson and RMSE computed across all 75 with no per-plate centering.",
           "- Neighbor-step calculations only cross adjacent plotted samples within the",
           "  same group and plate, and are NOT a validated acquisition-order trend.",
           "- No raw PG normalization or HELP calibration is performed here.",
           "- Selecting on this dataset then claiming calibrated improvement on the",
           "  same dataset would be circular without a holdout or independent QC test.",
           "",
           "## Files",
           "",
           "- 00_all96_ranking_by_overall_Pearson.csv",
           "- 01_all96_ranking_by_centered_RMSE.csv",
           "- 02_all96_ranking_by_neighbor_curve_changes.csv",
           "- 03_sample_curve_values.csv",
           "- summary.json"
    ]
    (OUT/"README.md").write_text("\n".join(lines)+"\n")
    print("RUN_SUMMARY",json.dumps(statlist,ensure_ascii=False))
    print("TOP_PEARSON",json.dumps(by_r.head(5)[[
       "HELP_sequence","ownPG_Pearson_r","ownPG_centered_RMSE_log2",
       "ownPG_SD_ratio","neighbor_step_pearson_r"]].to_dict("records"),ensure_ascii=False))
    print("TOP_RMSE",json.dumps(by_rmse.head(5)[[
       "HELP_sequence","ownPG_Pearson_r","ownPG_centered_RMSE_log2",
       "ownPG_SD_ratio","neighbor_step_pearson_r"]].to_dict("records"),ensure_ascii=False))
    print("PASS 96 complete HELP x 75 raw PG sample medians, two independent ranking metrics")


if __name__=="__main__":
    main()
