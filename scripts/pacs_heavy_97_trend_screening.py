#!/usr/bin/env python3
"""Discover 75/75-detected HELP Heavy peptides tracking raw PG median fluctuations.

PACS main 75 heterogeneous samples, raw pg_matrix.tsv and original H sheet.
For analysis only: candidates ranked by QC per-plate-demeaned PG/Heavy
fluctuations in 15 pooled QC, with descriptive full75 and HC/S diagnostics.
No correction is applied; no existing files are overwritten.
"""
from __future__ import annotations
import hashlib, json, math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from filter_samples import readbook, tab
from yeast_protein_qc import GROUPS, COLORS
import pacs_two_reference_plate_trial as old

OUT=old.ROOT/"results"/"pacs_help97_heavy_trend_vs_raw_PG_75samples"
PLATES=("1","2","3")
GROUP_COLOR=dict(zip(GROUPS,COLORS))
TOP=10


def raw_inputs():
    pg,names,info,logPG,_,_=old.load()
    shared,sheets=readbook(old.HELP)
    if len(sheets)!=3:raise ValueError("HELP workbook expected L/H/ratio_LH")
    header=tab(sheets[1],shared)
    if len(header)!=98:raise ValueError("Need exactly 97 Heavy HELP sequences")
    cols={str(v).strip():col for col,v in header[0].items() if col!="A"}
    if len(cols)!=75 or set(cols)!=set(names):raise ValueError("Heavy sample IDs mismatch PG IDs")
    sequences=[str(row.get("A","")).strip() for row in header[1:]]
    if len(sequences)!=97 or len(set(sequences))!=97 or not all(sequences):
        raise ValueError("HELP peptide rows not unique")
    H=np.full((97,75),np.nan)
    for k,row in enumerate(header[1:]):
        for j,sid in enumerate(names):
            try:
                val=float(row.get(cols[sid],"nan"))
                if np.isfinite(val) and val>0:H[k,j]=np.log2(val)
            except (TypeError,ValueError):
                pass
    full=np.isfinite(H).all(axis=1)
    if not np.any(full):raise ValueError("No Heavy HELP detected 75/75")
    groups=np.array([info[s]["group"] for s in names],dtype=object)
    plates=np.array([info[s]["plate"] for s in names],dtype=object)
    return pg,names,info,logPG,sequences,H,full,groups,plates


def median_reference_pg(raw,group,plate):
    """Return raw PG median for own detected set and fixed protein intersections.

    Fixed_75 PG determined before looking at Heavy. Fixed_15QC PG is valid in
    all 15 pooled QCs; applied to all 75 where available, with a QC-complete
    vector for QC correlation. No imputation.
    """
    valid=np.isfinite(raw)
    qc=group=="pooled QC"
    study=np.isin(group,("HC","S"))
    full75=valid.all(axis=1)
    fullQC=valid[:,qc].all(axis=1)
    fullStudy=valid[:,study].all(axis=1)
    if full75.sum()<50 or fullQC.sum()<500 or fullStudy.sum()<50:
        raise ValueError(f"Fixed PG set sizes unexpected: 75={full75.sum()} 15QC={fullQC.sum()} 54Study={fullStudy.sum()}")
    def get(mask):
        sub=raw[mask]
        with np.errstate(all="ignore"):
            return np.nanmedian(sub,axis=0)
    return dict(
        all_detected_PG=np.nanmedian(raw,axis=0),
        common_75_PG=np.median(raw[full75],axis=0),
        QC_common_15_PG=get(fullQC),
        Study_common_54_PG=get(fullStudy)
    ),dict(common_75=full75,QC_common_15=fullQC,Study_common_54=fullStudy)


def metrics(x,y,scope,seq,basis):
    keep=np.isfinite(x)&np.isfinite(y)
    a,b=x[keep],y[keep]
    row={"HELP_sequence":seq,"scope":scope,"PG_basis":basis,"n_pairs":int(len(a))}
    if len(a)<4 or min(np.std(a),np.std(b))<1e-12:
        return {**row,"pearson_r":np.nan,"pearson_p":np.nan,
                "spearman_rho":np.nan,"slope_PG_per_HELP":np.nan,
                "HELP_SD_log2":np.nan,"PG_SD_log2":np.nan,
                "deviation_rmse_log2":np.nan,
                "deviation_mae_log2":np.nan}
    ra=pearsonr(a,b)
    rank=spearmanr(a,b)
    # Median-center both signals (peptide/PG differ in absolute baseline).
    delta_a=a-np.median(a)
    delta_b=b-np.median(b)
    residual=delta_a-delta_b
    slope=float(np.sum((a-a.mean())*(b-b.mean()))/np.sum((a-a.mean())**2))
    return {**row,
        "pearson_r":float(ra.statistic),"pearson_p":float(ra.pvalue),
        "spearman_rho":float(rank.statistic),
        "slope_PG_per_HELP":slope,
        "HELP_SD_log2":float(a.std(ddof=1)),
        "PG_SD_log2":float(b.std(ddof=1)),
        "deviation_rmse_log2":float(np.sqrt(np.mean(residual**2))),
        "deviation_mae_log2":float(np.mean(np.abs(residual)))}


def qc_plate_residual(x,plates,qc):
    """Remove each plate's median in QC, so common plate offsets cannot inflate r."""
    z=x[qc].copy()
    p=plates[qc]
    for plate in PLATES:
        ix=(p==plate)
        if np.isfinite(z[ix]).sum()<3:raise ValueError(f"Insufficient QC {plate}")
        z[ix]-=np.nanmedian(z[ix])
    return z


def save(fig,path,dpi=185):
    fig.savefig(path,dpi=dpi,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    if path.stat().st_size<24000:raise AssertionError(f"Suspicious image {path.name}")


def axis_groups(ax,ids,groups,plates):
    for boundary in range(1,len(ids)):
        u,v=ids[boundary-1],ids[boundary]
        if groups[u]!=groups[v]:
            ax.axvline(boundary+.5,color="#41566A",linewidth=1.1,alpha=.65)
        elif plates[u]!=plates[v]:
            ax.axvline(boundary+.5,color="#8795A2",linewidth=.7,linestyle=":",alpha=.7)
    ax.grid(axis="y",alpha=.13)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def sample_order(names,groups,plates):
    return np.array(sorted(range(len(names)),key=lambda j:(
        GROUPS.index(groups[j]),int(plates[j]),names[j])),int)


def boxed(mat,order,groups,plate,ylim=None):
    fig,ax=plt.subplots(figsize=(24,8.2))
    fig.subplots_adjust(left=.055,right=.985,bottom=.22,top=.86)
    vecs=[mat[:,j][np.isfinite(mat[:,j])] for j in order]
    bp=ax.boxplot(vecs,positions=np.arange(1,len(order)+1),widths=.64,
                  patch_artist=True,showfliers=False,manage_ticks=False,
                  boxprops={"edgecolor":"#617788","linewidth":.7},
                  medianprops={"color":"#18394A","linewidth":1.7},
                  whiskerprops={"color":"#718494","linewidth":.8},
                  capprops={"color":"#718494","linewidth":.8})
    for patch,j in zip(bp["boxes"],order):
        patch.set_facecolor(GROUP_COLOR[groups[j]]);patch.set_alpha(.52)
    median=np.array([float(np.median(x)) for x in vecs])
    ax.scatter(np.arange(1,len(order)+1),median,color="#253A4C",s=10,zorder=4)
    if ylim:ax.set_ylim(*ylim)
    ax.set_xlim(.4,len(order)+.6)
    axis_groups(ax,order,groups,plate)
    ax.set_xticks(range(1,len(order)+1),[str(i) for i in range(1,len(order)+1)],
                  rotation=90,fontsize=6.8)
    handles=[Patch(facecolor=GROUP_COLOR[g],alpha=.6,label=g) for g in GROUPS]
    fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(.5,.025),
               ncol=4,frameon=False,fontsize=10)
    return fig,ax,median


def box_figures(raw,Hfull,names,groups,plates,order,pg_medians):
    fig,ax,med=boxed(raw,order,groups,plates)
    grand=float(np.median(pg_medians["all_detected_PG"]))
    ax.axhline(grand,c="#CB6D3D",ls="--",lw=1.6)
    ax.set_title("Raw PG: per-sample log2 protein-intensity distributions (all detected PG)",
                 pad=13,fontsize=15,fontweight="bold")
    ax.set_ylabel("Protein-group log2 intensity")
    ax.text(.01,.99,f"Overall median of the 75 per-sample PG medians = {grand:.4f}",
            transform=ax.transAxes,fontsize=11,va="top")
    save(fig,OUT/"01_raw_PG_75sample_boxplots_with_median.png")

    fig,ax,med=boxed(Hfull,order,groups,plates)
    overall=float(np.median(np.median(Hfull,axis=0)))
    ax.axhline(overall,c="#CB6D3D",ls="--",lw=1.6)
    ax.set_title(f"HELP Heavy: all {len(Hfull)} sequences detected in 100% of 75 samples",
                 fontsize=15,fontweight="bold",pad=13)
    ax.set_ylabel("Original HELP Heavy log2 intensity")
    ax.text(.01,.99,f"Median of 75 per-sample Heavy medians = {overall:.4f}",
            transform=ax.transAxes,fontsize=11,va="top")
    save(fig,OUT/"02_full_detection_Heavy_HELP_75sample_boxplots.png")


def median_lines(pg,Hfull,names,groups,plates,order):
    # Separately center the two samplewise medians before overlay, not unit variance.
    allpg=pg["all_detected_PG"]
    fixedpg=pg["common_75_PG"]
    helpmed=np.median(Hfull,axis=0)
    series=[
        ("PG median, own detected PG",allpg,"#263F54","-"),
        ("PG median, fixed 75-common PG",fixedpg,"#449379","-"),
        ("Heavy HELP box median (100% detected)",helpmed,"#CC7945","--")
    ]
    fig,ax=plt.subplots(figsize=(24,6.8))
    fig.subplots_adjust(left=.06,right=.985,bottom=.25,top=.87)
    x=np.arange(1,76)
    for label,values,col,ls in series:
        centered=values-np.median(values)
        ax.plot(x,centered[order],label=label,c=col,lw=1.5,marker=".",
                markersize=4,linestyle=ls,alpha=.92)
    axis_groups(ax,order,groups,plates)
    ax.axhline(0,c="#A1AAB3",lw=.85)
    ax.set_xticks(x,[str(i) for i in range(1,76)],fontsize=6.8,rotation=90)
    ax.set_xlabel("Sample index; see sample_order_and_medians.csv")
    ax.set_ylabel("Log2 median deviation from cross-sample median")
    ax.set_title("Median fluctuation shapes on one log2 deviation scale (no rescaling by SD)",
                 fontsize=14,fontweight="bold")
    handles,labels=ax.get_legend_handles_labels()
    fig.legend(handles,labels,loc="lower center",bbox_to_anchor=(.5,.025),
               frameon=False,ncol=3,fontsize=10)
    save(fig,OUT/"03_centered_median_trends_PG_vs_all_full_Heavy_HELP.png")


def trend_panels(pg,H,seqs,topseq,names,groups,plates,order):
    # Display best QC-within-plate candidates but compare to BOTH PG definitions.
    fig,axs=plt.subplots(5,2,figsize=(23,19.0),sharex=True)
    fig.subplots_adjust(left=.065,right=.985,bottom=.10,top=.92,hspace=.46,wspace=.2)
    series=[
       ("Sample own detected PG median", pg["all_detected_PG"],"#293D50"),
       ("Fixed 75-common PG median",pg["common_75_PG"],"#4A967F")
    ]
    x=np.arange(1,len(order)+1)
    for idx,seq in enumerate(topseq):
        heavy=H[seqs.index(seq)]
        target=heavy-np.median(heavy)
        for k,(label,pgy,col) in enumerate(series):
            ax=axs[idx,k]
            pgdelta=pgy-np.median(pgy)
            ax.plot(x,pgdelta[order],c=col,lw=1.45,marker=".",ms=3.1,
                    label="PG median deviation")
            ax.plot(x,target[order],c="#CA7946",lw=1.5,alpha=.85,
                    label=seq+" Heavy deviation")
            ax.axhline(0,c="#ABB7C1",lw=.7)
            axis_groups(ax,order,groups,plates)
            ax.set_title(f"{seq} against {label}",fontsize=10.8,pad=7)
            ax.set_ylabel("Centered log2 deviation",fontsize=9)
    for ax in axs[-1]:
        ax.set_xticks(x,[str(i) for i in x],rotation=90,fontsize=5.6)
        ax.set_xlabel("75 sample index")
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="lower center",bbox_to_anchor=(.5,.02),
               frameon=False,ncol=2,fontsize=10)
    fig.suptitle("Top 5 HELP from QC within-plate residual screening | unscaled log2 PG and Heavy deviations",
                 fontsize=14,fontweight="bold")
    save(fig,OUT/"04_top5_HELP_vs_PG_raw_trend_overlay.png",dpi=170)


def QC_overlay(pg,H,seqs,candidates,names,groups,plates,order):
    qidx=np.array([j for j in order if groups[j]=="pooled QC"],int)
    if len(qidx)!=15:raise ValueError("Need all 15 pooled QC")
    x=np.arange(1,16)
    # Anchor each plate to its OWN median: remove plate effects.
    pgq=qc_plate_residual(pg["QC_common_15_PG"],plates,groups=="pooled QC")
    # qc_plate_residual is in the ORIGINAL QC order; build an index mapping.
    raw_q_idx=np.where(groups=="pooled QC")[0]
    order_map={j:i for i,j in enumerate(raw_q_idx)}
    ordered=[order_map[j] for j in qidx]
    fig,axs=plt.subplots(5,1,figsize=(16,16),sharex=True)
    fig.subplots_adjust(left=.08,right=.985,bottom=.13,top=.92,hspace=.43)
    for k,(ax,seq) in enumerate(zip(axs,candidates)):
        hv=H[seqs.index(seq)]
        hh=qc_plate_residual(hv,plates,groups=="pooled QC")
        ax.plot(x,pgq[ordered],c="#293E52",lw=1.5,marker="o",ms=4,
                label="PG median (15-QC common PG), within-plate centered")
        ax.plot(x,hh[ordered],c="#D08343",lw=1.5,marker="o",ms=4,
                label=seq+" Heavy, within-plate centered")
        ax.axhline(0,c="#A2AFB8",lw=.85)
        ax.axvline(5.5,c="#73879A",lw=.8)
        ax.axvline(10.5,c="#73879A",lw=.8)
        ax.set_title(seq+" · only pooled QC technical-reference samples",pad=10,fontsize=11)
        ax.set_ylabel("Within-plate log2 deviation")
        ax.grid(axis="y",alpha=.15)
    axs[-1].set_xticks(x,[f"QC{j}" for j in range(1,16)],rotation=65,fontsize=8)
    handles,labels=axs[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="lower center",bbox_to_anchor=(.5,.02),
               ncol=2,frameon=False,fontsize=10)
    fig.suptitle("Pooled QC technical trends · plate-mean shifts removed by median centering",
                 fontweight="bold",fontsize=13)
    save(fig,OUT/"05_top5_HELP_vs_fixedQC_PG_within_plate_trends.png",dpi=165)


def heatmap_rank(table):
    view=table.head(25)
    cols=("overall75_allPG_r","overall75_fixedPG_r","pooledQC_fixedPG_r",
          "QC_within_plate_fixedPG_r","Study54_allPG_r","HC27_allPG_r","S27_allPG_r")
    values=view[list(cols)].to_numpy(float)
    fig,ax=plt.subplots(figsize=(13.2,11.2))
    fig.subplots_adjust(left=.20,right=.88,bottom=.2,top=.91)
    im=ax.imshow(values,vmin=-1,vmax=1,cmap="RdBu_r",aspect="auto")
    ax.set_yticks(range(len(view)),view.HELP_sequence,fontsize=8.5)
    ax.set_xticks(range(len(cols)),[
         "Overall 75 · own PG","Overall 75 · shared PG","15 pooled QC · shared PG",
         "15 QC · within-plate","54 Study · own PG","27 HC · own PG","27 S · own PG"
    ],rotation=41,ha="right",fontsize=8.5)
    for i in range(len(view)):
        for j in range(len(cols)):
            v=values[i,j]
            if np.isfinite(v):
                ax.text(j,i,f"{v:+.2f}",ha="center",va="center",
                        color="white" if abs(v)>.65 else "#243444",fontsize=8.1)
    fig.colorbar(im,ax=ax,fraction=.045,pad=.045,label="Pearson r")
    ax.set_title("Top 25 QC-trend candidate HELP | other scopes as diagnostics",
                 fontsize=13,pad=13)
    save(fig,OUT/"06_top25_correlation_heatmap.png",dpi=195)


def main():
    if OUT.exists():
        raise FileExistsError(f"Cannot replace existing trend-screening results: {OUT}")
    source_pg,names,info,raw,seqs,H,complete,groups,plates=raw_inputs()
    means,pgmask=median_reference_pg(raw,groups,plates)
    order=sample_order(names,groups,plates)
    qmask=groups=="pooled QC"
    study=np.isin(groups,("HC","S"))
    assert {g:int(sum(groups==g)) for g in GROUPS}=={
        "Internal-QC":6,"pooled QC":15,"HC":27,"S":27}
    assert {p:int(sum(plates==p)) for p in PLATES}=={"1":27,"2":27,"3":21}
    positive=np.isfinite(H)
    details=[];ranked=[]
    qcaudit={}

    for i,seq in enumerate(seqs):
        detected=int(positive[i].sum())
        item={"HELP_sequence":seq,"n_positive_of_75":detected,
              "detection_rate":detected/75.,"full_detected_75":detected==75}
        if detected:
            vals=np.exp2(H[i,positive[i]])
            item["Heavy_CV_pct_linear"]=float(np.std(vals,ddof=1)/np.mean(vals)*100) if len(vals)>1 else np.nan
            item["Heavy_median_log2"]=float(np.median(H[i,positive[i]]))
        else:
            item["Heavy_CV_pct_linear"]=np.nan
            item["Heavy_median_log2"]=np.nan
        details.append(item)
        if detected!=75:continue
        track=H[i]
        records=[]
        for scope,mask in (
            ("overall75",np.ones(75,bool)),
            ("pooledQC15",qmask),("Study54",study),
            ("HC27",groups=="HC"),("S27",groups=="S"),
            ("InternalQC6",groups=="Internal-QC"),
            ("plate1_all27",plates=="1"),
            ("plate2_all27",plates=="2"),
            ("plate3_all21",plates=="3")):
            for basis in ("all_detected_PG","common_75_PG",
                          "QC_common_15_PG","Study_common_54_PG"):
                if scope=="overall75" and basis in ("QC_common_15_PG","Study_common_54_PG"):continue
                if scope in ("plate1_all27","plate2_all27","plate3_all21") and basis in ("QC_common_15_PG","Study_common_54_PG"):continue
                records.append(metrics(track[mask],means[basis][mask],scope,seq,basis))
        qctrack=qc_plate_residual(track,plates,qmask)
        qcPG=qc_plate_residual(means["QC_common_15_PG"],plates,qmask)
        allqcpg=qc_plate_residual(means["all_detected_PG"],plates,qmask)
        records.append(metrics(qctrack,qcPG,"pooledQC15_within_plate_demeaned",seq,"QC_common_15_PG"))
        records.append(metrics(qctrack,allqcpg,"pooledQC15_within_plate_demeaned",seq,"all_detected_PG"))
        details_for_seq={ (r["scope"],r["PG_basis"]):r for r in records }
        def val(s,b,field):return details_for_seq[(s,b)][field]
        item.update({
            "overall75_allPG_r":val("overall75","all_detected_PG","pearson_r"),
            "overall75_fixedPG_r":val("overall75","common_75_PG","pearson_r"),
            "pooledQC_fixedPG_r":val("pooledQC15","QC_common_15_PG","pearson_r"),
            "QC_within_plate_fixedPG_r":val("pooledQC15_within_plate_demeaned","QC_common_15_PG","pearson_r"),
            "QC_within_plate_allPG_r":val("pooledQC15_within_plate_demeaned","all_detected_PG","pearson_r"),
            "QC_within_plate_spearman":val("pooledQC15_within_plate_demeaned","QC_common_15_PG","spearman_rho"),
            "QC_within_plate_slope_PG_per_HELP":val("pooledQC15_within_plate_demeaned","QC_common_15_PG","slope_PG_per_HELP"),
            "QC_within_plate_HELP_SD_log2":val("pooledQC15_within_plate_demeaned","QC_common_15_PG","HELP_SD_log2"),
            "QC_within_plate_PG_SD_log2":val("pooledQC15_within_plate_demeaned","QC_common_15_PG","PG_SD_log2"),
            "QC_within_plate_deviation_RMSE_log2":val("pooledQC15_within_plate_demeaned","QC_common_15_PG","deviation_rmse_log2"),
            "QC_within_plate_deviation_MAE_log2":val("pooledQC15_within_plate_demeaned","QC_common_15_PG","deviation_mae_log2"),
            "Study54_allPG_r":val("Study54","all_detected_PG","pearson_r"),
            "Study54_commonPG_r":val("Study54","Study_common_54_PG","pearson_r"),
            "HC27_allPG_r":val("HC27","all_detected_PG","pearson_r"),
            "S27_allPG_r":val("S27","all_detected_PG","pearson_r"),
            "InternalQC6_allPG_r":val("InternalQC6","all_detected_PG","pearson_r"),
            "plate1_allPG_r":val("plate1_all27","all_detected_PG","pearson_r"),
            "plate2_allPG_r":val("plate2_all27","all_detected_PG","pearson_r"),
            "plate3_allPG_r":val("plate3_all21","all_detected_PG","pearson_r"),
        })
        item["QC_direction_positive"]=bool(item["QC_within_plate_fixedPG_r"]>0)
        item["QC_slope_around_one"]=bool(.5<=item["QC_within_plate_slope_PG_per_HELP"]<=2.)
        item["QC_screening_eligible"]=bool(item["QC_direction_positive"] and item["QC_slope_around_one"])
        ranked.append(item)
        qcaudit[seq]=records

    if not ranked:raise AssertionError("No 75/75 HELP to rank")
    df=pd.DataFrame(ranked)
    # Rank by shape matching QC within-plate after removing constant plate shifts.
    # Don't make large cross-biological-sample correlations the primary criterion.
    # Require positive trend and approximate unit-scale relation where possible;
    # fall back to direction-only only when no peptide fulfills both.
    has_scale=bool(df["QC_screening_eligible"].any())
    if has_scale:pool=df[df.QC_screening_eligible]
    else:pool=df[df.QC_direction_positive]
    if pool.empty:pool=df
    pool=pool.sort_values(
         ["QC_within_plate_deviation_RMSE_log2","QC_within_plate_fixedPG_r","Heavy_CV_pct_linear"],
         ascending=[True,False,True],kind="mergesort")
    other=df[~df.HELP_sequence.isin(pool.HELP_sequence)].sort_values(
         ["QC_within_plate_fixedPG_r","QC_within_plate_deviation_RMSE_log2"],
         ascending=[False,True],kind="mergesort")
    ranked_df=pd.concat([pool,other],ignore_index=True)
    ranked_df.insert(0,"exploratory_QC_trend_rank",np.arange(1,len(ranked_df)+1))
    top5=ranked_df.HELP_sequence.head(min(5,len(ranked_df))).tolist()
    top10=ranked_df.HELP_sequence.head(min(TOP,len(ranked_df))).tolist()
    long=[]
    for seq in ranked_df.HELP_sequence:
        long.extend(qcaudit[seq])
    samples=[]
    helpmed=np.nanmedian(H[complete],axis=0)
    for k,j in enumerate(order):
        row={"sample_order_index":k+1,"sample_id":names[j],"group":groups[j],
            "plate":plates[j],"run":info[names[j]]["run"],
            "n_raw_PG_detected":int(np.isfinite(raw[:,j]).sum()),
            "n_original_HELP_Heavy_detected":int(positive[:,j].sum()),
            "n_75_complete_HELP":int(np.sum(complete)),
            "HELP_complete_Heavy_median_log2":float(helpmed[j])}
        for b,vector in means.items():
            row[b+"_log2_median"]=float(vector[j]) if np.isfinite(vector[j]) else ""
            row[b+"_relative_all75_median"]=float(vector[j]-np.nanmedian(vector)) if np.isfinite(vector[j]) else ""
        for seq in top10:
            v=H[seqs.index(seq),j]
            row[seq+"_log2_Heavy"]=float(v)
            row[seq+"_relative_75median_log2"]=float(v-np.median(H[seqs.index(seq)]))
        samples.append(row)

    OUT.mkdir(parents=True)
    pd.DataFrame(details).to_csv(OUT/"00_all97_HELP_detection_and_CV.csv",index=False)
    ranked_df.to_csv(OUT/"01_full_detect_HELP_QC_trend_ranking.csv",index=False)
    pd.DataFrame(long).to_csv(OUT/"02_all_spearman_pearson_slope_RMSE_by_scope.csv",index=False)
    pd.DataFrame(samples).to_csv(OUT/"03_sample_order_and_medians.csv",index=False)
    for label,mask in pgmask.items():
        # Same protein annotations and row numbers, no changed source matrix.
        rows=source_pg.iloc[np.where(mask)[0]][old.ANN].copy()
        rows.insert(0,"source_matrix_row_number_1based_including_header",
                    np.where(mask)[0]+2)
        rows.to_csv(OUT/f"04_{label}_PG_list.csv",index=False)
    box_figures(raw,H[complete],names,groups,plates,order,means)
    median_lines(means,H[complete],names,groups,plates,order)
    trend_panels(means,H,seqs,top5,names,groups,plates,order)
    QC_overlay(means,H,seqs,top5,names,groups,plates,order)
    heatmap_rank(ranked_df)

    qc=groups=="pooled QC"
    qcv=means["QC_common_15_PG"][qc]
    qcset=[]
    for p in PLATES:
        sel=qc&(plates==p)
        qcset.append({"plate":p,"n_QC":int(sel.sum()),
            "PG_fixed15_median_SD_log2":float(np.std(means["QC_common_15_PG"][sel],ddof=1)),
            "PG_all_detected_median_SD_log2":float(np.std(means["all_detected_PG"][sel],ddof=1))})
    pd.DataFrame(qcset).to_csv(OUT/"05_QC_PG_median_technical_SD_by_plate.csv",index=False)
    manifest={
        "input_files":{"PG":str(old.PG.name),"HELP":str(old.HELP.name),
                       "metadata":[str(old.META_QC.name),str(old.META_STUDY.name)]},
        "input_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in (old.PG,old.HELP,old.META_QC,old.META_STUDY)},
        "PACS_samples":75,"raw_PG_count":6675,
        "full_detected_HELP_count":int(complete.sum()),
        "HELP_total":97,
        "full_detected_sequence_order":[s for s,ok in zip(seqs,complete) if ok],
        "PG_fixed_sets":{"common75":int(pgmask["common_75"].sum()),
                         "QC_common15":int(pgmask["QC_common_15"].sum()),
                         "Study_common54":int(pgmask["Study_common_54"].sum())},
        "primary_screening":"QC within-plate median-centered log2 HELP versus fixed 15-QC PG medians; positive r and slope 0.5 to 2 then lowest unscaled RMSE",
        "min_slope_half_max_double_heuristic":True,
        "used_full_count":int(len(df)),
        "n_QC_scoring_eligible":int(df.QC_screening_eligible.sum()),
        "fallback_when_no_scale_candidates":not has_scale,
        "top5_exploratory":top5,"top10_exploratory":top10,
        "notes":["Screening on QC median trend is exploratory and reuses QC observations; not an independent validation.",
                 "Only 15 pooled QC (5/plate), so correlations are noisy.",
                 "Raw PG median and PG detection counts may be strongly coupled.",
                 "75 heterogenous HC/S/QC/Internal-QC means full-75 r is descriptive, not pure technical correlation.",
                 "Raw PG may already be upstream normalized; matching Heavy trend is not guaranteed to improve accuracy.",
                 "No calibration applied, no input data modified, no original experiment overwritten."],
    }
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    lines=[
        "# PACS: select 100%-detected HELP Heavy that tracks original PG median trends",
        "",
        "This is exploratory DISCOVERY, not a calibrated matrix. All references",
        "use original (uncorrected) PACS 75-sample PG and original 97 HELP Heavy H values.",
        "The three plates and all four groups are retained with original sample IDs.",
        "",
        "## Detection and fixed PG sets",
        "",
        f"- HELP H with full positive-finite detection 75/75: **{int(complete.sum())}/97**.",
        f"- PG common to all 75: **{int(pgmask['common_75'].sum())}** / 6675.",
        f"- PG common to all 15 pooled QC: **{int(pgmask['QC_common_15'].sum())}** / 6675.",
        f"- PG common to all 54 Study samples: **{int(pgmask['Study_common_54'].sum())}** / 6675.",
        "",
        "Each sample's PG log2 median is computed (a) from that sample's all",
        "detected PG, (b) a FIXED full75 common set. QC comparisons additionally",
        "use the much larger PG intersection fixed within the 15 pooled QC.",
        "All median shifts are defined by subtracting the MEDIAN across the",
        "corresponding sample population; no z-score normalization or imputation.",
        "",
        "## Ranking and safeguards",
        "",
        "1. Require one HELP Heavy sequence present in all 75 samples.",
        "2. Correlate each original Heavy H log2 trace against samplewise PG medians.",
        "3. For QC technical trend, remove each plate's QC median from both",
        "   HELP Heavy and PG (15-QC-shared PG) before Pearson, Spearman,",
        "   log2 slope, and unscaled residual RMSE. Thus plate offsets cannot",
        "   inflate within-plate technical tracking.",
        "4. Eligible candidates need positive QC within-plate Pearson r AND",
        "   slope PG vs Heavy between 0.5 and 2 (heuristic, not a validated cut).",
        "5. Sort eligible by lowest centered QC deviation RMSE; use other",
        "   correlations (overall 75, 54 Study, HC/S, per-plate) as diagnostics.",
        "If no eligible candidate exists, the script falls back to positive-r",
        "candidates, or all fully-detected sequences when none is positive.",
        "This is an explicitly exploratory same-data selection, NOT proof",
        "of superior spike-in calibration. p-values are unadjusted.",
        "",
        "## Top 10 QC within-plate trend candidates (exploratory)",
        "",
        "| Rank | HELP | QC residual Pearson r | QC residual slope PG/H | QC RMSE log2 | All75 own-PG r | All75 fixed-PG r |",
        "|---:|---|---:|---:|---:|---:|---:|"
    ]
    for x in ranked_df.head(10).itertuples():
        lines.append(
            f"| {x.exploratory_QC_trend_rank} | {x.HELP_sequence} | "
            f"{x.QC_within_plate_fixedPG_r:+.3f} | "
            f"{x.QC_within_plate_slope_PG_per_HELP:.3f} | "
            f"{x.QC_within_plate_deviation_RMSE_log2:.3f} | "
            f"{x.overall75_allPG_r:+.3f} | {x.overall75_fixedPG_r:+.3f} |")
    lines +=[
       "",
       "## Figures",
       "",
       "### All PG per-sample boxplots",
       "![Raw PG boxes](01_raw_PG_75sample_boxplots_with_median.png)",
       "",
       "### All 75/75-detected original HELP Heavy per-sample boxplots",
       "![Fully-detected HELP Heavy boxes](02_full_detection_Heavy_HELP_75sample_boxplots.png)",
       "",
       "### Compare sample median deviation trends",
       "![PG and Heavy sample medians](03_centered_median_trends_PG_vs_all_full_Heavy_HELP.png)",
       "",
       "### Top candidate sequence versus raw PG medians",
       "![Top 5 candidate traces](04_top5_HELP_vs_PG_raw_trend_overlay.png)",
       "",
       "### QC-only within-plate trend overlay: same pooled sample type",
       "![Top 5 QC traces](05_top5_HELP_vs_fixedQC_PG_within_plate_trends.png)",
       "",
       "### Top 25 correlations across comparison populations",
       "![QC and study correlation heatmap](06_top25_correlation_heatmap.png)",
       "",
       "## Important cautions",
       "",
       "- The 75 samples are heterogeneous, not technical replicates like the 39 NYD repeats.",
       "- HC and S intensity variation can reflect genuine biology, not instrument error.",
       "- The fixed all-75 PG set may be small due to different sample groups/missingness.",
       "- Strong trend similarity in the same dataset may simply learn already-normalized PG.",
       "- Selecting a HELP on QC then checking those same QC is NOT held-out validation.",
       "- Standard +fi factor uses QC reference median minus HELP sample log2 Heavy.",
       "  A raw-PG matching HELP with too little or too much amplitude can raise PG CV.",
       "- Validate any proposed factor on separate QC or an independent dataset",
       "  before promoting it to a technical-calibration method.",
       "",
       "## Data outputs",
       "",
       "- 00_all97_HELP_detection_and_CV.csv: 97 source HELP detection counts/Heavy CV.",
       "- 01_full_detect_HELP_QC_trend_ranking.csv: eligible candidates & diagnostics.",
       "- 02_all_spearman_pearson_slope_RMSE_by_scope.csv: per-peptide correlations.",
       "- 03_sample_order_and_medians.csv: 75 ordered sample indices and intensities.",
       "- 04_*_PG_list.csv: exact fixed 75/QC15/Study54 source protein identities.",
       "- 05_QC_PG_median_technical_SD_by_plate.csv: per-plate PG median stability.",
       "- run_manifest.json: source hashes, precise conditions and candidate names.",
       ""
    ]
    (OUT/"README.md").write_text("\n".join(lines),encoding="utf-8")
    print("DISCOVERY_RESULTS",json.dumps({
       "full_HELP":int(complete.sum()),"fixed_PG":manifest["PG_fixed_sets"],
       "n_eligible":manifest["n_QC_scoring_eligible"],"top5":top5,
       "top10":ranked_df.head(10)[[
          "HELP_sequence","QC_within_plate_fixedPG_r",
          "QC_within_plate_slope_PG_per_HELP",
          "QC_within_plate_deviation_RMSE_log2",
          "overall75_allPG_r","overall75_fixedPG_r",
          "pooledQC_fixedPG_r","Study54_allPG_r"]].to_dict("records")},ensure_ascii=False))
    print("PASS 97 Heavy HELP screening with 75-sample trace plot, QC ranking and source audits",flush=True)


if __name__=="__main__":
    main()
