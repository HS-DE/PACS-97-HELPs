#!/usr/bin/env python3
"""PACS Raw PG median-matching HELP screening: ONLY Study54 and pooledQC15.

Study = HC + S (54, never separately screened). pooledQC = 15 separately.
Internal-QC 6 excluded from ALL target references, detection filters and plots.
Use original H Heavy and original raw PG log2, no correction/imputation.
"""
from __future__ import annotations
import json
import hashlib
import math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import pacs_heavy_97_trend_screening as source

OUT=source.old.ROOT/"results"/"pacs_HELP_trend_Study54_pooledQC15_excludeInternalQC"
BATCHES=("1","2","3")
C={"PG":"#263F55","Heavy":"#C77946","Study":"#36799B","pooledQC":"#9470A6"}
GROUPS={
    "Study54":("HC","S"),
    "pooledQC15":("pooled QC",)
}


def met(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float)
    if len(x)!=len(y) or len(x)<4 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("Correlation input must be same length, finite, >=4")
    if np.std(x)<1e-12 or np.std(y)<1e-12:
        return {"pearson_r":np.nan,"spearman_rho":np.nan,
                "centered_RMSE_log2":np.nan,"centered_MAE_log2":np.nan,
                "Heavy_over_PG_SD":np.nan,"slope_PG_per_Heavy":np.nan}
    xc=x-np.median(x); yc=y-np.median(y)
    return {
        "pearson_r":float(pearsonr(x,y).statistic),
        "spearman_rho":float(spearmanr(x,y).statistic),
        "centered_RMSE_log2":float(np.sqrt(np.mean((xc-yc)**2))),
        "centered_MAE_log2":float(np.mean(abs(xc-yc))),
        "Heavy_over_PG_SD":float(np.std(x,ddof=1)/np.std(y,ddof=1)),
        "slope_PG_per_Heavy":float(np.dot(x-x.mean(),y-y.mean())/np.dot(x-x.mean(),x-x.mean()))
    }


def plate_residual(values,plates):
    z=np.array(values,float).copy()
    for p in BATCHES:
        mask=plates==p
        if np.count_nonzero(mask)!=5:raise ValueError("Expected 5 pooled QC in each plate")
        z[mask]-=np.median(z[mask])
    return z


def base_data():
    pg,sample_ids,info,pglog,peptides,H,_,group,plate=source.raw_inputs()
    if pglog.shape!=(6675,75) or H.shape!=(97,75):raise ValueError("Unexpected input")
    group_masks={k:np.isin(group,labels) for k,labels in GROUPS.items()}
    assert {k:int(m.sum()) for k,m in group_masks.items()}=={"Study54":54,"pooledQC15":15}
    assert np.sum(group=="Internal-QC")==6
    assert not np.any(group_masks["Study54"]&group_masks["pooledQC15"])
    assert int(np.sum(group_masks["Study54"]|group_masks["pooledQC15"]))==69
    return pg,sample_ids,info,pglog,peptides,H,group,plate,group_masks


def calc_group_data(PG,H,seqs,names,info,group,plate,mask):
    sample_ix=np.where(mask)[0]
    # Within each group sort only plate -> sample id, not HC vs S.
    ix=np.array(sorted(sample_ix,key=lambda j:(int(plate[j]),names[j])),dtype=int)
    p=plate[ix]
    filtered_PG=PG[:,ix]
    filtered_H=H[:,ix]
    fullPG=np.isfinite(filtered_PG).all(axis=1)
    fullHELP=np.isfinite(filtered_H).all(axis=1)
    n=int(len(ix))
    if n not in (54,15):raise ValueError("Unexpected subpopulation count")
    assert int(fullPG.sum())==({54:3965,15:5227}[n]),(n,int(fullPG.sum()))
    own=np.nanmedian(filtered_PG,axis=0)
    common=np.median(filtered_PG[fullPG],axis=0)
    assert np.isfinite(own).all() and np.isfinite(common).all()
    return {"index":ix,"plate":p,"H":filtered_H,"PG":filtered_PG,
            "fullPG":fullPG,"fullHELP":fullHELP,
            "ownPG":own,"fixedPG":common,"n":n}


def score_peptides(g,data,sequences):
    rows=[]
    d=data
    for i,seq in enumerate(sequences):
        if not d["fullHELP"][i]:continue
        H=d["H"][i]
        own=met(H,d["ownPG"])
        fixed=met(H,d["fixedPG"])
        out={"HELP_sequence":seq,"scope":g,"n_samples":d["n"],
             "Heavy_detected_n":int(np.isfinite(H).sum()),
             "PG_fixed_common_count":int(d["fullPG"].sum())}
        for col,stats in (("ownPG",own),("fixedPG",fixed)):
            out.update({col+"_"+key:val for key,val in stats.items()})
        if g=="pooledQC15":
            Hc=plate_residual(H,d["plate"])
            ownc=plate_residual(d["ownPG"],d["plate"])
            fixedc=plate_residual(d["fixedPG"],d["plate"])
            for name,target in (("ownPG",ownc),("fixedPG",fixedc)):
                out.update({"withinPlate_"+name+"_"+k:v for k,v in met(Hc,target).items()})
        rows.append(out)
    return pd.DataFrame(rows)


def section_boundaries(ax,plate):
    for k in range(1,len(plate)):
        if plate[k]!=plate[k-1]:
            ax.axvline(k+.5,c="#718697",ls="--",lw=.9,alpha=.7)
    ax.grid(axis="y",alpha=.15)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def figsave(fig,path):
    fig.savefig(path,dpi=178,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    if path.stat().st_size<25000:raise AssertionError("Unexpected small figure "+str(path))


def box_compare(g,d):
    # Different axes (Heavy and PG differ in absolute intensity) avoid misleading overlay.
    order_n=d["n"]
    fig,axs=plt.subplots(2,1,figsize=(max(13,order_n*.34),12.1),
                         sharex=True)
    fig.subplots_adjust(left=.075,right=.985,top=.90,bottom=.19,hspace=.37)
    pp=np.arange(1,order_n+1)
    for ax,values,title,color,main in (
        (axs[0],d["PG"],f"{g}: each sample's original PG intensity box",C["PG"],d["ownPG"]),
        (axs[1],d["H"][d["fullHELP"]],f"{g}: complete-detection Heavy HELP distribution",
         C["Heavy"],np.median(d["H"][d["fullHELP"]],axis=0)),
    ):
        cols=[values[:,j][np.isfinite(values[:,j])] for j in range(order_n)]
        bp=ax.boxplot(cols,positions=pp,widths=.6,showfliers=False,patch_artist=True,
                      manage_ticks=False,medianprops={"color":"#233E55","linewidth":1.5})
        for b in bp["boxes"]:
            b.set_facecolor(color);b.set_alpha(.32)
        ax.plot(pp,main,c=color,lw=1.35,marker=".",ms=5,label="Per-sample median")
        ax.axhline(np.median(main),c="#82929D",ls=":",lw=1.1)
        ax.set_title(title,loc="left",fontsize=12,pad=10,fontweight="bold")
        ax.set_ylabel("log2 raw intensity")
        section_boundaries(ax,d["plate"])
    axs[-1].set_xticks(pp,[str(x) for x in pp],rotation=90,fontsize=6.5)
    axs[-1].set_xlabel("Sample index, plate -> ID (see sample table)")
    handles,labels=axs[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="lower center",bbox_to_anchor=(.5,.015),
               ncol=2,frameon=False)
    fig.suptitle(f"Original PG vs all 100%-detected Heavy HELP | {g}",fontsize=15,fontweight="bold")
    figsave(fig,OUT/f"01_{g}_RawPG_fullHeavy_boxplots.png")


def overlay(g,d,seqs,ranks,metric,definition,png,label):
    selected=ranks.head(5).HELP_sequence.tolist()
    pg=d["ownPG"] if definition=="ownPG" else d["fixedPG"]
    y=pg-np.median(pg)
    x=np.arange(1,d["n"]+1)
    fig,axs=plt.subplots(5,1,figsize=(max(14,d["n"]*.35),17.4),sharex=True)
    fig.subplots_adjust(left=.075,right=.985,top=.90,bottom=.085,hspace=.42)
    for ax,seq in zip(axs,selected):
        idx=seqs.index(seq)
        orig=d["H"][idx]; hc=orig-np.median(orig)
        row=ranks[ranks.HELP_sequence==seq].iloc[0]
        ax.plot(x,y,c=C["PG"],marker=".",ms=5,lw=1.4,label="PG median deviation")
        ax.plot(x,hc,c=C["Heavy"],marker=".",ms=4,lw=1.4,label="Heavy HELP deviation")
        section_boundaries(ax,d["plate"])
        ax.axhline(0,c="#B0BEC9",lw=.75)
        ax.set_ylabel("log2 deviation")
        ax.set_title(f"{seq} · Pearson {row[definition+'_pearson_r']:+.3f} "
                     f"· RMSE {row[definition+'_centered_RMSE_log2']:.3f} "
                     f"· Heavy/PG SD {row[definition+'_Heavy_over_PG_SD']:.2f}",
                     loc="left",pad=9,fontsize=10.3)
    axs[-1].set_xticks(x,[str(k) for k in x],rotation=90,fontsize=7)
    axs[-1].set_xlabel("Sample position, sorted plate -> sample ID")
    h,l=axs[0].get_legend_handles_labels()
    fig.legend(h,l,loc="lower center",bbox_to_anchor=(.5,.014),ncol=2,frameon=False)
    fig.suptitle(f"{g}: top five by {label} against sample {definition} median",
                 fontsize=14,fontweight="bold")
    figsave(fig,OUT/png)


def pooled_within_overlay(d,seqs,ranking,png):
    pp=d["plate"]
    y=plate_residual(d["fixedPG"],pp)
    x=np.arange(1,16)
    fig,axs=plt.subplots(5,1,figsize=(14.5,17),sharex=True)
    fig.subplots_adjust(left=.09,right=.985,top=.91,bottom=.09,hspace=.40)
    for ax,seq in zip(axs,ranking.HELP_sequence.head(5)):
        H=d["H"][seqs.index(seq)]
        heavy=plate_residual(H,pp)
        r=ranking.set_index("HELP_sequence").loc[seq]
        ax.plot(x,y,lw=1.4,marker="o",ms=4,c=C["PG"],label="pooled QC fixed-PG median")
        ax.plot(x,heavy,lw=1.3,marker="o",ms=4,c=C["Heavy"],label="HELP Heavy")
        ax.axhline(0,c="#C0C8D0",lw=.7)
        section_boundaries(ax,pp)
        ax.set_ylabel("Within-plate log2 deviation")
        ax.set_title(f"{seq} · within-plate Pearson {r.withinPlate_fixedPG_pearson_r:+.3f}; "
                     f"RMSE {r.withinPlate_fixedPG_centered_RMSE_log2:.3f}",
                     fontsize=10.5,loc="left",pad=11)
    axs[-1].set_xticks(x,[str(i) for i in x])
    axs[-1].set_xlabel("Pooled QC sample index (five from each plate)")
    h,l=axs[0].get_legend_handles_labels()
    fig.legend(h,l,loc="lower center",bbox_to_anchor=(.5,.01),ncol=2,frameon=False)
    fig.suptitle("Only pooled QC 15 | fixed QC-common PG | plate-median centered",
                 fontsize=14,fontweight="bold")
    figsave(fig,OUT/png)


def main():
    if OUT.exists():raise FileExistsError(f"Do not overwrite historical run: {OUT}")
    pg,names,info,PG,seqs,H,groups,plates,masks=base_data()
    data={key:calc_group_data(PG,H,seqs,names,info,groups,plates,mask)
          for key,mask in masks.items()}
    scores={key:score_peptides(key,d,seqs) for key,d in data.items()}
    assert len(scores["Study54"])>=96 and len(scores["pooledQC15"])>=96
    OUT.mkdir(parents=True)
    figs=[]
    all_meta={}
    for group,d in data.items():
        group_rows=[]
        for k,j in enumerate(d["index"],1):
            row={"sample_order":k,"sample_id":names[j],
                 "research_group":group,"original_group":groups[j],
                 "plate":plates[j],
                 "PG_detected_count":int(np.isfinite(d["PG"][:,k-1]).sum()),
                 "HELP_Heavy_detected_count":int(np.isfinite(d["H"][:,k-1]).sum()),
                 "PG_all_detected_median_log2":float(d["ownPG"][k-1]),
                 "PG_fixed_group_common_median_log2":float(d["fixedPG"][k-1]),
                 "full_detected_HELP_median_log2":float(np.median(d["H"][d["fullHELP"],k-1]))}
            group_rows.append(row)
        pd.DataFrame(group_rows).to_csv(OUT/f"00_{group}_sample_order_and_medians.csv",index=False)
        pg_rows=pg.iloc[np.where(d["fullPG"])[0]][source.old.ANN].copy()
        pg_rows.insert(0,"source_pg_row_number_including_header",np.where(d["fullPG"])[0]+2)
        pg_rows.to_csv(OUT/f"00_{group}_common_PG_list.csv",index=False)
        stats=scores[group]
        overall_best_r=stats.sort_values("ownPG_pearson_r",ascending=False,kind="stable")
        overall_best_rmse=stats.sort_values("ownPG_centered_RMSE_log2",kind="stable")
        fixed_best_r=stats.sort_values("fixedPG_pearson_r",ascending=False,kind="stable")
        fixed_best_rmse=stats.sort_values("fixedPG_centered_RMSE_log2",kind="stable")
        stats.to_csv(OUT/f"01_{group}_completeHELP_all_metrics.csv",index=False)
        overall_best_r.to_csv(OUT/f"02_{group}_ranking_ownPG_Pearson.csv",index=False)
        overall_best_rmse.to_csv(OUT/f"03_{group}_ranking_ownPG_RMSE.csv",index=False)
        fixed_best_r.to_csv(OUT/f"04_{group}_ranking_fixedPG_Pearson.csv",index=False)
        fixed_best_rmse.to_csv(OUT/f"05_{group}_ranking_fixedPG_RMSE.csv",index=False)
        box_compare(group,d)
        figs.append(f"01_{group}_RawPG_fullHeavy_boxplots.png")
        for definition,metric,rank,label in (
            ("ownPG","pearson_r",overall_best_r,"highest Pearson"),
            ("ownPG","centered_RMSE_log2",overall_best_rmse,"lowest median-centered RMSE"),
            ("fixedPG","centered_RMSE_log2",fixed_best_rmse,"lowest fixed-PG RMSE")):
            name=f"06_{group}_{definition}_{metric}_top5_overlays.png"
            overlay(group,d,seqs,rank,metric,definition,name,label)
            figs.append(name)
        all_meta[group]={
            "sample_count":d["n"],"sample_IDs":[names[i] for i in d["index"]],
            "original_groups":list(GROUPS[group]),
            "PG_full_detection_proteins":int(d["fullPG"].sum()),
            "full_detection_HELP_count":int(d["fullHELP"].sum()),
            "full_detection_HELP_sequences":[seqs[i] for i in np.where(d["fullHELP"])[0]],
            "PG_median_log2_SD_own":float(np.std(d["ownPG"],ddof=1)),
            "PG_median_log2_SD_fixed":float(np.std(d["fixedPG"],ddof=1)),
            "best_by_own_r":overall_best_r.HELP_sequence.head(5).tolist(),
            "best_by_own_RMSE":overall_best_rmse.HELP_sequence.head(5).tolist(),
            "best_by_fixed_r":fixed_best_r.HELP_sequence.head(5).tolist(),
            "best_by_fixed_RMSE":fixed_best_rmse.HELP_sequence.head(5).tolist(),
            "max_own_r":float(overall_best_r.ownPG_pearson_r.max()),
            "max_fixed_r":float(fixed_best_r.fixedPG_pearson_r.max()),
            "count_own_r_ge_09":int((overall_best_r.ownPG_pearson_r>=.9).sum()),
            "count_fixed_r_ge_09":int((fixed_best_r.fixedPG_pearson_r>=.9).sum()),
        }

    qc=scores["pooledQC15"]
    qcp=qc.sort_values("withinPlate_fixedPG_pearson_r",ascending=False,kind="stable")
    qcr=qc.sort_values("withinPlate_fixedPG_centered_RMSE_log2",kind="stable")
    qcp.to_csv(OUT/"07_pooledQC15_withinPlate_fixedPG_Pearson_rank.csv",index=False)
    qcr.to_csv(OUT/"08_pooledQC15_withinPlate_fixedPG_RMSE_rank.csv",index=False)
    pooled_within_overlay(data["pooledQC15"],seqs,qcr,"09_pooledQC15_withinPlate_fixedPG_top5_RMSE.png")
    figs.append("09_pooledQC15_withinPlate_fixedPG_top5_RMSE.png")
    all_meta["pooledQC15"]["within_plate_top5_by_RMSE"]=qcr.HELP_sequence.head(5).tolist()
    all_meta["pooledQC15"]["within_plate_top5_by_r"]=qcp.HELP_sequence.head(5).tolist()
    all_meta["pooledQC15"]["within_plate_max_r"]=float(qcp.iloc[0].withinPlate_fixedPG_pearson_r)
    all_meta["pooledQC15"]["within_plate_min_RMSE"]=float(qcr.iloc[0].withinPlate_fixedPG_centered_RMSE_log2)

    manifest={
        "user_requested_scope":"Only Study54 (HC+S combined) and pooledQC15 separately; 6 Internal-QC excluded",
        "n_samples_researched":69,"n_excluded_InternalQC":6,
        "group_summaries":all_meta,
        "PG_matrix_raw":str(source.old.PG.name),
        "HELP_raw_H_sheet":str(source.old.HELP.name),
        "source_hashes_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                 (source.old.PG,source.old.HELP,source.old.META_QC,source.old.META_STUDY)},
        "definitions":{
             "ownPG":"For each sample, median(log2 positive raw PG); protein sets may vary",
             "fixedPG":"Group-specific PG intersection present in ALL 54 Study or ALL 15 pooled QC respectively",
             "Heavy_all_detected":"Separate full-detection test within each research group, not 75 samples",
             "correlation":"Pearson across samples IN THE GROUP only",
             "RMSE":"RMSE of two log2 curves after subtracting each group's own median, without amplitude scaling",
             "withinPlate_QC":"Remove median within each pooled QC plate from both signals before Pearson/RMSE",
             "plate_sort":"sample positions ordered plate -> sample ID; HC and S not separated or individually ranked",
             "no_data_corrections":"Only original raw data read; no PG/HELP normalization or missingness fill"
        },
        "generated_figures":figs
    }
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n")
    lines=[
      "# PACS Heavy HELP vs Raw protein sample-median trends | Study54 and pooledQC15",
      "",
      "**No Internal-QC.** All research is split into two separate and non-overlapping groups:",
      "**Study54 = HC27 + S27 combined**, and **pooledQC15 = 15 pooled QC**.",
      "Sample positions within each group are sorted plate then sample_id.",
      "",
      "## Fixed protein universes",
      "",
      f"- Study54 PG common to all 54 samples: {all_meta['Study54']['PG_full_detection_proteins']}.",
      f"- pooledQC15 PG common to all 15 QC samples: {all_meta['pooledQC15']['PG_full_detection_proteins']}.",
      f"- Fully detected Heavy HELP sequences in Study54: {all_meta['Study54']['full_detection_HELP_count']}.",
      f"- Fully detected Heavy HELP sequences in pooledQC15: {all_meta['pooledQC15']['full_detection_HELP_count']}.",
      "",
      "PG own-detected medians reflect varying detection sets. FixedPG means the",
      "SAME PG rows per group are used in every sample within that group.",
      "No z-score / regression scaling is applied to the matching RMSE.",
      "Only original, log2 positive intensities are used.",
      "",
      "## Best candidates for each group, using each sample's own-detected PG median",
      ""
    ]
    for group in GROUPS:
        lines += [f"### {group} (n={all_meta[group]['sample_count']})",
                  "",
                  "| Top rank | HELP | Pearson r | Centered log2 RMSE | SD Heavy / SD PG |",
                  "| ---: | --- | ---: | ---: | ---: |"]
        top=scores[group].sort_values("ownPG_centered_RMSE_log2").head(5)
        for k,row in enumerate(top.itertuples(),1):
            lines.append(f"| {k} | {row.HELP_sequence} | "
                         f"{row.ownPG_pearson_r:+.3f} | {row.ownPG_centered_RMSE_log2:.3f} | "
                         f"{row.ownPG_Heavy_over_PG_SD:.2f} |")
        lines+=["",f"![{group} boxes](01_{group}_RawPG_fullHeavy_boxplots.png)",
                "",f"![{group} own PG r](06_{group}_ownPG_pearson_r_top5_overlays.png)",
                "",f"![{group} own PG RMSE](06_{group}_ownPG_centered_RMSE_log2_top5_overlays.png)",
                "",f"![{group} fixed PG RMSE](06_{group}_fixedPG_centered_RMSE_log2_top5_overlays.png)",
                ""]
    lines += [
        "## pooledQC15 plate-centered technical-trend diagnostic",
        "",
        "Five pooled QC per plate. Both Heavy and PG curves have the own-plate",
        "median removed before comparing their technical fluctuations.",
        "",
        "![QC within plate](09_pooledQC15_withinPlate_fixedPG_top5_RMSE.png)",
        "",
        "## Cautions",
        "",
        "- HC and S are jointly analyzed as 54 Study samples; their real biology remains.",
        "- QC trends are more relevant to technical calibration, but only 15 QC exist.",
        "- PG may already be upstream-normalized; Heavy may be more noisy.",
        "- A HELP with high r but tiny Heavy/PG SD ratio does not match raw amplitude.",
        "- No selected candidate is claimed to be a validated calibration reference.",
        "- RMSE-based same-data selection must be tested on held-out QCs.",
        "",
        "## Audit files",
        "",
        "- 00_*_sample_order_and_medians.csv: exact included sample IDs",
        "- 00_*_common_PG_list.csv: source PG identity of fixed intersection",
        "- 01_*_completeHELP_all_metrics.csv: independent group-wise metrics",
        "- 02 to 05 group-specific ranking tables; 07,08 QC within-plate rankings",
        "- run_manifest.json includes exact source SHA256 and sample lists",
    ]
    (OUT/"README.md").write_text("\n".join(lines)+"\n")
    print("GROUP_RESULTS",json.dumps({k:{kk:v for kk,v in meta.items()
        if kk!="sample_IDs" and kk!="full_detection_HELP_sequences"}
       for k,meta in all_meta.items()},ensure_ascii=False))
    print("PASS 54 combined HC/S Study + 15 pooled QC; Internal-QC entirely excluded",flush=True)

if __name__=="__main__":
    main()
