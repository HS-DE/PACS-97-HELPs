#!/usr/bin/env python3
"""Scheme A exploratory screen. Heavy in all 69 study+pooled QC, pooled-QC
linear Heavy CV <10%, and 15 pooled QC median trend (Pearson, centered RMSE,
plate offset directions). Excludes Internal-QC and residual MAD.
No intensity correction or calibration is applied.
"""
from __future__ import annotations
import json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pacs_HELP_trend_study54_pooledQC15 as prev

OUT=prev.source.old.ROOT/"results"/"pacs_schemeA_QC15_Heavy_CV10_median_trend_screen"
PLATES=("1","2","3")
COLPG="#34485F"; COLH="#147D74"; EPS=.005

def cv(logs):
    x=np.exp2(np.asarray(logs,dtype=float))
    if not np.isfinite(x).all():return np.nan
    return float(100*np.std(x,ddof=1)/np.mean(x))

def plate_med(x,p):
    return np.array([np.median(x[p==z]) for z in PLATES])

def platediff(h,pg):
    # The three plate comparisons have only two degrees of freedom.
    rec={};total=0;important=0;n_imp=0
    for a,b,name in ((0,1,"2v1"),(0,2,"3v1"),(1,2,"3v2")):
        dh=float(h[b]-h[a]);dp=float(pg[b]-pg[a])
        same=bool(np.sign(dh)==np.sign(dp)) if dh!=0 and dp!=0 else False
        imp=bool(abs(dp)>=EPS)
        total+=int(same);important+=int(same and imp);n_imp+=int(imp)
        rec.update({"plate_"+name+"_Heavy_delta":dh,
                    "plate_"+name+"_PG_delta":dp,
                    "plate_"+name+"_same_direction":same,
                    "plate_"+name+"_PG_informative":imp})
    return {**rec,"plate_sign_same_of3":total,
            "plate_sign_same_informative":important,
            "plate_informative_count":n_imp,
            "plate_centered_RMSE_log2":float(np.sqrt(np.mean(
                ((h-np.median(h))-(pg-np.median(pg)))**2)))}

def figsave(fig,name):
    f=OUT/name;fig.savefig(f,dpi=185,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    assert f.stat().st_size>17000, (name,f.stat().st_size)

def boundaries(ax,plate):
    for j in range(1,len(plate)):
        if plate[j]!=plate[j-1]:ax.axvline(j+.5,c="#9BA9B5",ls="--",lw=.8)
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    ax.grid(axis="y",alpha=.15)

def overlay(d,seqs,rank,kind,criterion):
    n=min(5,len(rank))
    fig,axes=plt.subplots(n,1,figsize=(17,3.2*n+1.6),sharex=True)
    fig.subplots_adjust(left=.09,right=.985,top=.90,bottom=.09,hspace=.50)
    axes=np.atleast_1d(axes)
    x=np.arange(1,16);pg=d[kind]-np.median(d[kind])
    for ax,(_,row) in zip(axes,rank.head(n).iterrows()):
        seq=row.HELP_sequence;h=d["H"][seqs.index(seq)]
        h=h-np.median(h)
        ax.plot(x,pg,"-o",lw=1.5,ms=4.5,c=COLPG,label="QC PG median")
        ax.plot(x,h,"-D",lw=1.3,ms=3.7,c=COLH,label="Heavy "+seq)
        ax.axhline(0,c="#B8C6CE",lw=.7)
        boundaries(ax,d["plate"])
        ax.set_ylabel("Centered log2")
        ax.set_title(f"{seq} | CV {row.CV_QC15_linear_pct:.2f}% | "
                     f"r {row[kind+'_pearson_r']:+.3f} | "
                     f"RMSE {row[kind+'_centered_RMSE_log2']:.3f} | "
                     f"plate {int(row[kind+'_plate_sign_same_informative'])}/"
                     f"{int(row[kind+'_plate_informative_count'])}",
                     loc="left",fontsize=10.4,pad=8)
    axes[-1].set_xticks(x,[str(j) for j in x])
    axes[-1].set_xlabel("QC sample display index (plate -> sample ID, not acquisition time)")
    h,l=axes[0].get_legend_handles_labels()
    fig.legend(h,l,loc="lower center",bbox_to_anchor=(.5,.012),frameon=False,ncol=2)
    fig.suptitle(f"CV<10% Heavy: top {n} by {criterion} vs {kind} QC PG median",
                 fontsize=14,fontweight="bold")
    name=f"01_{kind}_{criterion}_top5_QC_curve.png"
    figsave(fig,name)
    return name

def plates_fig(d,seqs,rank,kind):
    n=min(8,len(rank));nr=int(np.ceil(n/4))
    fig,axes=plt.subplots(nr,4,figsize=(19,4.3*nr+1.5))
    fig.subplots_adjust(left=.065,right=.98,top=.85,bottom=.16,hspace=.48,wspace=.3)
    axlist=np.array(axes).ravel()
    ref=plate_med(d[kind],d["plate"]);ref-=np.median(ref)
    for i,(_,r) in enumerate(rank.head(n).iterrows()):
        ax=axlist[i]
        h=plate_med(d["H"][seqs.index(r.HELP_sequence)],d["plate"])
        h=h-np.median(h)
        ax.plot([1,2,3],ref,"-o",c=COLPG,lw=1.6,ms=6)
        ax.plot([1,2,3],h,"-D",c=COLH,lw=1.5,ms=5)
        ax.axhline(0,c="#C0C9D0",lw=.7)
        ax.set_title(f"{r.HELP_sequence} | CV {r.CV_QC15_linear_pct:.2f}%\n"
                     f"plate sign {int(r[kind+'_plate_sign_same_informative'])}/"
                     f"{int(r[kind+'_plate_informative_count'])}",fontsize=9.5)
        ax.set_xticks([1,2,3],["P1","P2","P3"])
        ax.set_ylabel("Plate median delta log2")
        ax.grid(axis="y",alpha=.17)
        ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    for ax in axlist[n:]:ax.axis("off")
    fig.suptitle(f"QC PG plate median vs Heavy | {kind} | best matching directions and magnitudes",
                 fontsize=13,fontweight="bold")
    fig.legend([Line2D([0],[0],c=COLPG,marker="o",lw=2),
                Line2D([0],[0],c=COLH,marker="D",lw=2)],
               ["QC PG plate median","HELP Heavy plate median"],
               ncol=2,loc="lower center",bbox_to_anchor=(.5,.028),frameon=False)
    name=f"02_{kind}_plate_directions_top8.png"
    figsave(fig,name)
    return name

def main():
    if OUT.exists():raise FileExistsError("Guarded output already exists: "+str(OUT))
    pg,names,info,PG,seqs,H,groups,plates,masks=prev.base_data()
    research=masks["Study54"]|masks["pooledQC15"]
    assert sum(research)==69 and sum(masks["pooledQC15"])==15
    d=prev.calc_group_data(PG,H,seqs,names,info,groups,plates,masks["pooledQC15"])
    assert len(seqs)==97 and d["n"]==15 and d["fullPG"].sum()==5227
    valid69=np.isfinite(H[:,research]).all(axis=1)
    rows=[]
    target={kind:plate_med(d[kind],d["plate"]) for kind in ("ownPG","fixedPG")}
    for i,seq in enumerate(seqs):
        h=d["H"][i];c=cv(h)
        ok=bool(valid69[i] and np.isfinite(c) and c<10)
        rec={"HELP_sequence":seq,"Heavy_detected_out_of_69":int(np.isfinite(H[i,research]).sum()),
             "full69_detected":bool(valid69[i]),
             "QC15_detected":int(np.isfinite(h).sum()),
             "CV_QC15_linear_pct":c,"pass_69_CV_lt10":ok,
             "QC_plate1_CV_linear_pct":cv(h[d["plate"]=="1"]),
             "QC_plate2_CV_linear_pct":cv(h[d["plate"]=="2"]),
             "QC_plate3_CV_linear_pct":cv(h[d["plate"]=="3"])}
        for kind in ("ownPG","fixedPG"):
            rec.update({kind+"_"+key:v for key,v in prev.met(h,d[kind]).items()})
            rec.update({kind+"_"+key:v for key,v in platediff(plate_med(h,d["plate"]),target[kind]).items()})
        rows.append(rec)
    allhelp=pd.DataFrame(rows)
    selected=allhelp[allhelp.pass_69_CV_lt10].copy()
    OUT.mkdir(parents=True)
    allhelp.to_csv(OUT/"00_all97_detection_CV_trend_metrics.csv",index=False)
    qc_rows=[]
    for k,j in enumerate(d["index"]):
        qc_rows.append({"QC_display_position":k+1,"sample_id":names[j],"plate":d["plate"][k],
                        "PG_detected_n":int(np.isfinite(d["PG"][:,k]).sum()),
                        "PG_own_sample_box_median_log2":d["ownPG"][k],
                        "PG_fixed5227_median_log2":d["fixedPG"][k]})
    pd.DataFrame(qc_rows).to_csv(OUT/"01_QC15_sample_medians.csv",index=False)
    plate_rows=[{"plate":pl,"n_QC":int(sum(d["plate"]==pl)),
                 "PG_own_median_log2":target["ownPG"][i],
                 "PG_fixed5227_median_log2":target["fixedPG"][i]}
                for i,pl in enumerate(PLATES)]
    pd.DataFrame(plate_rows).to_csv(OUT/"02_QC_plate_PG_medians.csv",index=False)
    meta={"n_original_HELP":97,"n_study":54,"n_pooled_QC":15,
          "excluded_Internal_QC":6,"n_full69":int(sum(valid69)),
          "n_pass_CV10":len(selected),
          "CV_definition":"100*sd(linear Heavy in pooledQC15,ddof=1)/mean(linear Heavy)",
          "CV_threshold_percent_strict_lt":10,
          "MAD_excluded":True,"normalized":False,"imputed":False,
          "ownPG":"Each sample box median of all its detected PG (variable protein set)",
          "fixedPG":"Each sample median of identical 5227 QC-common PG",
          "plate_direction":"sign of Heavy and PG plate-median differences; all 3 pairs reported, 2 independent",
          "sign_interesting_PG_difference_log2":EPS,
          "PG_plate_medians":plate_rows,
          "data_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (prev.source.old.PG,prev.source.old.HELP,
                                   prev.source.old.META_QC,prev.source.old.META_STUDY)}}
    figs=[]
    if not selected.empty:
        for kind in ("ownPG","fixedPG"):
            ranks={
              "pearson":selected.sort_values([kind+"_pearson_r",kind+"_centered_RMSE_log2"],
                      ascending=[False,True],kind="stable"),
              "RMSE":selected.sort_values([kind+"_centered_RMSE_log2",kind+"_pearson_r"],
                      ascending=[True,False],kind="stable"),
              "plate":selected.sort_values([kind+"_plate_sign_same_informative",
                    kind+"_plate_centered_RMSE_log2",kind+"_pearson_r"],
                      ascending=[False,True,False],kind="stable"),
            }
            for criterion,rank in ranks.items():
                rank.to_csv(OUT/f"03_{kind}_{criterion}_CV10_ranking.csv",index=False)
                meta[f"top5_{kind}_{criterion}"]=rank.HELP_sequence.head(5).tolist()
            figs.append(overlay(d,seqs,ranks["pearson"],kind,"pearson"))
            figs.append(overlay(d,seqs,ranks["RMSE"],kind,"RMSE"))
            figs.append(plates_fig(d,seqs,ranks["plate"],kind))
        selected.sort_values("ownPG_centered_RMSE_log2").to_csv(
            OUT/"04_all_CV10_peptides_metrics.csv",index=False)
    meta["generated_figures"]=figs
    (OUT/"run_manifest.json").write_text(json.dumps(meta,indent=2,ensure_ascii=False)+"\n")
    readme=["# Scheme A: pooled QC Heavy trend screening","",
      "- Heavy fully detected in all 69 non-Internal-QC samples.",
      "- QC linear Heavy CV must be strictly <10% across ALL 15 pooled QC samples.",
      "- No residual MAD, no normalization, no imputation.",
      "- Compare each CV-passed Heavy with 15 QC PG sample-box medians:",
      "  Pearson, unscaled median-centered log2 RMSE, plate-pair direction.",
      "- Also use fixed 5227 PG present in all 15 pooled QC for robustness.",
      "- Compare all 3 plate pairs but recognize only 2 are independent.",
      "- A PG plate contrast <0.005 log2 is flagged as weak evidence.",
      "",
      f"Original HELP: {len(seqs)}, full 69 detected: {int(valid69.sum())},",
      f"QC Heavy CV<10% passing: **{len(selected)}**.","",
      "## Reference QC PG plate medians","",
      "| Plate | per-sample box median log2 | fixed 5227 PG log2 |",
      "|---|---:|---:|"]
    for x in plate_rows:
        readme.append(f"| {x['plate']} | {x['PG_own_median_log2']:.5f} | {x['PG_fixed5227_median_log2']:.5f} |")
    if not selected.empty:
        for kind in ("ownPG","fixedPG"):
            rank=pd.read_csv(OUT/f"03_{kind}_RMSE_CV10_ranking.csv")
            readme+=["",f"## {kind}: top five by RMSE","",
                     "| HELP | CV QC15 % | Pearson | RMSE log2 | plate same/important |",
                     "|---|---:|---:|---:|---:|"]
            for r in rank.head(5).itertuples():
                readme.append(f"| {r.HELP_sequence} | {r.CV_QC15_linear_pct:.2f} | "
                              f"{getattr(r,kind+'_pearson_r'):+.3f} | "
                              f"{getattr(r,kind+'_centered_RMSE_log2'):.4f} | "
                              f"{int(getattr(r,kind+'_plate_sign_same_informative'))}/"
                              f"{int(getattr(r,kind+'_plate_informative_count'))} |")
        readme+=["","## Charts",""]
        for fig in figs:readme+=["![chart]("+fig+")",""]
    readme+=["","## Files",
      "- 00_all97_detection_CV_trend_metrics.csv (all 97 with gates and measures)",
      "- 01_QC15_sample_medians.csv (samples and actual PG reference)",
      "- 02_QC_plate_PG_medians.csv (plate-level PG reference)",
      "- 03_ownPG_*_CV10_ranking.csv / 03_fixedPG_*_CV10_ranking.csv",
      "- 04_all_CV10_peptides_metrics.csv (qualified only)",
      "- run_manifest.json (conditions and source hashes)","",
      "Reminder: HC+S used only for Heavy 69/69 detection gate. No HC/S PG median",
      "used to train QC trends. QC is identical material but 15 samples / 3 plates",
      "is too small to claim external validity. Primary and fixed QC PG medians",
      "can have opposing plate directions. No arbitrary Pearson cutoff introduced."]
    (OUT/"README.md").write_text("\n".join(readme)+"\n")
    print("RESULT",json.dumps({k:v for k,v in meta.items()
      if k in ("n_full69","n_pass_CV10","PG_plate_medians") or k.startswith("top5_")},
      ensure_ascii=False))
    if not selected.empty:
        for kind in ("ownPG","fixedPG"):
            for criterion in ("pearson","RMSE","plate"):
                t=pd.read_csv(OUT/f"03_{kind}_{criterion}_CV10_ranking.csv")
                print("TOP",kind,criterion,t.head(8)[[
                    "HELP_sequence","CV_QC15_linear_pct",
                    kind+"_pearson_r",kind+"_centered_RMSE_log2",
                    kind+"_plate_sign_same_of3",kind+"_plate_sign_same_informative",
                    kind+"_plate_informative_count",kind+"_plate_centered_RMSE_log2"
                ]].to_json(orient="records"))
    print("DONE")
if __name__=="__main__":main()
