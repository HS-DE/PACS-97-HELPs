#!/usr/bin/env python3
"""Study54 TLYSSSPR Heavy versus sample PG median, plate-color trends."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pacs_HELP_trend_study54_pooledQC15 as base

OUT=base.source.old.ROOT/"results"/"pacs_TLYSSSPR_Study54_plate_trends"
COL={"1":"#3075B6","2":"#D17C31","3":"#248878"}

def makefig(study,pg,heavy,kind,center):
    plates=study["plate"].to_numpy().astype(str)
    names=study["sample_id"].to_numpy()
    x=np.arange(1,55)
    pm=pg-np.median(pg) if center else pg
    hm=heavy-np.median(heavy) if center else heavy
    stats=base.met(heavy,pg)
    fig,ax=plt.subplots(figsize=(24,8.4))
    fig.subplots_adjust(left=.07,right=.99,top=.78,bottom=.32)
    if center:
        for p in ("1","2","3"):
            sel=plates==p
            ax.plot(x[sel],pm[sel],color=COL[p],marker="o",markersize=5,linestyle="-",linewidth=1.55)
            ax.plot(x[sel],hm[sel],color=COL[p],marker="D",markersize=4,linestyle="--",linewidth=1.55)
        ax.axhline(0,color="#A7B2BD",linewidth=.8)
        ax.set_ylabel("log2 deviation from Study54 median")
    else:
        # Original log2 levels differ greatly. Two y axes prevent invisibility.
        ax2=ax.twinx()
        for p in ("1","2","3"):
            sel=plates==p
            ax.plot(x[sel],pm[sel],color=COL[p],marker="o",markersize=5,linestyle="-",linewidth=1.55)
            ax2.plot(x[sel],hm[sel],color=COL[p],marker="D",markersize=4,linestyle="--",linewidth=1.55)
        ax.set_ylabel("Raw PG median log2")
        ax2.set_ylabel("TLYSSSPR Heavy log2")
    for i in range(1,54):
        if plates[i]!=plates[i-1]:ax.axvline(i+.5,color="#8496A7",linestyle=":",linewidth=1.15)
    for p in ("1","2","3"):
        ix=np.where(plates==p)[0]
        ax.text(np.mean(ix+1),1.04,"Plate "+p+" (n="+str(len(ix))+")",ha="center",
                color=COL[p],weight="bold",fontsize=10,transform=ax.get_xaxis_transform())
    ax.set_xlim(.5,54.5)
    ax.set_xticks(x,names,rotation=90,fontsize=6.7)
    ax.grid(axis="y",alpha=.15)
    ax.set_xlabel("54 Study samples, HC+S, ordered plate -> sample ID")
    ax.set_title("TLYSSSPR Heavy vs PG median | "+kind+
                 " | "+("median-centered" if center else "raw log2 with dual axes")+
                 "\nPearson r = {:+.4f} | median-centered RMSE = {:.4f} log2".format(
                     stats["pearson_r"],stats["centered_RMSE_log2"]),
                 loc="left",fontsize=15,weight="bold",pad=47)
    handles=[Line2D([0],[0],color=COL[p],lw=3,label="Plate "+p) for p in ("1","2","3")]
    handles.extend([Line2D([0],[0],color="#506278",marker="o",lw=1.5,label="PG median"),
                    Line2D([0],[0],color="#506278",marker="D",ls="--",lw=1.5,label="TLYSSSPR Heavy")])
    fig.legend(handles=handles,ncol=5,loc="lower center",bbox_to_anchor=(.5,.115),
               frameon=False,fontsize=10)
    fn=("02_" if center else "01_")+kind+("_centered.png" if center else "_raw_log2.png")
    fig.savefig(OUT/fn,dpi=180,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    return fn,stats

def main():
    if OUT.exists():raise FileExistsError(str(OUT))
    pg,names,info,logPG,seqs,H,groups,plates,masks=base.base_data()
    d=base.calc_group_data(logPG,H,seqs,names,info,groups,plates,masks["Study54"])
    assert d["n"]==54 and int(d["fullPG"].sum())==3965
    assert set(groups[d["index"]])=={"HC","S"}
    assert {p:int(sum(d["plate"]==p)) for p in ("1","2","3")}=={"1":20,"2":20,"3":14}
    heavy=d["H"][seqs.index("TLYSSSPR")]
    assert np.isfinite(heavy).all()
    study=pd.DataFrame({"sample_id":[names[j] for j in d["index"]],
                        "plate":d["plate"],
                        "HC_or_S":groups[d["index"]]})
    for key,col in [("PG_box_median",d["ownPG"]),("PG_fixed3965_median",d["fixedPG"])]:
        study[key+"_log2"]=col
        study[key+"_centered"]=col-np.median(col)
    study["TLYSSSPR_Heavy_log2"]=heavy
    study["Heavy_centered"]=heavy-np.median(heavy)
    OUT.mkdir(parents=True)
    study.to_csv(OUT/"Study54_TLYSSSPR_and_PG_values.csv",index=False)
    report={"n_Study":54,"n_HC":27,"n_S":27,"InternalQC_used":False,
            "pooledQC_used":False,"PG_fixed_count":3965,"figures":[],
            "metrics":{}}
    for name,col in [("PG_box_median",d["ownPG"]),("PG_fixed3965_median",d["fixedPG"])]:
        for centered in (False,True):
            fig,st=makefig(study,col,heavy,name,centered)
            report["figures"].append(fig)
            report["metrics"][name]=st
    (OUT/"summary.json").write_text(json.dumps(report,indent=2))
    readme=["# TLYSSSPR in Study54 (HC+S); plate colored","",
            "Original PG and Heavy from PACS. No Internal-QC, pooled QC, normalization, or missing-value fill.",
            "54 Study samples: Plate1=20, Plate2=20, Plate3=14. Solid circles=PG median, dashed diamonds=Heavy.",
            "",
            "| Reference | Study Pearson | Study centered RMSE (log2) |",
            "|---|---:|---:|"]
    for name,st in report["metrics"].items():
        readme.append("| "+name+" | {:+.4f} | {:.4f} |".format(st["pearson_r"],st["centered_RMSE_log2"]))
    for fig in report["figures"]:readme.extend(["","!["+fig+"]("+fig+")"])
    readme.extend(["","Full 54-sample values in Study54_TLYSSSPR_and_PG_values.csv."])
    (OUT/"README.md").write_text("\n".join(readme)+"\n")
    print("STUDY54_METRICS",json.dumps(report))
if __name__=="__main__":main()
