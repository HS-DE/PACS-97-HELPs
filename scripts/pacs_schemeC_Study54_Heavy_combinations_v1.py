#!/usr/bin/env python3
"""Scheme C, exploratory: match signed per-Study-sample Raw PG median deviations
by equally weighted log2 Heavy HELP combinations (1,2,3 sequences).

Study54 is HC27+S27 together; pooledQC15 is a separate diagnostic.
Internal-QC excluded completely. Prior QC rule: full69 and linear QC CV <10%.
Search both per-sample detected-PG and fixed Study54-common PG medians.
No calibration is applied to the original PG matrix; QC correction is an
unbiased-to-Study objective diagnostic, not evidence of external validation.
"""
from __future__ import annotations
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

import pacs_HELP_trend_study54_pooledQC15 as core

OUT=core.source.old.ROOT/"results"/"pacs_schemeC_Study54_Heavy_combinations_v1"
COL={"1":"#2D75B8","2":"#D17D38","3":"#29947F"}
CHUNK=2500


def cv_linear_by_protein(log2_matrix):
    v=np.exp2(log2_matrix)
    return 100*np.std(v,axis=1,ddof=1)/np.mean(v,axis=1)


def measured_for_matrix(h,ind,center=True):
    z=np.mean(h[np.array(ind)],axis=0)
    return z-np.median(z) if center else z


def batches(indices,kind):
    if kind==1:
        for x in indices:
            yield ((x,),)
    else:
        it=itertools.combinations(indices,kind)
        while True:
            b=list(itertools.islice(it,CHUNK))
            if not b:break
            yield b


def evaluate(pred_matrix,targets):
    # pred_matrix is N x Study54, already median-centered per combination.
    out={}
    for name,y in targets.items():
        y=y-np.median(y)
        dif=pred_matrix-y[None,:]
        rmse=np.sqrt(np.mean(dif*dif,axis=1))
        a=pred_matrix-pred_matrix.mean(axis=1,keepdims=True)
        b=y-y.mean()
        norm=np.sqrt(np.sum(a*a,axis=1)*np.sum(b*b))
        r=np.divide(np.sum(a*b[None,:],axis=1),norm,
                    out=np.full(len(norm),np.nan),where=norm>1e-12)
        yr=np.std(y,ddof=1)
        sd=np.std(pred_matrix,axis=1,ddof=1)/yr
        out[name]=(r,rmse,sd)
    return out


def search(seqs,allowed,H,targets):
    n=len(allowed)
    outcsv=OUT/"01_all_equal_weight_combinations_1to3.csv"
    rowcounter=0
    first=True
    for k in (1,2,3):
        for tup in batches(allowed,k):
            ind=np.array(tup,dtype=int)
            # shape nC x Study54; geometric mean in linear Heavy terms.
            p=np.mean(H[ind],axis=1)
            p=p-np.median(p,axis=1,keepdims=True)
            stats=evaluate(p,targets)
            rows=[]
            for i,t in enumerate(tup):
                row={"n_HELP":k,
                     "HELP_1":seqs[t[0]],
                     "HELP_2":seqs[t[1]] if k>=2 else "",
                     "HELP_3":seqs[t[2]] if k>=3 else ""}
                for name in targets:
                    r,rmse,ratio=stats[name]
                    row.update({name+"_Pearson":float(r[i]),
                                name+"_RMSE":float(rmse[i]),
                                name+"_SD_ratio":float(ratio[i])})
                rows.append(row)
            pd.DataFrame(rows).to_csv(outcsv,index=False,
                                      mode="w" if first else "a",header=first)
            first=False
            rowcounter+=len(rows)
        print("SEARCHED",k,rowcounter,flush=True)
    assert rowcounter==n+math.comb(n,2)+math.comb(n,3)
    return pd.read_csv(outcsv,keep_default_na=False)


def top_summaries(allrows,seqs,study,qc,Hstudy,Hqc):
    best=[]
    for pgkind in ("ownPG","fixedPG"):
        for size in (1,2,3):
            d=allrows[allrows.n_HELP==size]
            for goal,ascending in (("RMSE",True),("Pearson",False)):
                sortcol=pgkind+"_"+goal
                rank=d.sort_values(sortcol,ascending=ascending,kind="stable")
                rank.head(50).to_csv(
                    OUT/f"02_{pgkind}_top50_{size}HELP_{goal}.csv",index=False)
                best.append({"PG_target":pgkind,"size":size,
                             "ranking":goal,**rank.iloc[0].to_dict()})
    selected=pd.DataFrame(best)
    qcpg=qc["PG"][qc["fullPG"]]
    raw_cv=cv_linear_by_protein(qcpg)
    raw_median_cv=float(np.median(raw_cv))
    qc_orig=np.median(qc["ownPG"])
    for i,r in selected.iterrows():
        key=[r.HELP_1,r.HELP_2,r.HELP_3]
        idx=[seqs.index(v) for v in key if isinstance(v,str) and v]
        sig=measured_for_matrix(Hqc,idx)
        # Shift original pooled-QC common-PG log2 by heavy signed deviations.
        # Centering correction only changes an overall constant, not QC CV.
        corrected_cv=cv_linear_by_protein(qcpg-sig[None,:])
        target=qc["ownPG"]-np.median(qc["ownPG"])
        fixed=qc["fixedPG"]-np.median(qc["fixedPG"])
        sigc=sig-np.median(sig)
        selected.loc[i,"QC_before_median_proteinCV"]=raw_median_cv
        selected.loc[i,"QC_after_median_proteinCV"]=float(np.median(corrected_cv))
        selected.loc[i,"QC_delta_median_proteinCV"]=float(np.median(corrected_cv)-raw_median_cv)
        selected.loc[i,"QC_PGown_Pearson"]=float(np.corrcoef(sigc,target)[0,1])
        selected.loc[i,"QC_PGown_RMSE"]=float(np.sqrt(np.mean((sigc-target)**2)))
        selected.loc[i,"QC_PGfixed_Pearson"]=float(np.corrcoef(sigc,fixed)[0,1])
        selected.loc[i,"QC_PGfixed_RMSE"]=float(np.sqrt(np.mean((sigc-fixed)**2)))
    selected.to_csv(OUT/"03_winning_combinations_with_QC_diagnostics.csv",index=False)
    return selected


def make_plots(winners,seqs,study,qc,Hstudy,Hqc,names,groups):
    plots=[]
    plate=np.array(study["plate"],dtype=str)
    sample=[names[i] for i in study["index"]]
    x=np.arange(1,55)
    for targetname in ("ownPG","fixedPG"):
        tar=study[targetname]-np.median(study[targetname])
        fig,axes=plt.subplots(3,1,figsize=(23,15),sharex=True)
        fig.subplots_adjust(left=.075,right=.985,top=.88,bottom=.19,hspace=.41)
        for ax,size in zip(axes,(1,2,3)):
            r=winners[(winners.PG_target==targetname)&
                (winners["size"]==size)&(winners.ranking=="RMSE")].iloc[0]
            inds=[seqs.index(z) for z in (r.HELP_1,r.HELP_2,r.HELP_3) if z]
            fitted=measured_for_matrix(Hstudy,inds)
            desc=" + ".join(seqs[z] for z in inds)
            for p in ("1","2","3"):
                sel=plate==p
                ax.plot(x[sel],tar[sel],color=COL[p],marker="o",lw=1.6,
                        ms=4.3,linestyle="-")
                ax.plot(x[sel],fitted[sel],color=COL[p],marker="D",lw=1.4,
                        ms=3.6,linestyle="--")
            for k in range(1,54):
                if plate[k]!=plate[k-1]:ax.axvline(k+.5,color="#899AA8",ls=":",lw=1.2)
            ax.axhline(0,color="#9BAABA",lw=.75)
            ax.grid(axis="y",alpha=.15)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.set_ylabel("Signed log2 deviation")
            ax.set_title(f"{size} HELP: {desc}\nPearson {r[targetname+'_Pearson']:+.4f}"
                         f" | RMSE {r[targetname+'_RMSE']:.4f}"
                         f" | HELP/PG SD {r[targetname+'_SD_ratio']:.2f}",
                         loc="left",fontsize=10.8,pad=9)
        axes[-1].set_xticks(x,sample,rotation=90,fontsize=6.8)
        axes[-1].set_xlabel("Study54 samples sorted by plate -> sample ID")
        handles=[Line2D([0],[0],color=COL[p],lw=3,label="Plate "+p)
                 for p in ("1","2","3")]
        handles += [Line2D([0],[0],color="#526A7A",marker="o",lw=1.5,
                           label="PG sample median"),
                    Line2D([0],[0],color="#526A7A",marker="D",ls="--",lw=1.5,
                           label="HELP mean deviation")]
        fig.legend(handles=handles,ncol=5,loc="lower center",
                   bbox_to_anchor=(.5,.015),frameon=False)
        fig.suptitle(f"Scheme C: signed sample-median deviations | {targetname} | Study54",
                     fontsize=14,fontweight="bold")
        png=f"04_Study54_{targetname}_best1_2_3HELP_plate_colored.png"
        fig.savefig(OUT/png,dpi=180,bbox_inches="tight",facecolor="white")
        plt.close(fig)
        assert (OUT/png).stat().st_size>25000
        plots.append(png)
    return plots


def main():
    if OUT.exists():raise FileExistsError(str(OUT))
    _,names,info,pglog,seqs,H,groups,plates,masks=core.base_data()
    study=core.calc_group_data(pglog,H,seqs,names,info,groups,plates,masks["Study54"])
    qc=core.calc_group_data(pglog,H,seqs,names,info,groups,plates,masks["pooledQC15"])
    assert study["n"]==54 and qc["n"]==15
    assert int(study["fullPG"].sum())==3965
    assert int(qc["fullPG"].sum())==5227
    keep=masks["Study54"]|masks["pooledQC15"]
    assert keep.sum()==69 and np.isfinite(H[:,keep]).all()
    # Compute original raw linear QC CV; no data reduction beyond prior rule.
    qcCV=np.array([100*np.std(np.exp2(qc["H"][i]),ddof=1)/
                   np.mean(np.exp2(qc["H"][i])) for i in range(len(seqs))])
    allowed=np.where(qcCV<10)[0].tolist()
    assert len(allowed)==68 and len(seqs)==97
    Hs=study["H"]; Hq=qc["H"]
    targets={"ownPG":study["ownPG"],"fixedPG":study["fixedPG"]}
    OUT.mkdir(parents=True)
    rows=search(seqs,allowed,Hs,targets)
    expected=68+math.comb(68,2)+math.comb(68,3)
    assert len(rows)==expected
    winners=top_summaries(rows,seqs,study,qc,Hs,Hq)
    plots=make_plots(winners,seqs,study,qc,Hs,Hq,names,groups)
    qc_baseline=float(np.median(cv_linear_by_protein(qc["PG"][qc["fullPG"]])))
    summary={
        "n_Study54":54,"n_HC":27,"n_S":27,"n_pooledQC":15,
        "n_excluded_InternalQC":6,"CV_QC15_threshold_pct":10,
        "candidate_count":len(allowed),"candidate_names":[seqs[j] for j in allowed],
        "n_tested_combinations":expected,"combos_by_size":{"1":68,"2":math.comb(68,2),
                                      "3":math.comb(68,3)},
        "objective":"Minimize RMSE of signed per-sample log2 PG median deviations versus "
                    "equal weight mean log2 Heavy deviation, after independently median-centering "
                    "final combination on Study54; no amplitude rescaling.",
        "two_targets":"ownPG median of each sample's detected PG; fixedPG median of 3965 common Study PG",
        "primary_target":"ownPG",
        "pooled_QC_role":"QC linear Heavy CV gate and QC PG CV before/after exploratory diagnostic only",
        "QC_fixed_commonPG_median_CV_raw":qc_baseline,
        "training_validation_caveat":"Candidate search fit to Study raw PG sample medians; flattening "
                     "those same Study medians is circular and NOT a validation of batch-effect removal.",
        "plots":plots,
    }
    compact=[]
    for _,r in winners.iterrows():
        compact.append({"target":r.PG_target,"n_HELP":int(r["size"]),
                        "criterion":r.ranking,
                        "HELP":[z for z in (r.HELP_1,r.HELP_2,r.HELP_3) if z],
                        "study_Pearson":float(r[r.PG_target+"_Pearson"]),
                        "study_RMSE":float(r[r.PG_target+"_RMSE"]),
                        "QC_before_CV":float(r.QC_before_median_proteinCV),
                        "QC_after_CV":float(r.QC_after_median_proteinCV),
                        "QC_PGown_Pearson":float(r.QC_PGown_Pearson),
                        "QC_PGown_RMSE":float(r.QC_PGown_RMSE)})
    summary["winners"]=compact
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n")
    text=[
      "# Scheme C v1: sample signed-median deviations matched by HELP combinations",
      "",
      "Study54 (HC27+S27) only. Pooled QC15 only for CV filter and QC diagnostic.",
      "Internal-QC excluded. No original matrix changed.",
      "",
      "## Equations",
      "",
      "y[i] = median(log2 original detected-PG signals in Study sample i) - "
      "median_over_54(median log2 PG per Study sample).",
      "",
      "x[j,i] = log2 Heavy[j,i] - median_over_54(log2 Heavy[j,:]).",
      "",
      "predicted[i] = mean(x[j,i] for j in selected HELP), then subtract "
      "median_over_54(predicted).",
      "",
      "RMSE = sqrt(mean_over_54((predicted[i]-y[i])^2)), unscaled. "
      "Pearson and Heavy/PG SD ratio also reported.",
      "",
      "No negative weight, no free scaling, no intercept fitting. "
      "Equivalent to equal-weight combination in log2 (linear-scale geometric mean).",
      "",
      "## Search",
      "",
      f"- 97/97 detected 69/69, {len(allowed)} pass pooled QC15 raw Heavy linear CV<10%.",
      f"- Enumerated {math.comb(68,2)} pairs and {math.comb(68,3)} triples, "
      f"plus 68 singles, total {expected}.",
      "- Two independent PG targets: own-detected PG medians (main), "
      "fixed 3965 common Study PG (sensitivity).",
      "",
      "## Best by minimum Study RMSE",
      "",
      "| PG target | HELP count | HELP peptides | Study Pearson | Study RMSE (log2) | QC median PG CV raw -> after |",
      "|---|---:|---|---:|---:|---|"
    ]
    for x in compact:
        if x["criterion"]!="RMSE":continue
        text.append(f"| {x['target']} | {x['n_HELP']} | {' + '.join(x['HELP'])} | "
                    f"{x['study_Pearson']:+.4f} | {x['study_RMSE']:.4f} | "
                    f"{x['QC_before_CV']:.2f}% -> {x['QC_after_CV']:.2f}% |")
    text += ["","## Graphics",""]
    for png in plots:text+=["!["+png+"]("+png+")",""]
    text += [
      "## QC diagnostic definition",
      "",
      "For each winning combination, use its mean log2 Heavy signal within 15 QC, "
      "median-center across QC, subtract from all positive PG in the fixed "
      "5227-QC-common set; report median per-protein linear CV before/after.",
      "This is a diagnostic, not proof; optimization only targets Study PG medians.",
      "",
      "## Cautions",
      "",
      "- Perfectly flattening training PG sample medians does not demonstrate correction; "
      "the Study PG medians include true biological heterogeneity.",
      "- The own-detected PG target can be confounded by changing detected PG counts.",
      "- All 52k combination scores are in 01_all_equal_weight_combinations_1to3.csv.",
      "- Only three plates; leave-one-plate-out and external projects remain essential.",
      "- For independent technical benefit evaluate QC protein CV, plate dependence, "
      "and HC/S biological preservation; do not optimize PCA grouping alone.",
      "",
      "## Files",
      "",
      "- 01_all_equal_weight_combinations_1to3.csv",
      "- 02_*_top50_* .csv (independent rankings, each PG definition and size)",
      "- 03_winning_combinations_with_QC_diagnostics.csv",
      "- summary.json"
    ]
    (OUT/"README.md").write_text("\n".join(text)+"\n")
    print("COMBINATION_SUMMARY",json.dumps({
        "count":expected,"qc_raw":qc_baseline,"winners":compact,
        "plots":plots},ensure_ascii=False),flush=True)
    print("PASS complete",flush=True)

if __name__=="__main__":
    main()
