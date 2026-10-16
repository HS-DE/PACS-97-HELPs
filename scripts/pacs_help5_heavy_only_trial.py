#!/usr/bin/env python3
"""Five HELP Heavy-only normalization compared with matched five-HELP L/H.

Q_j = median[log2(Heavy_j)] across 15 pooled QC. For sample i,
log2 factor = mean_j(Q_j - log2(Heavy_ij)), i.e. geometric mean
of QC/clinical Heavy ratios. Matched L/H control uses the same formula.
Correct HC/S and pooled QC, retain original Internal-QC. A QC
leave-one-out sensitivity check excludes each QC from its own reference.
"""
from pathlib import Path
import hashlib,json,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from filter_samples import readbook,tab
from yeast_protein_qc import ROOT,GROUPS,COLORS
import pacs_two_reference_plate_trial as old

HELP=ROOT/"LH_ratio_tables_2026.10.13—ALL-HELP.xlsx"
OUT=ROOT/"results"/"help5_heavy_only_normalization"
SEQS=("FSVVYAK","STGSWSTLK","VTGWGNLK","VVEESELAR","VYFAGFPR")
UNIPROTS=("P02765","P00751","P00734","P02748","P07225")
GENES=("AHSG","CFB","F2","C9","PROS1")
COLORS_P=("#C84456","#2875AD","#8E63B5","#248D96","#BC7D2F")
LAYERS=("L","H","ratio_LH")
METHODS=("Raw","L/H control","Heavy only")

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def sheet_data(samples):
    ss,sheets=readbook(HELP)
    if len(sheets)!=3:raise ValueError("Expect L/H/ratio_LH sheets")
    result={}
    for label,sheet in zip(LAYERS,sheets):
        rows=tab(sheet,ss)
        header={v:k for k,v in rows[0].items() if k!="A"}
        byseq={r.get("A",""):r for r in rows[1:]}
        if len(rows)!=98 or len(byseq)!=97 or set(header)!=set(samples):
            raise ValueError(f"Invalid source layout for {label}")
        matrix=[]
        for seq in SEQS:
            if seq not in byseq:raise ValueError(f"Missing {seq}")
            v=np.array([float(byseq[seq].get(header[s],"nan")) for s in samples])
            if np.any(~np.isfinite(v)) or np.any(v<=0):
                raise ValueError(f"Missing or invalid {label} {seq}")
            matrix.append(np.log2(v))
        result[label]=np.array(matrix)
    error=result["L"]-result["H"]-result["ratio_LH"]
    if np.max(np.abs(error))>1e-10:raise ValueError("Source L/H and ratio inconsistent")
    return result,dict(zip(SEQS,np.max(np.abs(error),axis=1).tolist()))

def target_rows(pg):
    output={}
    for seq,accession,gene in zip(SEQS,UNIPROTS,GENES):
        matches=[]
        for i,r in pg.iterrows():
            tok=set(re.split(r"[;,|\s]+",r["Protein.Group"]))|set(re.split(r"[;,|\s]+",r["Protein.Ids"]))
            if accession in tok:matches.append(i)
        if len(matches)!=1:raise ValueError(f"{gene} mapped to {len(matches)} protein groups")
        output[seq]={"row":int(matches[0]),"accession":accession,"gene":gene}
    return output

def factors(signals,infos,samples):
    q=np.array([i for i,s in enumerate(samples) if infos[s]["group"]=="pooled QC"],int)
    study=np.array([i for i,s in enumerate(samples) if infos[s]["group"] in ("HC","S")],int)
    lm=np.array([i for i,s in enumerate(samples) if infos[s]["group"]=="Internal-QC"],int)
    assert (len(q),len(study),len(lm))==(15,54,6)
    out={};records=[];refs=[]
    for label in ("H","ratio_LH"):
        val=signals[label]
        qref=np.median(val[:,q],axis=1)
        per_peptide=qref[:,None]-val
        f=per_peptide.mean(axis=0)
        f[lm]=0
        qc_loo=np.full(len(samples),np.nan)
        for idx in q:
            peers=q[q!=idx]
            qc_loo[idx]=(np.median(val[:,peers],axis=1)-val[:,idx]).mean()
        out[label]={"factor":f,"qc_loo":qc_loo}
        for i,seq in enumerate(SEQS):
            refs.append({"source":label,"sequence":seq,
                         "pooled_QC_median_log2":float(qref[i])})
        for k,s in enumerate(samples):
            r={"source":label,"sample_id":s,"group":infos[s]["group"],
               "plate":infos[s]["plate"],"run":infos[s]["run"],
               "factor_log2":float(f[k]),"multiplicative_factor":float(2**f[k]),
               "qc_loo_factor_log2":qc_loo[k],"n_help":5}
            for j,seq in enumerate(SEQS):
                r["ref_log2_"+seq]=float(qref[j])
                r["signal_log2_"+seq]=float(val[j,k])
                r["individual_factor_log2_"+seq]=float(per_peptide[j,k])
            records.append(r)
    return out,records,refs,q,study,lm

def addf(raw,f):return raw+np.asarray(f)[None,:]

def savefig(fig,name):
    fig.savefig(OUT/name,dpi=180,bbox_inches="tight",facecolor="white")
    plt.close(fig)

def beautify(ax):
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    ax.grid(axis="y",alpha=.17);ax.set_axisbelow(True)

def sample_plot(stages,samples,info,targets):
    order=sorted(range(75),key=lambda i:(GROUPS.index(info[samples[i]]["group"]),
                                         int(info[samples[i]]["plate"]),samples[i]))
    fig,axes=plt.subplots(3,1,figsize=(23,17),sharex=True,layout="constrained")
    for ax,(method,x) in zip(axes,stages.items()):
        b=ax.boxplot([x[:,i][np.isfinite(x[:,i])] for i in order],
                     showfliers=False,positions=range(1,76),widths=.6,
                     patch_artist=True,manage_ticks=False)
        for box,i in zip(b["boxes"],order):
            box.set_facecolor(COLORS[GROUPS.index(info[samples[i]]["group"])])
            box.set_alpha(.48)
        for c,seq in zip(COLORS_P,SEQS):
            idx=targets[seq]["row"]
            ax.scatter(range(1,76),x[idx,order],s=25,color=c,label=seq,
                       edgecolor="white",linewidth=.3,zorder=5)
        for j in range(1,75):
            a=info[samples[order[j-1]]];z=info[samples[order[j]]]
            if a["group"]!=z["group"]:ax.axvline(j+.5,color="#77818D",lw=1.1)
            elif a["plate"]!=z["plate"]:ax.axvline(j+.5,color="#9AABB6",lw=.6,ls=":")
        ax.set_title(method,loc="left",fontweight="bold",fontsize=13)
        ax.set_ylabel("log2 protein intensity")
        ax.legend(frameon=False,loc="upper left",bbox_to_anchor=(1,1),ncol=1,fontsize=9)
        beautify(ax)
    axes[-1].set_xticks(range(1,76),[samples[i] for i in order],rotation=90,fontsize=6)
    fig.suptitle("Raw, L/H control and five-HELP Heavy-only · five source PG proteins highlighted",fontsize=15,fontweight="bold")
    savefig(fig,"01_all_sample_pg_boxes_highlight_five.png")

def cvplot(stats):
    fig,axes=plt.subplots(1,3,figsize=(16,6.6),sharey=True,layout="constrained")
    cvdat=dict(stats)
    for ax,name in zip(axes,METHODS):
        curves=cvdat[name]
        vec=[x[1][np.isfinite(x[1])] for x in curves]
        p=ax.violinplot(vec,showextrema=False)
        for art,c in zip(p["bodies"],COLORS):
            art.set_facecolor(c);art.set_alpha(.23)
        bx=ax.boxplot(vec,widths=.22,showfliers=False,patch_artist=True)
        for art,c in zip(bx["boxes"],COLORS):
            art.set_facecolor(c);art.set_alpha(.42)
        for n,(_,_,s) in enumerate(curves,1):
            ax.scatter(n,s["median_cv_pct"],s=35,color="#293744",zorder=6)
            ax.text(n,s["median_cv_pct"]+3,f'{s["median_cv_pct"]:.2f}%',
                    ha="center",fontsize=8.5)
        ax.set_xticks(range(1,5),["Internal-QC","pooled QC","HC","S"],
                      rotation=20,fontsize=9)
        ax.set_ylim(0,190)
        ax.set_title(name,loc="left",fontweight="bold")
        beautify(ax)
    axes[0].set_ylabel("Protein CV (%) · evaluated on linear intensities")
    fig.suptitle("Group protein CV distributions",fontsize=15,fontweight="bold")
    savefig(fig,"02_group_cv_violin_box.png")

def pcaplot(stages,samples,info):
    fig,axs=plt.subplots(2,3,figsize=(17,11),layout="constrained")
    rows=[]
    for row,groups in enumerate((["HC","S"],GROUPS)):
        _,_,_,_,mask=old.pca(stages["Raw"],samples,info,groups)
        for col,(name,matrix) in enumerate(stages.items()):
            r=old.scatter_pca(axs[row,col],matrix,samples,info,groups,mask,name+" "+("/".join(groups)))
            rows.append({"method":name,"view":"study" if row==0 else "all",**r})
    fig.suptitle("PCA · matched protein features for raw and corrected matrices",fontsize=15,fontweight="bold")
    savefig(fig,"03_pca_raw_ratio_heavy.png")
    return rows

def factorplot(records,samples,info):
    d=pd.DataFrame(records)
    fig,axs=plt.subplots(2,1,figsize=(19,9),sharex=True,layout="constrained")
    order=sorted(samples,key=lambda s:(GROUPS.index(info[s]["group"]),int(info[s]["plate"]),s))
    for ax,layer in zip(axs,("ratio_LH","H")):
        frame=d[d.source==layer].set_index("sample_id").loc[order]
        for g,col in zip(GROUPS,COLORS):
            ix=np.array([i for i,s in enumerate(order) if info[s]["group"]==g])
            ax.scatter(ix+1,frame.iloc[ix].factor_log2,s=30,color=col,alpha=.8,label=g)
        ax.axhline(0,color="#707C88",lw=.8)
        ax.set_ylabel("Log2 correction factor")
        ax.set_title("L/H control" if layer=="ratio_LH" else "Heavy only",
                     loc="left",fontweight="bold")
        beautify(ax)
    axs[0].legend(frameon=False,ncol=4)
    axs[-1].set_xticks(range(1,76),order,rotation=90,fontsize=6)
    fig.suptitle("Five-peptide mean log2 QC-to-sample signal offsets",fontsize=15,fontweight="bold")
    savefig(fig,"04_sample_factors_heavy_vs_ratio.png")

def within_plate_qc(log,info,samples,stage):
    result=[]
    for plate in ("1","2","3"):
        sub=log[:,[j for j,s in enumerate(samples) if info[s]["group"]=="pooled QC"
                 and info[s]["plate"]==plate]]
        if sub.shape[1]!=5:raise ValueError("Plate lacks five pooled QCs")
        linear=np.exp2(np.clip(sub,-100,100))
        linear[~np.isfinite(sub)]=np.nan
        with np.errstate(all="ignore"):
            n=np.isfinite(linear).sum(axis=1)
            mean=np.nanmean(linear,axis=1)
            sd=np.nanstd(linear,axis=1,ddof=1)
            cv=100*sd/mean
        cv[(n<4)|(mean<1e-8)|(~np.isfinite(cv))]=np.nan
        v=cv[np.isfinite(cv)]
        result.append({"method":stage,"plate":plate,"n_QC":5,
                       "n_protein_CV":len(v),"median_cv_pct":float(np.median(v))})
    return result

def qcplot(plate):
    df=pd.DataFrame(plate)
    fig,ax=plt.subplots(figsize=(10,5.6),layout="constrained")
    for name,col in zip(("Raw","L/H control","Heavy only","Heavy only QC LOO"),
                        ("#8792A1","#D29339","#259477","#8254A0")):
        s=df[df.method==name]
        ax.plot([int(x) for x in s.plate],s.median_cv_pct,marker="o",
                color=col,label=name,lw=1.6)
    ax.set_xticks([1,2,3],["Plate 1","Plate 2","Plate 3"])
    ax.set_ylabel("Within-Plate pooled-QC protein median CV (%)")
    ax.legend(frameon=False)
    ax.set_title("5 pooled QC per plate · independent-of-self LOO sensitivity",loc="left",
                 fontweight="bold")
    beautify(ax);savefig(fig,"05_pooled_qc_plate_cv.png")

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError("Existing results folder protected")
    pg,samples,info,raw,ratios,refs=old.load()
    source,arithmetic=sheet_data(samples)
    target=target_rows(pg)
    fc,rec,ref,QC,STUDY,LM=factors(source,info,samples)
    ratio=addf(raw,fc["ratio_LH"]["factor"])
    heavy=addf(raw,fc["H"]["factor"])
    loo=fc["H"]["factor"].copy();loo[QC]=fc["H"]["qc_loo"][QC]
    heavyloo=addf(raw,loo)
    assert np.array_equal(heavy[:,LM],raw[:,LM],equal_nan=True)
    assert np.array_equal(ratio[:,LM],raw[:,LM],equal_nan=True)
    stages={"Raw":raw,"L/H control":ratio,"Heavy only":heavy}
    OUT.mkdir(parents=True,exist_ok=False)
    cvstats={name:old.cv_by_group(mat,samples,info) for name,mat in stages.items()}
    cvstats["Heavy only QC LOO"]=old.cv_by_group(heavyloo,samples,info)
    summary=[]
    each=[]
    for name,r in cvstats.items():
        for g,cv,v in r:
            summary.append({"method":name,**v})
            if name!="Heavy only QC LOO":
                each.extend({"method":name,"group":g,"protein_row":i+1,
                             "cv_pct":v0} for i,v0 in enumerate(cv))
    pd.DataFrame(rec).to_csv(OUT/"sample_factors.csv",index=False)
    pd.DataFrame(ref).to_csv(OUT/"reference_five_peptide_QC_medians.csv",index=False)
    pd.DataFrame(summary).to_csv(OUT/"group_cv_summary.csv",index=False)
    pd.DataFrame(each).to_csv(OUT/"per_protein_cv.csv",index=False)
    qctable=[]
    for name,mat in (*stages.items(),("Heavy only QC LOO",heavyloo)):
        qctable+=within_plate_qc(mat,info,samples,name)
    pd.DataFrame(qctable).to_csv(OUT/"pooled_QC_plate_cv.csv",index=False)
    for name,mat,fname in (("Heavy only",heavy,"corrected_pg_log2_heavy_only.tsv"),
                          ("L/H control",ratio,"corrected_pg_log2_LH_control.tsv")):
        out=pg[old.ANN].copy()
        for j,sid in enumerate(samples):out[sid]=mat[:,j]
        out.to_csv(OUT/fname,sep="\t",index=False,na_rep="")
    sample_plot(stages,samples,info,target)
    cvplot(list(cvstats.items())[:3])
    pd.DataFrame(pcaplot(stages,samples,info)).to_csv(OUT/"pca_summary.csv",index=False)
    factorplot(rec,samples,info)
    qcplot(qctable)
    manifest={
      "sources":{p.name:sha(p) for p in (old.PG,HELP,old.META_QC,old.META_STUDY)},
      "proteins":len(pg),"samples":75,
      "five_peptides":list(SEQS),"source_proteins":target,
      "QC_reference_sample_count":15,"QC_LOO_others_count":14,
      "f_i":"arithmetic mean of median pooled-QC log2 Heavy minus sample log2 Heavy, over five peptides",
      "ratio_control":"same formula with L/H ratio values",
      "corrected_protein":"raw log2 PG + f_i for HC/S and pooled QC, Internal-QC unchanged",
      "CV":"sample sd (ddof=1)/mean on corrected linear-scale intensity",
      "log2_before_median":True,"no_imputation":True,
      "L_minus_H_minus_log2_ratio_max_abs":arithmetic
    }
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding="utf-8")
    groups=("HC","S","pooled QC","Internal-QC")
    lines=["# PACS Heavy-only five-HELP normalization","",
       "Raw Heavy intensity from 2026.10.13 workbook, not L/H.",
       "Reference: for each of the 5 peptides, median log2(Heavy) across 15 pooled QC.",
       "A sample log2 factor is the mean over all five individual differences",
       "between that peptide QC median and the sample log2 Heavy signal.",
       "Thus the original-scale multiplier is the geometric mean of",
       "the five H_QC/H_sample values.",
       "Matched L/H uses same formula only as a direct control.",
       "HC/S and pooled QC get their factor, while 6 Internal-QC are unchanged.",
       "No previous results or source workbooks are touched.","",
       "## Median per-protein CV (%) on corrected linear scale","",
       "| Method | HC | S | pooled QC | Internal-QC |",
       "| --- | ---: | ---: | ---: | ---: |"]
    for name in cvstats:
        d={r["group"]:r["median_cv_pct"] for r in summary if r["method"]==name}
        lines.append("| "+name+" | "+" | ".join(f"{d[g]:.2f}" for g in groups)+" |")
    lines.extend(["","## Visualizations","",
        "![75-sample boxplots](01_all_sample_pg_boxes_highlight_five.png)","",
        "![Group CV](02_group_cv_violin_box.png)","",
        "![PCA](03_pca_raw_ratio_heavy.png)","",
        "![Factors](04_sample_factors_heavy_vs_ratio.png)","",
        "![Within-Plate QC](05_pooled_qc_plate_cv.png)","",
        "## Outputs and caveats","",
        "sample_factors.csv stores 5 individual log2 contributions, 75 samples",
        "for each of Heavy and matched L/H control. Two complete 6675 x 75",
        "corrected log2 PG matrices, per-protein CVs, group summary, and PCA",
        "results are also available. Raw data unchanged.",
        "Leave-one-out pooled QC diagnostics exclude the evaluated QC sample",
        "from its own per-peptide median reference, but the QC pool still",
        "determines the reference population; not an independent held-out QC.",
        "HC/S group CV includes genuine biological variation and cannot",
        "alone prove technical calibration.",
        "Heavy spiked just before loading cannot account for earlier",
        "sample extraction, digestion or recovery differences.",
        "" ])
    (OUT/"README.md").write_text("\n".join(lines),encoding="utf-8")
    print("GROUP_CV",pd.DataFrame(summary)[["method","group","median_cv_pct"]].to_json(orient="records"),flush=True)
    print("PLATE_QC",pd.DataFrame(qctable).to_json(orient="records"),flush=True)
    print("VALIDATED: five Heavy, 75 samples, 3 methods, QC LOO, 5 PNG",flush=True)

if __name__=="__main__":
    main()
