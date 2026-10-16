#!/usr/bin/env python3
"""PACS joint HELP-heavy and yeast spike-in factor decomposition.

Five Heavy peptides added pre-injection, two yeast ENO1 protein-quant rows
from EN01.xlsx, and (secondary) two ENO1 yeast peptides from
sample-QC-InternalQC-HELP_yeast.xlsx.

Reference for EVERY feature: median original log2 intensity in 15 pooled QCs.
HELP h_i: geometric mean of five H_i/H_QC,j, in log2 (mean of 5 diffs).
Yeast protein sum: P00924 and P00925 rows summed on ORIGINAL LINEAR
intensity scale for each sample, then log2 and normalized to QC.
Yeast protein max: larger total row-sum protein as sensitivity analysis.
Yeast peptides: equal-weight mean of 2 peptide-wise log2 deviations, secondary.
Residual log2 r = log2 y - log2 h, after HELP correction.
Two-step HELP+yeast residual factor = -log2 h - log2 r = -log2 y
exactly identical to yeast-only, not an additional independent method.
Only HC/S and pooled QC normalized, Internal-QC untouched in all outputs.
"""
from pathlib import Path
import hashlib, json, math, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from filter_samples import readbook,tab
from yeast_protein_qc import ROOT,GROUPS,COLORS
import pacs_two_reference_plate_trial as base
import pacs_help5_heavy_only_trial as heavy

OUT=ROOT/"results"/"help5_heavy_yeast_decomposition"
PROTEIN=ROOT/"EN01.xlsx"
PEPTIDE=ROOT/"sample-QC-InternalQC-HELP_yeast.xlsx"
YEAST_IDS=("P00924_ENO1_YEAST","P00925_ENO1_YEAST")
YEAST_PEPS=("AADALLLK","VNQIGTLSESIK")
METHODS=("Raw","HELP-heavy only","Yeast protein sum","Yeast protein single","Yeast peptides (2)")
METHOD_COLORS=dict(zip(METHODS,("#8592A1","#3F7CB5","#2C986F","#C0803F","#9B6CB5")))
GRPC=dict(zip(GROUPS,COLORS))

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for t in iter(lambda:f.read(1<<20),b""):h.update(t)
    return h.hexdigest()

def linear_positive(x):
    try:
        v=float(x)
        return v if math.isfinite(v) and v>0 else np.nan
    except (ValueError,TypeError):
        return np.nan

def read_rows(path,expected,samples):
    ss,sheets=readbook(path)
    if len(sheets)!=1:raise ValueError(f"{path.name} should have 1 sheet")
    rows=tab(sheets[0],ss)
    headers={v:c for c,v in rows[0].items() if c!="A"}
    if len(headers)!=75 or set(headers)!=set(samples):
        raise ValueError(f"{path.name} samples mismatch")
    ids=[r.get("A","") for r in rows[1:]]
    if tuple(ids)!=tuple(expected):
        raise ValueError(f"Unexpected feature IDs {ids} in {path.name}")
    vals=np.array([[linear_positive(r.get(headers[s],"")) for s in samples]
                   for r in rows[1:]],float)
    return vals

def audited_inputs(names,info):
    P=read_rows(PROTEIN,YEAST_IDS,names)
    K=read_rows(PEPTIDE,YEAST_PEPS,names)
    study=np.array([i for i,s in enumerate(names) if info[s]["group"] in ("HC","S")])
    qc=np.array([i for i,s in enumerate(names) if info[s]["group"]=="pooled QC"])
    lm=np.array([i for i,s in enumerate(names) if info[s]["group"]=="Internal-QC"])
    assert (len(study),len(qc),len(lm))==(54,15,6)
    # Do not impute missing signals. Only 69 samples are corrected.
    active=np.r_[study,qc]
    for name,arr in ((PROTEIN.name,P),(PEPTIDE.name,K)):
        bad=np.argwhere(~np.isfinite(arr[:,active]))
        if len(bad):
            raise ValueError(f"{name} has {len(bad)} missing/nonpositive study/pooled QC signals: {bad[:15].tolist()}")
    selected=np.argmax(np.nansum(P,axis=1))
    # Main source: sum raw protein intensity, then log2; no double counting
    # of the two original protein-quant rows as independent proteins.
    combined=np.sum(P,axis=0)
    if np.any(~np.isfinite(P),axis=0).any():
        combined[~np.isfinite(P).all(axis=0)]=np.nan
    series={
        "Yeast protein sum":np.log2(combined)[None,:],
        "Yeast protein single":np.log2(P[selected:selected+1]),
        "Yeast peptides (2)":np.log2(K),
    }
    audit={
        "yeast_protein_source_rows":list(YEAST_IDS),
        "protein_row_sums_all75_linear":np.nansum(P,axis=1).tolist(),
        "selected_single_protein_index_zero_based":int(selected),
        "selected_single_protein_source_id":YEAST_IDS[selected],
        "yeast_peptide_source_rows":list(YEAST_PEPS),
        "missing_yeast_protein_by_original_row_all75":[int(np.isnan(a).sum()) for a in P],
        "missing_yeast_peptide_by_source_row_all75":[int(np.isnan(a).sum()) for a in K],
        "protein_sum":"sum both source rows in original linear intensity then log2",
        "note":"Source labels contain P00925_ENO1_YEAST; accession annotation identity is not independently verified"
    }
    return series,audit,study,qc,lm,P,K

def reference_shift(logvalues,qc,active):
    z=np.asarray(logvalues,float)
    reference=np.median(z[:,qc],axis=1)
    if not np.all(np.isfinite(reference)):
        raise ValueError("Invalid pooled QC reference")
    shifts=z-reference[:,None]
    signal=np.nanmean(shifts,axis=0)
    if not np.all(np.isfinite(signal[active])):
        raise ValueError("Invalid QC or study sample signal")
    # Each pooled QC's reference excludes that QC, used only in QC LOO sensitivity.
    loo=np.full(z.shape[1],np.nan)
    for j in qc:
        group=qc[qc!=j]
        assert len(group)==14
        loo[j]=np.mean(z[:,j]-np.median(z[:,group],axis=1))
    return signal,reference,loo

def data_and_factors(names,info):
    pg,sid,details,raw,_,_=base.load()
    if sid!=names or any(details[s]!=info[s] for s in sid):
        raise AssertionError("PG metadata independently reloaded inconsistently")
    sheets=heavy.sheet_data(names)[0]
    yeast,audit,study,qc,lm,protein,peptide=audited_inputs(names,info)
    active=np.r_[study,qc]
    logh,hq,loo_h=reference_shift(sheets["H"],qc,active)
    signals={"HELP-heavy only":logh}
    perref={"HELP-heavy only":hq}
    loo={"HELP-heavy only":loo_h}
    for method,values in yeast.items():
        shift,ref,loo_values=reference_shift(values,qc,active)
        signals[method]=shift
        perref[method]=ref
        loo[method]=loo_values
    offsets={}
    loo_offsets={}
    for method in METHODS:
        candidate=np.zeros(len(names)) if method=="Raw" else -signals[method].copy()
        candidate[lm]=0
        if not np.isfinite(candidate).all():
            raise ValueError(f"Invalid factor applied in {method}")
        offsets[method]=candidate
        loo_adj=candidate.copy()
        if method!="Raw":
            loo_adj[qc]=-loo[method][qc]
        loo_offsets[method]=loo_adj
    residual={}
    for yeastmethod in ("Yeast protein sum","Yeast protein single","Yeast peptides (2)"):
        residual[yeastmethod]=signals[yeastmethod]-logh
        lhs=offsets["HELP-heavy only"]-residual[yeastmethod]
        lhs[lm]=0
        if not np.allclose(lhs,offsets[yeastmethod],rtol=0,atol=1e-12):
            raise AssertionError("Algebraic HELP+yeast residual != yeast-only")
    frames=[]
    for i,sample in enumerate(names):
        row={
          "sample_id":sample,"group":info[sample]["group"],
          "plate":info[sample]["plate"],"run":info[sample]["run"],
          "applies_correction":info[sample]["group"]!="Internal-QC",
          "log2_h_HEAVY_relative_QC":logh[i],
          "HELP_only_factor_log2":offsets["HELP-heavy only"][i],
          "HELP_only_multiplier":np.exp2(offsets["HELP-heavy only"][i]),
          "Yeast_protein_sum_log2_y_relative_QC":signals["Yeast protein sum"][i],
          "Yeast_protein_single_log2_y_relative_QC":signals["Yeast protein single"][i],
          "Yeast_peptide_log2_y_relative_QC":signals["Yeast peptides (2)"][i],
          "sum_protein_residual_log2_r_y_div_h":residual["Yeast protein sum"][i],
          "sum_protein_r_y_div_h":np.exp2(residual["Yeast protein sum"][i]),
          "sum_protein_residual_factor_log2":-residual["Yeast protein sum"][i] if i not in lm else 0,
          "Yeast_sum_only_factor_log2":offsets["Yeast protein sum"][i],
          "Yeast_sum_only_multiplier":np.exp2(offsets["Yeast protein sum"][i]),
          "Yeast_single_only_factor_log2":offsets["Yeast protein single"][i],
          "Yeast_peptide_only_factor_log2":offsets["Yeast peptides (2)"][i],
          "HELP_LOO_factor_log2":loo_offsets["HELP-heavy only"][i] if i in qc else np.nan,
          "Yeast_sum_LOO_factor_log2":loo_offsets["Yeast protein sum"][i] if i in qc else np.nan,
          "Yeast_single_LOO_factor_log2":loo_offsets["Yeast protein single"][i] if i in qc else np.nan,
          "Yeast_peptide_LOO_factor_log2":loo_offsets["Yeast peptides (2)"][i] if i in qc else np.nan,
          **{f"yeast_protein_{x}_linear":protein[k,i] for k,x in enumerate(("P00924","P00925"))},
          "yeast_protein_sum_linear":np.sum(protein[:,i]) if np.isfinite(protein[:,i]).all() else np.nan,
          **{f"yeast_peptide_{s}_linear":peptide[k,i] for k,s in enumerate(YEAST_PEPS)},
        }
        for k,s in enumerate(heavy.SEQS):
            row["log2_HELP_H_"+s]=sheets["H"][k,i]
            row["QC_reference_HELP_H_"+s]=hq[k]
        frames.append(row)
    reference=[]
    for m,vals in perref.items():
        features=list(heavy.SEQS) if m=="HELP-heavy only" else (
           ["sum_P00924_and_P00925"] if m=="Yeast protein sum" else
           [audit["selected_single_protein_source_id"]] if m=="Yeast protein single"
           else list(YEAST_PEPS))
        if len(features)!=len(vals):raise ValueError("Mismatch between yeast feature labels and references")
        for f,v in zip(features,vals):
            reference.append({"method":m,"source_feature":f,"pooled_QC_median_log2":float(v)})
    return pg,raw,offsets,loo_offsets,residual,frames,reference,audit

def chart(fig,name):
    fig.savefig(OUT/name,bbox_inches="tight",dpi=160,facecolor="white")
    plt.close(fig)

def style(ax):
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    ax.grid(axis="y",alpha=.15);ax.set_axisbelow(True)

def ordered(names,info):
    return sorted(range(75),key=lambda i:(GROUPS.index(info[names[i]]["group"]),
                                         int(info[names[i]]["plate"]),names[i]))

def plot_factors(names,info,table):
    ordered_ids=[names[i] for i in ordered(names,info)]
    f=table.set_index("sample_id").loc[ordered_ids]
    fig,axs=plt.subplots(3,1,figsize=(20,12),sharex=True,layout="constrained")
    for ax,(field,title,color) in zip(axs,(
      ("log2_h_HEAVY_relative_QC","HELP Heavy log2 h (relative to QC)","#477BB4"),
      ("Yeast_protein_sum_log2_y_relative_QC","Yeast protein sum log2 y (relative to QC)","#379A6B"),
      ("sum_protein_residual_log2_r_y_div_h","Residual log2 r = log2 y - log2 h","#C78539")
    )):
        ax.scatter(np.arange(1,76),f[field],s=27,color=color,alpha=.85)
        ax.axhline(0,lw=.8,color="#738393")
        ax.set_ylabel(title)
        for k in range(1,75):
            s1,s2=ordered_ids[k-1],ordered_ids[k]
            if info[s1]["group"]!=info[s2]["group"]:
                ax.axvline(k+.5,color="#7A8794",lw=1)
        style(ax)
    axs[-1].set_xticks(range(1,76),ordered_ids,fontsize=6,rotation=90)
    fig.suptitle("HELP / Yeast stage decomposition on 75 original samples",fontsize=16,fontweight="bold")
    chart(fig,"01_help_yeast_residual_per_sample.png")

def plot_scatter(table):
    fig,ax=plt.subplots(figsize=(9,8),layout="constrained")
    groups=GROUPS[1:]
    for group in groups:
        subset=table[table.group==group]
        ax.scatter(subset["log2_h_HEAVY_relative_QC"],
                   subset["Yeast_protein_sum_log2_y_relative_QC"],
                   color=COLORS[GROUPS.index(group)],s=44,
                   label=f"{group} n={len(subset)}",edgecolors="white",linewidths=.4)
    xy=np.r_[table["log2_h_HEAVY_relative_QC"].dropna().to_numpy(),
             table["Yeast_protein_sum_log2_y_relative_QC"].dropna().to_numpy()]
    lim=max(.3,np.max(np.abs(xy[np.isfinite(xy)]))*1.07)
    ax.plot([-lim,lim],[-lim,lim],color="#8492A1",ls="--",lw=1,label="y=h (zero residual)")
    ax.axvline(0,c="#D5DADF",lw=.8);ax.axhline(0,c="#D5DADF",lw=.8)
    ax.set(xlim=(-lim,lim),ylim=(-lim,lim),
           xlabel="HELP Heavy relative log2 h",ylabel="Yeast protein sum relative log2 y")
    ax.set_title("Relative spike-in signals: post-spike vs full workflow",fontsize=14,fontweight="bold")
    ax.legend(frameon=False)
    style(ax);chart(fig,"02_yeast_vs_help_relative_QC_scatter.png")

def plot_cv(summaries):
    table=pd.DataFrame(summaries)
    fig,ax=plt.subplots(figsize=(15,6.3),layout="constrained")
    names=list(METHODS)
    x=np.arange(len(names))
    w=.19
    for j,g in enumerate(("HC","S","pooled QC")):
        subset=[float(table[(table.method==m)&(table.group==g)].median_cv_pct.iloc[0]) for m in names]
        bars=ax.bar(x+(j-1)*w,subset,width=w*.96,color=GRPC[g],alpha=.9,label=g)
        for b,v in zip(bars,subset):
            ax.text(b.get_x()+b.get_width()/2,v+.5,f"{v:.1f}",fontsize=8,
                    ha="center",va="bottom")
    ax.set_xticks(x,["Raw","HELP Heavy","Yeast sum","Yeast single","Yeast peptides"],fontsize=10)
    ax.set_ylim(0,max(s["median_cv_pct"] for s in summaries)*1.16)
    ax.set_ylabel("Group median protein CV (%)")
    ax.set_title("Raw vs HELP-only vs yeast-only (linear-scale protein CV)",loc="left",
                 fontsize=14,fontweight="bold")
    ax.legend(frameon=False)
    style(ax);chart(fig,"03_group_CV_comparison.png")

def plot_pca(stages,names,info):
    fig,axes=plt.subplots(2,5,figsize=(23,10.8),layout="constrained")
    info_out=[]
    for row,groups in enumerate((["HC","S"],GROUPS)):
        raw_ix,z,ex,n,keep=base.pca(stages["Raw"],names,info,groups)
        for col,(name,mat) in enumerate(stages.items()):
            s=base.scatter_pca(axes[row,col],mat,names,info,groups,keep,
                               name+"\n"+("HC/S" if row==0 else "All"))
            info_out.append({"method":name,"view":"study" if row==0 else "all",**s})
    fig.suptitle("PCA: matched detected protein set, same samples and transforms",fontsize=16,fontweight="bold")
    chart(fig,"04_PCA_methods.png")
    return info_out

def plot_boxes(stages,names,info):
    order=ordered(names,info)
    pos=np.arange(1,76)
    fig,axes=plt.subplots(5,1,figsize=(23,24),sharex=True,layout="constrained")
    for ax,(label,mat) in zip(axes,stages.items()):
        vec=[mat[:,i][np.isfinite(mat[:,i])] for i in order]
        box=ax.boxplot(vec,positions=pos,patch_artist=True,widths=.62,
                       showfliers=False,manage_ticks=False)
        for patch,i in zip(box["boxes"],order):
            patch.set_facecolor(GRPC[info[names[i]]["group"]])
            patch.set_alpha(.49)
        ax.set_title(label,loc="left",fontweight="bold")
        ax.set_ylabel("log2 PG intensity")
        ax.grid(axis="y",alpha=.15);ax.set_axisbelow(True)
        for j in range(1,75):
            x1,x2=info[names[order[j-1]]],info[names[order[j]]]
            if x1["group"]!=x2["group"]:
                ax.axvline(j+.5,lw=1.1,color="#81909A")
        style(ax)
    axes[-1].set_xticks(pos,[names[i] for i in order],fontsize=6,rotation=90)
    fig.suptitle("All 75 original PG samples · comparison of five normalization methods",fontsize=16,fontweight="bold")
    chart(fig,"05_sample_PG_boxplots_all_methods.png")

def qc_per_plate(stages,names,info):
    records=[]
    for method,mat in stages.items():
        for plate in ("1","2","3"):
            idx=[j for j,s in enumerate(names) if info[s]["group"]=="pooled QC" and info[s]["plate"]==plate]
            assert len(idx)==5
            x=mat[:,idx]
            linear=np.exp2(np.clip(x,-100,100))
            linear[~np.isfinite(x)]=np.nan
            with np.errstate(all="ignore"):
                n=np.isfinite(linear).sum(axis=1)
                means=np.nanmean(linear,axis=1)
                stds=np.nanstd(linear,axis=1,ddof=1)
                cv=100*stds/means
            cv[(n<4)|(means<1e-8)]=np.nan
            valid=cv[np.isfinite(cv)]
            records.append({"method":method,"plate":plate,"n_pooled_QC":5,
                            "n_protein_cv":len(valid),
                            "median_cv_pct":float(np.median(valid))})
    return records

def plot_qc_plate(records):
    df=pd.DataFrame(records)
    fig,ax=plt.subplots(figsize=(11,6.5),layout="constrained")
    for name in METHODS+("HELP-heavy only QC LOO","Yeast protein sum QC LOO","Yeast protein single QC LOO","Yeast peptides (2) QC LOO"):
        sub=df[df.method==name]
        if not len(sub):continue
        ax.plot([int(s) for s in sub.plate],sub.median_cv_pct,marker="o",
                lw=1.3,alpha=.82,color=METHOD_COLORS.get(name,"#697B8A"),
                ls="--" if "LOO" in name else "-",
                label=name)
    ax.set_xticks([1,2,3],["Plate 1","Plate 2","Plate 3"])
    ax.set_ylabel("Within-Plate pooled-QC median protein CV (%)")
    ax.set_title("Independent-of-self pooled QC check (LOO vs normal)",loc="left",fontweight="bold")
    ax.legend(frameon=False,fontsize=8,bbox_to_anchor=(1,1),loc="upper left")
    style(ax);chart(fig,"06_pooled_QC_plate_CV_and_LOO.png")

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError(f"Refusing to overwrite existing result directory {OUT}")
    info=base.get_samples()
    names=list(pd.read_csv(base.PG,sep="\t",nrows=0).columns[5:])
    pg,raw,offsets,loo_offsets,residual,records,reference,audit=data_and_factors(names,info)
    stages={name:raw+offsets[name][None,:] for name in METHODS}
    lm=np.array([i for i,s in enumerate(names) if info[s]["group"]=="Internal-QC"])
    if not all(np.array_equal(mat[:,lm],raw[:,lm],equal_nan=True) for mat in stages.values()):
        raise AssertionError("Internal-QC must not be corrected")
    # Ensure h then residual gives same result as direct yeast-only.
    sample_residual=np.asarray([r["sum_protein_residual_log2_r_y_div_h"] for r in records])
    combo=stages["HELP-heavy only"]-sample_residual[None,:]
    combo[:,lm]=raw[:,lm]
    if not np.allclose(combo,stages["Yeast protein sum"],rtol=0,atol=1e-12,equal_nan=True):
        raise AssertionError("Two-step residual must equal yeast-only")
    OUT.mkdir(parents=True,exist_ok=False)
    pd.DataFrame(records).to_csv(OUT/"sample_help_yeast_factors.csv",index=False)
    pd.DataFrame(reference).to_csv(OUT/"QC_reference_per_feature.csv",index=False)
    summary=[]
    detail=[]
    for method,mat in stages.items():
        cvs=base.cv_by_group(mat,names,info)
        for group,cv,stat in cvs:
            summary.append({"method":method,**stat})
            for i,val in enumerate(cv):
                detail.append({"method":method,"group":group,"protein_row_1based":i+1,"cv_pct":val})
        out=pg[base.ANN].copy()
        for j,s in enumerate(names):out[s]=mat[:,j]
        if method!="Raw":
            filename={
              "HELP-heavy only":"corrected_pg_log2_help_only.tsv",
              "Yeast protein sum":"corrected_pg_log2_yeast_sum.tsv",
              "Yeast protein single":"corrected_pg_log2_yeast_single.tsv",
              "Yeast peptides (2)":"corrected_pg_log2_yeast_peptide.tsv",
            }[method]
            out.to_csv(OUT/filename,sep="\t",index=False,na_rep="")
    loo_stages={}
    for name in METHODS[1:]:
        loo_stages[name+" QC LOO"]=raw+loo_offsets[name][None,:]
    loo_summary=[]
    for method,mat in loo_stages.items():
        for group,cv,stat in base.cv_by_group(mat,names,info):
            if group=="pooled QC":
                loo_summary.append({"method":method,**stat})
    pd.DataFrame(loo_summary).to_csv(OUT/"pooled_QC_LOO_group_cv.csv",index=False)
    pd.DataFrame(summary).to_csv(OUT/"all_group_protein_CV.csv",index=False)
    pd.DataFrame(detail).to_csv(OUT/"all_per_protein_CV.csv",index=False)
    plate=qc_per_plate(stages,names,info)
    plate+=qc_per_plate(loo_stages,names,info)
    pd.DataFrame(plate).to_csv(OUT/"pooled_QC_plate_CV.csv",index=False)
    plot_factors(names,info,pd.DataFrame(records))
    plot_scatter(pd.DataFrame(records))
    plot_cv(summary)
    pca_rows=plot_pca(stages,names,info)
    pd.DataFrame(pca_rows).to_csv(OUT/"PCA_summary.csv",index=False)
    plot_boxes(stages,names,info)
    plot_qc_plate(plate)
    manifest={
      "input_sha256":{p.name:sha(p) for p in
        (base.PG,base.HELP,base.META_QC,base.META_STUDY,PROTEIN,PEPTIDE)},
      "peptides_5_help":list(heavy.SEQS),
      "yeast_input_audit":audit,
      "n_pooled_QC_for_reference":15,
      "leave_one_out_n_peer_QC":14,
      "sample_counts":{"HC":27,"S":27,"pooled QC":15,"Internal-QC":6},
      "log2_h":"mean of five individual log2(H_sample/H_median_pooledQC)",
      "log2_y_yeast_peptides":"mean of two peptide-specific log2(Y_sample/Y_median_pooledQC)",
      "log2_y_protein_sum":"log2(sum two original linear protein signals) minus median in 15 pooled QC",
      "log2_r":"log2_y-log2_h",
      "two_step":"-log2_h - (log2_y-log2_h) == -log2_y, exactly same as yeast-only",
      "corrected_groups":["HC","S","pooled QC"],
      "unchanged_group":["Internal-QC"],
      "CV_evaluated_linear":True,
      "PCA_missing_value_median_imputation_only":True,
      "yeast_residual_is_not_proof_of_digestion_efficiency":True,
      "no_original_file_overwrites":True,
    }
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding="utf-8")
    def table_rows(cols,methods):
        lines=[]
        d=pd.DataFrame(summary)
        for method in methods:
            subset=d[d.method==method].set_index("group")
            lines.append("| "+method+" | "+" | ".join(f"{subset.loc[col,'median_cv_pct']:.2f}" for col in cols)+" |")
        return lines
    lines=[
      "# Five-HELP Heavy and Yeast full-protein / peptide spike-in comparison","",
      "**Independent 75-sample, 6675-protein normalization trial** using",
      "raw HELP Heavy from the October 13 file, EN01.xlsx protein signals,",
      "and sample-QC-InternalQC-HELP_yeast.xlsx yeast peptide signals.",
      "No source files or previous results were modified.","",
      "## Main calculations","",
      "For each of five fixed Heavy HELP, reference is its own median log2 intensity",
      "across 15 pooled QCs. HELP log2 relative signal h_i is the equally weighted",
      "mean of five (sample - QC reference) log2 differences.",
      "For yeast proteins, two original rows P00924_ENO1_YEAST and",
      "P00925_ENO1_YEAST are **summed in original linear intensity** per sample",
      "then log2-transformed and compared to the pooled-QC median.",
      "One selected row (larger all-sample linear row sum) and 2 peptide-wise yeast",
      "intensities are secondary robustness variants.",
      "For each method, protein log2 intensity receives negative reference signal.",
      "Pooled QC and HC/S get factors; Internal-QC remains unchanged.","",
      "## Stage decomposition","",
      "log2 h = mean HELP log2(sample/QC) ; HELP-only correction factor = -log2 h.",
      "log2 y = yeast log2(sample/QC) ; yeast-only correction = -log2 y.",
      "Residual after HELP: log2 r = log2 y - log2 h.",
      "Two-step correction HELP + residual = -log2 h - (log2 y - log2 h)",
      "= -log2 y, **mathematically exactly equal to yeast-only**.",
      "Therefore the combined two-stage pathway is reported as factor decomposition,",
      "not as a separate new corrected matrix.","",
      "## Median protein CV by group (%)","",
      "| Method | HC | S | pooled QC | Internal-QC |",
      "| --- | ---: | ---: | ---: | ---: |",
      *table_rows(["HC","S","pooled QC","Internal-QC"],METHODS),"",
      "### QC leave-one-out sensitivity","",
      "Each pooled QC is corrected using only the OTHER 14 pooled QCs",
      "to estimate each reference; not an independent held-out cohort.",
      "See pooled_QC_LOO_group_cv.csv and pooled_QC_plate_CV.csv.","",
      "## Figures","",
      "![HELP h, yeast y and residual r by sample](01_help_yeast_residual_per_sample.png)","",
      "![Yeast vs HELP relative QC](02_yeast_vs_help_relative_QC_scatter.png)","",
      "![HC/S and QC group CV comparisons](03_group_CV_comparison.png)","",
      "![PCA](04_PCA_methods.png)","",
      "![All per-sample PG boxplots](05_sample_PG_boxplots_all_methods.png)","",
      "![Pooled QC plate CV](06_pooled_QC_plate_CV_and_LOO.png)","",
      "## Source audit","",
      "Two original yeast protein rows are used without assuming their accession",
      "or repeated ENO1 protein naming is confirmed. The protein-sum main variant",
      "does not count the two rows as independent proteins.",
      "Check run_manifest.json for missing signals and the selected single-protein row.","",
      "## Interpretation caveats","",
      "Within-HC/S protein CV includes biological variance, not just technical noise.",
      "Yeast residual may reflect sample aliquoting or recovery, not solely digestion.",
      "The fixed-ng peptide aliquot vs fixed-volume step is unspecified;",
      "therefore exact source of pre-HELP variance cannot be identified.",
      "Same pooled QC reference is used for factors and internal QC summaries;",
      "LOO is a sensitivity analysis, not independent verification.",
      "No double multiplication of two unadjusted correction factors is performed.","",
      "## Files","",
      "sample_help_yeast_factors.csv: raw/relative H, yeast, r, corrected log2 factors,",
      "per-peptide H and per-row yeast signals, QC LOO references.",
      "QC_reference_per_feature.csv: per-reference pooled QC log2 median.",
      "Four complete corrected PG matrices in log2 (HELP, yeast sum, yeast single,",
      "yeast peptides), per-protein CV, PCA, batch QC CV tables and source hashes.",
    ]
    (OUT/"README.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print("BATCH_COUNTS",manifest["sample_counts"],flush=True)
    print("AUDIT",json.dumps(audit),flush=True)
    print("CV",pd.DataFrame(summary)[["method","group","median_cv_pct"]].to_json(orient="records"),flush=True)
    print("QC_LOO",pd.DataFrame(loo_summary)[["method","group","median_cv_pct"]].to_json(orient="records"),flush=True)
    print("PASS Six PNG, 4 corrected log2 matrices, 75 sample factors, H+Yeast factor decomposition",flush=True)
if __name__=="__main__":
    main()
