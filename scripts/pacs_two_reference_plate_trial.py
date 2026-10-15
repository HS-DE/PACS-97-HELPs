#!/usr/bin/env python3
"""PACS exploratory normalization with P00751/P07225 or P07225 alone.

10 methods: [two references / PROS1 only] x
 [a: plate HELP ratio, b: plate HELP ratio + sample HELP ratio,
  b2: sample HELP ratio only, c: plate PG protein, d: sample PG protein].

All factor calculations and matrices in log2. CV back-transforms to linear
units (the only statistically meaningful protein-intensity CV definition).
Input: 6675 PG groups x 75 samples, ratio_LH sheet, two original metadata.
Original input files and historical results are never modified.
"""
from __future__ import annotations
import hashlib, json, math, re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from filter_samples import readbook,tab
from yeast_protein_qc import ROOT, GROUPS, COLORS, metadata, excel_rows
from pacs_common import ANN

PG=ROOT/"pg_matrix.tsv"
HELP=ROOT/"LH_ratio_tables_2026.10.13—ALL-HELP.xlsx"
META_QC=ROOT/"all_QC_metadata.xlsx"
META_STUDY=ROOT/"all_sample_metadata.xlsx"
OUT=ROOT/"results"/"two_reference_plate_normalization_v1"
PLATES=("1","2","3")
ALL_REF={"STGSWSTLK":"P00751","VYFAGFPR":"P07225"}
SETS={
    "two_refs":("STGSWSTLK","VYFAGFPR"),
    "pros1_only":("VYFAGFPR",),
}
BRANCHES={
    "a_plate_ratio":"Plate-level correction from pooled QC HELP ratios",
    "b_plate_ratio_plus_sample":"Plate HELP ratio then study sample HELP ratio correction",
    "b2_sample_ratio_only":"Study sample HELP ratio correction directly from Raw",
    "c_plate_protein":"Plate-level correction from pooled QC PG protein intensities",
    "d_sample_protein_only":"Study sample PG protein correction directly from Raw",
}
STUDY=("HC","S")
CORRECT_PLATE=("HC","S","pooled QC")
PAL=dict(zip(GROUPS,COLORS))
MARK={"STGSWSTLK":"o","VYFAGFPR":"D"}
REFCOL={"STGSWSTLK":"#287BBB","VYFAGFPR":"#D47733"}
STAGE_TITLES={"raw":"Raw","corrected":"Corrected"}

def digest(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for v in iter(lambda:f.read(1<<20),b""):h.update(v)
    return h.hexdigest()

def valid_log2(a):
    x=np.asarray(a,dtype=float)
    z=np.full(x.shape,np.nan,dtype=float)
    ok=np.isfinite(x)&(x>0)
    z[ok]=np.log2(x[ok])
    return z

def get_samples():
    labels=metadata()
    info={}
    qc_ids=set()
    for file,expected in ((META_QC,{"Internal-QC","pooled QC"}),(META_STUDY,{"HC","S"})):
        for row in excel_rows(file):
            sid=row.get("sample_id","").strip()
            if not sid or sid not in labels:continue
            group=labels[sid]
            if group not in expected:raise ValueError(f"Wrong source metadata: {sid}")
            plate=str(row.get("Plate","")).strip()
            run=str(row.get("Run","")).strip()
            if plate not in PLATES or not run:
                raise ValueError(f"Invalid Plate or Run: {sid}, {plate}, {run}")
            if sid in info:raise ValueError(f"Duplicate sample: {sid}")
            info[sid]={"sample_id":sid,"group":group,"plate":plate,"run":run}
            if file==META_QC and group=="pooled QC":qc_ids.add(sid)
    if set(info)!=set(labels) or len(info)!=75:
        raise ValueError("All 75 sample IDs must match both metadata files")
    for p in PLATES:
        if sum(x["group"]=="pooled QC" and x["plate"]==p for x in info.values())!=5:
            raise ValueError(f"Plate {p} does not have exactly 5 pooled QCs")
    counts={g:sum(x["group"]==g for x in info.values()) for g in GROUPS}
    if counts!={"Internal-QC":6,"pooled QC":15,"HC":27,"S":27}:
        raise ValueError(str(counts))
    return info

def load():
    info=get_samples()
    pg=pd.read_csv(PG,sep="\t",encoding="utf-8-sig",dtype=str,keep_default_na=False,low_memory=False)
    names=pg.columns[5:].tolist()
    if list(pg.columns[:5])!=ANN or len(pg)!=6675 or len(names)!=75 or set(names)!=set(info):
        raise ValueError("Unexpected PG schema/count/mismatching 75 sample IDs")
    raw=pg[names].apply(pd.to_numeric,errors="coerce").to_numpy(dtype=float)
    logged=valid_log2(raw)
    seq_rows={}
    ss,sheets=readbook(HELP)
    if len(sheets)!=3:raise ValueError("Expected L/H/ratio_LH sheets")
    entries=tab(sheets[2],ss)
    headers={val:col for col,val in entries[0].items() if col!="A"}
    if len(headers)!=75 or set(headers)!=set(names):
        raise ValueError("HELP ratio sample IDs mismatch")
    byseq={r.get("A",""):r for r in entries[1:]}
    if len(byseq)!=97:raise ValueError("Expected 97 original HELP sequences")
    for seq in ALL_REF:
        if seq not in byseq:raise ValueError(f"Missing requested HELP {seq}")
        seq_rows[seq]=np.array([float(byseq[seq].get(headers[s],"nan")) for s in names])
        if not np.all(np.isfinite(seq_rows[seq])&(seq_rows[seq]>0)):
            raise ValueError(f"HELP {seq} ratio is missing for at least one sample")
    ratio={seq:np.log2(values) for seq,values in seq_rows.items()}
    refs={}
    for seq,accession in ALL_REF.items():
        match=[]
        for ix,r in pg.iterrows():
            tokens=set(re.split(r"[;,|\s]+",r["Protein.Group"]))|set(re.split(r"[;,|\s]+",r["Protein.Ids"]))
            if accession in tokens:match.append(ix)
        if len(match)!=1:raise ValueError(f"Expected unique PG row for {seq}/{accession}, got {len(match)}")
        row=match[0]
        if not np.all(np.isfinite(logged[row])):
            raise ValueError(f"Reference PG protein {accession} has missing/nonpositive measurements")
        refs[seq]={"id":accession,"row":row,"names":str(pg.iloc[row]["Protein.Names"]),
                   "genes":str(pg.iloc[row]["Genes"])}
    return pg,names,info,logged,ratio,refs

def indexed(names,info,predicate):
    return np.array([j for j,sid in enumerate(names) if predicate(info[sid])],dtype=int)

def build_factors(branch,seqs,names,info,rawlog,ratios,refs):
    n=len(names)
    idx_study=indexed(names,info,lambda x:x["group"] in STUDY)
    idx_plate=indexed(names,info,lambda x:x["group"] in CORRECT_PLATE)
    ratio_array=np.vstack([ratios[s] for s in seqs])
    protein_array=np.vstack([rawlog[refs[s]["row"],:] for s in seqs])
    f1=np.zeros(n)
    f2=np.zeros(n)
    plate_values={}
    reference_plate_global=None
    sample_reference_global=None
    if branch in ("a_plate_ratio","b_plate_ratio_plus_sample","c_plate_protein"):
        src=ratio_array if branch!="c_plate_protein" else protein_array
        for p in PLATES:
            q=indexed(names,info,lambda x:x["group"]=="pooled QC" and x["plate"]==p)
            if len(q)!=5:raise ValueError(f"Invalid QC count on Plate {p}")
            plate_values[p]=float(np.median(src[:,q].ravel()))
        reference_plate_global=float(np.median(list(plate_values.values())))
        for j in idx_plate:
            f1[j]=reference_plate_global-plate_values[info[names[j]]["plate"]]
    if branch=="b_plate_ratio_plus_sample":
        post=ratio_array+f1[None,:]
        study_medians=np.median(post[:,idx_study],axis=0)
        sample_reference_global=float(np.median(study_medians))
        f2[idx_study]=sample_reference_global-study_medians
    elif branch=="b2_sample_ratio_only":
        study_medians=np.median(ratio_array[:,idx_study],axis=0)
        sample_reference_global=float(np.median(study_medians))
        f2[idx_study]=sample_reference_global-study_medians
    elif branch=="d_sample_protein_only":
        study_medians=np.median(protein_array[:,idx_study],axis=0)
        sample_reference_global=float(np.median(study_medians))
        f2[idx_study]=sample_reference_global-study_medians
    # Internal-QC always has zero offsets.
    lm=indexed(names,info,lambda x:x["group"]=="Internal-QC")
    assert np.all(f1[lm]==0) and np.all(f2[lm]==0)
    # pooled QC receives plate correction but never sample-level factor.
    qc=indexed(names,info,lambda x:x["group"]=="pooled QC")
    assert np.all(f2[qc]==0)
    total=f1+f2
    corrected=rawlog+total[None,:]
    # Only plate-ratio methods modify HELP ratios. Second-step f2 is applied
    # to PG proteins exclusively, NOT to the corrected HELP ratio.
    adjusted_ratio=ratio_array+(f1[None,:] if branch in ("a_plate_ratio","b_plate_ratio_plus_sample") else 0)
    checks=[{
        "sample_id":sid,"group":info[sid]["group"],"plate":info[sid]["plate"],
        "run":info[sid]["run"],
        "raw_selected_HELP_log2_median":float(np.median(ratio_array[:,j])),
        "after_plate_selected_HELP_log2_median":float(np.median(adjusted_ratio[:,j])),
        "raw_selected_PG_log2_median":float(np.median(protein_array[:,j])),
        "f1_plate_log2":float(f1[j]),"f2_sample_log2":float(f2[j]),
        "total_factor_log2":float(total[j]),
        "corrected_selected_PG_log2_median":float(np.median(protein_array[:,j]+total[j])),
        **{f"ratio_raw_log2_{s}":float(ratio_array[k,j]) for k,s in enumerate(seqs)},
        **{f"ratio_after_plate_log2_{s}":float(adjusted_ratio[k,j]) for k,s in enumerate(seqs)},
    } for j,sid in enumerate(names)]
    return corrected,adjusted_ratio,checks,{
      "plate_reference_per_plate":plate_values,"plate_reference_global":reference_plate_global,
      "sample_reference_global":sample_reference_global,
      "n_plate_qc_per_plate":5,"n_study_sample_factors":len(idx_study),
      "ref_kind":("HELP_LH_ratio_log2" if branch.startswith(("a_","b_")) else
                  "PG_protein_log2"),
      "sample_factor_source":("HELP_LH_ratio_log2" if branch in ("b_plate_ratio_plus_sample","b2_sample_ratio_only") else
                              "PG_protein_log2" if branch=="d_sample_protein_only" else None),
    }

def cv_by_group(mat,names,info):
    # CV in original linear units, after undoing log2 solely for CV evaluation.
    results=[]
    for g in GROUPS:
        ix=indexed(names,info,lambda x:x["group"]==g)
        a=mat[:,ix]
        x=np.exp2(np.clip(a,-100,100))
        x[~np.isfinite(a)]=np.nan
        mask=np.isfinite(x)&(x>0)
        counts=mask.sum(axis=1)
        required=max(3,math.ceil(.7*len(ix)))
        cleaned=np.where(mask,x,0.0)
        means=np.divide(cleaned.sum(axis=1),counts,
                        out=np.full(len(a),np.nan),where=counts>0)
        deviations=np.where(mask,cleaned-means[:,None],0.0)
        sd=np.sqrt(np.divide(np.square(deviations).sum(axis=1),counts-1,
                    out=np.full(len(a),np.nan),where=counts>1))
        cv=100*sd/means
        cv[(counts<required)|(means<1e-8)|(~np.isfinite(cv))]=np.nan
        good=cv[np.isfinite(cv)]
        if len(good)==0:raise ValueError(f"No CV values for {g}")
        results.append((g,cv,{"group":g,"sample_n":len(ix),
            "n_proteins_with_cv":int(len(good)),"median_cv_pct":float(np.median(good)),
            "q1_cv_pct":float(np.percentile(good,25)),"q3_cv_pct":float(np.percentile(good,75))}))
    return results

def pca(mat,names,info,groups,common_feature_mask=None):
    ix=indexed(names,info,lambda r:r["group"] in groups)
    sub=mat[:,ix]
    if common_feature_mask is None:
        keep=np.sum(np.isfinite(sub),axis=1)>=math.ceil(.7*len(ix))
        keep&=np.nanstd(sub,axis=1)>1e-8
    else:keep=common_feature_mask
    work=sub[keep].copy()
    med=np.nanmedian(work,axis=1)
    work=np.where(np.isfinite(work),work,med[:,None])
    # SVD in sample space (~54 or 75) rather than across >6000 proteins.
    center=work.T-work.T.mean(axis=0)
    gram=center@center.T
    values,vec=np.linalg.eigh(gram)
    order=np.argsort(values)[::-1]
    values=np.maximum(values[order],0)
    vec=vec[:,order]
    coords=vec[:,:2]*np.sqrt(values[:2])[None,:]
    explained=values[:2]/max(values.sum(),1e-12)*100
    return ix,coords,explained,int(keep.sum()),keep

def scatter_pca(ax,mat,names,info,groups,keep,title):
    ix,z,ex,npro,mask=pca(mat,names,info,groups,keep)
    for g in groups:
        j=[k for k,i in enumerate(ix) if info[names[i]]["group"]==g]
        points=z[j,:]
        if len(points):
            ax.scatter(points[:,0],points[:,1],s=30 if len(ix)>60 else 46,
                       color=PAL[g],label=f"{g} (n={len(points)})",alpha=.83,
                       edgecolor="white",lw=.45)
    ax.set_xlabel(f"PC1 ({ex[0]:.1f}%)");ax.set_ylabel(f"PC2 ({ex[1]:.1f}%)")
    ax.set_title(title,fontsize=11,fontweight="semibold")
    ax.axhline(0,c="#C2CAD2",lw=.6);ax.axvline(0,c="#C2CAD2",lw=.6)
    ax.legend(frameon=False,fontsize=8,loc="best")
    return {"n_samples":len(ix),"n_proteins":npro,"PC1_pct":float(ex[0]),"PC2_pct":float(ex[1])}

def save(fig,path):
    fig.savefig(path,dpi=145,bbox_inches="tight",facecolor="white")
    plt.close(fig)

def beautify(ax):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.grid(axis="y",alpha=.16)
    ax.set_axisbelow(True)

def plot_sample_boxes(raw,corrected,names,info,refs,seqs,path):
    # Grouped samples. 75 boxes per panel, highlight actual P00751/P07225.
    ordered=sorted(range(len(names)),key=lambda j:(GROUPS.index(info[names[j]]["group"]),
                                                     int(info[names[j]]["plate"]),names[j]))
    fig,axs=plt.subplots(2,1,figsize=(22,13),sharex=True,layout="constrained")
    for ax,mat,title in zip(axs,(raw,corrected),("Raw log2 PG","Corrected log2 PG")):
        vectors=[mat[:,j][np.isfinite(mat[:,j])] for j in ordered]
        boxes=ax.boxplot(vectors,positions=np.arange(1,len(names)+1),widths=.63,
                         showfliers=False,patch_artist=True,manage_ticks=False)
        for b,j in zip(boxes["boxes"],ordered):
            b.set_facecolor(PAL[info[names[j]]["group"]]);b.set_alpha(.58)
        for seq in seqs:
            vals=mat[refs[seq]["row"],ordered]
            ax.scatter(np.arange(1,len(names)+1),vals,
                       s=25,color=REFCOL[seq],marker=MARK[seq],zorder=7,
                       label=f'{seq} / {refs[seq]["id"]}',edgecolor="white",lw=.3)
        ax.set_ylabel("log2 PG intensity")
        ax.set_title(title,loc="left",fontweight="bold")
        ax.legend(frameon=False,loc="upper left",bbox_to_anchor=(1.0,1),fontsize=9)
        beautify(ax)
    axs[-1].set_xticks(np.arange(1,len(names)+1),
                       [names[j] for j in ordered],rotation=90,fontsize=6)
    for ax in axs:
        prev=None
        for k,j in enumerate(ordered,1):
            g=info[names[j]]["group"]
            if prev is not None and g!=prev:ax.axvline(k-.5,color="#677B8D",lw=1.1)
            prev=g
    fig.suptitle("Per-sample protein boxplots · selected reference PG proteins highlighted",fontsize=16)
    save(fig,path)

def plot_cv(rawcv,newcv,path):
    fig,axes=plt.subplots(1,2,figsize=(15,6),sharey=True,layout="constrained")
    for ax,detail,title in zip(axes,(rawcv,newcv),("Raw","Corrected")):
        vectors=[r[1][np.isfinite(r[1])] for r in detail]
        violin=ax.violinplot(vectors,positions=range(1,5),showextrema=False,showmedians=False)
        for part,g in zip(violin["bodies"],GROUPS):
            part.set_facecolor(PAL[g]);part.set_alpha(.28)
        boxes=ax.boxplot(vectors,positions=range(1,5),widths=.19,
                         showfliers=False,patch_artist=True,manage_ticks=False)
        for b,g in zip(boxes["boxes"],GROUPS):
            b.set_facecolor(PAL[g]);b.set_alpha(.42)
        for k,(g,cv,stat) in enumerate(detail,1):
            ax.text(k,float(np.percentile(vectors[k-1],99.5))+6,
                    f'{stat["median_cv_pct"]:.2f}%',ha="center",fontsize=9,fontweight="bold")
        ax.set_ylim(0,max(100,min(250,max(np.nanpercentile(v,99.5) for v in vectors)+30)))
        ax.set_xticks(range(1,5),[g+"\n(n="+str(r[2]["sample_n"])+")" for g,r in zip(GROUPS,detail)],fontsize=9)
        ax.set_title(title,loc="left",fontweight="bold")
        ax.set_ylabel("Protein CV (%) on linear intensity")
        beautify(ax)
    fig.suptitle("Group-level per-protein CV · raw vs corrected",fontsize=15,fontweight="bold")
    save(fig,path)

def plot_pca(raw,corrected,names,info,path):
    fig,axes=plt.subplots(2,2,figsize=(14,11),layout="constrained")
    record=[]
    for row,groups in enumerate((["HC","S"],GROUPS)):
        ix,z,ex,n,keep=pca(raw,names,info,groups)
        for col,mat in enumerate((raw,corrected)):
            title=("Raw" if col==0 else "Corrected")+" · "+("HC / S" if row==0 else "All groups")
            record.append({"view":"study_only" if row==0 else "all_groups",
                           "stage":"raw" if col==0 else "corrected",
                           **scatter_pca(axes[row,col],mat,names,info,groups,keep,title)})
    fig.suptitle("PCA of log2 PG intensity · matched feature set per view",fontsize=16,fontweight="bold")
    save(fig,path)
    return record

def plot_factor(factors,seqs,branch,path):
    df=pd.DataFrame(factors)
    order=df.sort_values(["group","plate","sample_id"],kind="stable")
    fig,axes=plt.subplots(2,1,figsize=(15,8),layout="constrained")
    for ax,k in zip(axes,("f1_plate_log2","f2_sample_log2")):
        for g in GROUPS:
            sub=order[order.group==g]
            ax.scatter(np.arange(len(sub))+1+{"Internal-QC":0,"pooled QC":6,"HC":21,"S":48}[g],
                       sub[k],s=34,color=PAL[g],label=f"{g} (n={len(sub)})",alpha=.8)
        ax.axhline(0,c="#718093",lw=.8)
        ax.set_ylabel(k);beautify(ax)
    axes[0].legend(frameon=False,ncol=4,loc="upper center",bbox_to_anchor=(.5,1.16))
    axes[-1].set_xlabel("Samples sorted by group and Plate")
    fig.suptitle(f"Factors · {branch} · {len(seqs)} reference(s)",fontsize=15)
    save(fig,path)

def plot_plate_qc(refs,seqs,raw,new,names,info,path):
    qindex=indexed(names,info,lambda x:x["group"]=="pooled QC")
    fig,axs=plt.subplots(1,2,figsize=(13,5.8),layout="constrained")
    # Each QC sample PG global median; annotate each Plate's pooled QCs.
    for ax,mat,title in zip(axs,(raw,new),("Raw","Corrected")):
        for p in PLATES:
            ids=[i for i in qindex if info[names[i]]["plate"]==p]
            vals=np.nanmedian(mat[:,ids],axis=0)
            x=[int(p)]*len(vals)
            ax.scatter(np.array(x)+np.linspace(-.12,.12,len(vals)),vals,s=46,
                       color=PAL["pooled QC"],alpha=.8,edgecolor="white",lw=.35)
            ax.plot([int(p)-.22,int(p)+.22],[float(np.median(vals))]*2,
                    lw=2.2,color="#D07B2F")
        ax.set_xticks((1,2,3),[f"Plate {p}\n5 pooled QCs" for p in PLATES])
        ax.set_ylabel("Per-QC sample median log2 PG intensity")
        ax.set_title(title,loc="left",fontweight="bold")
        beautify(ax)
    fig.suptitle("Pooled QC plate-level PG global medians (descriptive check)",fontsize=14)
    save(fig,path)

def evaluate_case(group,branch,seqs,pg,names,info,raw,ratios,refs,baseline):
    folder=OUT/group/branch
    folder.mkdir(parents=True,exist_ok=False)
    log,adjratio,factors,notes=build_factors(branch,seqs,names,info,raw,ratios,refs)
    oldraw=raw[:,indexed(names,info,lambda x:x["group"]=="Internal-QC")]
    newlm=log[:,indexed(names,info,lambda x:x["group"]=="Internal-QC")]
    if not np.array_equal(oldraw,newlm,equal_nan=True):
        raise AssertionError("Never modify Internal-QC")
    deltas=np.array([f["total_factor_log2"] for f in factors])
    assert np.allclose(log,raw+deltas[None,:],atol=1e-10,equal_nan=True)
    output=pg[ANN].copy()
    for j,s in enumerate(names):output[s]=log[:,j]
    output.to_csv(folder/"corrected_pg_log2.tsv",sep="\t",index=False,na_rep="")
    pd.DataFrame(factors).to_csv(folder/"sample_factors.csv",index=False)
    ratios_out=pd.DataFrame({"sample_id":names,"group":[info[s]["group"] for s in names],
                             "plate":[info[s]["plate"] for s in names]})
    for j,s in enumerate(seqs):
        ratios_out["raw_log2_ratio_"+s]=ratios[s]
        ratios_out["post_plate_log2_ratio_"+s]=adjratio[j]
    ratios_out.to_csv(folder/"selected_HELP_ratio_log2.csv",index=False)
    details=cv_by_group(log,names,info)
    summary=[]
    cv_rows=[]
    for (g,cvs,r) in details:
        summary.append({"reference_set":group,"branch":branch,**r})
        for i,value in enumerate(cvs):
            cv_rows.append({"protein_row":i+1,"protein_accessions":pg.iloc[i]["Protein.Ids"],
                            "group":g,"cv_pct":value})
    pd.DataFrame(cv_rows).to_csv(folder/"protein_group_cv.csv",index=False)
    pd.DataFrame(summary).to_csv(folder/"group_cv_summary.csv",index=False)
    plot_sample_boxes(raw,log,names,info,refs,seqs,
                      folder/"01_sample_boxplots_targets_raw_vs_corrected.png")
    plot_cv(baseline,details,folder/"02_cv_violin_box_raw_vs_corrected.png")
    pca_summary=plot_pca(raw,log,names,info,folder/"03_pca_raw_vs_corrected.png")
    plot_factor(factors,seqs,branch,folder/"04_factor_offsets.png")
    plot_plate_qc(refs,seqs,raw,log,names,info,folder/"05_pooled_QC_plate_medians.png")
    pd.DataFrame(pca_summary).to_csv(folder/"pca_summary.csv",index=False)
    (folder/"factor_reference.json").write_text(
        json.dumps({"selected_sequences":seqs,"associated_proteins":[refs[s]["id"] for s in seqs],
                    "branch":branch,"description":BRANCHES[branch],**notes},
                   ensure_ascii=False,indent=2),encoding="utf-8")
    return summary,notes

def make_overview(summaries,rawsummary):
    df=pd.DataFrame(summaries)
    baseline=pd.DataFrame([x[2] for x in rawsummary])
    rawcv=baseline.set_index("group")["median_cv_pct"].to_dict()
    df["raw_median_cv_pct"]=df["group"].map(rawcv)
    df["delta_cv_percentage_points"]=df["median_cv_pct"]-df["raw_median_cv_pct"]
    df.to_csv(OUT/"all_methods_group_cv.csv",index=False)
    fig,ax=plt.subplots(figsize=(15,7),layout="constrained")
    stage_labels=[f"{g}/{b.split('_')[0]}" for g in SETS for b in BRANCHES]
    positions=np.arange(len(stage_labels))
    width=.28
    for shift,g in ((-.5,"HC"),(.5,"S")):
        y=[float(df[(df.reference_set==group)&(df.branch==method)&
                    (df.group==g)].median_cv_pct.iloc[0])
           for group in SETS for method in BRANCHES]
        ax.bar(positions+shift*width,y,width=.26,color=PAL[g],label=g,alpha=.83)
    ax.axhline(rawcv["HC"],c=PAL["HC"],ls="--",lw=1,label=f'HC Raw {rawcv["HC"]:.2f}%')
    ax.axhline(rawcv["S"],c=PAL["S"],ls="--",lw=1,label=f'S Raw {rawcv["S"]:.2f}%')
    ax.set_xticks(positions,stage_labels,rotation=35,ha="right",fontsize=9)
    ax.set_ylabel("Median per-protein CV (%) across Study samples")
    ax.set_title("Ten normalization branches · HC and S group CV compared with Raw",loc="left",fontweight="bold")
    ax.legend(frameon=False,ncol=4)
    beautify(ax)
    save(fig,OUT/"00_all_methods_study_cv_comparison.png")
    return df,rawcv

def report(df,rawcv,reference_log,pg,names,info,refs):
    counts=pd.Series([info[s]["plate"] for s in names]).value_counts().to_dict()
    matrix=pd.crosstab(pd.Series([info[s]["plate"] for s in names],name="Plate"),
                       pd.Series([info[s]["group"] for s in names],name="group")).to_dict()
    sha={p.name:digest(p) for p in (PG,HELP,META_QC,META_STUDY)}
    manifest={
        "input_sha256":sha,"protein_rows":len(pg),"samples":len(names),
        "plate_sample_counts":counts,"plate_group_counts":matrix,
        "source_refs":refs,"reference_sets":SETS,"branches":BRANCHES,
        "log2_for_factors_and_outputs":True,
        "CV":"100 * sample SD(linear intensities, ddof=1) / arithmetic mean(linear intensities)",
        "plate_QC_aggregation":"median of 5 QC x number-of-references log2 measurements within Plate; global median of 3 Plate medians",
        "study_sample_aggregation":"median of 1 or 2 log2 reference measurements within each Study sample; global median across 54 Study sample medians",
        "plate_affected_groups":["HC","S","pooled QC"],
        "sample_factor_affected_groups":["HC","S"],
        "Internal-QC":"Unchanged except universal output conversion from linear to log2",
        "b2_interpretation":"Direct Raw + sample HELP ratio offset, not plate-corrected source",
        "ratio_after_b":"Only f1 changes HELP ratios; f2 applied to protein matrix only",
        "per_target_independent_HC_S_biology":"Not protected by this exploratory model",
        "no_imputation":True
    }
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    summary=ref_summary=pg[ANN].copy()
    lines=[
      "# Two HELP references / two PG source proteins: exploratory normalization trial","",
      "**10 completed methods:** a, b, b2, c, d with either STGSWSTLK+VYFAGFPR",
      "(P00751+P07225) or only VYFAGFPR (P07225). No yeast or other 95 HELP peptides",
      "are used for factor estimation. Uncorrected original PG and ratio data preserved.","",
      "## Rules used in this reproducible experiment","",
      "- All input intensities and HELP L/H ratios are log2-transformed first;",
      "  nonpositive/missing source values are not imputed.",
      "- 3 Plates; 5 pooled QC per Plate; 54 Study (HC+S); Internal-QC 6.",
      "- Plate factor: within each plate median of the 5 pooled QC × 1 or 2",
      "  selected reference measurements, then median of the three plate medians;",
      "  f1_plate = median_of_plates - plate_median.",
      "- Plate factor is added to HC/S and pooled QC protein intensities.",
      "  For ratio-based a/b, it is ALSO added to the selected HELP ratios.",
      "  Other 95 HELP ratios are not changed.",
      "- For b, study factor is calculated from the stage-a adjusted HELP ratios;",
      "  for b2, the very same sample-level operation starts from raw HELP",
      "  ratios and raw PG intensity, WITHOUT f1 or plate corrections.",
      "- For d, study factor is calculated directly from raw PG protein intensities.",
      "- Study factor median_s = median of 54 study-sample reference medians;",
      "  f2_s = median_s - individual sample median; applied to HC/S PG only.",
      "- Internal-QC is **never corrected**; pooled QC receives ONLY f1",
      "  in plate-correction branches, and no sample-level f2.",
      "- For b, f2 adjusts PG intensities only, NOT the HELP ratios.",
      "- Every stage uses log2 shifts. CV is evaluated on the",
      "  **back-transformed linear corrected intensities** (standard CV).","",
      "## Reference PG mapping","",
      "| Peptide | UniProt accession | PG protein name | Genes | Data row |",
      "| --- | --- | --- | --- | ---: |"
    ]
    for seq,ref in refs.items():
        lines.append(f'| {seq} | {ref["id"]} | {ref["names"]} | {ref["genes"]} | {ref["row"]+1} |')
    lines += ["","## Ten methods and CV summary (%)","",
              f'Uncorrected HC median CV {rawcv["HC"]:.2f}%; S {rawcv["S"]:.2f}%.',
              "",
              "| Reference | Branch | HC median CV | S median CV | HC Δ vs Raw (pp) | S Δ vs Raw (pp) |",
              "| --- | --- | ---: | ---: | ---: | ---: |"]
    for group in SETS:
        for branch in BRANCHES:
            c={}
            for g in STUDY:
                r=df[(df.reference_set==group)&(df.branch==branch)&(df.group==g)].iloc[0]
                c[g]=r
            lines.append(f'| {group} | [{branch}]({group}/{branch}/) | '
                         f'{c["HC"].median_cv_pct:.2f}% | {c["S"].median_cv_pct:.2f}% | '
                         f'{c["HC"].delta_cv_percentage_points:+.2f} | '
                         f'{c["S"].delta_cv_percentage_points:+.2f} |')
    lines += ["","![Ten-method comparison](00_all_methods_study_cv_comparison.png)","",
      "## Files per method","",
      "Each method folder contains: corrected_pg_log2.tsv (6675 × 75),",
      "sample_factors.csv (75 records with f1 and f2), selected_HELP_ratio_log2.csv",
      "(selected peptides only, raw and after-plate values), protein_group_cv.csv,",
      "group_cv_summary.csv, pca_summary.csv, factor_reference.json, and",
      "five PNG figures: sample protein boxplots with both targets highlighted,",
      "per-protein CV violin/box comparison, PCA Raw vs Corrected (all/Study),",
      "sample factor offsets and pooled-QC plate median checks.","",
      "## Interpretation cautions","",
      "- Biological HC/S variation contributes to group CV, so lower/higher group",
      "  CV is not conclusive evidence of removing/preserving batch effects.",
      "- These two endogenous proteins may differ biologically; their raw variation",
      "  is not assumed to be purely technical. This is exploratory as requested.",
      "- Pooled-QC plate correction is not independent QC validation because the",
      "  same 15 pooled QC samples were used to derive the plate offsets.",
      "- This script does not select an optimal method or assert correction improves data.",
      "- Results produced independently of past methods; no originals or historical",
      "  result directories were altered.","",
      "## Audit metadata","",
      "See run_manifest.json for exact source hashes, reference use and definitions.",
    ]
    (OUT/"README.md").write_text("\n".join(lines)+"\n",encoding="utf-8")

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError(f"Existing result folder protected: {OUT}")
    pg,names,info,raw,ratios,refs=load()
    OUT.mkdir(parents=True,exist_ok=True)
    baseline=cv_by_group(raw,names,info)
    summaries=[];references=[]
    for group,seqs in SETS.items():
        for branch in BRANCHES:
            print("RUNNING",group,branch,flush=True)
            cv,meta=evaluate_case(group,branch,seqs,pg,names,info,raw,ratios,refs,baseline)
            summaries.extend(cv)
            references.append({"reference_set":group,"branch":branch,**meta})
    df,rawcv=make_overview(summaries,baseline)
    pd.DataFrame(references).to_json(OUT/"all_method_factor_references.json",
                                     orient="records",indent=2)
    report(df,rawcv,references,pg,names,info,refs)
    print("PASS: 10 completed combinations, 50 method figures, raw vs corrected CV/PCA; Internal-QC invariant",flush=True)

if __name__=="__main__":
    main()
