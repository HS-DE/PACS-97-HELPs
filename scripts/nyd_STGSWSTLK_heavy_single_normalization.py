#!/usr/bin/env python3
"""NYD 40 samples: STGSWSTLK single HELP Heavy reference normalization.

Exactly one STGSWSTLK Heavy intensity per sample; therefore the within-sample
median_hi is the log2 of that single intensity. Global median_h is the median
of these 40 log2 intensities; f_i = median_h - median_hi. Add f_i to every
finite-positive, log2-transformed protein-group intensity in that sample.

Input: NYD_data/40samples/{Heavy.xlsx,NYD_pg_matrix2.tsv,NYD_sampleinfo2.xlsx}
Output: independent new subfolder. Original Heavy/PG and previous results
untouched. PCA Raw fit ONCE, Corrected projected on identical Raw basis.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from filter_samples import readbook, tab
from nyd_pg_batch_boxplots import meta, BATCHES, COLORS, ANNOT

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"NYD_data"/"40samples"
PG=BASE/"NYD_pg_matrix2.tsv"
HELP=BASE/"Heavy.xlsx"
META=BASE/"NYD_sampleinfo2.xlsx"
OUT=BASE/"STGSWSTLK_Heavy_single_reference_normalization"
SEQ="STGSWSTLK"
BLUE="#2364E8"
YELLOW="#EABF28"
PROTEINS={"P00751":{"gene":"CFB","color":BLUE,"marker":"D"},
          "P07225":{"gene":"PROS1","color":YELLOW,"marker":"o"}}
BATCH_COUNTS={"T001":16,"T002":8,"T003":8,"T004":8}

def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for p in iter(lambda:f.read(1<<20),b""):h.update(p)
    return h.hexdigest()

def writecsv(filename,rows):
    if not rows:raise ValueError(f"Empty result {filename}")
    with (OUT/filename).open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def val(v):
    try:
        z=float(v)
        return z if math.isfinite(z) and z>0 else np.nan
    except (ValueError,TypeError):return np.nan

def log2valid(x):
    a=np.asarray(x,dtype=float)
    z=np.full(a.shape,np.nan)
    good=np.isfinite(a)&(a>0)
    z[good]=np.log2(a[good])
    return z

def inputs():
    md=meta()
    # Same order as established NYD 40-sample plots.
    records=[r for batch in BATCHES for r in md if r["batch"]==batch]
    names=[r["sample_id"] for r in records]
    assert len(names)==40 and len(set(names))==40
    if {b:sum(r["batch"]==b for r in records) for b in BATCHES}!=BATCH_COUNTS:
        raise ValueError("Unexpected batch counts")
    with PG.open("r",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        columns=reader.fieldnames
        if columns is None or list(columns[:5])!=list(ANNOT):
            raise ValueError("PG annotation fields don't match expected layout")
        if len(columns)!=45 or set(columns[5:])!=set(names):
            raise ValueError("PG sample columns don't match NYD metadata")
        annotations=[]
        intensities=[]
        for row in reader:
            annotations.append([row.get(field,"") for field in ANNOT])
            intensities.append([val(row.get(s)) for s in names])
    if len(annotations)!=4264:
        raise ValueError(f"Expected 4264 PG protein groups, got {len(annotations)}")
    raw=log2valid(intensities)
    if raw.shape!=(4264,40):raise ValueError("Unexpected PG matrix shape")
    pg_targets={}
    for accession,params in PROTEINS.items():
        rows=[]
        for i,anno in enumerate(annotations):
            tokens=set()
            for attr in anno[:2]:
                tokens.update(re.split(r"[;,|\s]+",attr))
            if accession in tokens:rows.append(i)
        if len(rows)!=1:raise ValueError(f"Protein {accession} mapped {len(rows)} times")
        ix=rows[0]
        if params["gene"] not in re.split(r"[;,|\s]+",annotations[ix][3]):
            raise ValueError(f"Protein {accession} gene annotation mismatch")
        if not np.isfinite(raw[ix]).all():
            raise ValueError(f"Protein {accession} missing in one or more samples")
        pg_targets[accession]=ix
    ss,sheets=readbook(HELP)
    if len(sheets)!=1:raise ValueError("Heavy.xlsx must have one sheet")
    rowdata=tab(sheets[0],ss)
    if len(rowdata)!=98:raise ValueError("Expected 97 HELP Heavy rows")
    column_ids={value:key for key,value in rowdata[0].items() if key!="A"}
    if len(column_ids)!=40 or set(column_ids)!=set(names):
        raise ValueError("Heavy sample names must match 40 metadata IDs")
    byseq={r.get("A",""):r for r in rowdata[1:]}
    if len(byseq)!=97 or SEQ not in byseq:
        raise ValueError("Expected exactly 97 unique HELP sequences with STGSWSTLK")
    heavy=np.array([val(byseq[SEQ].get(column_ids[s],"")) for s in names],float)
    if not np.isfinite(heavy).all():
        raise ValueError("STGSWSTLK Heavy must have finite positive values in all 40 samples")
    return records,names,annotations,raw,pg_targets,heavy

def apply_calibration(raw,heavy):
    # log2 FIRST. Only one STGSWSTLK observation per sample.
    median_hi=np.log2(heavy)
    median_h=float(np.median(median_hi))
    f=median_h-median_hi
    after=raw+f[None,:]
    # Aligns the STGSWSTLK Heavy signal theoretically, without modifying Heavy.xlsx.
    if not np.allclose(median_hi+f,median_h,rtol=0,atol=1e-12):
        raise AssertionError("STGSWSTLK Heavy did not align with global median")
    if not np.array_equal(np.isfinite(raw),np.isfinite(after)):
        raise AssertionError("Correction must preserve PG missingness")
    return median_hi,median_h,f,after

def boxplot(ax,matrix,records,names,target_rows,headline,full_x=True):
    pos=np.arange(1,41)
    valid=[x[np.isfinite(x)] for x in matrix.T]
    box=ax.boxplot(valid,positions=pos,widths=.63,
        patch_artist=True,showfliers=False,manage_ticks=False,
        boxprops=dict(color="#6B8191",linewidth=.75),
        whiskerprops=dict(color="#6B8191",linewidth=.75),
        capprops=dict(color="#6B8191",linewidth=.75),
        medianprops=dict(color="#18384C",linewidth=1.7))
    for patch,r in zip(box["boxes"],records):
        patch.set_facecolor(COLORS[r["batch"]]);patch.set_alpha(.44)
    for accession,params in PROTEINS.items():
        xoff=-.15 if accession=="P00751" else .15
        ax.scatter(pos+xoff,matrix[target_rows[accession]],s=49,
            marker=params["marker"],c=params["color"],label=f"{params['gene']} / {accession}",
            edgecolors="#263A45",linewidths=.55,zorder=8)
    start=1
    for batch in BATCHES:
        n=BATCH_COUNTS[batch]
        if start>1:ax.axvline(start-.5,lw=1.1,c="#718699")
        ax.text(start+(n-1)/2,1.017,f"{batch}  (n={n})",
            transform=ax.get_xaxis_transform(),color=COLORS[batch],
            fontsize=10.5,fontweight="bold",ha="center",clip_on=False)
        start+=n
    ax.set_xlim(.35,40.65)
    ax.set_xticks(pos,names if full_x else [""]*40,rotation=90,fontsize=6.7)
    ax.set_title(headline,loc="left",fontweight="bold",fontsize=13,pad=31)
    ax.set_ylabel("log2 protein-group intensity")
    ax.grid(axis="y",alpha=.16)
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)
    ax.legend(frameon=False,loc="upper right",ncol=2,fontsize=9)

def save(fig,name,dpi=195):
    path=OUT/name
    fig.savefig(path,dpi=dpi,facecolor="white",bbox_inches="tight")
    plt.close(fig)
    if path.stat().st_size<15000:raise ValueError("Suspiciously small figure: "+str(path))

def plot_pgs(raw,corrected,records,names,targets):
    fig,ax=plt.subplots(2,1,figsize=(20,14.5),sharex=True,sharey=True,layout="constrained")
    boxplot(ax[0],raw,records,names,targets,"Raw · 4264 PG protein groups",False)
    boxplot(ax[1],corrected,records,names,targets,
            "STGSWSTLK Heavy calibrated · raw log2 protein + sample factor")
    fig.suptitle("NYD 40-sample PG distributions · blue CFB/P00751 · yellow PROS1/P07225",
                 fontweight="bold",fontsize=16)
    save(fig,"01_PG_40samples_Raw_vs_corrected_highlight_CFB_PROS1.png")

def plot_heavy(median_hi,median_h,f,records,names):
    fig,(ax,bx)=plt.subplots(2,1,figsize=(19,9.8),sharex=True,layout="constrained")
    x=np.arange(1,41)
    for batch in BATCHES:
        mask=np.array([r["batch"]==batch for r in records])
        ax.scatter(x[mask],median_hi[mask],c=COLORS[batch],s=58,
                   edgecolors="white",linewidths=.5,label=batch)
        bx.bar(x[mask],f[mask],width=.77,color=COLORS[batch],alpha=.83)
    ax.scatter(x,median_hi+f,c="#223849",marker="_",s=65,zorder=5,
               label="Aligned STGSWSTLK H (diagnostic)")
    ax.axhline(median_h,color=BLUE,lw=1.5,ls="--",alpha=.9,
               label=f"median_h = {median_h:.4f}")
    ax.set_ylabel("log2 STGSWSTLK Heavy")
    ax.set_title("Original Heavy reference and theoretical alignment",loc="left",fontweight="bold")
    ax.legend(frameon=False,ncol=3,fontsize=9)
    bx.axhline(0,color="#394F62",lw=.8)
    bx.set_ylabel("fi = median_h − median_hi (log2)")
    bx.set_title("Applied per-sample additive protein offset",loc="left",fontweight="bold")
    for axis in (ax,bx):
        axis.grid(axis="y",alpha=.15);axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False);axis.spines["right"].set_visible(False)
        for cut in (16.5,24.5,32.5):axis.axvline(cut,c="#A1ABB6",lw=.9)
    bx.set_xticks(x,names,rotation=90,fontsize=7)
    fig.suptitle("STGSWSTLK Heavy-only reference · one measured Heavy value per sample",
                 fontweight="bold",fontsize=15)
    save(fig,"02_STGSWSTLK_Heavy_per_sample_offsets.png")

def stats_cv(raw,corrected,records,annotations):
    cv_by_batch=[]
    per_protein=[]
    arrays={}
    for stage,log in (("Raw",raw),("Corrected",corrected)):
        for batch in BATCHES:
            idx=np.array([j for j,r in enumerate(records) if r["batch"]==batch])
            n=int(len(idx))
            threshold=max(3,math.ceil(.70*n))
            x=log[:,idx]
            # Back-transform to linear BEFORE mean, SD, CV.
            linear=np.exp2(np.clip(x,-100,100))
            linear[~np.isfinite(x)]=np.nan
            present=np.isfinite(linear)&(linear>0)
            counts=np.sum(present,axis=1)
            cleaned=np.where(present,linear,0.)
            means=np.divide(np.sum(cleaned,axis=1),counts,
                            out=np.full(raw.shape[0],np.nan),where=counts>0)
            deviations=np.where(present,linear-means[:,None],0.0)
            sd=np.sqrt(np.divide(np.sum(deviations**2,axis=1),counts-1,
                      out=np.full(raw.shape[0],np.nan),where=counts>1))
            cvs=100*sd/means
            cvs[(counts<threshold)|(~np.isfinite(cvs))|(means<=1e-8)]=np.nan
            valid=cvs[np.isfinite(cvs)]
            if len(valid)<100:raise ValueError("Too few proteins for CV in "+batch)
            cv_by_batch.append({"stage":stage,"batch":batch,
                "n_samples":n,"minimum_positive_required":threshold,
                "n_proteins_with_valid_cv":len(valid),
                "median_protein_cv_pct":float(np.median(valid)),
                "q1_protein_cv_pct":float(np.percentile(valid,25)),
                "q3_protein_cv_pct":float(np.percentile(valid,75)),
                "mean_protein_cv_pct":float(np.mean(valid))})
            for k,v in enumerate(cvs):
                per_protein.append({"stage":stage,"batch":batch,
                    "protein_row_1based":k+1,"Protein.Group":annotations[k][0],
                    "Protein.Ids":annotations[k][1],
                    "n_positive_samples":int(counts[k]),
                    "cv_pct":float(v)})
            arrays[(stage,batch)]=valid
    return cv_by_batch,per_protein,arrays

def plot_cv(summary,cvvalues):
    fig,(ax,bx)=plt.subplots(1,2,figsize=(16,6.4),layout="constrained")
    raw=[next(r["median_protein_cv_pct"] for r in summary if r["stage"]=="Raw" and r["batch"]==b)
         for b in BATCHES]
    post=[next(r["median_protein_cv_pct"] for r in summary if r["stage"]=="Corrected" and r["batch"]==b)
          for b in BATCHES]
    xx=np.arange(4)
    ax.bar(xx-.19,raw,width=.36,color="#8293A2",label="Raw")
    ax.bar(xx+.19,post,width=.36,color=BLUE,label="Corrected")
    for pos,val in zip(list(xx-.19)+list(xx+.19),raw+post):
        ax.text(pos,val+.6,f"{val:.2f}%",ha="center",fontsize=9)
    ax.set_xticks(xx,[f"{b}\nn={BATCH_COUNTS[b]}" for b in BATCHES])
    ax.set_ylabel("Median per-protein CV (%) · original linear scale")
    ax.set_ylim(0,max(raw+post)*1.28);ax.legend(frameon=False)
    ax.set_title("Median protein CV by batch",loc="left",fontweight="bold")
    a=[cvvalues[("Raw",b)] for b in BATCHES]
    c=[cvvalues[("Corrected",b)] for b in BATCHES]
    # Explicitly avoid truncating any data in the statistics: view only up to
    # pooled 99th percentile, since long-tail CV can obscure the body.
    limit=float(np.percentile(np.concatenate(a+c),99))*1.09
    for k,(stage,dataset,shift,alpha) in enumerate(
        (("Raw",a,-.18,.23),("Corrected",c,.18,.45))):
        v=bx.violinplot(dataset,positions=np.arange(4)+1+shift,widths=.33,
                        showmedians=False,showextrema=False)
        for art,b in zip(v["bodies"],BATCHES):
            art.set_facecolor("#84909B" if stage=="Raw" else COLORS[b])
            art.set_alpha(alpha)
        bplot=bx.boxplot(dataset,positions=np.arange(4)+1+shift,widths=.12,
                         patch_artist=True,showfliers=False,manage_ticks=False)
        for art,b in zip(bplot["boxes"],BATCHES):
            art.set_facecolor("#84909B" if stage=="Raw" else COLORS[b])
            art.set_alpha(.63)
    bx.set_xticks(np.arange(1,5),list(BATCHES))
    bx.set_ylabel("Individual protein CV (%)")
    bx.set_ylim(0,limit)
    bx.set_title("CV distributions · Raw (left) / Corrected (right)\nDisplay capped at 99th percentile",
                 loc="left",fontweight="bold")
    for axis in (ax,bx):
        axis.grid(axis="y",alpha=.16);axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False);axis.spines["right"].set_visible(False)
    fig.suptitle("Raw vs single-HELP Heavy normalization: four NYD batches",fontweight="bold",fontsize=15)
    save(fig,"03_batch_protein_CV_raw_vs_corrected.png")

def raw_fit_pca(raw,corrected,names,records):
    # Feature criterion defined only on Raw, identical for both views.
    min_n=math.ceil(.70*raw.shape[1])
    selected=(np.isfinite(raw).sum(axis=1)>=min_n)&(np.nanstd(raw,axis=1)>1e-8)
    if np.sum(selected)<100:raise ValueError("Too few raw-fit PCA features")
    raw_features=raw[selected]
    corr_features=corrected[selected]
    med=np.nanmedian(raw_features,axis=1)
    rawfill=np.where(np.isfinite(raw_features),raw_features,med[:,None]).T
    corrfill=np.where(np.isfinite(corr_features),corr_features,med[:,None]).T
    raw_mean=rawfill.mean(axis=0)
    x=rawfill-raw_mean
    # Eigen decomposition of Raw's 40x40 Gram. Corrected never fit independently.
    gram=x@x.T
    eigenvalues,u=np.linalg.eigh(gram)
    order=np.argsort(eigenvalues)[::-1]
    eigenvalues=np.maximum(0,eigenvalues[order])
    u=u[:,order]
    if eigenvalues[1]<=1e-10:raise ValueError("Insufficient Raw PCA rank")
    components=(x.T@u[:,:2])/np.sqrt(eigenvalues[:2])[None,:]
    # Deterministic sign orientation from RAW-only loadings.
    for k in range(2):
        loc=int(np.argmax(np.abs(components[:,k])))
        if components[loc,k]<0:components[:,k]*=-1
    raw_score=x@components
    corrected_score=(corrfill-raw_mean)@components
    assert np.allclose(raw_score,raw_score.mean(axis=0)+raw_score-raw_score.mean(axis=0))
    explained=eigenvalues[:2]/max(eigenvalues.sum(),1e-10)*100
    scores=[]
    for j,n in enumerate(names):
        scores.append({"sample_id":n,"batch":records[j]["batch"],
           "Raw_PC1":float(raw_score[j,0]),"Raw_PC2":float(raw_score[j,1]),
           "Corrected_projected_PC1":float(corrected_score[j,0]),
           "Corrected_projected_PC2":float(corrected_score[j,1]),
           "PC1_delta":float(corrected_score[j,0]-raw_score[j,0]),
           "PC2_delta":float(corrected_score[j,1]-raw_score[j,1]),
           "projected_movement_distance":float(np.linalg.norm(corrected_score[j]-raw_score[j]))})
    audit={"fitted_on":"Raw PG log2 only",
        "corrected_data":"Projected onto original Raw components; not a second PCA fit",
        "feature_selection":"Raw detection >=70% of 40 and Raw SD>1e-8",
        "n_selected_proteins":int(np.sum(selected)),
        "missing_imputation":"Each selected protein's Raw median, identical for Raw and Corrected",
        "centering":"Each selected protein's Raw mean after Raw-only missing imputation",
        "PC1_raw_explained_variance_pct":float(explained[0]),
        "PC2_raw_explained_variance_pct":float(explained[1]),
        "PC_axis_orientation":"Raw loading with largest absolute magnitude is positive"}
    return scores,audit

def plot_pca(scores,audit):
    x=np.array([[r["Raw_PC1"],r["Raw_PC2"]] for r in scores])
    y=np.array([[r["Corrected_projected_PC1"],r["Corrected_projected_PC2"]] for r in scores])
    all_=np.vstack([x,y])
    xmin,xmax=all_[:,0].min(),all_[:,0].max()
    ymin,ymax=all_[:,1].min(),all_[:,1].max()
    pad_x=max((xmax-xmin)*.12,1.0)
    pad_y=max((ymax-ymin)*.12,1.0)
    fig,axes=plt.subplots(1,3,figsize=(18,6.6),sharex=True,sharey=True,layout="constrained")
    stages=[("Raw · fitted PCA",x),("Corrected · projected on Raw axes",y)]
    for ax,(title,coords) in zip(axes[:2],stages):
        for batch in BATCHES:
            selected=np.array([r["batch"]==batch for r in scores])
            points=coords[selected]
            ax.scatter(points[:,0],points[:,1],s=45,c=COLORS[batch],
                       edgecolors="white",linewidths=.5,alpha=.87,
                       label=f"{batch} (n={BATCH_COUNTS[batch]})")
        ax.set_title(title,loc="left",fontweight="bold")
        ax.legend(frameon=False,fontsize=8)
    ax=axes[2]
    for batch in BATCHES:
        selected=np.array([r["batch"]==batch for r in scores])
        for a,b in zip(x[selected],y[selected]):
            ax.plot([a[0],b[0]],[a[1],b[1]],c=COLORS[batch],alpha=.35,lw=.75)
        ax.scatter(x[selected,0],x[selected,1],marker="o",s=18,
                   c=COLORS[batch],alpha=.45)
        ax.scatter(y[selected,0],y[selected,1],marker="x",s=24,
                   c=COLORS[batch],alpha=.95)
    ax.set_title("Each sample: Raw ○ → Corrected ×",loc="left",fontweight="bold")
    for axis in axes:
        axis.set_xlim(xmin-pad_x,xmax+pad_x)
        axis.set_ylim(ymin-pad_y,ymax+pad_y)
        axis.set_xlabel(f"Raw PC1 ({audit['PC1_raw_explained_variance_pct']:.1f}%)")
        axis.set_ylabel(f"Raw PC2 ({audit['PC2_raw_explained_variance_pct']:.1f}%)")
        axis.axhline(0,c="#BCC5CC",lw=.7)
        axis.axvline(0,c="#BCC5CC",lw=.7)
        axis.grid(alpha=.13)
        axis.spines["top"].set_visible(False);axis.spines["right"].set_visible(False)
        axis.set_aspect("equal",adjustable="box")
    fig.suptitle("PCA: one RAW-fitted coordinate system for both data matrices",
                 fontsize=16,fontweight="bold")
    save(fig,"04_PCA_Raw_fitted_common_coordinates.png")

def plot_target_proteins(raw,after,records,targets):
    fig,axes=plt.subplots(2,2,figsize=(12.8,9.5),sharex=True,layout="constrained")
    for row,(acc,p) in enumerate(PROTEINS.items()):
        idx=targets[acc]
        lo=min(raw[idx].min(),after[idx].min())
        hi=max(raw[idx].max(),after[idx].max())
        for col,(stage,mat) in enumerate((("Raw",raw),("Corrected",after))):
            ax=axes[row,col]
            groups=[mat[idx,[j for j,r in enumerate(records) if r["batch"]==b]]
                    for b in BATCHES]
            bp=ax.boxplot(groups,positions=range(1,5),widths=.5,patch_artist=True,
                medianprops=dict(c="#273A48",lw=1.6),showfliers=False)
            for art,b in zip(bp["boxes"],BATCHES):
                art.set_facecolor(COLORS[b]);art.set_alpha(.44)
            for k,(b,group) in enumerate(zip(BATCHES,groups),1):
                ax.scatter(k+np.linspace(-.16,.16,len(group)),group,s=31,
                           c=p["color"],edgecolors="#3B454F",linewidths=.5,alpha=.83)
            ax.set_ylim(lo-(hi-lo)*.18,hi+(hi-lo)*.18)
            ax.set_xticks(range(1,5),list(BATCHES))
            ax.set_title(f"{p['gene']} / {acc} · {stage}",loc="left",fontweight="bold")
            ax.set_ylabel("log2 protein intensity")
            ax.grid(axis="y",alpha=.14)
            ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    fig.suptitle("Selected PG proteins by batch · original vs adjusted",fontweight="bold",fontsize=15)
    save(fig,"05_CFB_PROS1_by_batch_Raw_vs_corrected.png")

def plot_sample_medians(raw,after,names,records):
    r=np.nanmedian(raw,axis=0)
    a=np.nanmedian(after,axis=0)
    fig,ax=plt.subplots(figsize=(19,6.4),layout="constrained")
    x=np.arange(1,41)
    ax.plot(x,r,lw=1.1,c="#647787",alpha=.65,label="Raw sample PG median")
    ax.plot(x,a,lw=1.15,c=BLUE,alpha=.78,label="Corrected sample PG median")
    for j,record in enumerate(records):
        c=COLORS[record["batch"]]
        ax.plot([x[j],x[j]],[r[j],a[j]],c=c,alpha=.56,lw=1.2)
        ax.scatter([x[j],x[j]],[r[j],a[j]],c=["#66798A",c],
                   edgecolors="white",lw=.35,s=25,zorder=4)
    for cut in (16.5,24.5,32.5):ax.axvline(cut,c="#99A5B2",lw=.8)
    ax.set_xticks(x,names,fontsize=6.8,rotation=90)
    ax.set_ylabel("Median across detected log2 PG intensities")
    ax.set_title("Sample protein-distribution median: before and after Heavy-only shifts",
                 loc="left",fontweight="bold",fontsize=14)
    ax.legend(frameon=False)
    ax.grid(axis="y",alpha=.16);ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save(fig,"06_sample_PG_medians_Raw_vs_corrected.png")

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError("Existing trial directory protected: "+str(OUT))
    records,names,annotations,raw,targets,heavy=inputs()
    median_hi,median_h,factors,corrected=apply_calibration(raw,heavy)
    OUT.mkdir(parents=True,exist_ok=False)
    # Retain all annotation columns and sample IDs. Positive log2 intensities,
    # missing/nonpositive remain blank rather than being imputed.
    with (OUT/"corrected_NYD_pg_matrix2_log2.tsv").open("w",encoding="utf-8",newline="") as fl:
        writer=csv.writer(fl,delimiter="\t")
        writer.writerow(list(ANNOT)+names)
        for i,anno in enumerate(annotations):
            writer.writerow(anno+[float(v) if np.isfinite(v) else "" for v in corrected[i]])
    samples=[]
    for j,s in enumerate(names):
        r={"sample_id":s,"batch":records[j]["batch"],
           "Heavy_STGSWSTLK_original_linear":float(heavy[j]),
           "median_hi_log2_single_value":float(median_hi[j]),
           "median_h_log2_all_40":median_h,
           "fi_log2":float(factors[j]),
           "multiplier_linear":float(np.exp2(factors[j])),
           "STGSWSTLK_Heavy_after_factor_diagnostic_log2":float(median_hi[j]+factors[j]),
           "PG_median_log2_Raw":float(np.nanmedian(raw[:,j])),
           "PG_median_log2_Corrected":float(np.nanmedian(corrected[:,j])),
           "PG_detected_positive_n":int(np.isfinite(raw[:,j]).sum())}
        for acc,params in PROTEINS.items():
            i=targets[acc]
            r[f"{params['gene']}_{acc}_Raw_PG_log2"]=float(raw[i,j])
            r[f"{params['gene']}_{acc}_Corrected_PG_log2"]=float(corrected[i,j])
        samples.append(r)
    writecsv("sample_STGSWSTLK_Heavy_factors.csv",samples)
    summary,details,values=stats_cv(raw,corrected,records,annotations)
    writecsv("batch_protein_CV_summary.csv",summary)
    writecsv("per_protein_batch_CV_long.csv",details)
    pca_scores,pca_audit=raw_fit_pca(raw,corrected,names,records)
    writecsv("PCA_Raw_fitted_sample_coordinates.csv",pca_scores)
    plot_pgs(raw,corrected,records,names,targets)
    plot_heavy(median_hi,median_h,factors,records,names)
    plot_cv(summary,values)
    plot_pca(pca_scores,pca_audit)
    plot_target_proteins(raw,corrected,records,targets)
    plot_sample_medians(raw,corrected,names,records)
    cv={s:{b:next(r for r in summary if r["stage"]==s and r["batch"]==b) for b in BATCHES}
        for s in ("Raw","Corrected")}
    manifest={
        "sources":{path.name:sha256(path) for path in (PG,HELP,META)},
        "n_samples":40,"n_protein_groups":4264,"n_HELP_peptides_in_source":97,
        "one_selected_peptide":SEQ,
        "single_peptide_one_Heavy_measurement_per_sample":True,
        "40_sample_global_median_of_log2_STGSWSTLK_Heavy":median_h,
        "f_i":"median_h - median_hi, median_hi=log2 Heavy_STGSWSTLK_i (singleton)",
        "original_scale_factor":"2**(median_h - log2(H_i))",
        "applied_to":"All 4264 PG proteins in all 40 samples, not the HELP data file",
        "normalization_on_log2":True,"no_imputation":True,
        "raw_original_not_overwritten":True,
        "valid_PG_cells_Raw":int(np.isfinite(raw).sum()),
        "valid_PG_cells_Corrected":int(np.isfinite(corrected).sum()),
        "calibrated_Heavy_median_log2_range":float(np.ptp(median_hi+factors)),
        "batch_sizes":BATCH_COUNTS,"PCA":pca_audit,
        "PG_targets":{acc:{"gene":par["gene"],"data_row_1based":targets[acc]+1} for acc,par in PROTEINS.items()},
        "plots":[
            "01_PG_40samples_Raw_vs_corrected_highlight_CFB_PROS1.png",
            "02_STGSWSTLK_Heavy_per_sample_offsets.png",
            "03_batch_protein_CV_raw_vs_corrected.png",
            "04_PCA_Raw_fitted_common_coordinates.png",
            "05_CFB_PROS1_by_batch_Raw_vs_corrected.png",
            "06_sample_PG_medians_Raw_vs_corrected.png",
        ]}
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    lines=[
        "# NYD 40-sample single STGSWSTLK Heavy reference normalization","",
        "Independent, exploratory new trial. Original NYD PG and Heavy.xlsx unchanged.",
        "40 samples from metadata: T001=16, T002/T003/T004=8 each; 4264 PG rows.",
        "The ONE STGSWSTLK Heavy signal per sample is valid in all 40 samples.",
        "No pooled QC exists in this 40-sample metadata; the median reference is",
        "computed directly across the 40 original samples.","",
        "## Exact recipe","",
        "1. Log2-transform original positive PG and STGSWSTLK Heavy intensities.",
        "2. median_hi = median of the singleton STGSWSTLK Heavy observation in sample i",
        "   = log2(Heavy_i).",
        "3. median_h = median(median_hi) across ALL 40 samples.",
        "4. fi = median_h - median_hi.",
        "5. Corrected PG_i = raw PG log2_i + fi for every detected protein.",
        "6. Do NOT alter other HELP peptides, Heavy.xlsx, or original PG input.",
        "Nonpositive/unobserved PG values are missing in log2 outputs (not imputed).",
        "Correction multipliers on original linear intensities are 2**fi.",
        "The corrected-Heavy trace on the factor graph is a THEORETICAL check;",
        "the original Heavy workbook is intentionally unchanged.","",
        "## Median protein-group CV within each batch (%)","",
        "CV is computed on INVERSE-log2-transformed linear protein intensity",
        "with sample SD (ddof=1), requiring at least max(3,ceil(70% of batch n))",
        "positive measurements for each protein. HC/S disease labels are not used here.",
        "These are within-batch biological-sample CV values, NOT pure technical QC.","",
        "| Batch | n | Raw protein median CV | Corrected protein median CV | Delta (pp) |",
        "| --- | ---: | ---: | ---: | ---: |"]
    for b in BATCHES:
        before=cv["Raw"][b]["median_protein_cv_pct"]
        after=cv["Corrected"][b]["median_protein_cv_pct"]
        lines.append(f"| {b} | {BATCH_COUNTS[b]} | {before:.2f}% | {after:.2f}% | {after-before:+.2f} |")
    lines+=["","## Visualization 1: raw and corrected all-protein sample boxplots","",
            "CFB / P00751 is a blue diamond; PROS1 / P07225 is a yellow circle.",
            "![40 PG samples raw vs corrected](01_PG_40samples_Raw_vs_corrected_highlight_CFB_PROS1.png)","",
            "## Visualization 2: STGSWSTLK Heavy reference and fi factors","",
            "![Heavy reference factors](02_STGSWSTLK_Heavy_per_sample_offsets.png)","",
            "## Visualization 3: per-batch protein CV","",
            "![Raw and corrected CV](03_batch_protein_CV_raw_vs_corrected.png)","",
            "## Visualization 4: Raw-fitted PCA on common coordinates","",
            "Raw fit ONCE; corrected samples are projected onto identical Raw",
            "features, missing-value medians, centering and PCA loadings.",
            "PC1/PC2 explained percentages are from Raw only. Axis limits match",
            "across panels, and the third panel shows each sample movement.",
            "![Common-coordinates PCA](04_PCA_Raw_fitted_common_coordinates.png)","",
            "## Visualization 5: highlighted target proteins, Raw vs Corrected","",
            "![Target protein intensity](05_CFB_PROS1_by_batch_Raw_vs_corrected.png)","",
            "## Visualization 6: sample-level all-protein medians","",
            "![PG median shift](06_sample_PG_medians_Raw_vs_corrected.png)","",
            "## Exported files","",
            "- corrected_NYD_pg_matrix2_log2.tsv: complete 4264 x 40 corrected PG matrix.",
            "- sample_STGSWSTLK_Heavy_factors.csv: Heavy, medians, fi, linear multiplier,",
            "  before/after PG sample median and CFB/PROS1 PG intensities.",
            "- batch_protein_CV_summary.csv and per_protein_batch_CV_long.csv.",
            "- PCA_Raw_fitted_sample_coordinates.csv: shared PCA scores and movement.",
            "- run_manifest.json: exact original SHA256 and all methods.","",
            "## Caution","",
            "By using one Heavy peptide's deviation from an overall 40-sample",
            "median, the procedure assumes its signal changes largely reflect",
            "technical changes shared by other proteins. It can increase variation",
            "if this peptide has measurement noise or non-shared fluctuations.",
            "Lower across-clinical-sample CV does not prove better calibration.",
            ""]
    (OUT/"README.md").write_text("\n".join(lines),encoding="utf-8")
    print("GLOBAL_STGSWSTLK_MEDIAN_LOG2",median_h,flush=True)
    print("FACTORS_RANGE_LOG2",float(factors.min()),float(factors.max()),flush=True)
    print("CV_SUMMARY",json.dumps([{k:r[k] for k in ("stage","batch","median_protein_cv_pct","n_proteins_with_valid_cv")} for r in summary]),flush=True)
    print("PCA",json.dumps(pca_audit),flush=True)
    print("PASS: all 40 original STGSWSTLK Heavy -> fi, corrected PG, 6 PNG, one shared-frame PCA",flush=True)

if __name__=="__main__":
    main()
