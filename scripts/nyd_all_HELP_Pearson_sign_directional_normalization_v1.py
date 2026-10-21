#!/usr/bin/env python3
"""NYD exploratory 40-sample HELP-specific correlation-sign normalization.

For every 40/40 positive Heavy HELP:
  h_i = log2(one Heavy observation for sample i)
  median_h = median_i(h_i)
  fi = median_h - h_i
  r_g = Pearson across SAME 40 samples of raw log2 PG_g and h_i
  corrected_gi = raw_gi + sign(r_g)*fi, if Pearson defined;
  undefined/zero correlation PG rows are left unchanged.
Sign is determined once per PG over samples, NOT once per sample (undefined).
Full per-peptide corrected matrices are stored compressed, never overwrite raw.
Separate 39-sample exclusion sensitivity is independently refit and summarized.
"""
from __future__ import annotations
import csv, gzip, hashlib, json, math, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from filter_samples import readbook,tab
from nyd_pg_batch_boxplots import meta, pg_matrix,PG, META, BATCHES,COLORS, ANNOT
import nyd_STGSWSTLK_heavy_exclude_T003_outlier as older

BASE=Path(__file__).resolve().parents[1]/"NYD_data"/"40samples"
HEAVY=BASE/"Heavy.xlsx"
OUT=BASE/"HELP97_Pearson_sign_directional_correction_all_detected_v1"
MIN_PEARSON_OBS=3
MIN_CV_FRACTION=.7
EXCLUDE=older.EXCLUDED_ID

def load():
    records=[r for batch in BATCHES for r in meta() if r["batch"]==batch]
    names=[r["sample_id"] for r in records]
    if len(names)!=40 or len(set(names))!=40:
        raise ValueError("Expected 40 unique NYD samples")
    raw,nproteins,npositive=pg_matrix(records)
    assert raw.shape==(4264,40)
    annotation=pd.read_csv(PG,sep="\t",encoding="utf-8-sig",usecols=list(ANNOT),
                           dtype=str,keep_default_na=False)
    assert len(annotation)==4264
    shared,sheets=readbook(HEAVY)
    if len(sheets)!=1:raise ValueError("One Heavy sheet expected")
    rows=tab(sheets[0],shared)
    if len(rows)!=98:raise ValueError("Expected 97 Heavy sequences")
    header={str(v).strip():c for c,v in rows[0].items() if c!="A"}
    assert len(header)==40 and set(header)==set(names)
    seqs=[str(row.get("A","")).strip() for row in rows[1:]]
    if len(seqs)!=97 or len(set(seqs))!=97:raise ValueError("Invalid HELP sequences")
    heavy=np.array([[older.val(row.get(header[s],"")) for s in names]
                    for row in rows[1:]],float)
    assert heavy.shape==(97,40)
    return records,names,raw,annotation,seqs,heavy

def protein_pearson(logPG,h):
    """Pearson per PG over the sample dimension, pairwise positive detection."""
    present=np.isfinite(logPG)&np.isfinite(h)[None,:]
    n=present.sum(axis=1)
    x=np.where(present,logPG,0.)
    z=np.where(present,h[None,:],0.)
    den=np.maximum(n,1)
    mx=np.sum(x,axis=1)/den
    my=np.sum(z,axis=1)/den
    xc=np.where(present,logPG-mx[:,None],0.)
    yc=np.where(present,h[None,:]-my[:,None],0.)
    cross=np.sum(xc*yc,axis=1)
    vx=np.sum(xc*xc,axis=1)
    vy=np.sum(yc*yc,axis=1)
    norm=np.sqrt(vx*vy)
    good=(n>=MIN_PEARSON_OBS)&(norm>1e-10)
    pearson=np.full(logPG.shape[0],np.nan)
    pearson[good]=np.clip(cross[good]/norm[good],-1,1)
    direction=np.zeros(logPG.shape[0],np.int8)
    direction[pearson>0]=1
    direction[pearson<0]=-1
    return pearson,direction,n

def cv_array(m):
    n=np.isfinite(m).sum(axis=1)
    linear=np.where(np.isfinite(m),np.exp2(m),0.)
    summ=linear.sum(axis=1)
    mu=summ/np.maximum(n,1)
    sq=np.square(linear-mu[:,None])
    sq=np.where(np.isfinite(m),sq,0.).sum(axis=1)
    st=np.sqrt(sq/np.maximum(n-1,1))
    result=np.full(len(n),np.nan)
    good=(n>=2)&(mu>0)
    result[good]=100*st[good]/mu[good]
    return result,n

def summaryCV(m,records):
    out={}
    for b in ["ALL"]+list(BATCHES):
        ind=[i for i,r in enumerate(records) if b=="ALL" or r["batch"]==b]
        cv,n=cv_array(m[:,ind])
        minimum=max(3,math.ceil(len(ind)*MIN_CV_FRACTION))
        selected=cv[(n>=minimum)&np.isfinite(cv)]
        out[b]={"n_samples":len(ind),"n_proteins":len(selected),
                "median_CV_pct":float(np.median(selected)) if len(selected) else float("nan")}
    return out

def raw_pca_basis(raw):
    valid=np.isfinite(raw)
    n=valid.sum(axis=1)
    selected=(n>=math.ceil(raw.shape[1]*.7))&(
        np.nanvar(np.where(valid,raw,np.nan),axis=1)>1e-8)
    if selected.sum()<200:raise AssertionError("Too few PCA proteins")
    x=raw[selected]
    med=np.nanmedian(x,axis=1)
    imputed=np.where(np.isfinite(x),x,med[:,None])
    mat=imputed.T
    center=mat.mean(axis=0)
    u,s,vt=np.linalg.svd(mat-center,full_matrices=False)
    components=vt[:2].T
    score=(mat-center)@components
    variance=s*s
    explained=100*variance[:2]/variance.sum()
    return selected,med,center,components,score,explained

def project_pca(after,basis):
    selected,med,center,comp,rawscore,explained=basis
    x=after[selected]
    complete=np.where(np.isfinite(x),x,med[:,None])
    return (complete.T-center)@comp

def frame(ax):
    ax.grid(axis="y",alpha=.14)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def add_batch_labels(ax,records):
    for b in BATCHES:
        ids=np.array([i for i,r in enumerate(records) if r["batch"]==b])
        ax.text(float(ids.mean()+1),1.025,b+" ("+str(len(ids))+")",
                transform=ax.get_xaxis_transform(),ha="center",va="bottom",
                fontsize=8.5,color=COLORS[b],weight="bold",clip_on=False)
        last=int(ids[-1]+1)
        if last<len(records):
            ax.axvline(last+.5,color="#899DAD",ls=":",lw=1)

def chart_box(seq,raw,after,records,path):
    x=np.arange(1,len(records)+1)
    samplelabels=[r["sample_id"][-6:] for r in records]
    allm=np.concatenate([raw[np.isfinite(raw)],after[np.isfinite(after)]])
    lim=(np.quantile(allm,.004)-1,np.quantile(allm,.996)+1)
    fig,axes=plt.subplots(2,1,figsize=(19,11.3),sharex=True)
    fig.subplots_adjust(left=.065,right=.985,top=.86,bottom=.14,hspace=.32)
    for ax,m,label in zip(axes,[raw,after],["Original PG log2","Pearson-sign HELP corrected PG log2"]):
        box=ax.boxplot([a[np.isfinite(a)] for a in m.T],positions=x,widths=.65,
                       patch_artist=True,showfliers=False,manage_ticks=False,
                       medianprops={"color":"#33485B","lw":1.4})
        for patch,r in zip(box["boxes"],records):
            patch.set_facecolor(COLORS[r["batch"]]);patch.set_alpha(.58)
        add_batch_labels(ax,records);ax.set_ylim(*lim);ax.set_ylabel("log2 PG intensity")
        ax.set_title(label,loc="left",fontsize=11,pad=24);frame(ax)
    axes[-1].set_xticks(x,samplelabels,rotation=90,fontsize=6.7)
    fig.suptitle("NYD "+seq+" | all 40 samples | original / directional corrected PG",
                 fontsize=14,weight="bold")
    fig.savefig(path,dpi=115,bbox_inches="tight",facecolor="white")
    plt.close(fig)

def chart_factor(seq,raw,after,fi,records,path):
    x=np.arange(1,len(fi)+1)
    medraw=np.nanmedian(raw,axis=0)
    medafter=np.nanmedian(after,axis=0)
    fig,ax=plt.subplots(2,1,figsize=(18.5,9.4),sharex=True)
    fig.subplots_adjust(left=.07,right=.985,top=.84,bottom=.21,hspace=.36)
    for b in BATCHES:
        ix=np.array([k for k,r in enumerate(records) if r["batch"]==b])
        ax[0].plot(x[ix],medraw[ix],color=COLORS[b],marker="o",lw=1.5,ms=4.5)
        ax[0].plot(x[ix],medafter[ix],color=COLORS[b],marker="D",ls="--",lw=1.5,ms=4)
        ax[1].bar(x[ix],fi[ix],color=COLORS[b],alpha=.8)
    for a in ax:
        add_batch_labels(a,records);frame(a)
    ax[0].set_ylabel("PG sample median log2")
    ax[1].set_ylabel("fi = median(H) - H_i (log2)")
    ax[1].axhline(0,color="#90A2B1",lw=.8)
    ax[0].set_title("Sample PG median before / after correction",loc="left",fontsize=11,pad=26)
    ax[1].set_title("One-HELP signed sample correction factors",loc="left",fontsize=11,pad=26)
    ax[1].set_xticks(x,[r["sample_id"][-6:] for r in records],rotation=90,fontsize=6.8)
    handles=[Line2D([0],[0],color="#506175",lw=2,marker="o",label="Raw PG median"),
             Line2D([0],[0],color="#506175",lw=2,ls="--",marker="D",
                    label="Corrected PG median")]
    handles += [Line2D([0],[0],color=COLORS[b],lw=3,label=b) for b in BATCHES]
    fig.legend(handles=handles,ncol=6,loc="lower center",bbox_to_anchor=(.5,.04),frameon=False)
    fig.suptitle("NYD "+seq+" | Heavy fi and protein-median trajectories",
                 fontsize=14,weight="bold")
    fig.savefig(path,dpi=120,bbox_inches="tight",facecolor="white")
    plt.close(fig)

def chart_pca(seq,basis,after_scores,records,path):
    before=basis[4]; expl=basis[5]
    both=np.concatenate([before,after_scores],axis=0)
    low=both.min(axis=0); hi=both.max(axis=0)
    pad=np.maximum((hi-low)*.13,1.)
    fig,axes=plt.subplots(1,2,figsize=(14,6.2))
    fig.subplots_adjust(left=.08,right=.98,top=.85,bottom=.18,wspace=.23)
    for ax,points,title in zip(axes,[before,after_scores],
                    ["Raw PCA","Corrected, projected into Raw PCA space"]):
        for b in BATCHES:
            sel=np.array([r["batch"]==b for r in records])
            ax.scatter(points[sel,0],points[sel,1],s=39,color=COLORS[b],
                       edgecolor="white",lw=.5,label=b)
        ax.set_xlim(low[0]-pad[0],hi[0]+pad[0])
        ax.set_ylim(low[1]-pad[1],hi[1]+pad[1])
        ax.set_xlabel("Raw PC1 ({:.1f}%)".format(expl[0]))
        ax.set_ylabel("Raw PC2 ({:.1f}%)".format(expl[1]))
        ax.set_title(title,fontsize=10.5);frame(ax)
    axes[1].legend(frameon=False,ncol=4,fontsize=8,loc="upper center")
    fig.suptitle("NYD "+seq+" | one Raw-fitted PCA basis, batch colors",
                 fontweight="bold",fontsize=13)
    fig.savefig(path,dpi=125,bbox_inches="tight",facecolor="white")
    plt.close(fig)

def plot_overview(metrics,out):
    t=pd.DataFrame(metrics).sort_values("HELP_sequence")
    fig,axes=plt.subplots(1,2,figsize=(14,6.1))
    fig.subplots_adjust(left=.07,right=.99,top=.88,bottom=.15,wspace=.28)
    delta=t["ALL_delta_medianCV"].to_numpy(float)
    axes[0].hist(delta,bins=20,color="#427B96",edgecolor="white")
    axes[0].axvline(0,color="#C76A3F",lw=1.5,ls="--")
    axes[0].set_xlabel("After - before median PG CV (percentage points)")
    axes[0].set_ylabel("Number of Heavy peptides")
    axes[0].set_title("Overall PG CV change across full-40 HELP variants")
    axes[1].scatter(t["pearson_positive_PG_n"],delta,
                    c=t["pearson_negative_PG_n"],cmap="viridis",alpha=.8,s=40)
    axes[1].axhline(0,color="#C76A3F",lw=1.1,ls="--")
    axes[1].set_xlabel("PG with positive Pearson correlation")
    axes[1].set_ylabel("After - before median PG CV")
    axes[1].set_title("Correlation-sign splitting vs CV change")
    for a in axes:frame(a)
    fig.savefig(out/"00_overview_deltaCV_all_HELP.png",dpi=175,bbox_inches="tight")
    plt.close(fig)

def run(cohort,records,names,raw,annotations,seqs,H,save_all):
    n=len(names)
    assert n in (40,39)
    valid=np.isfinite(H)&(H>0)
    complete=valid.all(axis=1)
    candidates=np.flatnonzero(complete)
    if len(candidates)==0:raise RuntimeError("No fully detected Heavy peptides")
    basecv=summaryCV(raw,records)
    basis=raw_pca_basis(raw) if save_all else None
    original_names=["Raw_"+c for c in BATCHES]
    outdir=OUT/cohort
    outdir.mkdir(parents=True,exist_ok=False)
    peptide_results=[]
    factor_rows=[]
    pearson_rows=[]
    for counter,k in enumerate(candidates,1):
        seq=seqs[k]
        h=np.log2(H[k])
        ref=float(np.median(h))
        fi=ref-h
        assert np.all(np.isfinite(fi))
        pearson,sign,nobs=protein_pearson(raw,h)
        corrected=raw+sign[:,None]*fi[None,:]
        assert np.array_equal(np.isfinite(raw),np.isfinite(corrected))
        assert np.allclose(corrected-sign[:,None]*fi[None,:],raw,equal_nan=True,atol=1e-11)
        cv=summaryCV(corrected,records)
        metric={
            "HELP_sequence":seq,"n_samples":n,
            "Heavy_ref_median_log2":ref,
            "fi_min":float(fi.min()),"fi_max":float(fi.max()),
            "fi_std":float(fi.std(ddof=1)),
            "pearson_positive_PG_n":int((sign==1).sum()),
            "pearson_negative_PG_n":int((sign==-1).sum()),
            "pearson_undefined_or_zero_PG_n":int((sign==0).sum()),
            "pearson_median_valid":float(np.nanmedian(pearson)),
            "PG_median_sample_SD_Raw":float(np.nanstd(np.nanmedian(raw,axis=0),ddof=1)),
            "PG_median_sample_SD_corrected":float(np.nanstd(np.nanmedian(corrected,axis=0),ddof=1)),
        }
        for b in ["ALL"]+list(BATCHES):
            metric[b+"_raw_medianCV"]=basecv[b]["median_CV_pct"]
            metric[b+"_after_medianCV"]=cv[b]["median_CV_pct"]
            metric[b+"_delta_medianCV"]=cv[b]["median_CV_pct"]-basecv[b]["median_CV_pct"]
            metric[b+"_proteins_CV_n"]=cv[b]["n_proteins"]
        peptide_results.append(metric)
        for j,sid in enumerate(names):
            factor_rows.append({
                "HELP_sequence":seq,"sample_id":sid,"batch":records[j]["batch"],
                "median_hi_log2":float(h[j]),"median_h_log2":ref,
                "fi_log2":float(fi[j]),"linear_multiplier_if_positivePearson":float(np.exp2(fi[j])),
                "PG_median_raw":float(np.nanmedian(raw[:,j])),
                "PG_median_after":float(np.nanmedian(corrected[:,j]))
            })
        for i in range(len(sign)):
            pearson_rows.append({"HELP_sequence":seq,"protein_row_index":i+1,
                "Protein.Group":str(annotations.iloc[i]["Protein.Group"]),
                "Genes":str(annotations.iloc[i]["Genes"]),
                "n_pairwise_samples":int(nobs[i]),
                "Pearson_r":float(pearson[i]) if np.isfinite(pearson[i]) else np.nan,
                "correction_sign":int(sign[i])})
        if save_all:
            dest=outdir/"per_HELP"/seq
            dest.mkdir(parents=True,exist_ok=False)
            data=annotations.copy()
            for j,sid in enumerate(names):
                data[sid]=corrected[:,j]
            data.to_csv(dest/"corrected_PG_log2.tsv.gz",sep="\t",index=False,
                        na_rep="",compression={"method":"gzip","compresslevel":4})
            chart_box(seq,raw,corrected,records,dest/"01_Raw_vs_corrected_boxplots.png")
            chart_factor(seq,raw,corrected,fi,records,
                         dest/"02_Heavy_factor_and_PG_sample_medians.png")
            chart_pca(seq,basis,project_pca(corrected,basis),records,
                      dest/"03_Raw_vs_corrected_common_PCA.png")
        if counter%10==0 or counter==len(candidates):
            print("COHORT",cohort,"completed",counter,"of",len(candidates),flush=True)
    pd.DataFrame(peptide_results).sort_values("ALL_delta_medianCV").to_csv(
        outdir/"01_HELP_all_candidates_metrics_sorted_by_deltaCV.csv",index=False)
    pd.DataFrame(factor_rows).to_csv(
        outdir/"02_every_HELP_per_sample_factors.csv.gz",index=False,
        compression={"method":"gzip","compresslevel":4})
    pd.DataFrame(pearson_rows).to_csv(
        outdir/"03_every_HELP_per_protein_Pearson_signs.csv.gz",index=False,
        compression={"method":"gzip","compresslevel":4})
    if save_all:
        plot_overview(peptide_results,outdir)
    return {"n_samples":n,"full_detected_n":len(candidates),
            "full_detected_HELP":[seqs[i] for i in candidates],
            "batch_counts":{b:sum(r["batch"]==b for r in records) for b in BATCHES},
            "raw_CV":basecv,
            "n_helps_after_CV_improvement":sum(m["ALL_delta_medianCV"]<0 for m in peptide_results),
            "best_by_deltaCV":sorted(peptide_results,key=lambda x:x["ALL_delta_medianCV"])[:8],
            "per_HELP_full_outputs":bool(save_all)}

def main():
    if OUT.exists():raise FileExistsError("Result folder already exists: "+str(OUT))
    records,names,raw,ann,seqs,heavy=load()
    OUT.mkdir(parents=True)
    assert len(seqs)==97 and raw.shape==(4264,40)
    primary=run("full40",records,names,raw,ann,seqs,heavy,save_all=True)
    ix=[j for j,s in enumerate(names) if s!=EXCLUDE]
    assert len(ix)==39 and EXCLUDE in names
    records39=[records[j] for j in ix]
    names39=[names[j] for j in ix]
    sensitivity=run("exclude_T003_bad_sample_39",records39,names39,raw[:,ix],
                    ann,seqs,heavy[:,ix],save_all=False)
    summary={
        "experiment":"NYD 97 Heavy HELP Pearson-sign protein-wise normalization",
        "primary":primary,"sensitivity":sensitivity,"total_initial_HELP":97,
        "PG_count":4264,"Pearson_definition":"per protein across samples correlation with one peptide Heavy (log2)",
        "Pearson_min_pairwise":MIN_PEARSON_OBS,
        "Pearson_undefined":"unchanged protein across all samples",
        "factor":"fi=median_samples(log2 Heavy_j)-log2 Heavy_ji; one intensity per sample",
        "direction":"PG_corrected=PG_raw+sign(Pearson_gj)*fi",
        "sample_order":"batch T001,T002,T003,T004; metadata order within batch",
        "CV_definition":"linear PG sample SD (ddof=1)/mean; positive complete per batch >=70%, min 3",
        "CV_caveat":"no true identical-replicate QC; clinical within-batch protein CV includes biology",
        "PCA":"Raw fit with >=70% detected PG; Raw-only medians impute PCA, fixed basis projection",
        "full40_output":"corrected 4264x40 PG log2 matrix gz + 3 figures per fully detected peptide",
        "sensitivity39":"independent cohort, independently refit Pearson and Heavy median; summary only",
        "source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (PG,HEAVY,META)},
    }
    (OUT/"run_manifest.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    readme=[
        "# NYD: Every fully detected Heavy HELP with protein-specific Pearson-sign correction",
        "",
        "## Protocol",
        "",
        "Primary: all **40 original NYD samples**; each Heavy peptide must have 40/40 detected.",
        "Separately rerun after excluding known T003 outlier MAP26092326E01 (39 samples);",
        "do not silently mix candidate peptide sets, Pearson references, or CV results.",
        "",
        "For one peptide j and sample i:",
        "",
        "- h_ji = log2 Heavy intensity from the single peptide measurement.",
        "- ref_j = median over samples of h_ji.",
        "- fi = ref_j - h_ji (signed; distance not absolute value).",
        "- For each protein g, Pearson r_gj is computed across sample axis",
        "  between its original log2 PG signal and log2 Heavy signal j, pairwise available.",
        "- Positive Pearson: add fi to original log2 PG_gi.",
        "- Negative Pearson: subtract fi from original log2 PG_gi.",
        "- Undefined or exact zero: keep original PG_gi unchanged (no invented direction).",
        "",
        "IMPORTANT: a sample only has one Heavy reading for one peptide, so",
        "a per-sample Pearson is undefined. Signs attach to protein groups,",
        "not to individual samples. Missing original protein values stay missing.",
        "",
        "## Results",
        "",
        f"- Primary fully detected: {primary['full_detected_n']}/97 Heavy HELP.",
        f"- Sensitivity 39 fully detected: {sensitivity['full_detected_n']}/97 Heavy HELP.",
        f"- Primary HELP variants improving overall within-cohort median PG CV: "
        f"{primary['n_helps_after_CV_improvement']}/{primary['full_detected_n']}.",
        "",
        "## Main ranking (full 40)",
        "",
        "| HELP | raw median CV | corrected median CV | change (percentage points) | positive proteins | negative proteins |",
        "|---|---:|---:|---:|---:|---:|"
    ]
    for x in primary["best_by_deltaCV"][:8]:
        readme.append("| "+x["HELP_sequence"]+" | {:.2f} | {:.2f} | {:+.2f} | {} | {} |".format(
            x["ALL_raw_medianCV"],x["ALL_after_medianCV"],x["ALL_delta_medianCV"],
            x["pearson_positive_PG_n"],x["pearson_negative_PG_n"]))
    readme += [
        "",
        "## Outputs",
        "",
        "- full40/01_HELP_all_candidates_metrics_sorted_by_deltaCV.csv",
        "- full40/02_every_HELP_per_sample_factors.csv.gz",
        "- full40/03_every_HELP_per_protein_Pearson_signs.csv.gz",
        "- full40/per_HELP/<sequence>/corrected_PG_log2.tsv.gz",
        "- full40/per_HELP/<sequence>/01_Raw_vs_corrected_boxplots.png",
        "- full40/per_HELP/<sequence>/02_Heavy_factor_and_PG_sample_medians.png",
        "- full40/per_HELP/<sequence>/03_Raw_vs_corrected_common_PCA.png",
        "- full40/00_overview_deltaCV_all_HELP.png",
        "- exclude_T003_bad_sample_39/01_HELP_all_candidates_metrics_sorted_by_deltaCV.csv",
        "- exclude_T003_bad_sample_39/02_every_HELP_per_sample_factors.csv.gz",
        "- exclude_T003_bad_sample_39/03_every_HELP_per_protein_Pearson_signs.csv.gz",
        "",
        "## Limitations",
        "",
        "- This method chooses correction direction from raw protein-Heavy correlations",
        "  **on the same clinical samples**; that can encode biology and is not",
        "  independent evidence of removal of technical batch effects.",
        "- CV over 40 clinical individuals is not replicate-only technical CV.",
        "- Corrected sample-to-sample Pearson correlations need not remain unchanged:",
        "  signs vary by protein so transformation is NOT one constant per sample.",
        "- No QC or calibration reference beyond clinical samples is available.",
        "- No source files or earlier results overwritten.",
    ]
    (OUT/"README.md").write_text("\n".join(readme)+"\n",encoding="utf-8")
    print("FINISHED_NYD",json.dumps({
        "full40_HelpCount":primary["full_detected_n"],
        "full39_HelpCount":sensitivity["full_detected_n"],
        "full40_improved":primary["n_helps_after_CV_improvement"],
        "best8":primary["best_by_deltaCV"]
    },ensure_ascii=False),flush=True)
if __name__=="__main__":main()
