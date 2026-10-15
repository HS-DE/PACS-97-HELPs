#!/usr/bin/env python3
"""40-sample median intensity normalization, independent PG and PR matrices.

Based on the median-normalization *principle* reported by Tsantilas et al.
(2026, MacCoss group), not a claim that the exact authors' code is reproduced.
Every sample gets a multiplicative scale equal to 2**(reference_log2_median -
its own log2_median), with reference = median of the 40 sample log2 medians.
PG and PR matrices are independently normalized, using only non-yeast features.
There are no reference QC samples; yeast is assessed independently, NOT used to
calculate normalization factors. Original files remain untouched.
"""
import csv,hashlib,json,math,re,statistics
from collections import Counter
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nyd_pg_batch_boxplots import meta, PG, META, GROUPS, BATCHES, COLORS
from nyd_yeast_batch_cv import load as load_yeast

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"NYD_data"/"40samples"
PR=BASE/"NYD_pr_matrix2.tsv"
OUT=BASE/"median_normalization_40samples"
ANNOT_FIELDS=("Protein.Group","Protein.Ids","Protein.Names","Genes")
YEAST_PATTERN=re.compile(r"(?:^|[;,\s|])P00924(?:-\d+)?(?:$|[;,\s|])|ENO[12]_YEAST",re.I)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()

def numeric(val):
    try:
        x=float(val)
        return x if math.isfinite(x) and x>0 else np.nan
    except (TypeError,ValueError):return np.nan

def writecsv(path,data):
    with (OUT/path).open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]))
        w.writeheader();w.writerows(data)

def load_matrix(path,sids):
    with path.open("r",encoding="utf-8-sig",newline="") as f:
        reader=csv.reader(f,delimiter="\t")
        header=next(reader)
        if len(header)!=len(set(header)):raise ValueError("Duplicate headers "+path.name)
        if not set(sids).issubset(header) or len(sids)!=40:
            raise ValueError("Not all 40 samples present in "+path.name)
        indices=[header.index(x) for x in sids]
        other=[x for x in header if x not in sids]
        if not "Protein.Group" in other:
            raise ValueError("No Protein.Group annotation "+path.name)
        original=[]
        for line,r in enumerate(reader,start=2):
            if not r or all(not v for v in r):continue
            if len(r)!=len(header):
                raise ValueError(f"{path.name} malformed line {line}: {len(r)} instead of {len(header)}")
            original.append(r)
    if len(original)<500:raise ValueError("Suspiciously small protein/precursor matrix")
    intensity=np.array([[numeric(row[k]) for k in indices] for row in original],dtype=float)
    ann_indices=[header.index(k) for k in ANNOT_FIELDS if k in header]
    excluded=np.asarray([any(YEAST_PATTERN.search(row[k] or "")
                              for k in ann_indices) for row in original])
    if (~excluded).sum()<500:raise ValueError("Too few nonyeast features")
    log=np.full_like(intensity,np.nan)
    ok=np.isfinite(intensity)&(intensity>0)
    log[ok]=np.log2(intensity[ok])
    return header,original,indices,log,excluded

def calc_factors(log,excluded,meta_rows):
    ref_input=log[~excluded,:]
    med=np.nanmedian(ref_input,axis=0)
    count=np.isfinite(ref_input).sum(axis=0)
    if not np.all(np.isfinite(med)) or np.min(count)<50:
        raise ValueError("Insufficient measured nonyeast protein/precursor features")
    ref=float(np.median(med))
    delta=ref-med
    factors=np.exp2(delta)
    if not np.all(np.isfinite(factors)):
        raise ValueError("Nonfinite normalization factors")
    rows=[]
    for k,m in enumerate(meta_rows):
        rows.append(dict(
            sample_id=m["sample_id"],batch=m["batch"],
            n_features_for_factor=int(count[k]),
            raw_median_log2=float(med[k]),
            target_median_log2=ref,
            additive_log2_offset=float(delta[k]),
            multiplicative_factor=float(factors[k]),
            expected_median_normalized_log2=float(med[k]+delta[k])))
    return rows,factors,med,ref

def write_normalized(source_header,original,colidx,factors,filename):
    file=OUT/filename
    with file.open("w",encoding="utf-8",newline="") as dst:
        wr=csv.writer(dst,delimiter="\t",lineterminator="\n")
        wr.writerow(source_header)
        for row in original:
            copied=row.copy()
            for k,index in enumerate(colidx):
                value=numeric(row[index])
                if np.isfinite(value):
                    copied[index]=format(value*factors[k],".15g")
                # Keep original empty/zero/nonpositive/nonfinite entries exactly.
            wr.writerow(copied)
    return file

def matrix_stats(log,factors,source,excluded,meta_rows):
    normalized=log+np.log2(factors)[None,:]
    groups={b:[i for i,r in enumerate(meta_rows) if r["batch"]==b] for b in BATCHES}
    sample=[]
    for k,r in enumerate(meta_rows):
        observed=log[:,k][np.isfinite(log[:,k])]
        updated=normalized[:,k][np.isfinite(normalized[:,k])]
        sample.append(dict(matrix=source,sample_id=r["sample_id"],batch=r["batch"],
           positive_feature_n=len(observed),total_feature_n=len(log),
           raw_median_log2=float(np.median(observed)),
           corrected_median_log2=float(np.median(updated)),
           raw_q1_log2=float(np.percentile(observed,25)),
           raw_q3_log2=float(np.percentile(observed,75)),
           corrected_q1_log2=float(np.percentile(updated,25)),
           corrected_q3_log2=float(np.percentile(updated,75))))
    batch=[]
    for b,idx in groups.items():
        for stage,mat in (("raw",log),("median_normalized",normalized)):
            pooled=mat[:,idx].ravel()
            pooled=pooled[np.isfinite(pooled)]
            batch.append(dict(matrix=source,batch=b,stage=stage,
                n_samples=len(idx),n_positive_feature_sample_values=len(pooled),
                median_all_log2=float(np.median(pooled)),
                q1_all_log2=float(np.percentile(pooled,25)),
                q3_all_log2=float(np.percentile(pooled,75)),
                median_sample_medians_log2=float(np.median([
                    np.nanmedian(mat[:,k]) for k in idx]))))
    cvs=[]
    for b,idx in groups.items():
        n=len(idx);threshold=max(3,math.ceil(.70*n))
        for stage,mat in (("raw",log),("median_normalized",normalized)):
            # Linear intensity CV. Every row is a protein group (PG) or
            # precursor feature (PR); differences may include biology.
            linear=np.exp2(mat[~excluded][:,idx])
            good=np.isfinite(linear)&(linear>0)
            valid=good.sum(axis=1)>=threshold
            count=good.sum(axis=1)
            means=np.divide(np.nansum(linear,axis=1),np.maximum(count,1))
            delta=np.where(good,linear-means[:,None],0)
            ss=np.sum(delta**2,axis=1)
            sd=np.sqrt(np.divide(ss,np.maximum(count-1,1)))
            ok=valid&(means>0)
            value=np.full(len(means),np.nan)
            value[ok]=100*sd[ok]/means[ok]
            cv=value[np.isfinite(value)]
            if len(cv)<500:raise ValueError("Too few features for within-batch CV summary")
            cvs.append(dict(matrix=source,batch=b,stage=stage,
                min_positive_n=threshold,n_features_cv=len(cv),
                median_feature_cv_pct=float(np.median(cv)),
                q1_feature_cv_pct=float(np.percentile(cv,25)),
                q3_feature_cv_pct=float(np.percentile(cv,75))))
    return normalized,sample,batch,cvs

def style(ax):
    ax.grid(axis="y",alpha=.18)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)

def save(fig,name):
    fig.savefig(OUT/name,dpi=190,bbox_inches="tight",facecolor="white")
    plt.close(fig)

def per_sample_plot(raw,post,meta_rows,name):
    positions=np.arange(1,41)
    fig,axes=plt.subplots(2,1,figsize=(19,10.8),sharex=True,
                          gridspec_kw=dict(hspace=.24))
    allvals=np.concatenate([raw[np.isfinite(raw)],post[np.isfinite(post)]])
    ylo,yhi=np.percentile(allvals,[.25,99.75])
    for ax,mat,title in zip(axes,(raw,post),
                            ("Original signals","Median-normalized signals")):
        data=[mat[:,i][np.isfinite(mat[:,i])] for i in range(40)]
        boxes=ax.boxplot(data,positions=positions,widths=.63,patch_artist=True,
            showfliers=False,manage_ticks=False,
            medianprops=dict(color="#182F41",lw=1.55))
        for patch,r in zip(boxes["boxes"],meta_rows):
            patch.set_facecolor(COLORS[r["batch"]]);patch.set_alpha(.74)
        ax.set_ylim(ylo-.4,yhi+.45)
        ax.set_xlim(.4,40.6)
        ax.set_ylabel("log2 measured intensity",fontsize=11)
        ax.set_title(title,fontsize=14,fontweight="bold",loc="left")
        for edge in (16.5,24.5,32.5):
            ax.axvline(edge,color="#9CAAB7",lw=1)
        style(ax)
    axes[-1].set_xticks(positions,[r["sample_id"] for r in meta_rows],
                        rotation=90,fontsize=7)
    for b,(left,right) in zip(BATCHES,((1,16),(17,24),(25,32),(33,40))):
        axes[0].text((left+right)/2,1.035,f"{b} · n={GROUPS[b]}",
                  color=COLORS[b],transform=axes[0].get_xaxis_transform(),
                  ha="center",fontsize=11,fontweight="bold")
    fig.suptitle(f"NYD {name.upper()} · before vs after sample median normalization",
                 fontsize=16,fontweight="bold",y=.985)
    save(fig,f"01_{name}_40sample_boxplots_before_after.png")

def batch_medians_plot(pairs):
    fig,axes=plt.subplots(2,2,figsize=(13.5,10),layout="constrained")
    for ri,(name,raw,post) in enumerate(pairs):
        for ci,(title,mat) in enumerate((("Raw",raw),("Median normalized",post))):
            ax=axes[ri,ci]
            vals=[[float(np.nanmedian(mat[:,i])) for i,r in enumerate(ALL_META)
                   if r["batch"]==b] for b in BATCHES]
            boxes=ax.boxplot(vals,positions=range(1,5),patch_artist=True,
                             widths=.52,showfliers=False,manage_ticks=False)
            for patch,b in zip(boxes["boxes"],BATCHES):
                patch.set_facecolor(COLORS[b]);patch.set_alpha(.55)
            for j,(b,ys) in enumerate(zip(BATCHES,vals),1):
                ax.scatter(j+np.linspace(-.16,.16,len(ys)),ys,s=30,
                           color=COLORS[b],edgecolor="white",linewidth=.4,zorder=4)
            ax.set_title(name+" · "+title,loc="left",fontweight="bold")
            ax.set_xticks(range(1,5),[f"{b}\nn={GROUPS[b]}" for b in BATCHES])
            ax.set_ylabel("One log2 feature median per sample")
            style(ax)
    fig.suptitle("PG & PR · sample distribution medians by batch (alignment is by construction)",
                 fontsize=15,fontweight="bold")
    save(fig,"03_pg_pr_sample_medians_by_batch_before_after.png")

def yeast_check(info,pg_factor,pr_factor):
    rows,missing=load_yeast()
    if missing:raise ValueError("Yeast source has missing precursor/protein measurements")
    if {x["sample_id"] for x in rows}!={x["sample_id"] for x in info}:
        raise ValueError("Yeast sample IDs mismatch")
    byid={x["sample_id"]:x for x in rows}
    pf={r["sample_id"]:float(f) for r,f in zip(info,pg_factor)}
    rf={r["sample_id"]:float(f) for r,f in zip(info,pr_factor)}
    fields={"AADALLLK":("AADALLLK_sum_linear","PR"),
            "VNQIGTLSESIK":("VNQIGTLSESIK_sum_linear","PR"),
            "P00924_protein":("P00924_protein_linear","PG")}
    value_rows=[];out=[]
    for r in info:
        sid=r["sample_id"]
        x=byid[sid]
        for analyte,(field,source) in fields.items():
            original=float(x[field])
            factor=pf[sid] if source=="PG" else rf[sid]
            corrected=original*factor
            value_rows.append(dict(sample_id=sid,batch=r["batch"],
                yeast_analyte=analyte,factor_source=source,
                raw_intensity=original,applied_factor=factor,
                normalized_intensity=corrected))
    for b in ("All_40",*BATCHES):
        for analyte,(field,source) in fields.items():
            selected=[r for r in value_rows if r["yeast_analyte"]==analyte
                      and (b=="All_40" or r["batch"]==b)]
            v=np.array([r["raw_intensity"] for r in selected],float)
            u=np.array([r["normalized_intensity"] for r in selected],float)
            assert len(v)==(40 if b=="All_40" else GROUPS[b])
            out.append(dict(batch=b,analyte=analyte,
                yeast_factor_source=source,n_samples=len(v),
                raw_cv_pct=float(np.std(v,ddof=1)/np.mean(v)*100),
                normalized_cv_pct=float(np.std(u,ddof=1)/np.mean(u)*100),
                delta_cv_percentage_points=float(100*(np.std(u,ddof=1)/np.mean(u)-
                                                        np.std(v,ddof=1)/np.mean(v)))))
    return value_rows,out

def yeast_plot(rows):
    stages=("raw_cv_pct","normalized_cv_pct")
    fig,axes=plt.subplots(1,3,figsize=(17,6),sharey=True,layout="constrained")
    for ax,feat in zip(axes,("AADALLLK","VNQIGTLSESIK","P00924_protein")):
        data=[[next(r[col] for r in rows if r["batch"]==b and r["analyte"]==feat)
               for b in BATCHES] for col in stages]
        xx=np.arange(4)
        ax.bar(xx-.18,data[0],width=.35,color="#BBC6D1",label="Raw")
        ax.bar(xx+.18,data[1],width=.35,color="#3C82B6",label="Median normalized")
        for i,(raw,post) in enumerate(zip(*data)):
            ax.plot([i-.18,i+.18],[raw,post],color="#344B5E",alpha=.45,lw=1)
        ax.set_title(feat,fontsize=12,fontweight="bold")
        ax.set_xticks(xx,BATCHES)
        ax.set_ylabel("Within-batch yeast CV (%)" if feat=="AADALLLK" else "")
        ax.legend(frameon=False,fontsize=9)
        style(ax)
    fig.suptitle("Yeast QC reference · CV before/after (not used to calculate factors)",
                 fontsize=15,fontweight="bold")
    save(fig,"04_yeast_cv_before_after.png")

def cvplot(rows,kind):
    fig,ax=plt.subplots(figsize=(11,6),layout="constrained")
    width=.34;xx=np.arange(4)
    b0=[next(r for r in rows if r["matrix"]==kind and r["batch"]==b and r["stage"]=="raw")
        for b in BATCHES]
    b1=[next(r for r in rows if r["matrix"]==kind and r["batch"]==b and r["stage"]=="median_normalized")
        for b in BATCHES]
    for offset,data,color,label in ((-.18,b0,"#C2CAD4","Raw"),(.18,b1,"#397FB4","Normalized")):
        pos=xx+offset
        ax.bar(pos,[r["median_feature_cv_pct"] for r in data],width=.34,
               color=color,label=label,alpha=.92)
    ax.set_xticks(xx,[f"{b}\n(n={GROUPS[b]})" for b in BATCHES])
    ax.set_ylabel("Median feature CV (%) · original linear intensity")
    ax.set_title(f"{kind.upper()} within-batch feature CVs (descriptive)",
                 fontsize=14,fontweight="bold",loc="left")
    ax.legend(frameon=False)
    style(ax)
    save(fig,f"05_{kind}_batch_feature_cv_before_after.png")

def main():
    global ALL_META
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError("Refusing to overwrite prior results: "+str(OUT))
    ALL_META=meta()
    sids=[r["sample_id"] for r in ALL_META]
    output=[];all_samples=[];all_batches=[];all_cvs=[];factors={}
    arrays={}
    OUT.mkdir(parents=True,exist_ok=True)
    for name,path in (("pg",PG),("pr",PR)):
        header,rows,indices,log,excluded=load_matrix(path,sids)
        factor_rows,factor,median,ref=calc_factors(log,excluded,ALL_META)
        for r in factor_rows:r["matrix"]=name
        output.extend(factor_rows)
        normalized,per_sample,per_batch,per_cv=matrix_stats(log,factor,name,excluded,ALL_META)
        all_samples.extend(per_sample);all_batches.extend(per_batch);all_cvs.extend(per_cv)
        factors[name]=factor
        arrays[name]=(log,normalized)
        result=write_normalized(header,rows,indices,factor,
                f"NYD_{name}_matrix2_median_normalized.tsv")
        # Independent physical round-trip validation of every valid intensity
        verify_h,verify_rows,verify_indices,verify_log,_=load_matrix(result,sids)
        if verify_h!=header or len(verify_rows)!=len(rows):
            raise AssertionError("Normalized matrix schema/row count changed")
        orig_ok=np.isfinite(log)
        if not np.array_equal(np.isfinite(verify_log),orig_ok):
            raise AssertionError("Normalization changed detection/missingness")
        if not np.allclose(verify_log[orig_ok],normalized[orig_ok],rtol=0,atol=2e-10):
            raise AssertionError("Round-trip corrected intensity mismatch")
        print("NORMALIZED",name,"rows",len(rows),
              "excluded yeast features",int(excluded.sum()),
              "reference log2 median",ref,
              "factor range",float(factor.min()),float(factor.max()),flush=True)
        per_sample_plot(log,normalized,ALL_META,name)
    writecsv("median_normalization_sample_factors.csv",output)
    writecsv("sample_median_statistics.csv",all_samples)
    writecsv("batch_median_statistics.csv",all_batches)
    writecsv("batch_feature_cv_statistics.csv",all_cvs)
    batch_medians_plot([("PG",*arrays["pg"]),("PR",*arrays["pr"])])
    yeast_values,yeast_cv=yeast_check(ALL_META,factors["pg"],factors["pr"])
    writecsv("yeast_before_after_per_sample.csv",yeast_values)
    writecsv("yeast_cv_before_after.csv",yeast_cv)
    yeast_plot(yeast_cv)
    for name in ("pg","pr"):cvplot(all_cvs,name)
    manifest=dict(
        input_sources={"PG":str(PG.relative_to(ROOT)),"PR":str(PR.relative_to(ROOT)),
                       "META":str(META.relative_to(ROOT))},
        input_sha256={"PG":sha(PG),"PR":sha(PR),"META":sha(META),
                      "YEAST":sha(BASE/"YEAST.xlsx")},
        normalization="per-sample median of log2 finite-positive intensities, centered on median of all 40 sample medians",
        sample_count=40,batch_sizes=GROUPS,
        normalization_factors_derived_from="respective endogenous PG/PR matrix, excluding yeast P00924/ENO1_YEAST rows",
        yeast_used_for_normalization=False,
        sample_vs_batch_reference="global 40-sample median, NOT batch-specific",
        no_imputation=True,
        pg_pr_calibrated_independently=True,
        rounding_output="15 significant digits",
        figures=[
            "01_pg_40sample_boxplots_before_after.png",
            "01_pr_40sample_boxplots_before_after.png",
            "03_pg_pr_sample_medians_by_batch_before_after.png",
            "04_yeast_cv_before_after.png",
            "05_pg_batch_feature_cv_before_after.png",
            "05_pr_batch_feature_cv_before_after.png",
        ]
    )
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    ref=[
        "# NYD 40 samples: median normalization (MacCoss paper-inspired)","",
        "**Method:** per-sample median normalization performed directly on 40 biological samples.",
        "Paper reference: Tsantilas et al. (2026), *Circulating extracellular vesicles in",
        "plasma carry accessible molecular signatures of aging in mice*, MacCoss group.",
        "doi:10.64898/2026.07.10.737625.",
        "The paper reports median-normalized log intensities but does not specify its exact",
        "software-scale factor formula in the supplied text. Our explicit reproducible",
        "implementation is a paper-inspired choice, not a claim of exact software replication.","",
        "For each PG/PR matrix separately, m_s=median(log2(I_positive)),",
        "m_ref=median_s(m_s); I_corrected=I_raw*2^(m_ref-m_s).",
        "40 sample columns joined to NYD_sampleinfo2.xlsx by exact sample ID.",
        "All samples are treated individually, not separately calibrated by batch.",
        "The yeast P00924/ENO1_YEAST rows are excluded from computing factors",
        "when present, and used only as hold-out signal diagnostics.",
        "Raw nonpositive and missing entries remain exactly as observed;",
        "no imputation, no replacement of source files.","",
        "## Outputs","",
        "- NYD_pg_matrix2_median_normalized.tsv: corrected PG protein-group intensity matrix.",
        "- NYD_pr_matrix2_median_normalized.tsv: corrected PR feature intensity matrix.",
        "- median_normalization_sample_factors.csv: 80 sample/matrix factors (PG and PR separately).",
        "- sample_median_statistics.csv, batch_median_statistics.csv: distributions before/after.",
        "- batch_feature_cv_statistics.csv: descriptive within-batch feature CV changes.",
        "- yeast_before_after_per_sample.csv, yeast_cv_before_after.csv: independent yeast diagnostics.",
        "- run_manifest.json: original SHA256 and exact mathematical method.","",
        "## Corrected PG and PR boxplots by sample","",
        "![PG 40 samples before and after](01_pg_40sample_boxplots_before_after.png)","",
        "![PR 40 samples before and after](01_pr_40sample_boxplots_before_after.png)","",
        "## Batch distribution of sample-level medians","",
        "![PG PR batch median distributions](03_pg_pr_sample_medians_by_batch_before_after.png)","",
        "## Yeast holdout CV before and after","",
        "Peptide reference features use the **PR** correction factors while the one",
        "P00924 protein reference uses **PG** correction factors. Original peptide",
        "charge states are summed within each unique peptide as in prior yeast CV.",
        "CV is sample SD / mean * 100 calculated on original-linear signals.",
        "The yeast is not an input to the normalization factor.","",
        "![Yeast CV before and after](04_yeast_cv_before_after.png)","",
        "### Batch-specific yeast CV details","",
        "| Batch | Feature | Raw CV (%) | Post-normalization CV (%) | Delta pp |",
        "| --- | --- | ---: | ---: | ---: |"]
    for r in yeast_cv:
        if r["batch"]=="All_40":continue
        ref.append(f'| {r["batch"]} | {r["analyte"]} | {r["raw_cv_pct"]:.2f} | '
                   f'{r["normalized_cv_pct"]:.2f} | {r["delta_cv_percentage_points"]:+.2f} |')
    ref += ["","## Within-batch PG/PR feature CV comparisons","",
        "These represent variation of individual biological sample measurements, NOT pure",
        "technical precision. The batch sample-median contraction is guaranteed",
        "by construction and cannot alone prove quality improvement.","",
        "![PG batch feature CV](05_pg_batch_feature_cv_before_after.png)","",
        "![PR batch feature CV](05_pr_batch_feature_cv_before_after.png)","",
        "## Interpretation and limitations","",
        "Median normalization assumes most measured proteomic features have",
        "comparable central abundance across samples. Genuine global shifts in",
        "abundance, systematic missingness, and differing biological compositions",
        "can bias normalization; this is especially important for samples with",
        "unusually low detection rates. The earlier PG quality analysis identified",
        "a low-detection T003 sample (MAP26092326E01) that warrants particular review.",
        "Median normalization does not by itself prove removal of a batch effect.",
        "This analysis has no repeated pooled QC to independently establish",
        "batch correction quality; yeast is examined as an additional holdout.",
        ""]
    (OUT/"README.md").write_text("\n".join(ref),encoding="utf-8")
    print("YEAST_CHECK",json.dumps(
      [{"b":r["batch"],"feat":r["analyte"],"before":round(r["raw_cv_pct"],3),
        "after":round(r["normalized_cv_pct"],3)} for r in yeast_cv],
      ensure_ascii=False),flush=True)
    print("SUCCESS: calibrated 40 PG and PR sample columns, reproducible matrices and yeast holdout checks",flush=True)

if __name__=="__main__":
    main()
