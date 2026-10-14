#!/usr/bin/env python3
"""Trial 2 on reuploaded pg_matrix.tsv: Raw -> Factor1 -> Factor1+Factor2.

Method a and b both included. Show three stages TOGETHER in every chart,
with an all-groups and HC/S-only view for each method. PNG only.
All factors and PCA computed in log2; ordinary CV calculated on linear intensities.
Neither LM Internal-QC nor pooled QC is corrected. All source protein rows kept.
"""
import json
import math
import os
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pathlib import Path
import pacs_trial2 as old
from pacs_common import ROOT, PG, ANN, finite_positive, cv_by_feature
from yeast_protein_qc import GROUPS, COLORS, metadata

OUT = ROOT / "results"
old.OUT = OUT
# Trial 2 uses the reuploaded, metadata-filtered 2026-10-13 HELP workbook.
# Override only this experiment's HELP source; historical Trial 1 remains intact.
old.HELP = ROOT / "LH_ratio_tables_2026.10.13—ALL-HELP.xlsx"
METHODS = ("a", "b")
STAGES = ("raw", "factor1", "factor1_factor2")
TITLES = ("Raw", "Factor1", "Factor1 + Factor2")
VIEWS = {"all_groups": list(GROUPS), "study_only": ["HC", "S"]}
PALETTE = dict(zip(GROUPS, COLORS))
STUDY = {"HC", "S"}


def save_csv(path, frame):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(frame).to_csv(path, index=False)


def save_fig(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=210, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def log2_matrix(a):
    return np.log2(np.where(finite_positive(a), a, np.nan))


def select(samples, labels, groups):
    ids = [s for group in groups for s in samples if labels[s] == group]
    return ids, [samples.index(s) for s in ids]


def prepare_input():
    labels = metadata()
    df = pd.read_csv(PG, sep="\t", dtype=str, encoding="utf-8-sig",
                     keep_default_na=False, low_memory=False)
    if list(df.columns[:5]) != ANN:
        raise ValueError("Invalid protein annotation column names")
    old_columns = list(df.columns[5:])
    if len(set(old_columns)) != len(old_columns):
        raise ValueError("Duplicate source sample column")
    if not set(labels).issubset(set(old_columns)):
        raise ValueError("Metadata sample columns missing from pg_matrix: "
                         + str(sorted(set(labels) - set(old_columns))))
    retained = [x for x in old_columns if x in labels]
    if len(retained) != 75:
        raise ValueError("Expected 75 retained sample columns")
    if {g: sum(labels[s] == g for s in retained) for g in GROUPS} != {
        "Internal-QC": 6, "pooled QC": 15, "HC": 27, "S": 27
    }:
        raise ValueError("Study and QC group sizes differ from expected")
    data = df[ANN + retained].copy()
    data.to_csv(OUT / "pg_matrix_filtered.tsv", sep="\t", index=False)
    audit = {
        "original_pg_matrix": "pg_matrix.tsv",
        "input_protein_rows": int(len(df)),
        "output_protein_rows": int(len(data)),
        "input_sample_columns": len(old_columns),
        "retained_sample_columns": len(retained),
        "removed_sample_columns": [s for s in old_columns if s not in labels],
        "groups": {g: sum(labels[s] == g for s in retained) for g in GROUPS},
    }
    (OUT / "screening_report.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    raw = data[retained].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    if len(raw) < 100:
        raise ValueError("Unexpectedly small protein matrix")
    return data, labels, retained, raw, audit


def calculate_factors(samples, labels):
    fa, fb, ra, rb, b_minus_a = old.compute_yeast_factors(samples, labels)
    f2, help_medians, counts, ref_help = old.compute_help_factor2(samples, labels)
    factors, rows = {}, []
    for s in samples:
        is_study = labels[s] in STUDY
        factors[s] = {
            "a": fa[s] if is_study else 0,
            "b": fb[s] if is_study else 0,
            "factor2": f2[s] if is_study else 0,
        }
        rows.append({
            "sample_id": s, "group": labels[s],
            "applied_factor1_a": factors[s]["a"],
            "applied_factor1_b": factors[s]["b"],
            "applied_factor2": factors[s]["factor2"],
            "candidate_factor2": f2[s],
            "sample_HELP_mR_log2": help_medians[s],
            "reference_HELP_mR_log2": ref_help,
            "n_valid_HELP_ratios": counts[s],
        })
    save_csv(OUT / "sample_factors.csv", rows)
    return factors, {"ref_ye_a": ra, "ref_ye_b": rb,
                     "ref_ye_b_minus_a": b_minus_a, "ref_HELP_log2": ref_help}


def correction_stages(raw, samples, labels, factors, method):
    data = {"raw": raw.copy()}
    totals = {"raw": np.zeros(len(samples))}
    original_good = finite_positive(raw)
    is_study = np.asarray([labels[s] in STUDY for s in samples])
    valid_study = original_good & is_study[None, :]
    lograw = log2_matrix(raw)
    for stage in STAGES[1:]:
        o = np.asarray([
            0 if labels[s] not in STUDY else
            factors[s][method] + (factors[s]["factor2"] if stage == "factor1_factor2" else 0)
            for s in samples
        ])
        totals[stage] = o
        new = raw.copy()
        candidate = np.exp2(lograw + o[None, :])
        new[valid_study] = candidate[valid_study]
        if not np.array_equal(finite_positive(new), original_good):
            raise AssertionError("Protein detection mask changed")
        for j, s in enumerate(samples):
            if labels[s] not in STUDY and not np.allclose(
                new[:,j], raw[:,j], rtol=0, atol=0, equal_nan=True
            ):
                raise AssertionError("QC protein intensities changed")
        data[stage] = new
    return data, totals


def write_matrix(path, input_df, samples, labels, matrix, log2_scale=False):
    out = input_df[ANN].copy()
    values = log2_matrix(matrix) if log2_scale else matrix
    for j,s in enumerate(samples):
        if not log2_scale and labels[s] not in STUDY:
            out[s] = input_df[s]
        else:
            old = input_df[s].to_numpy()
            out[s] = [
                format(v, ".16g") if math.isfinite(v) and (log2_scale or v > 0)
                else ("" if log2_scale else old[i])
                for i,v in enumerate(values[:,j])
            ]
    out.to_csv(path, sep="\t", index=False)


def calculate_cv(mats, df, samples, labels, method, folder):
    per_group, per_protein, chart_data = [], [], {}
    for g in GROUPS:
        idx = [j for j,s in enumerate(samples) if labels[s] == g]
        values = []
        count = None
        for stage in STAGES:
            cv, n, min_n = cv_by_feature(mats[stage][:,idx])
            if count is not None and not np.array_equal(n,count):
                raise AssertionError("Detection counts changed")
            count = n
            values.append(cv)
        mask = np.logical_and.reduce([np.isfinite(v) for v in values])
        valid = [v[mask] for v in values]
        if not len(valid[0]):
            raise AssertionError("No eligible proteins for "+g)
        if g in ("Internal-QC","pooled QC"):
            for v in valid[1:]:
                if not np.allclose(valid[0],v,atol=0,rtol=0,equal_nan=True):
                    raise AssertionError("CV changed for uncorrected QC group "+g)
        chart_data[g] = valid
        for stage,v in zip(STAGES,valid):
            per_group.append({
                "method":method, "stage":stage, "group":g,
                "sample_n":len(idx), "cv_min_detected_n":min_n,
                "eligible_proteins":int(mask.sum()),
                "median_cv_pct":float(np.median(v)),
                "median_paired_change_pp":float(np.median(v-valid[0])),
                "fraction_proteins_improved_pct":float(np.mean(v<valid[0])*100)
            })
        per_protein.append(pd.DataFrame({
            "row_1based":np.arange(1,len(df)+1),
            "Protein.Group":df["Protein.Group"].to_numpy(),
            "Protein.Names":df["Protein.Names"].to_numpy(),
            "group":g,"detected_n":count,
            "cv_raw_pct":values[0],"cv_factor1_pct":values[1],
            "cv_factor1_factor2_pct":values[2],
        }))
    save_csv(folder/"group_cv_stages.csv",per_group)
    pd.concat(per_protein,ignore_index=True).to_csv(
        folder/"per_protein_cv_stages.csv",index=False)
    return per_group, chart_data


def plot_sample_boxes(mats,samples,labels,groups,folder):
    selected,idx=select(samples,labels,groups)
    arrs=[log2_matrix(mats[stage][:,idx]) for stage in STAGES]
    bottom=min(np.nanpercentile(a,2) for a in arrs)-.6
    top=max(np.nanpercentile(a,98) for a in arrs)+.6
    fig,axes=plt.subplots(3,1,sharex=True,figsize=(max(14,len(selected)*.22),12.6),
                          layout="constrained")
    for ax,a,title in zip(axes,arrs,TITLES):
        b=ax.boxplot([a[:,i][np.isfinite(a[:,i])] for i in range(a.shape[1])],
                     patch_artist=True,showfliers=False,widths=.6)
        for p,s in zip(b["boxes"],selected):
            p.set_facecolor(PALETTE[labels[s]])
            p.set_alpha(.55)
        ax.set_title(title);ax.set_ylabel("Protein log2 intensity")
        ax.set_ylim(bottom,top)
        n=0
        for g in groups[:-1]:
            n+=sum(labels[s]==g for s in selected)
            ax.axvline(n+.5,c="#8899aa",ls=":",lw=.8)
    axes[-1].set_xticks(range(1,len(selected)+1),selected,fontsize=6,rotation=90)
    fig.suptitle("Per-sample protein boxplots: Raw / Factor1 / Factor1+Factor2")
    save_fig(fig,folder/"01_sample_boxplots_3stages.png")


def plot_group_medians(mats,samples,labels,groups,view,folder):
    selected,idx=select(samples,labels,groups)
    values=[np.nanmedian(log2_matrix(mats[stage][:,idx]),axis=0) for stage in STAGES]
    bottom=min(np.min(v) for v in values)-.4
    top=max(np.max(v) for v in values)+.4
    fig,axes=plt.subplots(1,3,sharey=True,figsize=(max(13,len(groups)*3.6),5.8),
                          layout="constrained")
    rows=[]
    for ax,stage,v,title in zip(axes,STAGES,values,TITLES):
        bins=[[v[i] for i,s in enumerate(selected) if labels[s]==g] for g in groups]
        boxes=ax.boxplot(bins,patch_artist=True,showfliers=False,widths=.45)
        for p,g in zip(boxes["boxes"],groups):
            p.set_facecolor(PALETTE[g]);p.set_alpha(.4)
        for i,(g,ys) in enumerate(zip(groups,bins),1):
            ax.scatter(i+np.linspace(-.15,.15,len(ys)),ys,
                       s=17,color=PALETTE[g],alpha=.7)
        ax.set_xticks(range(1,len(groups)+1),groups,rotation=12)
        ax.set_ylim(bottom,top);ax.set_title(title)
        ax.set_ylabel("Per-sample median log2 protein intensity")
        rows.extend({
            "view":view,"stage":stage,"sample_id":s,
            "group":labels[s],"median_log2_intensity":float(x)
        } for s,x in zip(selected,v))
    fig.suptitle("Group comparison: Raw / Factor1 / Factor1+Factor2")
    save_fig(fig,folder/"02_group_boxes_3stages.png")
    return rows


def plot_cv(values,groups,folder):
    fig,axes=plt.subplots(1,3,sharey=True,figsize=(max(14,3.7*len(groups)),6.5),
                          layout="constrained")
    for ax,stage_index,title in zip(axes,range(3),TITLES):
        for i,g in enumerate(groups,1):
            original=values[g][stage_index]
            clipped=np.minimum(original,200)
            if len(clipped)>1 and np.ptp(clipped)>0:
                v=ax.violinplot([clipped],positions=[i],widths=.75,showextrema=False)
                for p in v["bodies"]:
                    p.set_facecolor(PALETTE[g]);p.set_alpha(.38)
            b=ax.boxplot([clipped],positions=[i],widths=.17,
                         patch_artist=True,showfliers=False,manage_ticks=False)
            for p in b["boxes"]:
                p.set_facecolor(PALETTE[g]);p.set_alpha(.65)
            rng=np.random.default_rng(20261012+i)
            ax.scatter(i+rng.uniform(-.2,.2,len(clipped)),clipped,s=3,
                       color=PALETTE[g],alpha=.045,rasterized=True,linewidth=0)
            ax.text(i,207,f"{np.median(original):.2f}%",ha="center",fontsize=8)
        ax.set_xticks(range(1,len(groups)+1),groups,rotation=12)
        ax.set_xlim(.4,len(groups)+.6);ax.set_ylim(0,222)
        ax.set_title(title);ax.set_ylabel("Protein CV of linear intensity (%)")
        ax.grid(axis="y",alpha=.16)
    fig.suptitle("Per-protein CV (violin + box + scatter): 3 sequential stages")
    axes[-1].text(.98,.015,">200% capped for display only",
                  transform=axes[-1].transAxes,fontsize=8,ha="right")
    save_fig(fig,folder/"04_cv_3stages.png")


def plot_cv_medians(groupcv,groups,folder):
    df=pd.DataFrame(groupcv)
    fig,ax=plt.subplots(figsize=(max(10,2.6*len(groups)),5.5),layout="constrained")
    x=np.arange(len(groups))
    for i,(stage,title) in enumerate(zip(STAGES,TITLES)):
        med=[df[(df.stage==stage)&(df.group==g)].median_cv_pct.iloc[0] for g in groups]
        b=ax.bar(x+(i-1)*.24,med,width=.23,label=title)
        ax.bar_label(b,fmt="%.1f",fontsize=8,padding=2)
    ax.set_xticks(x,groups);ax.set_ylabel("Median protein CV (%)")
    ax.set_title("CV median across Raw / Factor1 / Factor1+Factor2")
    ax.legend(frameon=False);ax.grid(axis="y",alpha=.15)
    save_fig(fig,folder/"05_cv_median_comparison.png")


def plot_offsets(samples,labels,groups,factors,method,folder):
    selected,_=select(samples,labels,groups)
    x=np.arange(len(selected))
    f1=np.array([factors[s][method] for s in selected])
    f2=np.array([factors[s]["factor2"] for s in selected])
    fig,axes=plt.subplots(2,1,figsize=(max(12,len(selected)*.23),8.2),
                          sharex=True,layout="constrained")
    axes[0].bar(x-.2,f1,width=.38,label="factor1")
    axes[0].bar(x+.2,f2,width=.38,label="factor2")
    axes[0].legend(frameon=False);axes[0].axhline(0,c="black",lw=.8)
    axes[0].set_ylabel("Individual log2 offsets")
    axes[1].bar(x,f1+f2,color=[PALETTE[labels[s]] for s in selected])
    axes[1].axhline(0,c="black",lw=.8)
    axes[1].set_ylabel("Combined log2 offset")
    axes[1].set_xticks(x,selected,rotation=90,fontsize=6)
    fig.suptitle("PACS Trial 2: yeast and HELP offsets (both QC groups zero)")
    save_fig(fig,folder/"06_factor_components.png")


def pca_three(mats,samples,labels,groups,view,folder):
    selected,idx=select(samples,labels,groups)
    min_n=math.ceil(.70*len(selected))
    common=finite_positive(mats["raw"][:,idx]).sum(axis=1)>=min_n
    logged=[]
    for stage in STAGES:
        a=log2_matrix(mats[stage][common][:,idx])
        m=np.nanmedian(a,axis=1)
        logged.append(np.where(np.isfinite(a),a,m[:,None]))
    varying=np.logical_and.reduce([np.std(a,axis=1)>1e-10 for a in logged])
    if varying.sum()<100:
        raise ValueError("Insufficient PCA protein features")
    fig,axes=plt.subplots(1,3,figsize=(20,6.2),layout="constrained")
    summary,scores=[],[]
    for ax,stage,title,a in zip(axes,STAGES,TITLES,logged):
        model=PCA(n_components=2,svd_solver="full")
        z=model.fit_transform(a[varying].T)
        explained=100*model.explained_variance_ratio_
        summary.append({
            "view":view,"stage":stage,"n_samples":len(selected),
            "n_pca_proteins":int(varying.sum()),
            "PC1_pct":float(explained[0]),"PC2_pct":float(explained[1])
        })
        scores.extend({
            "view":view,"stage":stage,"sample_id":s,"group":labels[s],
            "PC1":float(pt[0]),"PC2":float(pt[1]),
            "PC1_explained_pct":float(explained[0]),
            "PC2_explained_pct":float(explained[1])
        } for s,pt in zip(selected,z))
        for g in groups:
            p=np.array([v for s,v in zip(selected,z) if labels[s]==g])
            ax.scatter(p[:,0],p[:,1],color=PALETTE[g],s=42,
                       alpha=.8,edgecolors="white",linewidth=.5,
                       label=f"{g} (n={len(p)})")
        ax.set_xlabel(f"PC1 ({explained[0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({explained[1]:.1f}%)")
        ax.axhline(0,c="#adb8c5",lw=.6)
        ax.axvline(0,c="#adb8c5",lw=.6)
        ax.set_title(title);ax.legend(frameon=False,fontsize=8)
    fig.suptitle(f"PCA 3 stages (same {varying.sum()} proteins within view)")
    save_fig(fig,folder/"03_pca_3stages.png")
    return scores,summary


def make_report(info,refs,cv,pca):
    lines=[
        "# PACS Trial 2, new uploaded protein matrix (stepwise)", "",
        f"Previously generated results/ has been archived at {os.environ.get('PACS_RESULTS_ARCHIVE', 'results_vN')}/.",
        "The new results/ is an independent Trial 2 rerun on the filtered pg_matrix.tsv.",
        "The root pg_matrix.tsv and older archive files are not changed.", "",
        "## Input provenance", "",
        f"- Protein source: pg_matrix.tsv; sha256={hashlib.sha256(PG.read_bytes()).hexdigest()}.",
        f"- HELP source: {old.HELP.name}; sha256={hashlib.sha256(old.HELP.read_bytes()).hexdigest()}.",
        f"- Yeast source: {old.YEAST.name}; sha256={hashlib.sha256(old.YEAST.read_bytes()).hexdigest()}.",
        "Both peptide workbooks were metadata-filtered to 75 samples before this run.",
        "Previous trial outputs and filtered source workbooks were not modified.",
        "",
        "## Sample screening", "",
        f"- New protein rows: {info['input_protein_rows']}.",
        f"- Original sample columns: {info['input_sample_columns']}.",
        f"- Kept: {info['retained_sample_columns']} (Internal-QC 6, pooled QC 15, HC 27, S 27).",
        f"- Removed sample columns: {len(info['removed_sample_columns'])}.",
        "- No imputation, filtering of protein rows, normalization or transformation of input raw data.",
        "",
        "## Definitions (log2 except CV)", "",
        "For each Study sample with yeast peptide log2 intensities y1 and y2:",
        "- factor1a = median_study((y1+y2)/2) - (y1+y2)/2.",
        "- factor1b = mean(median_study(y1), median_study(y2)) - (y1+y2)/2.",
        "- factor2 = median_of_pooled_QC_sample_medians(log2 HELP L/H) - sample_median(log2 HELP L/H).",
        "- Raw: original protein intensity.",
        "- After factor1: log2 corrected protein = log2 raw protein + factor1.",
        "- After factor1+factor2: log2 corrected protein = log2 raw + factor1 + factor2.",
        "- Only HC and S are corrected. Internal-QC and pooled QC are entirely unchanged.",
        "- CV uses unlogged positive intensity, SD(ddof=1)/arithmetic mean*100%,",
        "  detection >= max(3,ceil(0.7*group n)), minimum mean 1e-8.",
        "- PCA and boxplots use log2 intensities. PCA-only per-protein median imputation;",
        "  each view selects >=70% detected proteins, same set across stages, centered not z-scored.",
        "- PCA is separately fitted at each stage, so PC coordinates are not aligned between stages.",
        "",
        "## References", "",
        f"- Yeast ref a: {refs['ref_ye_a']:.12g}.",
        f"- Yeast ref b: {refs['ref_ye_b']:.12g}.",
        f"- Common factor1 b minus a: {refs['ref_ye_b_minus_a']:.12g}.",
        f"- Pooled QC HELP reference log2: {refs['ref_HELP_log2']:.12g}.",
        "",
        "## Per-group median protein CV (%)", "",
        "| Method | Stage | Internal-QC | pooled QC | HC | S |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    df=pd.DataFrame(cv)
    for method in METHODS:
        for stage,title in zip(STAGES,TITLES):
            vals=[df[(df.method==method)&(df.stage==stage)&(df.group==g)]
                  .median_cv_pct.iloc[0] for g in GROUPS]
            lines.append("| "+method+" | "+title+" | "
                         +" | ".join(f"{x:.2f}%" for x in vals)+" |")
    lines+=["", "## Joint figures: Raw / Factor1 / Factor1+Factor2", "",
            "Two yeast reference methods x two plotting views = 4 figure sets.",
            "Each set contains six PNG charts (24 PNGs total), no SVG.", ""]
    chartfiles=[
        ("01_sample_boxplots_3stages.png","Sample boxplots"),
        ("02_group_boxes_3stages.png","Group sample-median boxplots"),
        ("03_pca_3stages.png","PCA"),
        ("04_cv_3stages.png","Per-protein CV violin/box/scatter"),
        ("05_cv_median_comparison.png","Median CV changes"),
        ("06_factor_components.png","Factors"),
    ]
    for method in METHODS:
        for view in VIEWS:
            lines.extend(["### Method "+method+" / "+view,""])
            for name,desc in chartfiles:
                lines.extend(["!["+desc+"](method_"+method+"/figures/"+view+"/"+name+")",""])
    lines+=["## PCA explained variance", "",
            "| Method | View | Stage | Proteins | PC1 | PC2 |",
            "| --- | --- | --- | ---: | ---: | ---: |"]
    for row in pca:
        lines.append(f"| {row['method']} | {row['view']} | {row['stage']} | "
                     f"{row['n_pca_proteins']} | {row['PC1_pct']:.2f}% | "
                     f"{row['PC2_pct']:.2f}% |")
    lines+=["", "## Data products", "",
            "- pg_matrix_filtered.tsv: screened source, unmodified protein intensities.",
            "- sample_factors.csv: per-sample a/b factor1, factor2, and HELP ratio statistics.",
            "- yeast_factor_details.csv: peptide-level log2 intensities and both reference methods.",
            "- all_stage_cv_summary.csv and all_stage_pca_summary.csv: method x stage x group/view summaries.",
            "- method_a/ and method_b/: linear and log2 factor1/factor1+factor2 corrected matrices;",
            "  stage_offsets.csv, per_protein_cv_stages.csv, group_cv_stages.csv,",
            "  pca_scores.csv, pca_summary.csv, sample_log2_medians.csv and PNG figures.",
            "",
            "## Caution", "",
            "CV for both QC groups must remain unchanged by design, and cannot validate",
            "normalization effectiveness. HC/S CV includes biological variability.",
            "Methods a/b differ by a common log2 offset and are expected to have",
            "identical within-group CV; their all-groups PCA can differ because QC remains fixed.",
            "PCA explained variance increasing does not imply better normalization.",""]
    (OUT/"README.md").write_text("\n".join(lines),encoding="utf-8")


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    df,labels,samples,raw,audit=prepare_input()
    factors,refs=calculate_factors(samples,labels)
    cv_rows=[];pca_rows=[]
    for method in METHODS:
        folder=OUT/("method_"+method)
        folder.mkdir(parents=True,exist_ok=True)
        mats,offsets=correction_stages(raw,samples,labels,factors,method)
        for stage in STAGES[1:]:
            write_matrix(folder/("pg_matrix_"+stage+"_linear.tsv"),
                         df,samples,labels,mats[stage])
            write_matrix(folder/("pg_matrix_"+stage+"_log2.tsv"),
                         df,samples,labels,mats[stage],log2_scale=True)
        save_csv(folder/"stage_offsets.csv",[
            {"sample_id":s,"group":labels[s],"stage":stage,
             "applied_log2_offset":float(offset),
             "applied_linear_multiplier":float(np.exp2(offset))}
            for stage in STAGES for s,offset in zip(samples,offsets[stage])
        ])
        cvsummary,cvdata=calculate_cv(mats,df,samples,labels,method,folder)
        cv_rows.extend(cvsummary)
        ps=[];pca_summary=[];medians=[]
        for view,groups in VIEWS.items():
            figures=folder/"figures"/view
            plot_sample_boxes(mats,samples,labels,groups,figures)
            medians+=plot_group_medians(mats,samples,labels,groups,view,figures)
            scores,pinfo=pca_three(mats,samples,labels,groups,view,figures)
            ps+=scores
            pca_summary+=pinfo
            pca_rows.extend([{"method":method,**r} for r in pinfo])
            plot_cv(cvdata,groups,figures)
            plot_cv_medians(cvsummary,groups,figures)
            plot_offsets(samples,labels,groups,factors,method,figures)
        save_csv(folder/"pca_scores.csv",ps)
        save_csv(folder/"pca_summary.csv",pca_summary)
        save_csv(folder/"sample_log2_medians.csv",medians)
    save_csv(OUT/"all_stage_cv_summary.csv",cv_rows)
    save_csv(OUT/"all_stage_pca_summary.csv",pca_rows)
    cvdf=pd.DataFrame(cv_rows)
    for stage in STAGES[1:]:
        for g in GROUPS:
            vals=[cvdf[(cvdf.method==m)&(cvdf.stage==stage)&(cvdf.group==g)]
                  .median_cv_pct.iloc[0] for m in METHODS]
            if not np.allclose(vals[0],vals[1],rtol=0,atol=1e-8):
                raise AssertionError("a and b protein CVs should match")
    make_report(audit,refs,cv_rows,pca_rows)
    print("SUCCESS: stepwise PACS Trial 2, metadata screening and joint charts")
    print("protein rows",len(df),"sample columns",len(samples),
          "removed columns",len(audit["removed_sample_columns"]))
    print("references",refs)
    print(cvdf.pivot(index=["method","stage"],columns="group",
                     values="median_cv_pct").to_string())

if __name__=="__main__":
    main()
