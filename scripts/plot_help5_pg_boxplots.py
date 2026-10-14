#!/usr/bin/env python3
"""Visualize five HELP-source proteins against all raw PG intensities by metadata group."""
from pathlib import Path
import json
import math
import hashlib
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from yeast_protein_qc import ROOT, GROUPS, COLORS, metadata

PG=ROOT/"pg_matrix.tsv"
OUT=ROOT/"results"/"pg_help5_group_boxplots"
ANNOT=["Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description"]
TARGETS=[
    ("P02765","FETUA_HUMAN","AHSG","FSVVYAK"),
    ("P00751","CFAB_HUMAN","CFB","STGSWSTLK"),
    ("P00734","THRB_HUMAN","F2","VTGWGNLK"),
    ("P02748","CO9_HUMAN","C9","VVEESELAR"),
    ("P07225","PROS_HUMAN","PROS1","VYFAGFPR"),
]
G_COL=dict(zip(GROUPS,COLORS))
T_COL=dict(AHSG="#CC5148",CFB="#297BB9",F2="#8652AE",C9="#238A69",PROS1="#D58B28")
MARK=dict(AHSG="D",CFB="o",F2="^",C9="s",PROS1="P")

def tokens(s):
    return {v for v in re.split(r"[;,|\s]+",str(s)) if v}

def load():
    meta=metadata()
    pg=pd.read_csv(PG,sep="\t",dtype=str,keep_default_na=False,encoding="utf-8-sig",low_memory=False)
    if pg.columns[:5].tolist()!=ANNOT or len(pg)!=6675:raise ValueError("Unexpected PG schema or row count")
    samples=pg.columns[5:].tolist()
    if len(samples)!=75 or len(set(samples))!=75 or set(samples)!=set(meta):
        raise ValueError("PG and metadata sample IDs do not match (75 required)")
    groups={g:sum(meta[s]==g for s in samples) for g in GROUPS}
    if groups!={"Internal-QC":6,"pooled QC":15,"HC":27,"S":27}:
        raise ValueError(f"Unexpected groups {groups}")
    raw=pg[samples].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    valid=np.isfinite(raw)&(raw>0)
    logged=np.full(raw.shape,np.nan)
    logged[valid]=np.log2(raw[valid])
    return pg,meta,samples,raw,logged,groups

def map_targets(pg):
    result=[]
    # Strict exact accession token matching; no substring matching.
    for accession,protein,gene,peptide in TARGETS:
        selected=[]
        for i in range(len(pg)):
            row=pg.iloc[i]
            if accession in tokens(row["Protein.Group"])|tokens(row["Protein.Ids"]):
                selected.append(i)
        if len(selected)!=1:
            raise ValueError(f"Target {gene} / {accession} mapped to {len(selected)} PG rows: {selected}")
        i=selected[0]
        r=pg.iloc[i]
        result.append({
            "accession":accession,"protein_name_expected":protein,"gene":gene,
            "help_peptide":peptide,"pg_data_row_1based":i+1,"pg_tsv_line_1based":i+2,
            **{k:r[k] for k in ANNOT},
            "protein_name_match":protein in tokens(r["Protein.Names"]),
            "gene_name_match":gene in tokens(r["Genes"]),
        })
    return result

def save(fig,name):
    fig.savefig(OUT/name,dpi=210,facecolor="white",bbox_inches="tight")
    plt.close(fig)

def box(ax,data,at,colors,width=0.62):
    prepared=[np.asarray(d)[np.isfinite(d)] for d in data]
    boxes=ax.boxplot(prepared,positions=at,widths=width,patch_artist=True,
                     showfliers=False,manage_ticks=False,
                     medianprops=dict(color="#233342",lw=1.7),
                     boxprops=dict(color="#7B8A97",lw=.8),
                     whiskerprops=dict(color="#7B8A97",lw=.8),
                     capprops=dict(color="#7B8A97",lw=.8))
    for patch,color in zip(boxes["boxes"],colors):
        patch.set_facecolor(color)
        patch.set_alpha(.77)

def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y",alpha=.16,color="#687A8C")
    ax.set_axisbelow(True)

def group_overview(log,meta,samples,targets):
    idx=[[i for i,s in enumerate(samples) if meta[s]==g] for g in GROUPS]
    data=[log[:,j].ravel() for j in idx]
    fig,ax=plt.subplots(figsize=(12.7,7.3),layout="constrained")
    box(ax,data,range(1,5),["#CED6DE"]*4,.75)
    offs=np.linspace(-.25,.25,5)
    for gi,subset in enumerate(idx,1):
        for delta,t in zip(offs,targets):
            values=log[t["pg_data_row_1based"]-1,subset]
            values=values[np.isfinite(values)]
            if not len(values):continue
            low,mid,high=np.percentile(values,[25,50,75])
            color=T_COL[t["gene"]]
            ax.plot([gi+delta]*2,[low,high],color=color,lw=2.5,zorder=5)
            ax.scatter(gi+delta,mid,color=color,s=92,zorder=6,
                       marker=MARK[t["gene"]],edgecolors="white",linewidths=.7)
    ax.set_xticks(range(1,5),[f"{g}\n(n={len(j)})" for g,j in zip(GROUPS,idx)])
    ax.set_xlim(.5,4.5)
    ax.set_ylabel("log2 raw PG protein intensity (positive measurements)")
    ax.set_title("Raw PG distributions by sample group\n"
                 "All proteins: grey boxplots  |  HELP-source proteins: median + IQR",
                 fontsize=14,fontweight="bold",loc="left")
    legend=[Patch(facecolor="#CED6DE",label="All PG protein × sample measurements")]
    legend.extend(Line2D([],[],marker=MARK[t["gene"]],color=T_COL[t["gene"]],
                markersize=9,label=f'{t["gene"]} / {t["accession"]}')
                for t in targets)
    ax.legend(handles=legend,loc="upper left",bbox_to_anchor=(1.01,1),frameon=False)
    style(ax)
    save(fig,"01_group_background_boxplots_highlight_5proteins.png")

def target_panels(log,meta,samples,targets):
    fig,axarr=plt.subplots(2,3,figsize=(16,10),sharey=True,layout="constrained")
    allvals=np.concatenate([log[t["pg_data_row_1based"]-1] for t in targets])
    allvals=allvals[np.isfinite(allvals)]
    low,high=np.percentile(allvals,[0.5,99.5])
    pad=max(1.0,.14*(high-low))
    for panel,t in zip(axarr.flat,targets):
        vals=[log[t["pg_data_row_1based"]-1,
                  [i for i,s in enumerate(samples) if meta[s]==g]] for g in GROUPS]
        box(panel,vals,range(1,5),[G_COL[g] for g in GROUPS],.55)
        for ix,v in enumerate(vals,1):
            v=v[np.isfinite(v)]
            x=ix+np.linspace(-.16,.16,len(v))
            panel.scatter(x,v,color=T_COL[t["gene"]],s=16,zorder=5,
                          edgecolors="white",linewidths=.25,alpha=.83)
            panel.text(ix,low-pad*.77,f"{len(v)}/{len(vals[ix-1])}",
                       ha="center",fontsize=8,color="#5B6876")
        panel.set_xticks(range(1,5),["Internal\nQC","pooled\nQC","HC","S"])
        panel.set_ylim(low-pad,high+pad)
        panel.set_xlim(.55,4.45)
        panel.set_title(f'{t["gene"]} · {t["accession"]}\nHELP: {t["help_peptide"]}',
                        fontsize=12,fontweight="semibold",loc="left")
        style(panel)
    axarr.flat[5].axis("off")
    fig.supylabel("log2 raw PG protein intensity")
    fig.suptitle("Five HELP-source PG proteins across 4 groups (boxes + raw sample points)",
                 fontsize=15,fontweight="bold")
    save(fig,"02_five_protein_group_boxplots.png")

def sample_panels(log,meta,samples,targets):
    ordered=[s for g in GROUPS for s in samples if meta[s]==g]
    col={s:i for i,s in enumerate(samples)}
    positions=np.arange(1,76)
    fig,ax=plt.subplots(figsize=(21,8.6),layout="constrained")
    box(ax,[log[:,col[s]] for s in ordered],positions,
        [G_COL[meta[s]] for s in ordered],.62)
    for t in targets:
        vals=log[t["pg_data_row_1based"]-1,[col[s] for s in ordered]]
        valid=np.isfinite(vals)
        ax.scatter(positions[valid],vals[valid],color=T_COL[t["gene"]],
                   s=27,marker=MARK[t["gene"]],label=t["gene"],
                   edgecolors="white",linewidths=.3,zorder=7)
    start=1
    for g in GROUPS:
        n=sum(meta[s]==g for s in ordered)
        if start>1:ax.axvline(start-.5,lw=1,color="#81909F",alpha=.75)
        ax.text(start+(n-1)/2,1.025,f"{g}  (n={n})",
                ha="center",va="bottom",color=G_COL[g],fontweight="bold",
                transform=ax.get_xaxis_transform(),fontsize=12)
        start+=n
    ax.set_xticks(positions,ordered,rotation=90,fontsize=6)
    ax.set_xlim(.3,75.7)
    ax.set_ylabel("log2 raw PG intensity")
    ax.set_title("Raw per-sample PG boxplots with five HELP-related proteins highlighted",
                 fontsize=15,fontweight="bold",pad=30)
    ax.legend(loc="upper left",bbox_to_anchor=(1.005,1),frameon=False)
    style(ax)
    save(fig,"03_per_sample_boxplots_highlight_5proteins.png")

def output_tables(raw,log,meta,samples,targets):
    audit,groupstats,each=[],[],[]
    grp={g:[i for i,s in enumerate(samples) if meta[s]==g] for g in GROUPS}
    for t in targets:
        row=t["pg_data_row_1based"]-1
        audit.append({**t,"positive_detected_75":int(np.isfinite(log[row]).sum())})
        for i,s in enumerate(samples):
            value=raw[row,i]
            each.append({
                "sample_id":s,"group":meta[s],"gene":t["gene"],
                "accession":t["accession"],"help_peptide":t["help_peptide"],
                "raw_linear_intensity":value if np.isfinite(value) else np.nan,
                "raw_log2_intensity":log[row,i],
                "finite_positive":bool(np.isfinite(log[row,i])),
            })
        for g,indices in grp.items():
            values=raw[row,indices]
            measured=values[np.isfinite(values)&(values>0)]
            lv=np.log2(measured)
            groupstats.append({
                "gene":t["gene"],"accession":t["accession"],
                "help_peptide":t["help_peptide"],"group":g,
                "sample_n":len(indices),"positive_detected_n":len(measured),
                "not_detected_n":len(indices)-len(measured),
                "detection_rate_pct":100*len(measured)/len(indices),
                "linear_intensity_median":np.median(measured) if len(measured) else np.nan,
                "log2_median":np.median(lv) if len(lv) else np.nan,
                "log2_q1":np.percentile(lv,25) if len(lv) else np.nan,
                "log2_q3":np.percentile(lv,75) if len(lv) else np.nan,
                "cv_pct_original_linear":100*np.std(measured,ddof=1)/np.mean(measured)
                    if len(measured)>1 and np.mean(measured)>0 else np.nan,
            })
    pd.DataFrame(audit).to_csv(OUT/"matched_pg_protein_rows.csv",index=False)
    pd.DataFrame(groupstats).to_csv(OUT/"target_group_statistics.csv",index=False)
    pd.DataFrame(each).to_csv(OUT/"five_protein_sample_intensities.csv",index=False)
    return audit,groupstats

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError("Existing audit outputs protected: refusing overwrite")
    pg,meta,samples,raw,log,groups=load()
    targets=map_targets(pg)
    OUT.mkdir(parents=True,exist_ok=True)
    audit,stats=output_tables(raw,log,meta,samples,targets)
    group_overview(log,meta,samples,targets)
    target_panels(log,meta,samples,targets)
    sample_panels(log,meta,samples,targets)
    diag=[]
    for g in GROUPS:
        cols=[i for i,s in enumerate(samples) if meta[s]==g]
        values=log[:,cols].ravel()
        values=values[np.isfinite(values)]
        diag.append(dict(group=g,samples=groups[g],
                         positive_protein_x_sample_values=len(values),
                         possible_protein_x_sample_values=len(pg)*len(cols),
                         log2_median=np.median(values),
                         log2_q1=np.percentile(values,25),log2_q3=np.percentile(values,75)))
    pd.DataFrame(diag).to_csv(OUT/"all_pg_background_group_stats.csv",index=False)
    digest=hashlib.sha256(PG.read_bytes()).hexdigest()
    (OUT/"run_manifest.json").write_text(json.dumps({
        "pg_source":"pg_matrix.tsv","pg_sha256":digest,"proteins":len(pg),
        "samples":len(samples),"group_counts":groups,"matched":audit,
        "no_correction":True,"log2_display_only":True,"no_imputation":True,
    },ensure_ascii=False,indent=2),encoding="utf-8")
    lines=[
        "# Five HELP-source proteins · grouped raw PG boxplots","",
        "This is an exploratory check of the original, UNCORRECTED protein-level pg_matrix.tsv.",
        "Protein intensities are log2-transformed only for plotting.",
        "Missing/nonfinite/nonpositive values are excluded from plotting without imputation.","",
        "## Matching and dataset","",
        "- Original protein matrix: 6675 protein groups × 75 sample columns.",
        "- Group sizes: Internal-QC 6, pooled QC 15, HC 27, S 27.",
        "- Target identity is matched exactly by accession token in Protein.Group or Protein.Ids.",
        "- Each supplied accession must map to exactly one PG row (otherwise abort).",
        "- Expected gene and protein names are independently audited; see matched_pg_protein_rows.csv.",
        "- Peptide charge or heavy precursor entries are not extra rows in protein-group data.","",
        "| Accession | Protein name | Gene | HELP peptide | PG data row | Name matched | Gene matched | Detected samples |",
        "| --- | --- | --- | --- | ---: | --- | --- | ---: |"
    ]
    for t in audit:
        lines.append(f'| {t["accession"]} | {t["protein_name_expected"]} | {t["gene"]} | '
            f'{t["help_peptide"]} | {t["pg_data_row_1based"]} | '
            f'{"Yes" if t["protein_name_match"] else "No - check"} | '
            f'{"Yes" if t["gene_name_match"] else "No - check"} | '
            f'{t["positive_detected_75"]}/75 |')
    lines.extend([
        "","## Figures","",
        "### 1. All protein × sample intensities pooled by group, with five proteins highlighted","",
        "Grey boxplots: all valid raw PG protein × sample intensities pooled within each group.",
        "Five colored shapes: target protein group median; vertical colored strokes: target IQR.",
        "Do not interpret the grey background as across-sample CV or per-protein variability.","",
        "![Whole-proteome background and five proteins](01_group_background_boxplots_highlight_5proteins.png)","",
        "### 2. Four group boxplots separately for each target protein","",
        "Each target panel compares samples within the four groups; colored dots are measured samples.",
        "Numbers under the boxes: detected / all sample count. Common y-axis across panels.","",
        "![Five separate target protein boxplots](02_five_protein_group_boxplots.png)","",
        "### 3. All 75 sample-level PG intensity boxplots, with five proteins overlaid","",
        "Each box represents all valid raw PG protein intensities in one sample.",
        "Colored protein markers are the five matched proteins from that sample.","",
        "![Individual sample boxplots](03_per_sample_boxplots_highlight_5proteins.png)","",
        "## Target protein statistics by group","",
        "| Gene | Group | Detected/Total | Median log2 | Q1–Q3 log2 | CV (linear) |",
        "| --- | --- | ---: | ---: | ---: | ---: |"
    ])
    for st in stats:
        fmt=lambda x: f"{x:.2f}" if np.isfinite(x) else "NA"
        lines.append(f'| {st["gene"]} | {st["group"]} | {st["positive_detected_n"]}/{st["sample_n"]} | '
            f'{fmt(st["log2_median"])} | {fmt(st["log2_q1"])}–{fmt(st["log2_q3"])} | '
            f'{fmt(st["cv_pct_original_linear"])}% |')
    lines.extend([
        "","## Data outputs","",
        "- matched_pg_protein_rows.csv: strict accession matches and source annotations.",
        "- five_protein_sample_intensities.csv: 5 × 75 sample-level raw/log2 values.",
        "- target_group_statistics.csv: each target's group distributions and CV.",
        "- all_pg_background_group_stats.csv: pooled all-protein values per group.",
        "- run_manifest.json: source SHA-256, group counts, original data provenance.",
        "- Source data and prior correction results were not modified.",""])
    (OUT/"README.md").write_text("\n".join(lines),encoding="utf-8")
    print("MATCH:",[(t["gene"],t["Protein.Group"],t["protein_name_match"],t["gene_name_match"],t["positive_detected_75"]) for t in audit])
    print("PASS: 5 matched proteins, 4 group comparisons, 75 raw samples, 3 PNG figures.")

if __name__=="__main__":
    main()
