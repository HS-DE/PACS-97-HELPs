#!/usr/bin/env python3
"""Raw 97-HELP L / H / L:H boxplots: 75 samples, 2 peptides highlighted.

Read-only original 2026.10.13 HELP workbook and QC/Study metadata.
No f1/f2 normalization. All values log2-transformed for display.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from filter_samples import readbook, tab
from pacs_two_reference_plate_trial import HELP, get_samples
from yeast_protein_qc import GROUPS, COLORS, ROOT

OUT=ROOT/"results"/"help97_l_h_ratio_sample_boxplots"
LAYERS=("L","H","ratio_LH")
FOCUS=("STGSWSTLK","VYFAGFPR")
COLOR={"STGSWSTLK":"#C94552","VYFAGFPR":"#168DA0"}
MARK={"STGSWSTLK":"D","VYFAGFPR":"o"}
G_COLOR=dict(zip(GROUPS,COLORS))
FIGS={
  "L":"01_all_HELP_L_75sample_boxplots_highlight_two.png",
  "H":"02_all_HELP_H_75sample_boxplots_highlight_two.png",
  "ratio_LH":"03_all_HELP_ratio_75sample_boxplots_highlight_two.png"
}

def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()

def log2positive(array):
    a=np.asarray(array,dtype=float)
    out=np.full(a.shape,np.nan)
    good=np.isfinite(a)&(a>0)
    out[good]=np.log2(a[good])
    return out

def get_data(info):
    ss,sheets=readbook(HELP)
    if len(sheets)!=3:
        raise ValueError("Expected L, H and ratio_LH worksheets")
    data={}
    samples=None
    sequences=None
    for label,sh in zip(LAYERS,sheets):
        rows=tab(sh,ss)
        if len(rows)!=98:
            raise ValueError(f"{label}: need exactly 97 HELP peptide rows")
        header={v:k for k,v in rows[0].items() if k!="A"}
        if len(header)!=75 or set(header)!=set(info):
            raise ValueError(f"{label}: sample metadata mismatch")
        seq=[r.get("A","").strip() for r in rows[1:]]
        if len(set(seq))!=97 or not all(seq):
            raise ValueError(f"{label}: duplicate or empty peptide sequence")
        if sequences is None:sequences=seq
        if seq!=sequences:
            raise ValueError(f"{label}: sequences differ from sheet L")
        if samples is None:samples=list(header)
        linear=np.full((97,75),np.nan)
        for i,row in enumerate(rows[1:]):
            for j,sid in enumerate(samples):
                try:
                    val=float(row.get(header[sid],"nan"))
                    linear[i,j]=val if np.isfinite(val) else np.nan
                except (ValueError,TypeError):
                    pass
        logged=log2positive(linear)
        for focus in FOCUS:
            if not np.isfinite(logged[sequences.index(focus)]).all():
                raise ValueError(f"{label}: missing/nonpositive focus {focus}")
        data[label]={"linear":linear,"log2":logged}
    return sequences,samples,data

def sort_samples(info):
    names=sorted(info,key=lambda s:(
        GROUPS.index(info[s]["group"]),int(info[s]["plate"]),s))
    assert len(names)==75 and len(set(names))==75
    return names

def annotate_groups(ax,names,info):
    for g in GROUPS:
        index=[i+1 for i,sid in enumerate(names) if info[sid]["group"]==g]
        left,right=min(index),max(index)
        if left>1:ax.axvline(left-.5,color="#617083",lw=1.15,alpha=.7)
        ax.text((left+right)/2,1.022,f"{g} (n={len(index)})",
                transform=ax.get_xaxis_transform(),ha="center",va="bottom",
                fontsize=11,fontweight="bold",color=G_COLOR[g],clip_on=False)
        for k in range(left+1,right+1):
            if info[names[k-2]]["plate"]!=info[names[k-1]]["plate"]:
                ax.axvline(k-.5,color="#8797A6",lw=.75,ls=":",alpha=.65)

def plot_axis(ax,layer,names,source_names,sequences,data,info):
    ids=[source_names.index(s) for s in names]
    data_log=data[layer]["log2"][:,ids]
    vectors=[v[np.isfinite(v)] for v in data_log.T]
    if any(len(v)==0 for v in vectors):
        raise ValueError(f"Empty sample distribution in {layer}")
    pos=np.arange(1,76)
    bp=ax.boxplot(vectors,positions=pos,widths=.62,
                  showfliers=False,patch_artist=True,manage_ticks=False,
                  medianprops={"color":"#243646","linewidth":1.5},
                  boxprops={"color":"#7B8896","linewidth":.7},
                  whiskerprops={"color":"#7B8896","linewidth":.7},
                  capprops={"color":"#7B8896","linewidth":.7})
    for box,sid in zip(bp["boxes"],names):
        box.set_facecolor(G_COLOR[info[sid]["group"]])
        box.set_alpha(.48)
    for peptide,xoff in zip(FOCUS,(-.14,.14)):
        line=data_log[sequences.index(peptide)]
        ax.scatter(pos+xoff,line,s=35,marker=MARK[peptide],
                   color=COLOR[peptide],edgecolor="white",lw=.45,zorder=7,
                   label=f"{peptide} ({'P00751/CFB' if peptide=='STGSWSTLK' else 'P07225/PROS1'})")
    annotate_groups(ax,names,info)
    if layer=="ratio_LH":ax.axhline(0,c="#5D697A",ls="--",lw=.8,alpha=.6)
    ax.set_xlim(.35,75.65)
    ax.set_xticks(pos,names,rotation=90,fontsize=6)
    ax.set_ylabel({"L":"log2 Light intensity","H":"log2 Heavy intensity",
                   "ratio_LH":"log2 L/H ratio"}[layer],fontsize=11)
    ax.set_title({"L":"Light (L)","H":"Heavy (H)",
                  "ratio_LH":"L/H ratio"}[layer]+
                 " · all 97 HELP peptides / sample",
                 loc="left",fontsize=14,fontweight="bold",pad=37)
    ax.legend(loc="upper right",frameon=False,ncol=2,fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y",alpha=.15)
    ax.set_axisbelow(True)

def save(fig,path):
    fig.savefig(path,bbox_inches="tight",dpi=205,facecolor="white")
    plt.close(fig)
    if path.stat().st_size<12000:
        raise ValueError(f"Suspicious PNG file size for {path}")

def make_plots(names,source_names,seqs,data,info):
    for layer in LAYERS:
        fig,ax=plt.subplots(figsize=(23,8.1),layout="constrained")
        plot_axis(ax,layer,names,source_names,seqs,data,info)
        save(fig,OUT/FIGS[layer])
    fig,axes=plt.subplots(3,1,figsize=(23,19.5),layout="constrained")
    for ax,layer in zip(axes,LAYERS):
        plot_axis(ax,layer,names,source_names,seqs,data,info)
        if layer!="ratio_LH":
            ax.set_xticklabels([])
            ax.tick_params(axis="x",length=0)
    fig.suptitle("Raw HELP: Light / Heavy / L:H · 75 samples · 97 peptides",
                 fontsize=16,fontweight="bold")
    save(fig,OUT/"00_combined_three_layer_overview.png")

def make_tables(names,source_names,seqs,data,info):
    distribution=[]
    focus_rows=[]
    ix={sid:source_names.index(sid) for sid in names}
    for layer in LAYERS:
        log=data[layer]["log2"]
        linear=data[layer]["linear"]
        for sid in names:
            j=ix[sid]
            a=log[:,j]
            good=a[np.isfinite(a)]
            if not len(good):raise ValueError("Empty box")
            q1,median,q3=np.percentile(good,[25,50,75])
            distribution.append({
                "sheet":layer,"sample_id":sid,"group":info[sid]["group"],
                "plate":info[sid]["plate"],"run":info[sid]["run"],
                "valid_of_97":len(good),"q1_log2":float(q1),
                "median_log2":float(median),"q3_log2":float(q3),
                "min_log2":float(good.min()),"max_log2":float(good.max())})
            for peptide in FOCUS:
                i=seqs.index(peptide)
                y=float(log[i,j])
                focus_rows.append({
                    "sheet":layer,"sample_id":sid,"group":info[sid]["group"],
                    "plate":info[sid]["plate"],"run":info[sid]["run"],
                    "peptide":peptide,
                    "protein":"P00751/CFB" if peptide=="STGSWSTLK" else "P07225/PROS1",
                    "raw_value":float(linear[i,j]),"log2_value":y,
                    "sample_all_HELP_median_log2":float(median),
                    "difference_to_box_median_log2":y-float(median),
                    "percentile_among_valid_HELP":float(np.mean(good<=y)*100),
                    "within_box_IQR":bool(q1<=y<=q3),
                    "n_valid_HELP":len(good)})
    pd.DataFrame(distribution).to_csv(OUT/"all_HELP_sample_box_summaries.csv",index=False)
    pd.DataFrame(focus_rows).to_csv(OUT/"two_HELP_positions_by_sample.csv",index=False)
    assert len(distribution)==225 and len(focus_rows)==450
    return distribution,focus_rows

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError("Independent audit folder exists; no overwriting previous results")
    info=get_samples()
    seqs,source_names,data=get_data(info)
    names=sort_samples(info)
    OUT.mkdir(parents=True,exist_ok=False)
    make_tables(names,source_names,seqs,data,info)
    make_plots(names,source_names,seqs,data,info)
    checks={}
    for layer in LAYERS:
        x=data[layer]["log2"]
        checks[layer]={"valid":int(np.isfinite(x).sum()),
                       "nonpositive_missing_or_invalid":int((~np.isfinite(x)).sum()),
                       "total":int(x.size)}
    discrepancy={}
    for seq in FOCUS:
        i=seqs.index(seq)
        r=data["L"]["log2"][i]-data["H"]["log2"][i]-data["ratio_LH"]["log2"][i]
        if not np.all(np.isfinite(r)) or not np.allclose(r,0,atol=1e-10):
            raise ValueError(f"Original L/H/ratio discrepancy in {seq}")
        discrepancy[seq]=float(np.max(np.abs(r)))
    manifest={
        "source":"LH_ratio_tables_2026.10.13—ALL-HELP.xlsx",
        "source_sha256":digest(HELP),
        "sample_count":75,"HELP_count":97,
        "group_counts":{g:sum(info[s]["group"]==g for s in names) for g in GROUPS},
        "plots":list(FIGS.values())+["00_combined_three_layer_overview.png"],
        "sheet_valid_counts":checks,
        "selected_peptides":list(FOCUS),
        "ratio_L_equals_ratio_L_H_for_selected_max_abs_log2_error":discrepancy,
        "log2_display_only":True,"no_normalization":True,
        "no_imputation":True,
        "sample_sort":"group -> plate -> sample id",
        "marker_percentile":"fraction of positive valid HELP values <= selected HELP value * 100",
    }
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding="utf-8")
    readme=[
        "# 97 HELP peptides: Light / Heavy / L:H 75-sample boxplots","",
        "All boxplots show RAW 2026.10.13 HELP source workbook values.",
        "Exactly 97 HELP peptide sequences constitute each per-sample boxplot",
        "(excluding only nonpositive / invalid measurements, with no imputation).",
        "All values shown on log2 scale. No Plate/sample correction or factor1/factor2.","",
        "Two highlighted peptides in every sample:",
        "- Red diamond: STGSWSTLK (P00751/CFB).",
        "- Blue circle: VYFAGFPR (P07225/PROS1).","",
        "## Figures","",
        "- [Light L](01_all_HELP_L_75sample_boxplots_highlight_two.png)",
        "- [Heavy H](02_all_HELP_H_75sample_boxplots_highlight_two.png)",
        "- [L/H ratio](03_all_HELP_ratio_75sample_boxplots_highlight_two.png)",
        "- [Three-sheet overview](00_combined_three_layer_overview.png)","",
        "![Overview](00_combined_three_layer_overview.png)","",
        "Boxes are sample distributions across up to 97 HELP peptide signals.",
        "Colored points are individual highlighted peptide values on that same",
        "log2 axis. Minor horizontal offsets are for visibility only.",
        "Samples sorted Internal-QC / pooled QC / HC / S, then Plate and sample ID.",
        "Solid vertical boundaries denote groups; dotted lines denote Plates.","",
        "## Source checks","",
        "| Source sheet | Finite positive / 7275 | Missing or nonpositive |",
        "| --- | ---: | ---: |"]
    for layer in LAYERS:
        d=checks[layer]
        readme.append(f'| {layer} | {d["valid"]} | {d["nonpositive_missing_or_invalid"]} |')
    readme+=["","All 75 columns match metadata. 97 peptides uniquely identified on every sheet.",
              "Both focused peptides have valid values for all 75 samples in all sheets.",
              "For these two sequences, log2(L/H) matches log2(L)-log2(H),",
              "up to floating-point precision.","",
              "## Audit tables","",
              "- all_HELP_sample_box_summaries.csv: 225 observations, 75 samples x 3 sheets.",
              "- two_HELP_positions_by_sample.csv: 450 observations, 75 samples x 3 sheets x 2 peptides,",
              "  with source values, log2 values, sample median and percentile positions.",
              "- run_manifest.json: source hash, checks, and plotting definitions.","",
              "No existing experimental results or original spreadsheets were modified.",""]
    (OUT/"README.md").write_text("\n".join(readme),encoding="utf-8")
    print("VALID_COUNTS",checks,flush=True)
    print("SOURCE_RATIO_AGREEMENT",discrepancy,flush=True)
    print("PASS: 3 standalone PNGs + overview; 75 samples x 97 HELP; 450 highlighted positions",flush=True)

if __name__=="__main__":
    main()
