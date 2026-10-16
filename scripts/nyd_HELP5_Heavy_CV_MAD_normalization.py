#!/usr/bin/env python3
"""NYD: choose 5 Heavy HELP by full detection, low CV and low relative MAD;
normalize 39 NYD PG samples by the *opposite* of their mean Heavy log2 shift.

Confirmed T003 poor-quality sample MAP26092326E01 is excluded BEFORE all
ranking, reference and protein statistics, consistent with previous trials.
Each HELP reference is the median of its own log2 Heavy intensities across 39
retained clinical samples (NOT a pooled QC reference; none present).

Ranking: full finite-positive detection 39/39, relative MAD (raw linear
MAD divided by median, %) <= median relative MAD of fully detected peptides;
then choose lowest 5 linear-scale CV (sample SD ddof=1 / mean, %) with
relative MAD and peptide sequence as deterministic tiebreakers.
Report raw MAD, relative MAD, log2 MAD and complete 97-HELP ranking.

For each sample i, h_i = mean_j(log2(H_ij) - median_39(log2 H_j)).
The correction factor is f_i = -h_i; log2 PG_corrected = log2 PG_raw + f_i.
Raw-only PCA fitting; corrected projects on same PC directions, shared bounds;
PCA Raw and Corrected are TWO SEPARATE standalone charts.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import nyd_STGSWSTLK_heavy_exclude_T003_outlier as prev
import nyd_STGSWSTLK_heavy_batch_internal_trial as batch_plot
from filter_samples import readbook, tab
from nyd_pg_batch_boxplots import BATCHES, COLORS, ANNOT

OUT=prev.BASE/"HELP5_Heavy_CV_MAD_normalization_without_T003_outlier"
SINGLE_PLUS=prev.BASE/"STGSWSTLK_Heavy_single_reference_normalization_without_T003_outlier"
SINGLE_MINUS=prev.BASE/"STGSWSTLK_Heavy_subtract_factor_39samples"/"global_subtract_fi"
REQUIRED=5
BLUE="#2364E8"
YELLOW="#EABF28"
PALETTE=("#2166AC","#9B59B6","#24A189","#E28C3C","#D14B6E")
TARGETS={"STGSWSTLK":{"gene":"CFB","color":BLUE,"marker":"D"},
         "VYFAGFPR":{"gene":"PROS1","color":YELLOW,"marker":"o"}}

def csvwrite(path,rows):
    if not rows:raise ValueError("Empty CSV rows")
    with path.open("w",encoding="utf-8-sig",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def heavy_97(sample_ids):
    ss,sheets=readbook(prev.HELP)
    if len(sheets)!=1:raise ValueError("NYD Heavy workbook needs one sheet")
    rows=tab(sheets[0],ss)
    if len(rows)!=98:raise ValueError("NYD Heavy source must contain 97 HELP peptide rows")
    header={str(s).strip():col for col,s in rows[0].items() if col!="A"}
    if len(header)!=40 or len(set(header))!=40 or set(sample_ids)-set(header):
        raise ValueError("Heavy workbook missing expected study sample IDs")
    sequences=[str(row.get("A","")).strip() for row in rows[1:]]
    if len(set(sequences))!=97 or any(not s for s in sequences):
        raise ValueError("97 nonempty unique HELP peptide sequences required")
    a=np.array([[prev.val(row.get(header[s],"")) for s in sample_ids]
                for row in rows[1:]],dtype=float)
    if a.shape!=(97,39):raise AssertionError("Unexpected Heavy matrix shape")
    if not all(x in sequences for x in TARGETS):
        raise ValueError("Prior two target HELP peptides missing")
    return sequences,a

def peptide_metrics(sequences,heavy):
    rows=[]
    for seq,signals in zip(sequences,heavy):
        positive=signals[np.isfinite(signals)&(signals>0)]
        rate=len(positive)/len(signals)
        if positive.size:
            med=float(np.median(positive))
            absolute_mad=float(np.median(np.abs(positive-med)))
            relative_mad=float(100*absolute_mad/med) if med>0 else np.nan
            log=np.log2(positive)
            log_mad=float(np.median(np.abs(log-np.median(log))))
            cv=float(100*np.std(positive,ddof=1)/np.mean(positive)) if len(positive)>1 else np.nan
        else:
            med=absolute_mad=relative_mad=log_mad=cv=np.nan
        rows.append({
            "HELP_sequence":seq,
            "n_detected_positive":int(len(positive)),
            "n_total_samples":int(heavy.shape[1]),
            "detection_rate":rate,
            "CV_pct_original_linear_Heavy":cv,
            "MAD_original_linear_Heavy":absolute_mad,
            "relative_MAD_pct":relative_mad,
            "MAD_log2_Heavy":log_mad,
            "median_original_linear_Heavy":med
        })
    complete=[r for r in rows if r["n_detected_positive"]==39
              and np.isfinite(r["CV_pct_original_linear_Heavy"])
              and np.isfinite(r["relative_MAD_pct"])]
    if len(complete)<10:
        raise ValueError(f"Too few fully detected peptides: {len(complete)}")
    threshold=float(np.median([r["relative_MAD_pct"] for r in complete]))
    for row in rows:
        row["full_detection_39_of_39"]=row["n_detected_positive"]==39
        row["relative_MAD_guardrail_pct"]=threshold
        row["passes_low_MAD_guardrail"]=bool(
            row["full_detection_39_of_39"] and
            np.isfinite(row["relative_MAD_pct"]) and row["relative_MAD_pct"]<=threshold)
    allowed=sorted([r for r in complete if r["passes_low_MAD_guardrail"]],
                   key=lambda r:(r["CV_pct_original_linear_Heavy"],
                                 r["relative_MAD_pct"],r["HELP_sequence"]))
    if len(allowed)<REQUIRED:
        raise ValueError(f"Only {len(allowed)} fully detected low-MAD peptides")
    chosen=allowed[:REQUIRED]
    ids=[r["HELP_sequence"] for r in chosen]
    cv_ranked=sorted(complete,key=lambda r:(
        r["CV_pct_original_linear_Heavy"],r["relative_MAD_pct"],r["HELP_sequence"]))
    mad_ranked=sorted(complete,key=lambda r:(
        r["relative_MAD_pct"],r["CV_pct_original_linear_Heavy"],r["HELP_sequence"]))
    rc={r["HELP_sequence"]:k+1 for k,r in enumerate(cv_ranked)}
    rm={r["HELP_sequence"]:k+1 for k,r in enumerate(mad_ranked)}
    for r in rows:
        r["CV_rank_among_100pct_detected"]=rc.get(r["HELP_sequence"],"")
        r["relative_MAD_rank_among_100pct_detected"]=rm.get(r["HELP_sequence"],"")
        r["selected_top5"]=r["HELP_sequence"] in ids
        r["selected_rank_1_to_5"]=ids.index(r["HELP_sequence"])+1 if r["selected_top5"] else ""
    rows.sort(key=lambda r:(
        not r["full_detection_39_of_39"],
        r["CV_pct_original_linear_Heavy"] if np.isfinite(r["CV_pct_original_linear_Heavy"]) else float("inf"),
        r["HELP_sequence"]))
    return rows,chosen,threshold,len(complete)

def reference_and_factor(sequences,heavy,chosen):
    chosen_names=[r["HELP_sequence"] for r in chosen]
    chosen_idx=[sequences.index(s) for s in chosen_names]
    sample_log2=np.log2(heavy[chosen_idx,:])
    if not np.isfinite(sample_log2).all():
        raise ValueError("Chosen reference peptides must be valid in all 39 samples")
    q=np.median(sample_log2,axis=1)
    single_deviation=sample_log2-q[:,None]
    h=single_deviation.mean(axis=0)
    correction=-h
    assert np.allclose(correction,
        (q[:,None]-sample_log2).mean(axis=0),atol=1e-13)
    references=[{
        "rank":k+1,
        "HELP_sequence":s,
        "n_detected":39,
        "CV_pct_original_linear_Heavy":float(chosen[k]["CV_pct_original_linear_Heavy"]),
        "MAD_original_linear_Heavy":float(chosen[k]["MAD_original_linear_Heavy"]),
        "relative_MAD_pct":float(chosen[k]["relative_MAD_pct"]),
        "MAD_log2_Heavy":float(chosen[k]["MAD_log2_Heavy"]),
        "median_ref_log2_Heavy":float(q[k]),
        "median_ref_linear_Heavy":float(np.exp2(q[k])),
    } for k,s in enumerate(chosen_names)]
    return chosen_names,sample_log2,single_deviation,h,correction,references

def save_selection_plot(metrics,selected,threshold):
    vals=[r for r in metrics if r["full_detection_39_of_39"]]
    colors=["#2D8390" if r["passes_low_MAD_guardrail"] else "#B7AFA3"
            for r in vals]
    fig,ax=plt.subplots(figsize=(10.6,7.5),layout="constrained")
    ax.scatter([r["relative_MAD_pct"] for r in vals],
               [r["CV_pct_original_linear_Heavy"] for r in vals],
               color=colors,s=38,alpha=.72,edgecolors="white",linewidths=.4,
               label="100%-detected HELP")
    ax.axvline(threshold,c="#BA8041",ls="--",lw=1.4,
               label=f"Relative MAD guardrail = {threshold:.2f}%")
    for j,r in enumerate(selected):
        x,y=r["relative_MAD_pct"],r["CV_pct_original_linear_Heavy"]
        ax.scatter(x,y,color=PALETTE[j],s=155,marker="D",
                   edgecolors="#263747",lw=.9,zorder=7)
        ax.annotate(f"{j+1}. {r['HELP_sequence']}",xy=(x,y),xytext=(7,7+j%2*9),
                    textcoords="offset points",fontsize=9,fontweight="bold",
                    color=PALETTE[j],zorder=10)
    ax.set_xlabel("Relative MAD = raw Heavy MAD / median Heavy (%)")
    ax.set_ylabel("CV = SD / mean of raw linear Heavy (%)")
    ax.set_title("97 NYD HELP peptides · 39-sample CV/MAD selection",
                 loc="left",fontweight="bold",fontsize=15)
    ax.legend(frameon=False,fontsize=9)
    ax.grid(alpha=.14)
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    prev.save(fig,"01_all_HELP_CV_vs_relative_MAD_selected5.png")

def plot_help_boxes(sequences,heavy,names,records,selected):
    log=np.full(heavy.shape,np.nan)
    good=np.isfinite(heavy)&(heavy>0)
    log[good]=np.log2(heavy[good])
    fig,ax=plt.subplots(figsize=(20,8.7),layout="constrained")
    pos=np.arange(1,len(names)+1)
    b=ax.boxplot([v[np.isfinite(v)] for v in log.T],
                 positions=pos,widths=.65,patch_artist=True,
                 showfliers=False,manage_ticks=False,
                 medianprops={"color":"#283A46","lw":1.6})
    for patch,r in zip(b["boxes"],records):
        patch.set_facecolor(COLORS[r["batch"]]);patch.set_alpha(.46)
    for k,row in enumerate(selected):
        i=sequences.index(row["HELP_sequence"])
        ax.scatter(pos+(-.22+k*.11),log[i],marker="o",c=PALETTE[k],
                   s=32,edgecolors="white",lw=.4,zorder=6,
                   label=f"Selected {k+1}: {row['HELP_sequence']}")
    # Preserve the familiar two markers from earlier figures as separate
    # targets if not already among the five selected peptides.
    selected_ids={r["HELP_sequence"] for r in selected}
    for seq,spec in TARGETS.items():
        if seq in selected_ids:continue
        i=sequences.index(seq)
        ax.scatter(pos,log[i],s=26,color=spec["color"],marker=spec["marker"],
                   edgecolors="#273645",lw=.4,alpha=.88,
                   label=f"Reference marker: {seq}",zorder=5)
    for cut in (16.5,24.5,31.5):ax.axvline(cut,color="#8D9EAB",lw=1)
    ax.set_xticks(pos,names,rotation=90,fontsize=6.5)
    ax.set_ylabel("log2 raw Heavy intensity · valid 97 HELP peptides per sample")
    ax.set_title("All 97 HELP Heavy measurements per sample; selected five overlaid",
                 loc="left",fontweight="bold",fontsize=14)
    ax.legend(frameon=False,ncol=4,fontsize=8,loc="upper right")
    ax.grid(axis="y",alpha=.14);ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    prev.save(fig,"02_97_HELP_Heavy_39sample_boxes_selected5.png")

def plot_h_factors(samples,chosen):
    x=np.arange(1,len(samples)+1)
    h=np.array([r["h_mean_5HELP_log2_relative_all_samples"] for r in samples])
    correction=np.array([r["factor_minus_h_log2"] for r in samples])
    fig,axs=plt.subplots(2,1,figsize=(18,9.6),sharex=True,layout="constrained")
    for b in BATCHES:
        ix=np.array([j for j,r in enumerate(samples) if r["batch"]==b])
        axs[0].scatter(x[ix],h[ix],color=COLORS[b],s=55,
                       edgecolors="white",lw=.5,label=b)
        axs[1].bar(x[ix],correction[ix],color=COLORS[b],width=.77,alpha=.78)
    axs[0].axhline(0,color="#607887",lw=.8)
    axs[1].axhline(0,color="#607887",lw=.8)
    axs[0].set_ylabel("h = mean_j [log2(H_ij) − ref_j]")
    axs[1].set_ylabel("Applied fi = −h (log2)")
    axs[0].set_title("HELP relative offset h · retained for later Yeast comparison",
                     loc="left",fontweight="bold")
    axs[1].set_title("Add −h to every measured log2 protein-group intensity",
                     loc="left",fontweight="bold")
    axs[0].legend(frameon=False,ncol=4)
    for ax in axs:
        for cut in (16.5,24.5,31.5):ax.axvline(cut,c="#95A5B2",lw=.8)
        ax.grid(axis="y",alpha=.14);ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    axs[1].set_xticks(x,[r["sample_id"] for r in samples],rotation=90,fontsize=6.5)
    fig.suptitle("Five equally weighted HELP Heavy references · "+", ".join(chosen),
                 fontsize=14,fontweight="bold")
    prev.save(fig,"04_sample_HELP_h_and_negative_correction.png")

def cv_comparison_chart(comparison):
    fig,ax=plt.subplots(figsize=(13.9,6.5),layout="constrained")
    xx=np.arange(4)
    methods=[
        ("Raw","#8494A2",-.30),
        ("STGSWSTLK +fi","#CB9A43",-.10),
        ("STGSWSTLK −fi","#AC6BA9",.10),
        ("Selected 5 HELP −h","#268F9C",.30)
    ]
    for name,color,offset in methods:
        vals=[float(row[name]) for row in comparison]
        ax.bar(xx+offset,vals,width=.19,color=color,alpha=.95,label=name)
    ax.set_xticks(xx,[f"{b}\nn={prev.BATCH_COUNTS[b]}" for b in BATCHES])
    ax.set_ylim(0,max(row[m[0]] for row in comparison for m in methods)*1.22)
    ax.set_ylabel("Per-protein median within-batch CV (%) · linear scale")
    ax.set_title("One Heavy vs five selected Heavy references · exact same 39 samples",
                 loc="left",fontweight="bold",fontsize=14)
    ax.legend(frameon=False,ncol=4,fontsize=9)
    ax.grid(axis="y",alpha=.15);ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    prev.save(fig,"05_Raw_one_HELP_vs_five_HELP_batch_CV.png")

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError(f"Will not overwrite existing independent results: {OUT}")
    for refdir in (SINGLE_PLUS,SINGLE_MINUS):
        if not (refdir/"batch_protein_CV_summary.csv").exists():
            raise FileNotFoundError(f"Earlier comparable 39-sample trial missing: {refdir}")
    records,names,annotations,raw,protein_targets,original_single_h,original_40_ref=prev.inputs()
    assert len(names)==39 and prev.EXCLUDED_ID not in names
    sequences,heavy=heavy_97(names)
    if not np.allclose(heavy[sequences.index("STGSWSTLK")],original_single_h,rtol=0,atol=1e-10):
        raise AssertionError("Heavy source conflicts with earlier single STGSWSTLK extraction")
    metrics,chosen,threshold,ncomplete=peptide_metrics(sequences,heavy)
    selected_names,selected_log2,individual,h,fi,references=reference_and_factor(sequences,heavy,chosen)
    after=raw+fi[None,:]
    assert np.array_equal(np.isfinite(after),np.isfinite(raw))
    assert np.allclose(after-fi[None,:],raw,equal_nan=True,atol=1e-12)
    OUT.mkdir(parents=True,exist_ok=False)
    prev.OUT=OUT
    prev.writecsv("HELP97_Heavy_CV_MAD_detection_ranking.csv",metrics)
    prev.writecsv("selected_5_HELP_reference_medians.csv",references)
    details=[]
    for j,sid in enumerate(names):
        item={
          "sample_id":sid,"batch":records[j]["batch"],
          "h_mean_5HELP_log2_relative_all_samples":float(h[j]),
          "factor_minus_h_log2":float(fi[j]),
          "multiplier_linear_2_power_minus_h":float(np.exp2(fi[j])),
          "n_selected_HELP_used":5,
          "Raw_PG_sample_median_log2":float(np.nanmedian(raw[:,j])),
          "Corrected_PG_sample_median_log2":float(np.nanmedian(after[:,j])),
          "n_PG_positive":int(np.isfinite(raw[:,j]).sum()),
        }
        for k,seq in enumerate(selected_names):
            item[f"{seq}_Heavy_linear"]=float(np.exp2(selected_log2[k,j]))
            item[f"{seq}_Heavy_log2"]=float(selected_log2[k,j])
            item[f"{seq}_reference_log2"]=float(references[k]["median_ref_log2_Heavy"])
            item[f"{seq}_delta_sample_minus_reference_log2"]=float(individual[k,j])
        for acc,params in prev.PROTEINS.items():
            i=protein_targets[acc]
            item[f"{params['gene']}_{acc}_Raw_log2"]=float(raw[i,j])
            item[f"{params['gene']}_{acc}_Corrected_log2"]=float(after[i,j])
        details.append(item)
    prev.writecsv("sample_five_HELP_h_and_negative_factors.csv",details)
    with (OUT/"corrected_NYD_pg_matrix2_log2.tsv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.writer(f,delimiter="\t")
        writer.writerow(list(ANNOT)+names)
        for i,ann in enumerate(annotations):
            writer.writerow(list(ann)+[float(z) if np.isfinite(z) else "" for z in after[i]])
    cv,by_protein,values=prev.stats_cv(raw,after,records,annotations)
    prev.writecsv("batch_protein_CV_summary.csv",cv)
    prev.writecsv("per_protein_batch_CV_long.csv",by_protein)
    scores,pca_meta=prev.raw_fit_pca(raw,after,names,records)
    prev.writecsv("PCA_Raw_fitted_sample_coordinates.csv",scores)
    # Separated Raw PCA / Corrected PCA on exactly one Raw PCA frame.
    limits=batch_plot.separate_common_pca(scores,pca_meta)
    selected_rows=[dict(x) for x in chosen]
    save_selection_plot(metrics,selected_rows,threshold)
    plot_help_boxes(sequences,heavy,names,records,selected_rows)
    prev.plot_pgs(raw,after,records,names,protein_targets)
    plot_h_factors(details,selected_names)
    prev.plot_cv(cv,values)
    prev.plot_target_proteins(raw,after,records,protein_targets)
    prev.plot_sample_medians(raw,after,names,records)
    previous_plus=pd.read_csv(SINGLE_PLUS/"batch_protein_CV_summary.csv").set_index(["stage","batch"])
    previous_minus=pd.read_csv(SINGLE_MINUS/"batch_protein_CV_summary.csv").set_index(["stage","batch"])
    comp=[]
    for b in BATCHES:
        rawcv=next(r["median_protein_cv_pct"] for r in cv if r["stage"]=="Raw" and r["batch"]==b)
        v=next(r["median_protein_cv_pct"] for r in cv if r["stage"]=="Corrected" and r["batch"]==b)
        p=float(previous_plus.loc[("Corrected",b),"median_protein_cv_pct"])
        m=float(previous_minus.loc[("Corrected",b),"median_protein_cv_pct"])
        if not np.isclose(rawcv,float(previous_plus.loc[("Raw",b),"median_protein_cv_pct"]),rtol=0,atol=1e-10):
            raise AssertionError("Raw baseline changed")
        comp.append({"batch":b,"n":prev.BATCH_COUNTS[b],
            "Raw":rawcv,"STGSWSTLK +fi":p,"STGSWSTLK −fi":m,
            "Selected 5 HELP −h":v,"Selected_5_HELP_minus_Raw_pp":v-rawcv})
    prev.writecsv("Raw_single_HELP_vs_5HELP_batch_CV_comparison.csv",comp)
    cv_comparison_chart(comp)
    images=[
       "01_all_HELP_CV_vs_relative_MAD_selected5.png",
       "02_97_HELP_Heavy_39sample_boxes_selected5.png",
       "01_PG_39samples_Raw_vs_corrected_highlight_CFB_PROS1.png",
       "04_sample_HELP_h_and_negative_correction.png",
       "03_batch_protein_CV_raw_vs_corrected.png",
       "05_Raw_one_HELP_vs_five_HELP_batch_CV.png",
       "04a_PCA_Raw.png","04b_PCA_Corrected.png",
       "05_CFB_PROS1_by_batch_Raw_vs_corrected.png",
       "06_sample_PG_medians_Raw_vs_corrected.png"
    ]
    manifest={
      "source_sha256":{x.name:prev.sha256(x) for x in (prev.PG,prev.HELP,prev.META)},
      "n_original_samples":40,"n_retained_samples":39,
      "excluded_sample_id":prev.EXCLUDED_ID,
      "excluded_before_CV_MAD_detection_and_reference":True,
      "number_HELP_in_Heavy":97,"n_fully_detected_HELP":ncomplete,
      "selection_rule":"among 39/39 detected, relative MAD% <= median relative MAD% of full-detection peptides; select 5 lowest original-linear-scale CV%",
      "relative_MAD_guardrail_pct":threshold,
      "relative_MAD_definition":"100 * median(abs(original_linear_H - median(original_linear_H))) / median(original_linear_H)",
      "MAD_original_linear_definition":"median(abs(H - median(H)))",
      "MAD_log2_definition":"median(abs(log2 H - median(log2 H)))",
      "CV_pct_definition":"100 * sample_standard_deviation(H,ddof=1) / mean(H)",
      "selected_HELP_in_CV_rank_order":selected_names,
      "n_selected_HELP":5,
      "reference":"Per HELP median(log2 Heavy) across 39 retained clinical samples; not pooled QC",
      "h":"mean of the five [log2 Heavy of sample - its own 39-sample reference] values",
      "fi":"-h","corrected_PG":"raw PG log2 + (-h)",
      "no_extra_weights":True,"equal_weight_per_HELP":0.2,
      "no_imputation_or_factor_capping":True,
      "no_source_or_previous_results_modified":True,
      "n_protein_groups":4264,
      "valid_original_PG_cells":int(np.isfinite(raw).sum()),
      "valid_corrected_PG_cells":int(np.isfinite(after).sum()),
      "PCA":pca_meta,"PCA_same_axis_limits":limits,
      "PCA_display":"Raw and Corrected separate PNGs using one Raw fit",
      "CV_on_linear_PG":True,"figures":images
    }
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    lines=[
       "# NYD: five HELP Heavy references selected by CV and MAD",
       "",
       "Input is the 97 HELP rows of NYD_data/40samples/Heavy.xlsx;",
       "exclude known problematic T003 MAP26092326E01 before ranking, factor, CV and PCA.",
       "Retained 39 samples: T001=16, T002=8, T003=7, T004=8.",
       "No independent pooled QC samples exist here; all 39 samples provide the reference.",
       "",
       "## Automatic selection (no manual choice)",
       "",
       "1. For EVERY one of 97 HELP sequences, compute detection rate, CV on",
       "   original linear Heavy intensity, raw MAD, normalized relative MAD",
       "   and MAD on log2 Heavy intensity, across all 39 retained samples.",
       "2. Require detection rate exactly 100% (39/39 finite positive).",
       "3. Low-MAD guardrail: relative MAD no greater than the median relative",
       "   MAD among all peptides meeting step 2. This controls for differences",
       "   in raw intensity scale. Guardrail threshold: "+f"{threshold:.4f}%.",
       "4. From the low-MAD candidate set, choose FIVE lowest CV% peptides",
       "   (relative MAD ascending, sequence alphabetic for ties).",
       "5. Full screening results for all 97 in HELP97_Heavy_CV_MAD_detection_ranking.csv.",
       "",
       "## Selected five HELP peptides",
       "",
       "| Rank | HELP sequence | CV (%) | Raw linear MAD | Relative MAD (%) | log2 MAD | log2 reference |",
       "| ---: | --- | ---: | ---: | ---: | ---: | ---: |"]
    for r in references:
        lines.append(f"| {r['rank']} | {r['HELP_sequence']} | {r['CV_pct_original_linear_Heavy']:.3f} | {r['MAD_original_linear_Heavy']:.3f} | {r['relative_MAD_pct']:.3f} | {r['MAD_log2_Heavy']:.5f} | {r['median_ref_log2_Heavy']:.6f} |")
    lines+=["",
       "## Heavy-only factor formula",
       "",
       "For each chosen peptide j, reference_j = median of its log2 Heavy",
       "among all 39 retained samples (same peptide; NEVER pooled QC).",
       "For each sample i and each selected peptide j:",
       "delta_ij = log2(Heavy_ij) - reference_j.",
       "h_i = arithmetic mean of five delta_ij (each weight 1/5).",
       "fi = -h_i, and corrected log2 PG_i = original log2 PG_i + fi.",
       "Store h and five individual delta contributions in sample_five_HELP_h_and_negative_factors.csv",
       "so future Yeast relative shifts can be compared.",
       "No other HELP sequence contributes to any protein correction.",
       "",
       "## Same-39-sample within-batch median protein CV (%)",
       "",
       "| Batch | n | Raw | Single STGSWSTLK +fi | Single STGSWSTLK −fi | Five HELP −h |",
       "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for r in comp:
        lines.append(f"| {r['batch']} | {r['n']} | {r['Raw']:.2f} | {r['STGSWSTLK +fi']:.2f} | {r['STGSWSTLK −fi']:.2f} | {r['Selected 5 HELP −h']:.2f} |")
    lines += [
       "",
       "Batch CV is calculated on original linear protein intensity",
       "with sample SD / mean (ddof=1), requiring at least 70% detected samples.",
       "This CV includes biological differences between clinical samples and",
       "is not an independent technical QC measure.",
       "",
       "## Selection and distribution figures",
       "",
       "![97 HELP CV / relative MAD selection](01_all_HELP_CV_vs_relative_MAD_selected5.png)",
       "",
       "![97 HELP per-sample Heavy distributions and selected five](02_97_HELP_Heavy_39sample_boxes_selected5.png)",
       "",
       "## PG sample boxplots (blue CFB, yellow PROS1)",
       "",
       "![Raw vs corrected PG](01_PG_39samples_Raw_vs_corrected_highlight_CFB_PROS1.png)",
       "",
       "## Sample h and −h correction",
       "",
       "![Per-sample Heavy h and correction](04_sample_HELP_h_and_negative_correction.png)",
       "",
       "## CV comparisons",
       "",
       "![Raw/Corr CV](03_batch_protein_CV_raw_vs_corrected.png)",
       "",
       "![Raw/single/five HELP CV](05_Raw_one_HELP_vs_five_HELP_batch_CV.png)",
       "",
       "## Separate PCA charts, COMMON Raw-fitted basis",
       "",
       "One fit to Raw PG only; selected features, Raw imputation medians, Raw",
       "centering and eigenvectors identical for corrected projection.",
       "Same plotted x/y limits, no arrows or combined Raw/Corrected plot.",
       "",
       "### Raw PCA",
       "",
       "![Raw PCA](04a_PCA_Raw.png)",
       "",
       "### Corrected PCA",
       "",
       "![Corrected PCA](04b_PCA_Corrected.png)",
       "",
       "## Other PG validation plots",
       "",
       "![CFB and PROS1 by batch](05_CFB_PROS1_by_batch_Raw_vs_corrected.png)",
       "",
       "![Sample PG medians](06_sample_PG_medians_Raw_vs_corrected.png)",
       "",
       "## Reproducibility",
       "",
       "CSV: full screening 97 HELP; five peptide reference medians;",
       "39 sample offsets and all five individual contributions;",
       "Raw-vs-corrected batch CV and per-PG per-batch CV;",
       "Raw-only-fit PCA scores; full corrected 4264 x 39 log2 PG matrix.",
       "Manifest includes source SHA256 and exact selection threshold.",
       "No source or previous result was changed.",
       ""
    ]
    (OUT/"README.md").write_text("\n".join(lines),encoding="utf-8")
    print("ALL_HELP_SCREENING",json.dumps({"n_full_detection":ncomplete,"low_MAD_threshold_pct":threshold}),flush=True)
    print("SELECTED_5_HELP",json.dumps(references),flush=True)
    print("FACTOR_H_RANGE",json.dumps({"h_min":float(h.min()),"h_max":float(h.max()),"h_median":float(np.median(h))}),flush=True)
    print("BATCH_CV",json.dumps(comp),flush=True)
    print("PASS: 97 HELP screening, 5 selected, 39 h and -h factors, 4264x39 PG, 10 PNG, two common-frame PCA",flush=True)

if __name__=="__main__":
    main()
