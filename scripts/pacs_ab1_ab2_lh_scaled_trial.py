#!/usr/bin/env python3
"""Explore b-based ab1 Light and ab2 Heavy rescaled sample correction.

Use only STGSWSTLK and VYFAGFPR in L/H/ratio_LH. Read original PG matrix,
metadata, and old two-reference b's factor formula from independent source.
All inputs are log2 BEFORE medians. Scaling is quotient of log2 medians.
No prior results, source Excel or Internal-QC measurements are modified.
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pacs_two_reference_plate_trial as base
from filter_samples import readbook, tab

ROOT=base.ROOT
OUT=ROOT/"results"/"two_reference_plate_normalization_ab1_ab2"
SEQS=("STGSWSTLK","VYFAGFPR")
VARIANTS={"ab1":"L","ab2":"H"}

def lh_data(sample_ids):
    ss,worksheets=readbook(base.HELP)
    if len(worksheets)!=3:
        raise ValueError("Expected three worksheets: L, H, ratio_LH")
    output={}
    for name,sheet in zip(("L","H","ratio_LH"),worksheets):
        records=tab(sheet,ss)
        headers={sid:c for c,sid in records[0].items() if c!="A"}
        if len(headers)!=75 or set(headers)!=set(sample_ids):
            raise ValueError(f"{name} sample IDs inconsistent with PG")
        byseq={r.get("A",""):r for r in records[1:]}
        if len(records)!=98 or len(byseq)!=97:
            raise ValueError(f"{name} needs 97 distinct HELP rows")
        arr=[]
        for seq in SEQS:
            if seq not in byseq:
                raise ValueError(f"Missing {name} peptide {seq}")
            values=np.array([float(byseq[seq].get(headers[s],"nan")) for s in sample_ids])
            if not np.all(np.isfinite(values)&(values>0)):
                raise ValueError(f"Nonpositive/missing {name} intensities for {seq}")
            arr.append(np.log2(values))
        output[name]=np.vstack(arr)
    return output

def compute(sample_ids,info,raw,ratios,refs,tabulated):
    # Inherit earlier b, including original plate f1 and adjusted ratio.
    b_pg,after_ratio,factor_records,b_ref=base.build_factors(
        "b_plate_ratio_plus_sample",SEQS,sample_ids,info,raw,ratios,refs)
    frame=pd.DataFrame(factor_records).set_index("sample_id").loc[sample_ids].reset_index()
    study=base.indexed(sample_ids,info,lambda r:r["group"] in base.STUDY)
    qc=base.indexed(sample_ids,info,lambda r:r["group"]=="pooled QC")
    internal=base.indexed(sample_ids,info,lambda r:r["group"]=="Internal-QC")
    assert (len(study),len(qc),len(internal))==(54,15,6)
    f1=frame["f1_plate_log2"].to_numpy(float)
    f2i1=frame["f2_sample_log2"].to_numpy(float)
    stage_a=np.vstack([ratios[s] for s in SEQS])+f1[None,:]
    assert np.allclose(stage_a,after_ratio)
    m_si=np.median(stage_a[:,study],axis=0)
    m_s=float(np.median(m_si))
    assert np.allclose(f2i1[study],m_s-m_si,atol=1e-12)
    assert np.allclose(f2i1[np.r_[qc,internal]],0)
    outputs={}
    for branch,scale_source in VARIANTS.items():
        component=tabulated[scale_source]
        per_sample=np.median(component,axis=0)
        m_global=float(np.median(per_sample[study]))
        if not np.isfinite(m_global) or abs(m_global)<1e-8:
            raise ValueError(f"{branch} global Light/Heavy log2 median invalid")
        denominator=per_sample/m_global
        if np.any(~np.isfinite(denominator[study])) or np.any(denominator[study]<=0):
            raise ValueError(f"{branch} nonpositive Light/Heavy log2 normalization factor")
        f2=np.zeros(len(sample_ids))
        f2[study]=f2i1[study]/denominator[study]
        corrected=raw+(f1+f2)[None,:]
        if not np.array_equal(raw[:,internal],corrected[:,internal],equal_nan=True):
            raise AssertionError("Internal-QC inadvertently changed")
        if not np.allclose(corrected[:,qc],raw[:,qc]+f1[None,qc],equal_nan=True):
            raise AssertionError("pooled QC must only receive f1")
        table=frame[["sample_id","group","plate","run"]].copy()
        table["median_si_after_a_log2"]=np.median(stage_a,axis=0)
        table["median_s_54_study"]=m_s
        table["f1_plate_log2"]=f1
        table["f2i_1_unscaled_log2"]=f2i1
        table["source"]=("Light" if branch=="ab1" else "Heavy")
        table["median_li_or_hi_log2"]=per_sample
        table["median_l_or_h_global_study_log2"]=m_global
        table["f2l_or_f2h_scale"]=denominator
        table["f2_final_scaled_log2"]=f2
        table["total_factor_log2"]=f1+f2
        for i,seq in enumerate(SEQS):
            table["raw_log2_ratio_"+seq]=ratios[seq]
            table["after_a_log2_ratio_"+seq]=stage_a[i]
            table["raw_log2_L_"+seq]=tabulated["L"][i]
            table["raw_log2_H_"+seq]=tabulated["H"][i]
        info_out={
            "variant":branch,
            "factor1_plate_reference":b_ref["plate_reference_per_plate"],
            "factor1_global_reference":b_ref["plate_reference_global"],
            "median_s_after_a_study_log2":m_s,
            "median_l_or_h_study_global_log2":m_global,
            "scale_min_54_study":float(denominator[study].min()),
            "scale_median_54_study":float(np.median(denominator[study])),
            "scale_max_54_study":float(denominator[study].max()),
            "f2_unscaled_min":float(f2i1[study].min()),
            "f2_unscaled_max":float(f2i1[study].max()),
            "f2_scaled_min":float(f2[study].min()),
            "f2_scaled_max":float(f2[study].max()),
        }
        outputs[branch]=(corrected,table,info_out)
    return outputs,study

def save_method(branch,corrected,table,source_info,pg,names,info,raw,refs,rawcv):
    folder=OUT/branch
    folder.mkdir(parents=True,exist_ok=False)
    mat=pg[base.ANN].copy()
    for j,sid in enumerate(names):
        mat[sid]=corrected[:,j]
    mat.to_csv(folder/"corrected_pg_log2.tsv",sep="\t",index=False,na_rep="")
    table.to_csv(folder/"sample_factors_and_LH.csv",index=False)
    groupcv=base.cv_by_group(corrected,names,info)
    cv_summary=[{"variant":branch,**r} for _,_,r in groupcv]
    pd.DataFrame(cv_summary).to_csv(folder/"group_cv_summary.csv",index=False)
    rows=[]
    for g,values,summary in groupcv:
        rows.extend({"protein_data_row_1based":i+1,"protein_ids":pg.iloc[i]["Protein.Ids"],
                     "group":g,"cv_pct":v} for i,v in enumerate(values))
    pd.DataFrame(rows).to_csv(folder/"protein_group_cv.csv",index=False)
    base.plot_sample_boxes(raw,corrected,names,info,refs,SEQS,
                           folder/"01_sample_boxplots_targets_raw_vs_corrected.png")
    base.plot_cv(rawcv,groupcv,folder/"02_cv_violin_box_raw_vs_corrected.png")
    pca=base.plot_pca(raw,corrected,names,info,
                      folder/"03_pca_raw_vs_corrected.png")
    pd.DataFrame(pca).to_csv(folder/"pca_summary.csv",index=False)
    legacy=table.rename(columns={"f2_final_scaled_log2":"f2_sample_log2"})
    base.plot_factor(legacy.to_dict("records"),SEQS,branch,
                     folder/"04_factor_offsets.png")
    base.plot_plate_qc(refs,SEQS,raw,corrected,names,info,
                       folder/"05_pooled_QC_plate_medians.png")
    (folder/"factor_reference.json").write_text(
        json.dumps(source_info,indent=2,ensure_ascii=False),encoding="utf-8")
    return cv_summary

def figures_joint(outputs,study,rawcv,cvresults):
    fig,axs=plt.subplots(2,1,figsize=(16.5,10),sharex=True,layout="constrained")
    x=np.arange(len(study))
    old=outputs["ab1"][1]["f2i_1_unscaled_log2"].to_numpy()[study]
    for i,(branch,unit) in enumerate(VARIANTS.items()):
        frame=outputs[branch][1]
        scale=frame["f2l_or_f2h_scale"].to_numpy()[study]
        corrected=frame["f2_final_scaled_log2"].to_numpy()[study]
        ax=axs[i]
        ax.plot(x,old,color="#8B949E",label="original b: f2i_1",lw=1.4)
        ax.scatter(x,corrected,color="#3176B4" if branch=="ab1" else "#D68435",
                   s=26,label=f"{branch} scaled f2")
        ax.axhline(0,color="#A8B0BB",lw=.7)
        ax.set_ylabel("log2 additive factor f2")
        ax.legend(frameon=False,loc="upper left")
        base.beautify(ax)
        ax2=ax.twinx()
        ax2.plot(x,scale,color="#8B63A5",alpha=.55,ls="--",lw=1.2)
        ax2.set_ylabel(f"log2 {unit} median / study median",color="#7B5297")
        ax.set_title(f"{branch}: original b vs scaled sample offset",
                     loc="left",fontweight="bold")
    axs[-1].set_xticks(x,outputs["ab1"][1].iloc[study]["sample_id"],fontsize=6,rotation=90)
    fig.suptitle("54 HC/S samples · scaled vs original HELP sample factor",
                 fontsize=15,fontweight="bold")
    base.save(fig,OUT/"00_f2_factor_scaling_comparison.png")
    fig,ax=plt.subplots(figsize=(10.8,6.8),layout="constrained")
    group_values={s["variant"]:{x["group"]:x["median_cv_pct"] for x in cvresults
                  if x["variant"]==s["variant"]} for s in
                  [{"variant":"ab1"},{"variant":"ab2"}]}
    rawmed={g:r["median_cv_pct"] for g,_,r in rawcv}
    for j,g in enumerate(("HC","S")):
        values=[rawmed[g],group_values["ab1"][g],group_values["ab2"][g]]
        xx=np.arange(3)+(j-.5)*.31
        ax.bar(xx,values,width=.30,color=base.PAL[g],alpha=.85,label=g)
        for a,v in zip(xx,values):
            ax.text(a,v+.3,f"{v:.2f}%",ha="center",fontsize=9)
    ax.set_xticks(range(3),["Raw","ab1 Light","ab2 Heavy"])
    ax.set_ylim(0,max(rawmed.values())*1.7)
    ax.set_ylabel("Median protein CV on original linear intensities (%)")
    ax.legend(frameon=False)
    ax.set_title("HC/S protein CV · new scaling experiments",loc="left",
                 fontsize=14,fontweight="bold")
    base.beautify(ax)
    base.save(fig,OUT/"00_ab1_ab2_study_cv.png")

def write_readme(outputs,cvresults,rawcv,source_consistency,pg,names):
    raw={g:r["median_cv_pct"] for g,_,r in rawcv}
    prior={}
    for branch in ("a_plate_ratio","b_plate_ratio_plus_sample"):
        p=base.OUT/"two_refs"/branch/"group_cv_summary.csv"
        prior[branch]=pd.read_csv(p).set_index("group")["median_cv_pct"].to_dict()
    summary={b:{r["group"]:r["median_cv_pct"] for r in cvresults if r["variant"]==b}
             for b in VARIANTS}
    manifest={
        "input_sha256":{p.name:base.digest(p) for p in
                        (base.PG,base.HELP,base.META_QC,base.META_STUDY)},
        "protein_rows":len(pg),"samples":len(names),
        "reference_peptides":list(SEQS),
        "reference_proteins":["P00751","P07225"],
        "log2_conversion_before_all_medians":True,
        "Light_Heavy_scale":"median of log2 source intensity / study median of per-sample log2 medians",
        "factor2_scaled":"f2i1 / ratio_of_log2_medians",
        "sample_scaling_medians_from_HC_S_only":54,
        "plate_reference_f1":"same as prior two_refs b, 5 pooled QC per Plate",
        "correction_groups":{"f1":["HC","S","pooled QC"],"f2":["HC","S"],"Internal-QC":[]},
        "ab1":outputs["ab1"][2],"ab2":outputs["ab2"][2],
        "ratio_L_over_H_log2_consistency":source_consistency,
        "cv_evaluation":"sample SD / mean on back-transformed linear intensities",
        "no_imputation":True,
    }
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding="utf-8")
    lines=[
       "# PACS ab1 and ab2 · Light/Heavy-scaled b correction","",
       "Two new branches, both using STGSWSTLK (P00751 / CFB) and VYFAGFPR (P07225 / PROS1).",
       "The original Plate-level a factor1 and all historical results are unchanged.","",
       "## Formula and scope","",
       "Every raw L/H ratio, Light intensity, Heavy intensity and PG intensity",
       "is first converted to log2. We retain original b's pooled-QC Plate f1.",
       "For the 54 Study samples: median_si = median of two log2 selected HELP ratios",
       "after a's Plate correction; median_s = median across these 54 median_si;",
       "f2i1 = median_s - median_si.","",
       "- ab1: median_li = median of the two log2 Light values in sample i;",
       "  median_l = median of 54 Study median_li; f2l_i = median_li / median_l;",
       "  final f2_i = f2i1 / f2l_i.",
       "- ab2: median_hi = median of the two log2 Heavy values in sample i;",
       "  median_h = median of 54 Study median_hi; f2h_i = median_hi / median_h;",
       "  final f2_i = f2i1 / f2h_i.",
       "- HC and S: corrected_log2_PG = raw_log2_PG + f1_plate + f2_i.",
       "- pooled QC: raw_log2_PG + f1_plate, no f2; Internal-QC: unchanged.",
       "- f2 is NOT further added to the HELP ratios. All other 95 HELP peptides unused.",
       "- CV uses back-transformed original linear intensities, not log2 intensities.","",
       "The requested scaling is specifically the quotient of log2 medians.",
       "This can be greater or smaller than 1 and does not guarantee shrinkage.","",
       "## Median protein CV (%)","",
       "| Method | HC | S | pooled QC | Internal-QC |",
       "| --- | ---: | ---: | ---: | ---: |",
    ]
    for label,values in [("Raw",raw),("prior a",prior["a_plate_ratio"]),
                         ("prior b",prior["b_plate_ratio_plus_sample"]),
                         ("ab1 Light",summary["ab1"]),("ab2 Heavy",summary["ab2"])]:
        lines.append("| "+label+" | "+" | ".join(f"{values[g]:.2f}" for g in
                    ("HC","S","pooled QC","Internal-QC"))+" |")
    lines+=["","![Study CV comparison](00_ab1_ab2_study_cv.png)","",
            "![Scaled versus unscaled f2](00_f2_factor_scaling_comparison.png)","",
            "## Output files per branch","",
            "Each of ab1 and ab2 includes the full 6675x75 log2 corrected PG matrix;",
            "a per-sample CSV with f1, unscaled f2i1, log2 Light/Heavy medians,",
            "global Light/Heavy medians, scale factor, final f2 and selected L/H",
            "values; CV summaries, per-protein CV, PCA, reference metadata and five",
            "plots: sample boxplots with targets, CV violin/box, PCA,",
            "factor offsets, and pooled-QC Plate medians.","",
            "## L/H ratio consistency diagnostics","",
            "| Peptide | max(abs(log2 Light - log2 Heavy - log2 ratio)) |",
            "| --- | ---: |"]
    for x in source_consistency:
        lines.append(f'| {x["peptide"]} | {x["max_abs_difference_log2"]:.8f} |')
    lines+=["","No original files or older result directories are overwritten.",
           "This is an exploratory trial; HC/S CV may contain biological variation."]
    (OUT/"README.md").write_text("\n".join(lines)+"\n",encoding="utf-8")

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError("Results directory already exists; no overwrite allowed")
    pg,names,info,raw,ratios,refs=base.load()
    if len(pg)!=6675 or len(names)!=75:raise ValueError("Unexpected input size")
    sheet=lh_data(names)
    source_consistency=[]
    for i,seq in enumerate(SEQS):
        dif=sheet["L"][i]-sheet["H"][i]-sheet["ratio_LH"][i]
        source_consistency.append({"peptide":seq,
                                   "max_abs_difference_log2":float(np.abs(dif).max()),
                                   "median_abs_difference_log2":float(np.median(np.abs(dif)))})
    if not np.allclose(sheet["ratio_LH"],np.vstack([ratios[s] for s in SEQS]),atol=1e-12):
        raise ValueError("Selected ratio workbook and original b source do not match")
    methods,study=compute(names,info,raw,ratios,refs,sheet)
    OUT.mkdir(parents=True,exist_ok=False)
    baseline=base.cv_by_group(raw,names,info)
    summaries=[]
    for branch,(matrix,frame,details) in methods.items():
        print("START",branch,flush=True)
        summaries.extend(save_method(branch,matrix,frame,details,
                                     pg,names,info,raw,refs,baseline))
    pd.DataFrame(summaries).to_csv(OUT/"all_branches_group_cv.csv",index=False)
    figures_joint(methods,study,baseline,summaries)
    write_readme(methods,summaries,baseline,source_consistency,pg,names)
    print("PASS: ab1 ab2, 2 full log2 matrices, 10 method PNG + 2 comparison PNG",flush=True)
    print("CV RESULTS",summaries,flush=True)
    print("SCALING",methods["ab1"][2],methods["ab2"][2],flush=True)
    print("RATIO CONCORDANCE",source_consistency,flush=True)

if __name__=="__main__":
    main()
