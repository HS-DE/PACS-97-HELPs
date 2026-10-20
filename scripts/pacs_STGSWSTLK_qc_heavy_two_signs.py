#!/usr/bin/env python3
"""PACS 75-sample STGSWSTLK Heavy calibration against pooled QC references.

A: median log2 Heavy across 15 pooled QC samples.
B: median log2 Heavy across the five pooled QC samples on the same plate.
For each: test +fi and -fi, fi = reference - log2(H_sample).
All 75 original sample columns (HC/S/pooled QC/Internal-QC) are corrected.
No NYD batch-specific signs are transferred to the three PACS plates.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from filter_samples import readbook, tab
from yeast_protein_qc import ROOT, GROUPS, COLORS
import pacs_two_reference_plate_trial as old
import pacs_help5_heavy_only_trial as prior

OUT = ROOT / "results" / "stgswstlk_heavy_qc_samplewise_sign_comparison"
PEPTIDE = "STGSWSTLK"
PLATES = ("1", "2", "3")
REFERENCES = ("overall_reference", "within_plate_reference")
DIRECTIONS = ("plus_fi", "minus_fi")
STAGES = ("Raw",) + tuple(f"{ref}/{direction}" for ref in REFERENCES for direction in DIRECTIONS)
GROUP_COLORS = dict(zip(GROUPS, COLORS))


def original():
    info = old.get_samples()
    pg = pd.read_csv(old.PG, sep="\t", encoding="utf-8-sig",
                     dtype=str, keep_default_na=False, low_memory=False)
    names = list(pg.columns[5:])
    if list(pg.columns[:5]) != old.ANN or pg.shape != (6675, 80):
        raise ValueError(f"Unexpected input PG matrix layout {pg.shape}")
    if len(names) != 75 or set(names) != set(info) or len(set(names)) != 75:
        raise ValueError("PACS metadata and raw PG IDs mismatch")
    raw = old.valid_log2(pg[names].apply(pd.to_numeric, errors="coerce").to_numpy(float))
    if not any(np.isfinite(raw).ravel()):
        raise ValueError("Empty PG log2 matrix")

    shared, sheets = readbook(prior.HELP)
    if len(sheets) != 3:
        raise ValueError("Expected L/H/ratio_LH sheets")
    rows = tab(sheets[1], shared)  # exactly the raw Heavy sheet
    header = {str(v).strip(): col for col, v in rows[0].items() if col != "A"}
    by_sequence = {str(row.get("A", "")).strip(): row for row in rows[1:]}
    if len(rows) != 98 or len(by_sequence) != 97 or set(header) != set(names):
        raise ValueError("Unexpected 97 HELP raw Heavy layout")
    if PEPTIDE not in by_sequence:
        raise ValueError("STGSWSTLK absent from Heavy sheet")
    heavy = np.array([float(by_sequence[PEPTIDE].get(header[sid], "nan"))
                      for sid in names], dtype=float)
    if not np.all(np.isfinite(heavy) & (heavy > 0)):
        raise ValueError("Invalid STGSWSTLK Heavy in at least one PACS sample")
    pooled = np.array([j for j, sid in enumerate(names)
                       if info[sid]["group"] == "pooled QC"], dtype=int)
    internal = np.array([j for j, sid in enumerate(names)
                         if info[sid]["group"] == "Internal-QC"], dtype=int)
    assert len(pooled) == 15 and len(internal) == 6
    assert {g: sum(info[s]["group"] == g for s in names) for g in GROUPS} == {
        "Internal-QC": 6, "pooled QC": 15, "HC": 27, "S": 27}
    for p in PLATES:
        assert sum(info[names[j]]["plate"] == p for j in pooled) == 5
    return pg, names, info, raw, heavy, pooled, internal


def write_matrix(pg, names, matrix, target):
    output = pg[old.ANN].copy()
    for j, sid in enumerate(names):
        output[sid] = matrix[:, j]
    output.to_csv(target, sep="\t", index=False, na_rep="")


def refs_and_factors(names, info, heavy, pooled):
    h = np.log2(heavy)
    ref_global = float(np.median(h[pooled]))
    ref_plate = {}
    for p in PLATES:
        ids = [j for j in pooled if info[names[j]]["plate"] == p]
        ref_plate[p] = float(np.median(h[ids]))
    fs = {
        "overall_reference": np.full(len(names), ref_global) - h,
        "within_plate_reference": np.array(
            [ref_plate[info[sid]["plate"]] for sid in names]) - h,
    }
    # LOO is a descriptive QC sensitivity check only; corrections use exact medians.
    qc_loo = {}
    for ref in REFERENCES:
        off = np.full(len(names), np.nan)
        for j in pooled:
            own_plate = info[names[j]]["plate"]
            peers = [k for k in pooled
                     if k != j and (ref == "overall_reference" or
                                    info[names[k]]["plate"] == own_plate)]
            if len(peers) != (14 if ref == "overall_reference" else 4):
                raise AssertionError("Invalid pooled-QC leave-one-out reference")
            off[j] = float(np.median(h[peers]) - h[j])
        qc_loo[ref] = off
    return h, ref_global, ref_plate, fs, qc_loo


def qc_plate_cv(log_mat, names, info, stage):
    # Use 4 out of 5 detected in each plate, consistent with historic QC plots.
    result = []
    for p in PLATES:
        idx = [j for j,s in enumerate(names)
               if info[s]["group"] == "pooled QC" and info[s]["plate"] == p]
        if len(idx) != 5:
            raise AssertionError("Five pooled QC replicates required per plate")
        x = log_mat[:, idx]
        linear = np.exp2(np.clip(x, -100, 100))
        linear[~np.isfinite(x)] = np.nan
        valid = np.isfinite(linear) & (linear > 0)
        count = valid.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.nansum(linear, axis=1) / count
            sd = np.sqrt(np.nansum(np.where(valid, (linear-mean[:, None])**2, 0),
                                   axis=1)/(count-1))
            cv = sd/mean*100
        cv[(count<4)|(mean<=1e-8)|(~np.isfinite(cv))] = np.nan
        vals = cv[np.isfinite(cv)]
        result.append({"method":stage, "plate":p,"n_pooled_QC":5,
                       "n_proteins_with_CV":len(vals),
                       "median_protein_CV_pct":float(np.median(vals))})
    return result


def sample_boxes(raw, corrected, names, info, out, headline):
    order = sorted(range(len(names)),key=lambda j:(
        GROUPS.index(info[names[j]]["group"]),
        int(info[names[j]]["plate"]), names[j]))
    fig, axes = plt.subplots(2,1,figsize=(21,14),sharex=True,sharey=True,layout="constrained")
    for ax,mat,title in zip(axes,(raw,corrected),("Original protein intensities","Corrected protein intensities")):
        bp=ax.boxplot([mat[:,j][np.isfinite(mat[:,j])] for j in order],
                      positions=np.arange(1,76), widths=.62,
                      showfliers=False, patch_artist=True, manage_ticks=False,
                      medianprops={"color":"#283849","linewidth":1.6})
        for patch,j in zip(bp["boxes"],order):
            patch.set_facecolor(GROUP_COLORS[info[names[j]]["group"]])
            patch.set_alpha(.50)
        for j in range(1,75):
            x,y=info[names[order[j-1]]],info[names[order[j]]]
            if x["group"]!=y["group"]:ax.axvline(j+.5,c="#6E7F90",lw=1.05)
            elif x["plate"]!=y["plate"]:ax.axvline(j+.5,c="#9EABBA",lw=.8,ls=":")
        ax.set_ylabel("PG log2 intensity")
        ax.set_title(title, loc="left", fontweight="bold", fontsize=13)
        ax.grid(axis="y",alpha=.15)
    axes[1].set_xticks(np.arange(1,76),[names[j] for j in order],
                       rotation=90,fontsize=5.5)
    fig.suptitle(headline+" | 75 samples, 6,675 PG",fontsize=15,fontweight="bold")
    fig.savefig(out,dpi=170,bbox_inches="tight",facecolor="white")
    plt.close(fig)


def reference_plot(h, names, info, q, ref_g, ref_p):
    fig, ax = plt.subplots(figsize=(13,5.5),layout="constrained")
    for p in PLATES:
        idx=[j for j in q if info[names[j]]["plate"]==p]
        ax.scatter([int(p)]*len(idx),h[idx],s=45,alpha=.65,label=f"Plate {p} pooled QC")
        ax.hlines(ref_p[p],int(p)-.28,int(p)+.28,color="#2C465A",linewidth=3)
    ax.axhline(ref_g,ls="--",lw=1.4,c="#BF6E49",label=f"Global 15-QC median {ref_g:.4f}")
    ax.set_xticks([1,2,3],["Plate 1","Plate 2","Plate 3"])
    ax.set_ylabel("Original STGSWSTLK Heavy log2 intensity")
    ax.set_title("Only pooled QC samples define Heavy reference",fontweight="bold")
    handles,labels=ax.get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.5,-.045),
               ncol=4,frameon=False,fontsize=9)
    ax.grid(alpha=.15)
    fig.savefig(OUT/"01_QC_Heavy_references.png",dpi=170,bbox_inches="tight")
    plt.close(fig)


def cv_figure(table, qc_table):
    fig, axes = plt.subplots(1,2,figsize=(16,5.8),layout="constrained")
    stages=list(STAGES)
    gg=["pooled QC","Internal-QC","HC","S"]
    x=np.arange(4)
    width=.15
    colors=["#8C98A8","#298EB5","#C77854","#438E74","#8B64B0"]
    for i,stage in enumerate(stages):
        values=[float(table[(table.method==stage)&(table.group==g)].iloc[0].median_cv_pct) for g in gg]
        axes[0].bar(x+(i-2)*width, values, width*.94,label=stage,color=colors[i])
    axes[0].set_xticks(x,gg,rotation=17)
    axes[0].set_ylabel("Median protein CV (%) | unlogged scale")
    axes[0].set_title("Across sample groups (HC/S include biology)")
    for i,stage in enumerate(stages):
        values=[float(qc_table[(qc_table.method==stage)&(qc_table.plate==p)].iloc[0].median_protein_CV_pct) for p in PLATES]
        axes[1].plot([1,2,3],values,marker="o",lw=1.7,color=colors[i],label=stage)
    axes[1].set_xticks([1,2,3],["Plate 1","Plate 2","Plate 3"])
    axes[1].set_ylabel("Within-plate pooled QC median protein CV (%)")
    axes[1].set_title("5 pooled QC per plate; includes reference reuse")
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.5,-.035),
               ncol=3,frameon=False,fontsize=8)
    fig.suptitle("PACS 75-sample STGSWSTLK Heavy normalization | 4 methods + Raw",
                 fontweight="bold")
    fig.savefig(OUT/"02_CV_group_and_QC_plate_comparison.png",dpi=170,bbox_inches="tight")
    plt.close(fig)


def pca_chart(stages, names, info, groups, filename):
    idx=[j for j,s in enumerate(names) if info[s]["group"] in groups]
    raw=stages["Raw"][:,idx]
    good = (np.isfinite(raw).sum(axis=1)>=math.ceil(.7*len(idx))) & (
        np.nanstd(raw,axis=1)>1e-8)
    source=raw[good]
    median=np.nanmedian(source,axis=1)
    source=np.where(np.isfinite(source),source,median[:,None]).T
    mean=source.mean(axis=0)
    centered=source-mean
    gram=centered@centered.T
    eig, vectors=np.linalg.eigh(gram)
    ix=np.argsort(eig)[::-1]
    eig=np.maximum(eig[ix],0)
    comps=(centered.T@vectors[:,ix[:2]])/np.sqrt(eig[:2])[None,:]
    for k in range(2):
        pivot=int(np.argmax(abs(comps[:,k])))
        if comps[pivot,k]<0:comps[:,k]*=-1
    explained=eig[:2]/eig.sum()*100
    scores={}
    for stage,mat in stages.items():
        x=mat[good][:,idx]
        x=np.where(np.isfinite(x),x,median[:,None]).T
        scores[stage]=(x-mean)@comps
    allxy=np.vstack(list(scores.values()))
    vmin=allxy.min(axis=0);vmax=allxy.max(axis=0)
    pad=np.maximum((vmax-vmin)*.08,1.)
    fig,axs=plt.subplots(1,5,figsize=(21,4.75),sharex=True,sharey=True,layout="constrained")
    for ax,stage in zip(axs,STAGES):
        coords=scores[stage]
        for group in groups:
            positions=np.array([info[names[j]]["group"]==group for j in idx])
            ax.scatter(coords[positions,0],coords[positions,1],s=17,
                       alpha=.78,color=GROUP_COLORS[group],label=group)
        ax.set_title(stage.replace("/","\n"),fontsize=10)
        ax.set_xlim(vmin[0]-pad[0],vmax[0]+pad[0])
        ax.set_ylim(vmin[1]-pad[1],vmax[1]+pad[1])
        ax.set_xlabel(f"Raw PC1 {explained[0]:.1f}%")
        ax.grid(alpha=.12)
    axs[0].set_ylabel(f"Raw PC2 {explained[1]:.1f}%")
    handles,labels=axs[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.5,-.04),
               ncol=4,frameon=False,fontsize=8)
    fig.suptitle(f"Raw-fitted PCA; corrected projected onto same components | {'/'.join(groups)}",
                 fontsize=13,fontweight="bold")
    fig.savefig(OUT/filename,dpi=160,bbox_inches="tight")
    plt.close(fig)
    return {"view":"/".join(groups),"n_samples":len(idx),"n_selected_PG":int(good.sum()),
            "raw_PC1_explained_pct":float(explained[0]),"raw_PC2_explained_pct":float(explained[1])}


def main():
    if OUT.exists():
        raise FileExistsError(f"Preserving original output; path already exists: {OUT}")
    pg,names,info,raw,heavy,pooled,internal=original()
    h,ref_global,ref_plate,factors,loo=refs_and_factors(names,info,heavy,pooled)
    matrix={"Raw":raw}
    for ref in REFERENCES:
        for direction in DIRECTIONS:
            signed = 1 if direction=="plus_fi" else -1
            result = raw + signed*factors[ref][None,:]
            if not np.array_equal(np.isfinite(result), np.isfinite(raw)):
                raise AssertionError("Missingness changed")
            if np.nanmax(np.abs(result-raw-signed*factors[ref][None,:]))>1e-10:
                raise AssertionError("Incorrect per-sample shift")
            matrix[f"{ref}/{direction}"]=result

    # References differ by a plate-specific constant, so same-plate QC CV must match.
    # Preserve source files and old results: all outputs go to a new directory.
    OUT.mkdir(parents=True)
    pd.DataFrame([{"reference_type":"overall_reference","plate":"all_15",
                   "pooled_QC_count":15,"reference_STGSWSTLK_Heavy_log2":ref_global}]
                 +[{"reference_type":"within_plate_reference","plate":p,"pooled_QC_count":5,
                    "reference_STGSWSTLK_Heavy_log2":ref_plate[p]} for p in PLATES]).to_csv(
        OUT/"QC_reference_values.csv",index=False)

    detail=[]
    for j,sid in enumerate(names):
        record={"sample_id":sid,"group":info[sid]["group"],
                "plate":info[sid]["plate"],"run":info[sid]["run"],
                "original_STGSWSTLK_Heavy_linear":float(heavy[j]),
                "original_STGSWSTLK_Heavy_log2":float(h[j]),
                "raw_PG_median_log2":float(np.nanmedian(raw[:,j]))}
        for ref in REFERENCES:
            target=ref_global if ref=="overall_reference" else ref_plate[info[sid]["plate"]]
            f=float(factors[ref][j])
            record[ref+"_reference_log2"]=float(target)
            record[ref+"_fi"]=f
            record[ref+"_QC_LOO_fi"]=float(loo[ref][j]) if np.isfinite(loo[ref][j]) else np.nan
            for direction in DIRECTIONS:
                sign=1 if direction=="plus_fi" else -1
                prefix=ref+"_"+direction
                record[prefix+"_actual_shift_log2"]=sign*f
                record[prefix+"_multiplier_linear"]=float(np.exp2(sign*f))
                record[prefix+"_PG_median_log2"]=float(np.nanmedian(matrix[ref+"/"+direction][:,j]))
        detail.append(record)
    pd.DataFrame(detail).to_csv(OUT/"sample_75_factors_and_medians.csv",index=False,na_rep="")

    summaries=[];pooled_plates=[];per_protein=[]
    for stage,mat in matrix.items():
        for group, cv, summary in old.cv_by_group(mat,names,info):
            summaries.append({"method":stage,**summary})
            per_protein.extend({"method":stage,"group":group,
                                "PG_row_1based":int(i+1),"cv_pct":float(v) if np.isfinite(v) else np.nan}
                               for i,v in enumerate(cv))
        pooled_plates+=qc_plate_cv(mat,names,info,stage)
        if stage=="Raw":
            continue
        ref,dir_=stage.split("/")
        folder=OUT/ref
        folder.mkdir(exist_ok=True)
        write_matrix(pg,names,mat,folder/f"corrected_pg_log2_{dir_}.tsv")
        sample_boxes(raw,mat,names,info,folder/f"03_pg_boxes_raw_vs_{dir_}.png",
                     f"{ref} · {dir_}")
    pd.DataFrame(summaries).to_csv(OUT/"group_protein_CV_summary.csv",index=False)
    pd.DataFrame(pooled_plates).to_csv(OUT/"pooled_QC_within_plate_CV.csv",index=False)
    pd.DataFrame(per_protein).to_csv(OUT/"per_protein_CV_long.csv",index=False,na_rep="")
    q_loo=[]
    for ref in REFERENCES:
        for direction in DIRECTIONS:
            sign=1 if direction=="plus_fi" else -1
            loo_matrix=matrix[ref+"/"+direction].copy()
            loo_matrix[:,pooled]=raw[:,pooled]+sign*loo[ref][pooled][None,:]
            q_loo.extend(qc_plate_cv(loo_matrix,names,info,f"{ref}/{direction}/QC_LOO"))
    pd.DataFrame(q_loo).to_csv(OUT/"QC_leave_one_out_sensitivity.csv",index=False)

    s=pd.DataFrame(summaries)
    p=pd.DataFrame(pooled_plates)
    for direction in DIRECTIONS:
        a=p[p.method=="overall_reference/"+direction].sort_values("plate")
        b=p[p.method=="within_plate_reference/"+direction].sort_values("plate")
        if not np.allclose(a.median_protein_CV_pct,b.median_protein_CV_pct,atol=1e-9,rtol=0):
            raise AssertionError(f"Within-plate QC CV should not depend on reference type: {direction}")
        for plate in PLATES:
            ids=[j for j,n in enumerate(names) if info[n]["plate"]==plate]
            delta=(matrix["overall_reference/"+direction]-matrix["within_plate_reference/"+direction])[:,ids]
            expected=(1 if direction=="plus_fi" else -1)*(ref_global-ref_plate[plate])
            if np.nanmax(np.abs(delta-expected))>1e-10:
                raise AssertionError("Reference change not a plate-constant shift")

    reference_plot(h,names,info,pooled,ref_global,ref_plate)
    cv_figure(s,p)
    pca_details=[pca_chart(matrix,names,info,["HC","S"],"04_raw_fitted_PCA_54_study.png"),
                 pca_chart(matrix,names,info,GROUPS,"05_raw_fitted_PCA_all_75.png")]
    (OUT/"PCA_raw_fitted_details.json").write_text(json.dumps(pca_details,indent=2)+"\n")
    manifest={
        "original_data":{"PG":"pg_matrix.tsv","HELP":"LH_ratio_tables_2026.10.13—ALL-HELP.xlsx",
                         "metadata":["all_QC_metadata.xlsx","all_sample_metadata.xlsx"]},
        "sources_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (old.PG,prior.HELP,old.META_QC,old.META_STUDY)},
        "n_samples":75,"n_PG":6675,"groups":{"HC":27,"S":27,"pooled QC":15,"Internal-QC":6},
        "plates":{p:sum(info[sid]["plate"]==p for sid in names) for p in PLATES},
        "reference_source":"pooled QC only; never Internal-QC or Study",
        "reference_count_global":15,"reference_count_per_plate":5,
        "STGSWSTLK_only":True,"used_sheet":"H",
        "median_reference_overall_log2":ref_global,
        "median_reference_by_plate_log2":ref_plate,
        "factor":"reference_median - sample_log2_STGSWSTLK_Heavy",
        "application":"raw_pg_log2 + fi for plus_fi; raw_pg_log2 - fi for minus_fi",
        "corrected_groups":"all 75 including Internal-QC; latter never reference",
        "n_nan_original":int(np.size(raw)-np.isfinite(raw).sum()),
        "no_imputation_for_CV_or_export":True,
        "QC_LOO_diagnostics":"QC reference excludes evaluated QC from its own reference, separately for 15 or 5 QC group; sensitivity only",
        "PCA":"Raw-fitted projections (study-only and all-groups) with identical feature set, center and axes per view",
        "NYD_sign_pattern_not_transferred":True
    }
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+"\n")
    lines=["# PACS 75-sample STGSWSTLK Heavy correction from pooled QC","",
           "Unlike NYD, these are not identical biological repeats. We do not",
           "transfer the NYD T001–T004 sign pattern to the three PACS plates.",
           "Instead the same raw Heavy-derived factor is evaluated in both signs.","",
           "## Input and exactly applied formulas","",
           "- Original pg_matrix.tsv (6,675 PG x 75 samples) and H worksheet",
           "  of LH_ratio_tables_2026.10.13—ALL-HELP.xlsx.",
           "- STGSWSTLK only (no other HELP enters the reference).",
           "- Overall reference: median log2 Heavy of all 15 pooled QC.",
           "- Within-plate reference: median log2 Heavy of 5 pooled QC from the same plate.",
           "- fi(sample) = reference - log2(STGSWSTLK Heavy in that sample).",
           "- plus_fi = raw PG log2 + fi; minus_fi = raw PG log2 - fi.",
           "- Apply factor to all 75 samples individually (HC, S, pooled QC, Internal-QC).",
           "- Internal-QC does NOT define the reference. The 15 pooled QC do.",
           "- Use the full-precision reference values; preserve missingness.","",
           "## QC reference values (log2 Heavy)","",
           f"- Overall 15-QC reference: {ref_global:.12f}",
           *[f"- Plate {plate} 5-QC reference: {ref_plate[plate]:.12f}" for plate in PLATES],
           "",
           "## Median per-protein CV by sample group (%)","",
           "| Method | pooled QC | Internal-QC | HC | S |",
           "|---|---:|---:|---:|---:|"]
    for stage in STAGES:
        values={r["group"]:r["median_cv_pct"] for r in summaries if r["method"]==stage}
        lines.append("| "+stage+" | "+" | ".join(f"{values[g]:.2f}" for g in ("pooled QC","Internal-QC","HC","S"))+" |")
    lines += ["","## Within-plate pooled QC median protein CV (%)","",
              "| Method | Plate 1 | Plate 2 | Plate 3 |",
              "|---|---:|---:|---:|"]
    for stage in STAGES:
        values={r["plate"]:r["median_protein_CV_pct"] for r in pooled_plates if r["method"]==stage}
        lines.append("| "+stage+" | "+" | ".join(f"{values[p]:.2f}" for p in PLATES)+" |")
    lines += ["","![Heavy QC references](01_QC_Heavy_references.png)","",
              "![PG protein CV comparison](02_CV_group_and_QC_plate_comparison.png)","",
              "![Raw-fitted PCA 54 Study](04_raw_fitted_PCA_54_study.png)","",
              "![Raw-fitted PCA all 75](05_raw_fitted_PCA_all_75.png)","",
              "## Outputs","",
              "- Four full corrected PG log2 TSV matrices, under overall_reference/ and within_plate_reference/.",
              "- Each matrix has 6,675 PG and 75 sample columns plus the 5 PG annotations.",
              "- Four Raw-vs-corrected boxplots with matching axis limits.",
              "- sample_75_factors_and_medians.csv: actual reference, original Heavy and per-sample offsets.",
              "- group_protein_CV_summary.csv and per_protein_CV_long.csv.",
              "- pooled_QC_within_plate_CV.csv: CV within each plate's 5 pooled QC.",
              "- QC_leave_one_out_sensitivity.csv: pooled QC leave-one-out only, not held-out QC.",
              "- PCA_raw_fitted_details.json: Raw feature/axes/projection provenance.",
              "- run_manifest.json: exact source hashes, formulas and provenance.","",
              "## Interpretation","",
              "Pooled QC-based factor evaluation is partially circular because these QC",
              "samples are the ones used for reference estimation. Leave-one-out",
              "diagnostics remove own-sample reuse but are not independent external QCs.",
              "HC/S CV mixes biological and technical differences.",
              "Internal-QC may be a different preparation and is corrected as requested,",
              "not assumed to have the same Heavy baseline as pooled QC.",
              "Batch-internal QC CV is identical for global vs within-plate reference",
              "at the same sign because those factors differ only by a plate constant.",
              "Neither an increased PCA separation nor a flat sample boxplot proves calibration."]
    (OUT/"README.md").write_text("\n".join(lines)+"\n")
    print("REFERENCES",json.dumps({"global":ref_global,"plates":ref_plate}))
    print("GROUP_CV",s.to_json(orient="records"))
    print("POOLED_QC_PLATE_CV",p.to_json(orient="records"))
    print("PASS 75 samples, two QC references x two factor directions, four complete PG matrices",flush=True)


if __name__=="__main__":
    main()
