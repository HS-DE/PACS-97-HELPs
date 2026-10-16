#!/usr/bin/env python3
"""NYD 40-sample raw PG + 97 HELP Heavy per-sample boxplots.

Raw PG source: NYD_data/40samples/NYD_pg_matrix2.tsv (4264 x 40)
Raw HELP source: NYD_data/40samples/Heavy.xlsx (97 x 40)
Sample/batch source: NYD_data/40samples/NYD_sampleinfo2.xlsx

Each sample box uses ALL valid positive entries of the corresponding source.
Original values are log2 transformed only for plotting (no normalization).
Blue diamond: CFB/P00751 protein or STGSWSTLK Heavy HELP peptide.
Yellow circle: PROS1/P07225 protein or VYFAGFPR Heavy HELP peptide.
Independent results, source files and existing NYD analyses unchanged.

For any FUTURE before-after PCA: fit PCA on Raw alone, transform Corrected
using the Raw feature selection, imputation constants, centering and axes,
and keep matched PC limits. This plotting task does not calculate PCA.
"""
import hashlib
import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from nyd_pg_batch_boxplots import meta, BATCHES, COLORS

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"NYD_data"/"40samples"
PG=BASE/"NYD_pg_matrix2.tsv"
HEAVY=BASE/"Heavy.xlsx"
SAMPLEINFO=BASE/"NYD_sampleinfo2.xlsx"
OUT=BASE/"PG_and_HELP_Heavy_two_reference_highlight"
PG_ANNOT=("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description")
TARGETS=(
    {"peptide":"STGSWSTLK","accession":"P00751","gene":"CFB",
     "color":"#2364E8","marker":"D","text_color":"blue"},
    {"peptide":"VYFAGFPR","accession":"P07225","gene":"PROS1",
     "color":"#F2C94C","marker":"o","text_color":"yellow"},
)
IMAGE_FILES={
    "protein":"01_NYD_40_sample_PG_boxplots_highlight_CFB_PROS1.png",
    "heavy":"02_NYD_40_sample_HELP_Heavy_boxplots_highlight_two.png",
    "combined":"00_NYD_PG_and_HELP_Heavy_combined.png"
}

def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as handle:
        for piece in iter(lambda:handle.read(1<<20),b""):
            h.update(piece)
    return h.hexdigest()

def log2_positive(source):
    arr=np.asarray(source,dtype=float)
    out=np.full(arr.shape,np.nan)
    good=np.isfinite(arr)&(arr>0)
    out[good]=np.log2(arr[good])
    return out

def samples_metadata():
    records=meta()
    ids=[r["sample_id"] for r in records]
    assert len(ids)==40 and len(set(ids))==40
    grouped=[r for batch in BATCHES for r in records if r["batch"]==batch]
    names=[r["sample_id"] for r in grouped]
    if len(names)!=40 or set(names)!=set(ids):
        raise ValueError("Metadata mismatch")
    return grouped,names

def match_accession(row,accession):
    # Use complete protein accession tokens rather than unsafe substring matches.
    for col in ("Protein.Group","Protein.Ids"):
        if accession in set(re.split(r"[;|,\s]+",str(row[col]).strip())):
            return True
    return False

def pg_source(sample_ids):
    data=pd.read_csv(PG,sep="\t",encoding="utf-8-sig",
                     dtype=str,keep_default_na=False,low_memory=False)
    if tuple(data.columns[:5])!=PG_ANNOT:
        raise ValueError("PG annotation columns unexpectedly changed")
    cols=list(data.columns[5:])
    if len(cols)!=40 or len(set(cols))!=40 or set(cols)!=set(sample_ids):
        raise ValueError(f"PG samples do not match metadata: {set(sample_ids)^set(cols)}")
    if len(data)!=4264:
        raise ValueError(f"Expected 4264 PG rows, got {len(data)}")
    values=data[sample_ids].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    log=log2_positive(values)
    targets={}
    for target in TARGETS:
        mask=data.apply(lambda r:match_accession(r,target["accession"]),axis=1)
        index=np.flatnonzero(mask.to_numpy())
        if len(index)!=1:
            raise ValueError(f"{target['accession']}: expected unique PG group, found {index.tolist()}")
        idx=int(index[0])
        targets[target["peptide"]]=idx
        if not np.isfinite(log[idx]).all():
            raise ValueError(f"{target['accession']}: not quantified positively in all 40 samples")
        if target["gene"] not in set(re.split(r"[;,|\s]+",data.iloc[idx]["Genes"])):
            raise ValueError(f"Unexpected gene annotation {target['accession']}: {data.iloc[idx]['Genes']}")
    return data,values,log,targets

def heavy_source(sample_ids):
    frame=pd.read_excel(HEAVY,sheet_name=0,engine="openpyxl",dtype=object)
    if frame.shape!=(97,41):
        raise ValueError(f"Expected 97 HELP peptides x (name+40 samples), got {frame.shape}")
    name_column=frame.columns[0]
    if pd.notna(name_column) and str(name_column).strip() not in ("Unnamed: 0",""):
        raise ValueError(f"Unexpected first HELP peptide column {name_column!r}")
    cols=list(frame.columns[1:])
    if len(set(cols))!=40 or set(cols)!=set(sample_ids):
        raise ValueError(f"HELP Heavy sample IDs do not match metadata: {set(sample_ids)^set(cols)}")
    seq=[str(v).strip() for v in frame.iloc[:,0].tolist()]
    if len(set(seq))!=97 or any(not item or item=="nan" for item in seq):
        raise ValueError("HELP Heavy 97 stripped sequences are not unique")
    values=frame[sample_ids].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    log=log2_positive(values)
    targets={}
    for target in TARGETS:
        match=[i for i,pep in enumerate(seq) if pep==target["peptide"]]
        if len(match)!=1:
            raise ValueError(f"{target['peptide']}: expected unique HELP Heavy row, got {len(match)}")
        targets[target["peptide"]]=match[0]
        if not np.isfinite(log[match[0]]).all():
            raise ValueError(f"{target['peptide']}: invalid/missing Heavy in one or more 40 samples")
    return seq,values,log,targets

def figstyle(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y",color="#637587",alpha=.16)
    ax.set_axisbelow(True)

def show_box(ax,log,targets,names,metadata,kind,compact=False):
    pos=np.arange(1,41)
    vectors=[log[:,j][np.isfinite(log[:,j])] for j in range(40)]
    if any(len(v)==0 for v in vectors):
        raise ValueError("Some sample has no measurable features in its boxplot")
    bp=ax.boxplot(vectors,positions=pos,widths=.63,patch_artist=True,
                  showfliers=False,manage_ticks=False,
                  medianprops={"color":"#243648","linewidth":1.7},
                  boxprops={"color":"#748595","linewidth":.85},
                  whiskerprops={"color":"#748595","linewidth":.85},
                  capprops={"color":"#748595","linewidth":.85})
    for artist,record in zip(bp["boxes"],metadata):
        artist.set_facecolor(COLORS[record["batch"]])
        artist.set_alpha(.46)
    for t,offset in zip(TARGETS,(-.15,.15)):
        row=targets[t["peptide"]]
        y=log[row]
        label=(f"{t['gene']} / {t['accession']}" if kind=="protein"
               else f"{t['peptide']} / {t['gene']}")
        ax.scatter(pos+offset,y,s=55 if not compact else 40,
                   facecolors=t["color"],edgecolors="#273545",
                   linewidths=.7,marker=t["marker"],zorder=8,label=label)
    nstart=1
    for batch in BATCHES:
        count=sum(r["batch"]==batch for r in metadata)
        if nstart>1:ax.axvline(nstart-.5,color="#6D7C8B",lw=1.2)
        ax.text(nstart+(count-1)/2,1.018,f"{batch} · n={count}",
                transform=ax.get_xaxis_transform(),ha="center",va="bottom",
                fontweight="bold",color=COLORS[batch],fontsize=11,clip_on=False)
        nstart+=count
    ax.set_xlim(.35,40.65)
    ax.set_xticks(pos,names,rotation=90,fontsize=7 if not compact else 6)
    ax.set_ylabel("log2 PG protein-group intensity" if kind=="protein"
                  else "log2 HELP Heavy intensity",fontsize=11)
    ax.set_title(
        "NYD raw PG: all 4,264 proteins per sample" if kind=="protein"
        else "NYD raw Heavy: all 97 HELP peptides per sample",
        fontsize=14,fontweight="bold",pad=32,loc="left"
    )
    ax.legend(loc="upper right",frameon=False,ncol=2,fontsize=9)
    figstyle(ax)

def plots(pg_log,pgtargets,help_log,helptargets,names,records):
    for kind,mat,targets in (
        ("protein",pg_log,pgtargets),("heavy",help_log,helptargets)):
        fig,ax=plt.subplots(figsize=(20,8.5),layout="constrained")
        show_box(ax,mat,targets,names,records,kind)
        fig.savefig(OUT/IMAGE_FILES[kind],dpi=220,facecolor="white",bbox_inches="tight")
        plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(20,15.3),layout="constrained",sharex=True)
    show_box(axes[0],pg_log,pgtargets,names,records,"protein",compact=True)
    show_box(axes[1],help_log,helptargets,names,records,"heavy",compact=True)
    axes[0].tick_params(axis="x",labelbottom=False,bottom=False)
    fig.suptitle("40 NYD samples · original PG proteins and 97 HELP Heavy peptides",
                 fontsize=17,fontweight="bold")
    fig.savefig(OUT/IMAGE_FILES["combined"],dpi=200,facecolor="white",bbox_inches="tight")
    plt.close(fig)
    for path in IMAGE_FILES.values():
        if (OUT/path).stat().st_size<18000:
            raise ValueError(f"Rendered figure suspiciously small: {path}")

def tables(pg_values,pg_log,pg_targets,help_values,help_log,help_targets,names,records):
    summary=[]
    focus=[]
    for typ,x,log,indices in [
        ("PG_protein",pg_values,pg_log,pg_targets),
        ("HELP_Heavy",help_values,help_log,help_targets),
    ]:
        nrows=x.shape[0]
        for j,sid in enumerate(names):
            good=log[:,j][np.isfinite(log[:,j])]
            q1,med,q3=np.percentile(good,[25,50,75])
            summary.append({
                "source":typ,"sample_id":sid,"batch":records[j]["batch"],
                "total_features":nrows,"finite_positive_features":len(good),
                "missing_invalid_or_nonpositive":nrows-len(good),
                "q1_log2":float(q1),"median_log2":float(med),
                "q3_log2":float(q3),
                "sample_min_log2":float(good.min()),"sample_max_log2":float(good.max()),
            })
            for t in TARGETS:
                i=indices[t["peptide"]]
                v=float(log[i,j])
                focus.append({
                    "source":typ,"sample_id":sid,"batch":records[j]["batch"],
                    "peptide":t["peptide"],"protein":t["gene"],
                    "accession":t["accession"],"highlight_color":t["text_color"],
                    "source_feature_index_zero_based":int(i),
                    "raw_linear_intensity":float(x[i,j]),
                    "log2_intensity":v,
                    "all_feature_sample_median_log2":float(med),
                    "log2_difference_from_sample_median":v-med,
                    "within_sample_IQR":bool(q1<=v<=q3),
                    "percentile_among_valid_features":float(100*np.mean(good<=v)),
                    "background_n_valid":len(good)
                })
    sm=pd.DataFrame(summary)
    tdf=pd.DataFrame(focus)
    sm.to_csv(OUT/"sample_distribution_summary.csv",index=False)
    tdf.to_csv(OUT/"highlighted_two_targets_positions.csv",index=False)
    if len(sm)!=80 or len(tdf)!=160:
        raise AssertionError("Unexpected sample/target diagnostic table row count")
    return sm,tdf

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError("Refusing to overwrite previously generated NYD results")
    records,names=samples_metadata()
    pg,pg_values,pg_log,pg_targets=pg_source(names)
    sequences,help_values,help_log,help_targets=heavy_source(names)
    OUT.mkdir(parents=True,exist_ok=False)
    summary,focus=tables(pg_values,pg_log,pg_targets,help_values,
                         help_log,help_targets,names,records)
    plots(pg_log,pg_targets,help_log,help_targets,names,records)
    pg_matches=[{
        "peptide":t["peptide"],"accession":t["accession"],"gene":t["gene"],
        "row_1based_excluding_header":pg_targets[t["peptide"]]+1,
        "protein_group":pg.iloc[pg_targets[t["peptide"]]]["Protein.Group"],
        "Protein.Ids":pg.iloc[pg_targets[t["peptide"]]]["Protein.Ids"],
        "Genes":pg.iloc[pg_targets[t["peptide"]]]["Genes"],
    } for t in TARGETS]
    counts={}
    for source,count in (("PG_protein",4264),("HELP_Heavy",97)):
        s=summary[summary.source==source]
        total=int(s.total_features.sum())
        valid=int(s.finite_positive_features.sum())
        counts[source]={"total":total,"valid_positive":valid,
                        "missing_invalid_nonpositive":total-valid}
        if total!=count*40:raise AssertionError("Source size mismatch")
    manifest={
        "PG_source":str(PG.relative_to(ROOT)),
        "Heavy_source":str(HEAVY.relative_to(ROOT)),
        "metadata_source":str(SAMPLEINFO.relative_to(ROOT)),
        "source_sha256":{f.name:sha256(f) for f in (PG,HEAVY,SAMPLEINFO)},
        "sample_count":40,"batch_counts":{k:int(sum(s["batch"]==k for s in records)) for k in BATCHES},
        "protein_group_count":len(pg),"HELP_unique_peptide_count":len(sequences),
        "PG_exact_accession_matches":pg_matches,
        "reference_peptides":[t["peptide"] for t in TARGETS],
        "colors":{"STGSWSTLK/P00751":"blue #2364E8",
                  "VYFAGFPR/P07225":"yellow #F2C94C"},
        "original_raw_linear_positive_then_log2":True,
        "no_normalization":True,"no_imputation":True,
        "per_layer_validity":counts,
        "matching":"exact sample identifiers and exact protein accession tokens",
        "sample_order":"NYD_sampleinfo2.xlsx by batch T001/T002/T003/T004, original metadata order within batch",
        "PCA_policy_for_future_work":"Raw-fit PCA, fixed Raw axes/feature selection/centering; transform corrected values into same space",
        "figures":list(IMAGE_FILES.values()),
    }
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    md=[
        "# NYD 40 samples · raw PG + HELP Heavy highlighting two references","",
        "**The two inputs are original, unnormalized signals:**",
        "- PG background: NYD_pg_matrix2.tsv, 4,264 PG protein groups × 40 samples.",
        "- HELP background: Heavy.xlsx, 97 stripped HELP peptides × 40 samples.",
        "- Sample grouping: NYD_sampleinfo2.xlsx, T001 n=16, T002/T003/T004 n=8 each.",
        "- Finite positive raw values are converted to log2 for plotting, with no imputation.",
        "- Every sample box includes ALL available source features, NOT only two targets.",
        "- Blue diamond: STGSWSTLK → CFB / P00751.",
        "- Yellow circle: VYFAGFPR → PROS1 / P07225.",
        "- Two colored markers are drawn at the true source value for each sample.",
        "- Solid vertical lines distinguish batches; figures use identical sample ordering.","",
        "## Two individual PNGs","",
        "### 1. Original PG: 4,264 protein-group intensity boxplots (40 sample boxes)","",
        f"![PG protein boxplots]({IMAGE_FILES['protein']})","",
        "### 2. Heavy.xlsx: 97 HELP peptide intensity boxplots (40 sample boxes)","",
        f"![HELP Heavy boxplots]({IMAGE_FILES['heavy']})","",
        "### Combined overview","",
        f"![Combined original PG and HELP Heavy]({IMAGE_FILES['combined']})","",
        "## Data and target audits","",
        "| Source | Valid positive | Missing/nonpositive/invalid |",
        "| --- | ---: | ---: |",
    ]
    for typ in ("PG_protein","HELP_Heavy"):
        c=counts[typ]
        md.append(f"| {typ} | {c['valid_positive']}/{c['total']} | {c['missing_invalid_nonpositive']} |")
    md.extend(["","| HELP | PG protein group | 1-based PG row |",
               "| --- | --- | ---: |"])
    for t,match in zip(TARGETS,pg_matches):
        md.append(f"| {t['peptide']} | {t['gene']} / {t['accession']} | {match['row_1based_excluding_header']} |")
    md+=["","Both selected proteins and both selected HELP Heavy peptides have",
         "valid positive values in every one of the 40 samples.",
         "The two proteins are each uniquely matched to a PG row by exact accession.",
         "No samples or low-detection samples were removed.","",
         "## Audit exports","",
         "- sample_distribution_summary.csv: 80 rows (2 sources × 40 samples), all-feature medians and quartiles.",
         "- highlighted_two_targets_positions.csv: 160 rows (2 sources × 2 targets × 40 samples),",
         "  original raw value, log2 value, difference from sample box median, percentile, and IQR flag.",
         "- run_manifest.json: original input SHA256 and exact rules.","",
         "## Important PCA rule for future Raw vs Corrected plots","",
         "Fit PCA once on Raw and project Corrected into Raw's exact feature",
         "selection, imputation constants, centering, eigenvectors and axis limits.",
         "Do not separately fit Raw and Corrected PCA when interpreting individual",
         "sample movement. No PCA was calculated for this boxplot-only task.","",
         "Original matrices, Heavy workbook and previous analyses remain unchanged.",""]
    (OUT/"README.md").write_text("\n".join(md),encoding="utf-8")
    print("PASS: 40 samples, 4264 proteins, 97 HELP Heavy peptides, three boxplot PNGs",flush=True)
    print("MATCHED_PG_ROWS",pg_matches,flush=True)
    print("VALID_COUNTS",counts,flush=True)

if __name__=="__main__":
    main()
