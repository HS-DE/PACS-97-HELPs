#!/usr/bin/env python3
"""NYD yeast 40-sample batch CV, two charge-summed peptides plus one protein."""
from pathlib import Path
import csv, hashlib, json, math, statistics
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from filter_samples import readbook, tab

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/"NYD_data"/"40samples"
YEAST=DIR/"YEAST.xlsx"
META=DIR/"NYD_sampleinfo2.xlsx"
OUT=DIR/"YEAST_batch_CV_results"
BATCHES={"T001":16,"T002":8,"T003":8,"T004":8}
FEATURES=["AADALLLK","VNQIGTLSESIK","P00924_protein"]
SOURCE={"AADALLLK":("AADALLLK",1,2),"VNQIGTLSESIK":("VNQIGTLSESIK",2,3)}
COLORS={"AADALLLK":"#2C71B8","VNQIGTLSESIK":"#D88934","P00924_protein":"#269479"}
BCOLORS=["#4E79A7","#59A14F","#F28E2B","#B07AA1"]

def rows(path):
    ss,sheets=readbook(path)
    if len(sheets)!=1:raise ValueError(f"{path.name} expected one sheet")
    return tab(sheets[0],ss)

def positive(value):
    try:
        x=float(value)
        return x if math.isfinite(x) and x>0 else np.nan
    except (TypeError,ValueError):return np.nan

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write_csv(name,lines):
    with (OUT/name).open("w",newline="",encoding="utf-8-sig") as f:
        writer=csv.DictWriter(f,fieldnames=list(lines[0]))
        writer.writeheader();writer.writerows(lines)

def load():
    mr=rows(META)
    assert mr[0].get("A")=="样品编号" and mr[0].get("C")=="batch"
    metadata=[]
    for r in mr[1:]:
        sid=r.get("A","").strip()
        if not sid:continue
        batch=r.get("C","").strip()
        if sid!=r.get("B","").strip() or batch not in BATCHES:
            raise ValueError(f"Invalid metadata row {r}")
        if batch!=r.get("D") or batch!=r.get("E"):
            raise ValueError(f"Metadata batch and group disagree for {sid}")
        metadata.append({"sample_id":sid,"batch":batch})
    names=[r["sample_id"] for r in metadata]
    if len(names)!=40 or len(set(names))!=40 or Counter(r["batch"] for r in metadata)!=BATCHES:
        raise ValueError("40 sample IDs and batch sizes do not match expected values")
    mr=rows(YEAST)
    phs=[(i,r) for i,r in enumerate(mr) if r.get("A")=="Protein.Group" and r.get("B")=="Stripped.Sequence"]
    pgs=[(i,r) for i,r in enumerate(mr) if r.get("A")=="Protein.Group" and r.get("B") in names]
    if len(phs)!=1 or len(pgs)!=1:raise ValueError("Missing peptide/protein section headers")
    pi,ph=phs[0]; qi,protein_head=pgs[0]
    if qi<=pi:raise ValueError("Wrong section order")
    peptide_cols={v:c for c,v in ph.items() if c not in ("A","B","C","D")}
    protein_cols={v:c for c,v in protein_head.items() if c!="A"}
    if len(peptide_cols)!=40 or set(peptide_cols)!=set(names):
        raise ValueError("Peptide section sample IDs mismatch")
    if len(protein_cols)!=40 or set(protein_cols)!=set(names):
        raise ValueError("Protein section sample IDs mismatch")
    precursors=[r for r in mr[pi+1:qi] if r.get("B") and r.get("A")=="P00924"]
    keys={(r["B"],int(float(r["C"]))) for r in precursors}
    if len(precursors)!=4 or keys!={("AADALLLK",1),("AADALLLK",2),
                                      ("VNQIGTLSESIK",2),("VNQIGTLSESIK",3)}:
        raise ValueError(f"Unexpected precursor rows: {keys}")
    ptr={(r["B"],int(float(r["C"]))):r for r in precursors}
    proteins=[r for r in mr[qi+1:] if r.get("A")]
    if len(proteins)!=1 or proteins[0].get("A")!="P00924":
        raise ValueError("Protein part must contain exactly one P00924 row")
    protein=proteins[0]
    samples=[];missing=[]
    for r in metadata:
        sid=r["sample_id"]
        record=dict(r)
        for sequence,c1,c2 in SOURCE.values():
            values=[positive(ptr[(sequence,c)].get(peptide_cols[sid],"")) for c in (c1,c2)]
            for c,v in zip((c1,c2),values):
                record[sequence+"_z"+str(c)]=v
                if not np.isfinite(v):missing.append((sid,sequence,c))
            good=[v for v in values if np.isfinite(v)]
            record[sequence+"_sum_linear"]=float(sum(good)) if good else np.nan
            record[sequence+"_measured_charges"]=len(good)
        value=positive(protein.get(protein_cols[sid],""))
        record["P00924_protein_linear"]=value
        if not np.isfinite(value):missing.append((sid,"P00924_protein","protein"))
        samples.append(record)
    return samples,missing

def arr(records,feature):
    field="P00924_protein_linear" if feature=="P00924_protein" else feature+"_sum_linear"
    return np.array([r[field] for r in records],float)

def stats(samples):
    rows=[]
    for batch,n in BATCHES.items():
        rec=[r for r in samples if r["batch"]==batch]
        for feat in FEATURES:
            a=arr(rec,feat)
            a=a[np.isfinite(a)&(a>0)]
            count=len(a)
            need=max(3,math.ceil(.70*n))
            mean=float(np.mean(a)) if count else np.nan
            sd=float(np.std(a,ddof=1)) if count>1 else np.nan
            cv=float(sd/mean*100) if count>=need and mean>=1e-8 else np.nan
            rows.append({"batch":batch,"analyte":feat,
                "n_batch_samples":n,"n_positive_detected":count,
                "n_min_required":need,"detection_rate_pct":100*count/n,
                "mean_linear":mean,"sample_sd_linear_ddof1":sd,
                "CV_pct":cv,"median_linear":float(np.median(a)) if count else np.nan,
                "median_log2":float(np.median(np.log2(a))) if count else np.nan})
    if any(not np.isfinite(r["CV_pct"]) for r in rows):
        raise ValueError("One or more group CVs invalid due to inadequate data")
    return rows

def save(fig,name):
    fig.savefig(OUT/name,dpi=220,bbox_inches="tight",facecolor="white")
    plt.close(fig)

def style(ax):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.grid(axis="y",alpha=.19)
    ax.set_axisbelow(True)

def plots(samples,stat):
    cv={(r["batch"],r["analyte"]):r["CV_pct"] for r in stat}
    groups=list(BATCHES)
    # Chart 1: exact CVs by batch, grouped bars + labeled heatmap
    fig,(ax,bx)=plt.subplots(1,2,figsize=(15.6,6.5),layout="constrained")
    maximum=max(cv.values())
    pos=np.arange(4)
    for k,feat in enumerate(FEATURES):
        x=pos+(k-1)*.23
        y=[cv[b,feat] for b in groups]
        ax.bar(x,y,.215,color=COLORS[feat],label=feat,alpha=.87)
        for xx,yy in zip(x,y):
            ax.text(xx,yy+maximum*.018,f"{yy:.1f}",ha="center",fontsize=9)
    ax.set_xticks(pos,[f"{b}\nn={BATCHES[b]}" for b in groups])
    ax.set_ylim(0,maximum*1.28)
    ax.set_title("CV by analyte",loc="left",fontweight="bold")
    ax.set_ylabel("CV (%) · unlogged linear intensities")
    ax.legend(frameon=False,fontsize=9)
    style(ax)
    matrix=np.array([[cv[b,feat] for b in groups] for feat in FEATURES])
    image=bx.imshow(matrix,cmap="YlGnBu",aspect="auto",vmin=0,vmax=maximum*1.15)
    bx.set_xticks(range(4),groups)
    bx.set_yticks(range(3),["AADALLLK (z1+z2)","VNQIGTLSESIK (z2+z3)","P00924 protein"],fontsize=9)
    bx.set_title("CV heatmap",loc="left",fontweight="bold")
    for i in range(3):
        for j in range(4):
            yy=matrix[i,j]
            bx.text(j,i,f"{yy:.2f}%",ha="center",va="center",fontsize=10,
                    color="white" if yy>maximum*.58 else "#183040",fontweight="bold")
    fig.colorbar(image,ax=bx,fraction=.048,pad=.03,label="CV (%)")
    fig.suptitle("NYD YEAST | 40 samples across 4 batches",fontsize=16,fontweight="bold")
    save(fig,"01_batch_cv_bars_heatmap.png")
    # Chart 2: violin + box + scatter, 3 analyte CVs per batch
    fig,ax=plt.subplots(figsize=(11,7),layout="constrained")
    data=[[cv[b,f] for f in FEATURES] for b in groups]
    vio=ax.violinplot(data,positions=np.arange(1,5),widths=.7,
                      showmedians=False,showextrema=False)
    for artist,col in zip(vio["bodies"],BCOLORS):
        artist.set_facecolor(col);artist.set_alpha(.23);artist.set_edgecolor(col)
    boxes=ax.boxplot(data,positions=np.arange(1,5),widths=.22,patch_artist=True,
                     showfliers=False,manage_ticks=False,
                     medianprops={"color":"#263B4A","linewidth":2})
    for box,col in zip(boxes["boxes"],BCOLORS):
        box.set_facecolor(col);box.set_alpha(.35)
    for ix,b in enumerate(groups,1):
        for offset,feat in zip((-.13,0,.13),FEATURES):
            ax.scatter(ix+offset,cv[b,feat],s=115,c=COLORS[feat],
                       edgecolor="white",linewidth=.9,zorder=5)
        ax.text(ix,max(data[ix-1])+maximum*.035,
                f"Median={statistics.median(data[ix-1]):.2f}%",
                ha="center",fontsize=10,fontweight="bold")
    handles=[Line2D([],[],marker="o",color="none",markerfacecolor=COLORS[f],
                    markeredgecolor="white",markersize=11,label=f) for f in FEATURES]
    ax.legend(handles=handles,loc="upper left",bbox_to_anchor=(1,1),frameon=False)
    ax.set_xticks(range(1,5),[f"{b}\n(n={BATCHES[b]} samples)" for b in groups])
    ax.set_xlim(.5,4.5);ax.set_ylim(0,maximum*1.29)
    ax.set_ylabel("CV (%)")
    ax.set_title("Batch CV · violin / box / three measured analytes",
                 loc="left",fontweight="bold",fontsize=15)
    style(ax)
    fig.text(.07,.012,"Each batch has just 3 feature-level CV values; violin shape is descriptive, not a population estimate.",
             fontsize=9,color="#6B7685")
    save(fig,"02_batch_cv_violin_box_scatter.png")
    # Chart 3: sample-level raw-signal boxplots, log2 only on DISPLAY
    fig,axes=plt.subplots(1,3,figsize=(17,6),layout="constrained")
    for ax,feat in zip(axes,FEATURES):
        values=[np.log2(arr([r for r in samples if r["batch"]==b],feat)) for b in groups]
        boxes=ax.boxplot(values,positions=range(1,5),widths=.53,patch_artist=True,
                         showfliers=False,manage_ticks=False)
        for box,col in zip(boxes["boxes"],BCOLORS):
            box.set_facecolor(col);box.set_alpha(.58)
        for ix,(b,xx) in enumerate(zip(groups,values),1):
            ax.scatter(ix+np.linspace(-.17,.17,len(xx)),xx,s=28,c=BCOLORS[ix-1],
                       edgecolor="white",linewidth=.3,zorder=4)
        ax.set_title(feat,loc="left",fontweight="bold")
        ax.set_xticks(range(1,5),[f"{b}\nn={BATCHES[b]}" for b in groups])
        if feat==FEATURES[0]:ax.set_ylabel("log2 signal (only for display)")
        style(ax)
    fig.suptitle("YEAST intensities by batch · samples as points",fontsize=16,fontweight="bold")
    save(fig,"03_batch_log2_intensity_boxplots.png")

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError("Output exists. Refusing to overwrite previously generated results.")
    samples,missing=load()
    st=stats(samples)
    OUT.mkdir(parents=True,exist_ok=True)
    write_csv("yeast_per_sample_charge_sums.csv",samples)
    write_csv("yeast_batch_cv_long.csv",st)
    wide=[]
    for batch,n in BATCHES.items():
        row={"batch":batch,"n_samples":n}
        for f in FEATURES:row[f]=next(x["CV_pct"] for x in st if x["batch"]==batch and x["analyte"]==f)
        row["median_three_features_CV_pct"]=float(statistics.median(row[f] for f in FEATURES))
        wide.append(row)
    write_csv("yeast_batch_cv_wide.csv",wide)
    plots(samples,st)
    manifest={
        "sources":{"YEAST.xlsx":digest(YEAST),"NYD_sampleinfo2.xlsx":digest(META)},
        "batch_sizes":BATCHES,"unique_peptides":["AADALLLK","VNQIGTLSESIK"],
        "protein_accession":"P00924","precursors":{"AADALLLK":[1,2],"VNQIGTLSESIK":[2,3]},
        "sample_level_peptides":"sum of observed positive charge-state intensities on linear scale",
        "cv_definition":"100 * sample sd(ddof=1) / arithmetic mean (linear)",
        "cv_min_n":"max(3,ceil(70% of batch sample count))",
        "log2_only_for_display":True,"no_imputation":True,"no_normalization":True,
        "missing_precursor_or_protein_entries":missing}
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    lines=["# NYD YEAST: 40-sample batch CV and visualizations","",
        "Original workbooks were not changed. Data were joined by exact sample ID.",
        "The metadata has 16 T001 samples and 8 samples each for T002/T003/T004.",
        "The peptide portion has AADALLLK precursor charge 1 and 2, and",
        "VNQIGTLSESIK charge 2 and 3. For each sample, intensities are summed",
        "in linear units across the two charge rows for each stripped sequence.",
        "P00924 protein uses the single protein row without summation.","",
        "CV = 100 x sample SD (ddof=1) / arithmetic mean, on ORIGINAL LINEAR intensity.",
        "Only finite positive measurements contribute. No missing-value imputation.",
        "Valid CV requires >= max(3,ceil(0.7 x batch size)) observed values.",
        "Log2 is only used for visualizing original intensities.","",
        "## Batch CV (%)","",
        "| Batch | n | AADALLLK (z1+z2) | VNQIGTLSESIK (z2+z3) | P00924 protein | Median of 3 |",
        "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for r in wide:
        lines.append(f'| {r["batch"]} | {r["n_samples"]} | {r["AADALLLK"]:.2f}% | '
              f'{r["VNQIGTLSESIK"]:.2f}% | {r["P00924_protein"]:.2f}% | '
              f'{r["median_three_features_CV_pct"]:.2f}% |')
    lines.extend(["",f"Missing/invalid original values: **{len(missing)}**.","",
       "## Figures","",
       "### 1. CV grouped bars and heatmap","",
       "![Batch CV summary](01_batch_cv_bars_heatmap.png)","",
       "### 2. CV violin / boxplot / scatter","",
       "![Batch CV violin](02_batch_cv_violin_box_scatter.png)","",
       "Note: Each violin summarizes just three CV values per batch. Its smoothed",
       "shape is not a statistically robust distribution estimate.","",
       "### 3. Sample-level log2 intensity boxplots (four batches)","",
       "![Sample intensity boxplots](03_batch_log2_intensity_boxplots.png)","",
       "## Output CSVs","",
       "- yeast_per_sample_charge_sums.csv: original precursor values and each summed peptide.",
       "- yeast_batch_cv_long.csv: n, mean, sample SD, detection, CV and median by feature/batch.",
       "- yeast_batch_cv_wide.csv: four-batch comparison.",
       "- run_manifest.json: source SHA256, calculation policy and missing-value audit.","",
       "These are within-batch sample CVs and should not automatically be called",
       "pure technical CV unless the samples are known homogeneous technical QCs.",""])
    (OUT/"README.md").write_text("\n".join(lines),encoding="utf-8")
    print("BATCH_CV_RESULTS",json.dumps(wide))
    print("MISSING_INPUTS",missing)
    print("SUCCESS generated 12 CVs, 40-sample summary, 3 PNGs")
if __name__=="__main__":
    main()
