#!/usr/bin/env python3
"""NYD YEAST separate peptide/protein multivariate views; do not sum 2 peptides.

Input: same 40-sample YEAST.xlsx, metadata and independently checked precursor
sums used for previous batch-CV results. The 2 *different stripped peptides*
stay in separate features. The precursor charge states WITHIN each peptide
are added in original linear intensity, following prior agreed aggregation.

Peptides: true 40 x 2 PCA, 2x2 Pearson correlation between peptide features,
          2 x 40 standardized log2-intensity hierarchical sample clustering.
Protein:  40 x 1 first (and only) PCA coordinate (PC2 does not exist),
          mathematically trivial 1x1 Pearson self-correlation plus clearly
          labeled *Euclidean distance*, not fictitious pairwise sample
          Pearson correlations, and 1 x 40 sample hierarchical clustering.

Neither a batch effect nor a quality improvement can be established from only
one or two reference features. All inputs remain untouched.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from nyd_yeast_batch_cv import load as parse_source_yeast, digest

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "NYD_data" / "40samples"
OUT = BASE / "YEAST_PCA_correlation_clustering"
CV_INPUT = BASE / "YEAST_batch_CV_results"
PEPTIDES = ("AADALLLK", "VNQIGTLSESIK")
BATCHES = ("T001", "T002", "T003", "T004")
BATCH_COLORS = {
    "T001": "#4E79A7",
    "T002": "#F28E2B",
    "T003": "#59A14F",
    "T004": "#B07AA1",
}
BATCH_N = {"T001":16,"T002":8,"T003":8,"T004":8}
SAMPLE_COLUMNS = {
    "AADALLLK":"AADALLLK_sum_linear",
    "VNQIGTLSESIK":"VNQIGTLSESIK_sum_linear",
    "P00924":"P00924_protein_linear",
}
FILES = (
    "01_peptide_2feature_PCA.png",
    "02_peptide_feature_correlation.png",
    "03_peptide_sample_cluster_heatmap.png",
    "04_protein_1D_PCA.png",
    "05_protein_correlation_limit_and_distance.png",
    "06_protein_sample_cluster_heatmap.png",
)

def csvrows(filename):
    with (CV_INPUT / filename).open("r",encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def save_csv(name, rows):
    if not rows:raise ValueError("No CSV rows: "+name)
    with (OUT / name).open("w",encoding="utf-8-sig",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def save(fig,name):
    fig.savefig(OUT/name,dpi=220,facecolor="white",bbox_inches="tight")
    plt.close(fig)

def check_source():
    raw,invalid=parse_source_yeast()
    if invalid:raise ValueError("Yeast source had missing intensity entries: "+str(invalid))
    previous=csvrows("yeast_per_sample_charge_sums.csv")
    if len(raw)!=40 or len(previous)!=40:
        raise ValueError("Expected 40 source and prepared sample records")
    checks=("AADALLLK_z1","AADALLLK_z2","AADALLLK_sum_linear",
            "VNQIGTLSESIK_z2","VNQIGTLSESIK_z3","VNQIGTLSESIK_sum_linear",
            "P00924_protein_linear")
    by_id={r["sample_id"]:r for r in previous}
    if len(by_id)!=40 or set(by_id)!={r["sample_id"] for r in raw}:
        raise ValueError("Prepared samples no longer match original Excel")
    for r in raw:
        prev=by_id[r["sample_id"]]
        if r["batch"]!=prev["batch"]:
            raise AssertionError("Batch changed for "+r["sample_id"])
        for field in checks:
            if not np.isclose(float(r[field]),float(prev[field]),rtol=0,atol=1e-7):
                raise AssertionError("Original yeast sheet differs from prior charge-summed CSV: "+field)
        for sequence,charges in [("AADALLLK",("z1","z2")),
                                 ("VNQIGTLSESIK",("z2","z3"))]:
            expected=sum(float(r[sequence+"_"+c]) for c in charges)
            if not math.isclose(expected,float(r[sequence+"_sum_linear"]),rel_tol=1e-12):
                raise AssertionError("Peptide precursor sum mismatch")
    manifest=json.loads((CV_INPUT/"run_manifest.json").read_text(encoding="utf-8"))
    for file in ("YEAST.xlsx","NYD_sampleinfo2.xlsx"):
        if digest(BASE/file)!=manifest["sources"][file]:
            raise ValueError("Original workbook changed since prior CV analysis")
    # Immutable deterministic ordering: batches T001 -> T004, order from metadata
    raw.sort(key=lambda r:BATCHES.index(r["batch"]))
    if Counter(r["batch"] for r in raw)!=BATCH_N:
        raise ValueError("Unexpected batch counts")
    names=[r["sample_id"] for r in raw]
    if len(names)!=len(set(names)):
        raise ValueError("Duplicated samples")
    values=np.array([[float(r[SAMPLE_COLUMNS[c]]) for c in
        ("AADALLLK","VNQIGTLSESIK","P00924")] for r in raw],dtype=float)
    if values.shape!=(40,3) or np.any(~np.isfinite(values)) or np.any(values<=0):
        raise ValueError("Unexpected nonpositive yeast intensities")
    return raw,np.log2(values),manifest

def scaled(matrix):
    # Each analyte independently centered and divided by its sample SD.
    # No missing data were observed in the 40 x 3 inputs.
    means=matrix.mean(axis=0)
    sample_sd=matrix.std(axis=0,ddof=1)
    if np.any(sample_sd<=0):
        raise ValueError("Feature with no variance cannot be standardized")
    result=(matrix-means)/sample_sd
    return result,means,sample_sd

def draw_batch_legend(ax,loc="best"):
    handles=[Line2D([],[],ls="",marker="o",markersize=9,
                    markerfacecolor=BATCH_COLORS[b],markeredgecolor="white",
                    label=f"{b} (n={BATCH_N[b]})") for b in BATCHES]
    return ax.legend(handles=handles,loc=loc,frameon=False,fontsize=10)

def cosmetics(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(color="#778899",alpha=.15)
    ax.set_axisbelow(True)

def make_peptide_pca(raw,z):
    if z.shape!=(40,2):raise ValueError("Expected 40x2 peptide matrix")
    pca=PCA(n_components=2,svd_solver="full")
    scores=pca.fit_transform(z)
    variance=pca.explained_variance_ratio_
    fig,ax=plt.subplots(figsize=(10,8),layout="constrained")
    for b in BATCHES:
        idx=np.array([r["batch"]==b for r in raw])
        ax.scatter(scores[idx,0],scores[idx,1],s=64,alpha=.89,
                   color=BATCH_COLORS[b],label=f"{b} (n={BATCH_N[b]})",
                   edgecolors="white",linewidths=.75)
    ax.axhline(0,color="#8294A5",lw=.8,alpha=.5)
    ax.axvline(0,color="#8294A5",lw=.8,alpha=.5)
    ax.set_xlabel(f"PC1 ({variance[0]*100:.2f}% variance)",fontsize=12)
    ax.set_ylabel(f"PC2 ({variance[1]*100:.2f}% variance)",fontsize=12)
    ax.set_title("Yeast peptides: two independent peptide features\n"
                 "PCA of log2 intensities, standardized per peptide",
                 fontsize=14,fontweight="bold",loc="left")
    ax.legend(frameon=False,loc="upper left",bbox_to_anchor=(1.02,1))
    cosmetics(ax)
    save(fig,FILES[0])
    rows=[{"sample_id":r["sample_id"],"batch":r["batch"],
           "PC1":float(scores[i,0]),"PC2":float(scores[i,1])}
          for i,r in enumerate(raw)]
    loads=[{"peptide":name,
            "PC1_loading":float(pca.components_[0,i]),
            "PC2_loading":float(pca.components_[1,i])}
           for i,name in enumerate(PEPTIDES)]
    return rows,loads,variance.tolist()

def make_peptide_correlation(raw,log):
    if log.shape!=(40,2):
        raise ValueError("Expected 40x2 log2 peptide feature input")
    corr=np.corrcoef(log,rowvar=False)
    if corr.shape!=(2,2) or not np.all(np.isfinite(corr)):
        raise ValueError("Peptide correlation invalid")
    fig,axes=plt.subplots(1,5,figsize=(17.5,4.25),layout="constrained")
    sets=[("All samples",np.ones(40,dtype=bool))]
    sets.extend((b,np.array([r["batch"]==b for r in raw])) for b in BATCHES)
    records=[]
    for ax,(label,mask) in zip(axes,sets):
        c=np.corrcoef(log[mask],rowvar=False)
        heat=ax.imshow(c,cmap="RdBu_r",vmin=-1,vmax=1)
        ax.set_xticks([0,1],PEPTIDES,rotation=43,ha="right",fontsize=9)
        ax.set_yticks([0,1],PEPTIDES,fontsize=9)
        ax.set_title(f"{label}\n(n={int(mask.sum())})",fontsize=11,fontweight="bold")
        for i in range(2):
            for j in range(2):
                val=float(c[i,j])
                ax.text(j,i,f"{val:+.2f}",ha="center",va="center",
                        fontsize=13,fontweight="bold",
                        color="white" if abs(val)>.65 else "#283846")
        records.append({"group":label,"n_samples":int(mask.sum()),
                        "Pearson_r_AADALLLK_vs_VNQIGTLSESIK":float(c[0,1])})
    fig.colorbar(heat,ax=axes,fraction=.022,pad=.01,label="Pearson r")
    fig.suptitle("Pearson correlation BETWEEN the two yeast peptides\n"
                 "(40 samples overall and within each batch; log2 intensities)",
                 fontsize=14,fontweight="bold")
    save(fig,FILES[1])
    return records

def reordered_cluster(features, raw, feature_labels, title, outfile, metric="euclidean"):
    """Hierarchical heatmap with dedicated, non-overlapping title/legend/colorbar.

    The dendrogram, batch annotation, and feature heatmap also share the
    *same x-coordinate system*: one leaf = 10 dendrogram coordinate units.
    This preserves true leaf-to-heatmap column alignment.
    """
    if features.ndim != 2 or features.shape[0] != 40:
        raise ValueError("Sample-feature orientation invalid")
    distances = pdist(features, metric=metric)
    if not len(distances) or np.any(~np.isfinite(distances)):
        raise ValueError("Invalid Euclidean distances")

    hierarchical = linkage(distances, method="average", optimal_ordering=True)
    ordered = leaves_list(hierarchical)
    order = [raw[int(i)]["sample_id"] for i in ordered]
    bnames = [raw[int(i)]["batch"] for i in ordered]
    matrix = features[ordered, :].T

    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    n_features = len(feature_labels)
    fig = plt.figure(figsize=(20, 10.2 if n_features == 2 else 9.2))
    # Legend occupies its own row at the TOP rather than beneath 40 long
    # vertically rotated sample IDs. The bottom 25% is kept free for labels.
    gs = fig.add_gridspec(
        4, 2,
        width_ratios=(38, 1.5),
        height_ratios=(0.25, 2.4, 0.19, 1.05 if n_features == 2 else 0.65),
        left=0.095, right=0.90, top=0.895, bottom=0.255,
        wspace=0.17, hspace=0.16,
    )

    # 1: Dedicated legend row: cannot collide with x-axis sample labels.
    ax_legend = fig.add_subplot(gs[0, :])
    ax_legend.axis("off")
    legend = ax_legend.legend(
        handles=[Patch(facecolor=BATCH_COLORS[b],
                       label=f"{b} (n={BATCH_N[b]})") for b in BATCHES],
        ncol=4, loc="center", bbox_to_anchor=(0.5, 0.5),
        frameon=False, fontsize=11,
    )

    # 2: Dendrogram: 40 leaves at x=5,15,...,395; align shared axes exactly.
    ax_dend = fig.add_subplot(gs[1, 0])
    plotted = dendrogram(
        hierarchical, ax=ax_dend, no_labels=True, color_threshold=0,
        above_threshold_color="#64768B",
    )
    if list(plotted["leaves"]) != [int(x) for x in ordered]:
        raise AssertionError("Dendrogram order differs from heatmap leaf order")
    ax_dend.set_xlim(0, 400)
    ax_dend.set_xticks([])
    ax_dend.set_ylabel("Euclidean\nlinkage height", fontsize=10)
    ax_dend.spines["top"].set_visible(False)
    ax_dend.spines["right"].set_visible(False)
    ax_dend.spines["bottom"].set_visible(False)

    # 3: Per-leaf batch annotation on the identical 0..400 horizontal scale.
    ax_batch = fig.add_subplot(gs[2, 0], sharex=ax_dend)
    batch_row = np.array([[BATCHES.index(b) for b in bnames]])
    ax_batch.imshow(
        batch_row, cmap=ListedColormap([BATCH_COLORS[b] for b in BATCHES]),
        vmin=-0.5, vmax=3.5, interpolation="nearest",
        aspect="auto", extent=(0, 400, 1, 0),
    )
    ax_batch.set_xlim(0, 400)
    ax_batch.set_xticks([])
    ax_batch.set_yticks([])
    for sp in ax_batch.spines.values():
        sp.set_visible(False)

    # 4: Feature matrix: use dendrogram coordinates for exact x alignment.
    ax_heat = fig.add_subplot(gs[3, 0], sharex=ax_dend)
    limit = max(2.0, float(np.max(np.abs(matrix))) * 1.03)
    img = ax_heat.imshow(
        matrix, cmap="RdBu_r", vmin=-limit, vmax=limit,
        aspect="auto", interpolation="nearest",
        extent=(0, 400, n_features - 0.5, -0.5),
    )
    ax_heat.set_xlim(0, 400)
    ax_heat.set_yticks(np.arange(n_features), feature_labels, fontsize=11)
    ax_heat.set_xticks(np.arange(40) * 10 + 5)
    ax_heat.set_xticklabels(order, rotation=90, fontsize=7)
    ax_heat.tick_params(axis="x", length=0, pad=3)
    ax_heat.tick_params(axis="y", length=0, pad=7)

    # Dedicated colorbar column, independent of all three data axes.
    ax_cbar = fig.add_subplot(gs[3, 1])
    cbar = fig.colorbar(img, cax=ax_cbar)
    cbar.set_label("Feature z-score of log2 intensity", fontsize=10)

    fig.suptitle(title, fontsize=15, fontweight="bold", y=0.978)

    # Geometrical audit of the *rendered* layout before saving the PNG.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    legend_bounds = legend.get_window_extent(renderer)
    heat_bounds = ax_heat.get_tightbbox(renderer)
    if legend_bounds.overlaps(heat_bounds):
        raise AssertionError("Batch legend still overlaps the sample heatmap labels")
    fig_bounds = fig.bbox
    if (heat_bounds.x0 < -1 or heat_bounds.y0 < -1
            or heat_bounds.x1 > fig_bounds.x1 + 1):
        raise AssertionError("Heatmap and rotated labels exceed figure bounds")
    if ax_cbar.get_window_extent(renderer).overlaps(ax_heat.get_window_extent(renderer)):
        raise AssertionError("Colorbar intersects the heatmap body")

    save(fig, outfile)
    return [
        {"cluster_leaf_order": i + 1, "sample_id": sample_id, "batch": batch,
         "sample_linkage": "average",
         "sample_distance": "Euclidean (standardized log2 input)"}
        for i, (sample_id, batch) in enumerate(zip(order, bnames))
    ]


def make_protein_1d(raw,z):
    if z.shape!=(40,1):
        raise ValueError("P00924 is a single protein, expected 40 x 1")
    pca=PCA(n_components=1,svd_solver="full")
    scores=pca.fit_transform(z).ravel()
    fig,ax=plt.subplots(figsize=(11,6),layout="constrained")
    for row,g in enumerate(BATCHES):
        ids=[j for j,s in enumerate(raw) if s["batch"]==g]
        yy=np.full(len(ids),len(BATCHES)-row)
        # Display y categories, never claim the y-axis is PC2.
        ax.scatter(scores[ids],yy,s=78,color=BATCH_COLORS[g],
                   edgecolors="white",linewidth=.9,alpha=.92)
    ax.set_yticks([4,3,2,1],
                  [f"{b} (n={BATCH_N[b]})" for b in BATCHES])
    ax.set_ylim(.4,4.65)
    ax.set_xlabel("PC1 (100.0% of variance; single feature)",fontsize=12)
    ax.set_ylabel("Batch (categorical; NOT PC2)",fontsize=11)
    ax.set_title("Yeast protein P00924: only one PCA dimension is possible\n"
                 "PCA of standardized log2 protein intensity",
                 fontsize=14,fontweight="bold",loc="left")
    cosmetics(ax)
    save(fig,FILES[3])
    return [{"sample_id":r["sample_id"],"batch":r["batch"],
             "PC1":float(scores[i]),"PC2":"NOT_DEFINED_SINGLE_FEATURE"}
            for i,r in enumerate(raw)]

def make_protein_correlation_limit(raw,z):
    if z.shape!=(40,1):
        raise ValueError("Expected 40x1 protein feature")
    x=z.ravel()
    # sample correlation with only ONE feature is not defined; preserve that
    # distinction instead of outputting spurious 40x40 Pearson correlations.
    dist=squareform(pdist(z,metric="euclidean"))
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(15,7.2),
                               gridspec_kw={"width_ratios":[.75,2.1]},
                               layout="constrained")
    ax1.imshow([[1.0]],cmap="Blues",vmin=0,vmax=1)
    ax1.set_xticks([0],["P00924"])
    ax1.set_yticks([0],["P00924"])
    ax1.text(0,0,"1.00",ha="center",va="center",fontsize=20,
             fontweight="bold",color="white")
    ax1.set_title("Protein-feature Pearson correlation\n"
                  "Self-correlation only (trivial)",fontsize=12,fontweight="bold")
    ax1.set_xlabel("No second protein to correlate",fontsize=10)
    order=[i for b in BATCHES for i,r in enumerate(raw) if r["batch"]==b]
    dd=dist[np.ix_(order,order)]
    image=ax2.imshow(dd,cmap="magma_r",aspect="equal",interpolation="nearest")
    ax2.set_xticks(range(40),[raw[i]["sample_id"] for i in order],
                   rotation=90,fontsize=5)
    ax2.set_yticks(range(40),[raw[i]["sample_id"] for i in order],
                   fontsize=5)
    ax2.set_title("Sample-pair distance from P00924\n"
                  "Euclidean on standardized log2 signal (NOT Pearson r)",
                  fontsize=12,fontweight="bold")
    fig.colorbar(image,ax=ax2,fraction=.048,pad=.025,
                 label="Absolute standardized log2 difference")
    fig.suptitle("Why sample Pearson correlation is not estimable from ONE protein",
                 fontsize=15,fontweight="bold")
    save(fig,FILES[4])
    return [{"sample_id_a":raw[i]["sample_id"],
             "sample_id_b":raw[j]["sample_id"],
             "abs_zscore_difference":float(dist[i,j])}
            for i in range(40) for j in range(40)]

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError("Refusing to replace preexisting result directory: "+str(OUT))
    raw,log_data,source=check_source()
    if log_data.shape!=(40,3):
        raise ValueError("Wrong input shape")
    peptide_std,pmean,pstd=scaled(log_data[:,:2])
    protein_std,pmean1,pstd1=scaled(log_data[:,2:3])
    # Across-feature standardization makes 2 peptides contribute equally;
    # peptides remain TWO separate columns, never added to each other.
    OUT.mkdir(parents=True,exist_ok=True)

    pep_scores,loadings,ratio=make_peptide_pca(raw,peptide_std)
    corr=make_peptide_correlation(raw,log_data[:,:2])
    peptree=reordered_cluster(peptide_std,raw,PEPTIDES,
         "Two yeast peptide features: sample hierarchical clustering (average / Euclidean)",
         FILES[2])
    prot_scores=make_protein_1d(raw,protein_std)
    distance=make_protein_correlation_limit(raw,protein_std)
    prottree=reordered_cluster(protein_std,raw,["P00924"],
         "Yeast protein P00924: single-feature sample clustering (Euclidean)",
         FILES[5])
    save_csv("peptide_pca_scores.csv",pep_scores)
    save_csv("peptide_pca_loadings.csv",loadings)
    save_csv("peptide_feature_pearson_by_batch.csv",corr)
    save_csv("peptide_sample_cluster_order.csv",peptree)
    save_csv("protein_1D_pca_scores.csv",prot_scores)
    save_csv("protein_sample_euclidean_distance_long.csv",distance)
    save_csv("protein_sample_cluster_order.csv",prottree)
    feature_matrix=[
        {"sample_id":r["sample_id"],"batch":r["batch"],
         "AADALLLK_log2":float(log_data[i,0]),
         "VNQIGTLSESIK_log2":float(log_data[i,1]),
         "P00924_protein_log2":float(log_data[i,2]),
         "AADALLLK_zscore_log2":float(peptide_std[i,0]),
         "VNQIGTLSESIK_zscore_log2":float(peptide_std[i,1]),
         "P00924_protein_zscore_log2":float(protein_std[i,0]),
         }
         for i,r in enumerate(raw)
    ]
    save_csv("feature_matrix_40samples.csv",feature_matrix)
    manifest={
      "source":"NYD_data/40samples/YEAST.xlsx + NYD_sampleinfo2.xlsx",
      "prepared_intensities":"NYD_data/40samples/YEAST_batch_CV_results/yeast_per_sample_charge_sums.csv",
      "source_sha256":source["sources"],
      "peptide_features":list(PEPTIDES),
      "protein_features":["P00924"],
      "sample_count":40,
      "batch_counts":BATCH_N,
      "charge_states_aggregated_within_each_unique_peptide":{"AADALLLK":[1,2],"VNQIGTLSESIK":[2,3]},
      "different_peptides_summed_together":False,
      "data_preprocessing":"log2 intensity then per-feature centered+sampleSD-scaled (z-score)",
      "pca_peptide_explained_variance_ratio":ratio,
      "pca_protein_only_one_component":True,
      "sample_Pearson_correlation_with_one_protein":"undefined",
      "peptide_correlation_kind":"Pearson between two peptide features across samples (overall and per batch)",
      "protein_feature_correlation":"1x1 trivial self-correlation; companion sample heatmap shows Euclidean distance",
      "clustering":"average hierarchical linkage on Euclidean distance of standardized log2 features",
      "no_imputation":True,"no_normalization":True,
      "figures":list(FILES),
    }
    (OUT/"run_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    text=[
      "# NYD 40 samples: yeast peptide vs protein PCA, correlations and clustering","",
      "This is an exploratory analysis of YEAST.xlsx using NYD_sampleinfo2.xlsx",
      "(T001 n=16, T002/T003/T004 n=8). Original source Excel files and previous",
      "CV results are preserved without change.","",
      "## Exact features used","",
      "**Peptides (40 samples x 2 independent features):** AADALLLK, VNQIGTLSESIK.",
      "A peptide's own charge-state precursors are linearly summed first (1+2 and 2+3",
      "respectively), as in the previously verified charge-summed peptide CV.",
      "**The two different peptides are NEVER added together.**",
      "**Protein (40 samples x 1 feature):** P00924, directly from the protein row.",
      "For exploratory PCA and clustering, each independent feature is transformed to",
      "log2 and then centered/scaled to unit variance. The CV values previously",
      "reported remain on their original linear intensity scale.","",
      "## 1. Two-peptide PCA","",
      "True two-feature PCA; points are samples, colors are metadata batches.",
      f"Explained variance: PC1 **{ratio[0]*100:.2f}%**, PC2 **{ratio[1]*100:.2f}%**.","",
      "![Peptide PCA](01_peptide_2feature_PCA.png)","",
      "## 2. Peptide correlation heatmaps","",
      "Pearson correlations are computed **between the two unique peptide features**",
      "across all 40 samples and within each batch. This is a 2x2 FEATURE correlation",
      "matrix per group, not an invalid 40x40 sample correlation based on 2 values.",
      "The 2x2 diagonal is self-correlation (1.00).","",
      "![Peptide feature correlation](02_peptide_feature_correlation.png)","",
      "## 3. Peptide sample cluster heatmap","",
      "Each column is a sample; rows are the two independently standardized peptide",
      "signals. Hierarchical clustering uses sample Euclidean distances and average",
      "linkage. The top stripe shows original batch labels in clustered order.","",
      "![Peptide sample clustering](03_peptide_sample_cluster_heatmap.png)","",
      "## 4. Protein PCA: mathematical limitation","",
      "**With one protein feature only PC1 exists.** A conventional PC1-versus-PC2",
      "scatter would invent information: PC2 is undefined. The figure therefore plots",
      "true one-dimensional PC1 scores on the horizontal axis and **batch category**",
      "on the vertical axis (not PC2). PC1 accounts for all available variance.","",
      "![Protein one-dimensional PCA](04_protein_1D_PCA.png)","",
      "## 5. Protein correlation: mathematical limitation and substitute","",
      "With only one protein, the protein-feature Pearson matrix is 1x1 and 1.00",
      "(self-correlation), which tells nothing about between-sample similarity.",
      "A 40x40 **sample Pearson correlation is undefined with one feature**.",
      "To make this view useful, the right-hand heatmap instead shows the",
      "**Euclidean distance between samples** on the single standardized log2",
      "protein signal; its color is NOT Pearson correlation.","",
      "![Protein self-correlation and distance](05_protein_correlation_limit_and_distance.png)","",
      "## 6. Protein sample cluster heatmap","",
      "Hierarchical clustering of the 40 samples on one standardized log2 P00924",
      "signal using Euclidean distance and average linkage. With only one protein",
      "this clustering can be fragile and must not be taken as evidence that",
      "biological sample states or technical batches were reliably inferred.","",
      "![Protein sample clustering](06_protein_sample_cluster_heatmap.png)","",
      "## Exported audit tables","",
      "- feature_matrix_40samples.csv: all 40 sample IDs, batches, log2 and z-scored features.",
      "- peptide_pca_scores.csv and peptide_pca_loadings.csv.",
      "- peptide_feature_pearson_by_batch.csv: five 2x2 correlation summaries.",
      "- peptide_sample_cluster_order.csv and protein_sample_cluster_order.csv.",
      "- protein_1D_pca_scores.csv: PC2 explicitly shown as undefined.",
      "- protein_sample_euclidean_distance_long.csv: 40x40 distance matrix (long-form).",
      "- run_manifest.json: source hashes, feature declarations and method limitations.","",
      "**Interpretation warning:** Two peptides and one protein are insufficient",
      "to assess global proteomic batch effects. Evaluate batch separation",
      "descriptively, not as an instrument-quality conclusion.","",
    ]
    (OUT/"README.md").write_text("\n".join(text),encoding="utf-8")
    print("PASS: peptide PCs",ratio,"Pearson r overall",corr[0]["Pearson_r_AADALLLK_vs_VNQIGTLSESIK"],flush=True)
    print("PASS: protein only one PC and no pairwise Pearson; valid Euclidean clustering",flush=True)

if __name__=="__main__":
    main()
