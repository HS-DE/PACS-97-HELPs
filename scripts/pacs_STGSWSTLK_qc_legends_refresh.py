#!/usr/bin/env python3
"""Refresh only the four latest PACS STGSWSTLK QC normalization figures.

The original numerical tables and matrices are NEVER overwritten. Legends
are placed outside the plot axes (not floating over data), all PCA views share
the same Raw-fitted coordinates, and original image paths remain usable.
"""
import numpy as np
import pandas as pd

import pacs_STGSWSTLK_qc_heavy_two_signs as experiment
from yeast_protein_qc import GROUPS

def main():
    root=experiment.OUT
    if not (root/"run_manifest.json").is_file():
        raise FileNotFoundError("PACS STGSWSTLK experiment not completed")
    pg,names,info,raw,heavy,pooled,internal=experiment.original()
    h,g,p,refs,_=experiment.refs_and_factors(names,info,heavy,pooled)
    variants={"Raw":raw}
    for reference in experiment.REFERENCES:
        for direction in experiment.DIRECTIONS:
            label=f"{reference}/{direction}"
            adjusted=raw+(1 if direction=="plus_fi" else -1)*refs[reference][None,:]
            path=root/reference/f"corrected_pg_log2_{direction}.tsv"
            if not path.is_file():
                raise FileNotFoundError(path)
            # Validate regenerated plots represent exactly the saved matrix.
            table=pd.read_csv(path,sep="\t",low_memory=False)
            if table.shape!=(6675,80) or list(table.columns[5:])!=names:
                raise ValueError(f"Matrix misalignment {label}")
            matrix=table[names].to_numpy(float)
            if not np.array_equal(np.isfinite(matrix),np.isfinite(adjusted)):
                raise AssertionError(f"Matrix missingness changed {label}")
            if np.nanmax(abs(matrix-adjusted))>1e-9:
                raise AssertionError(f"Saved matrix does not match factors {label}")
            variants[label]=adjusted
    cv=pd.read_csv(root/"group_protein_CV_summary.csv")
    qc=pd.read_csv(root/"pooled_QC_within_plate_CV.csv")
    qc["plate"]=qc["plate"].astype(str)
    if len(cv)!=20 or len(qc)!=15:
        raise ValueError("Existing summary incomplete")
    experiment.reference_plot(h,names,info,pooled,g,p)
    experiment.cv_figure(cv,qc)
    for subset,filename in (
        (["HC","S"],"04_raw_fitted_PCA_54_study.png"),
        (GROUPS,"05_raw_fitted_PCA_all_75.png")
    ):
        experiment.pca_chart(variants,names,info,subset,filename)
    image_names=(
        "01_QC_Heavy_references.png",
        "02_CV_group_and_QC_plate_comparison.png",
        "04_raw_fitted_PCA_54_study.png",
        "05_raw_fitted_PCA_all_75.png")
    for name in image_names:
        if (root/name).stat().st_size<20000:
            raise AssertionError(f"Refreshed image suspiciously small: {name}")
    print("PASS: four PACS recent experiment figures refreshed with legends outside axes; all numerical results untouched.")

if __name__=="__main__":
    main()
