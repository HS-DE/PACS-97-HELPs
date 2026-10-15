#!/usr/bin/env python3
"""Delete and redraw ONLY two overlapping NYD YEAST hierarchical cluster PNGs.

This one-time repair verifies all previously published numerical results and
other figures remain bit-for-bit intact, and verifies newly computed sample
clustering order exactly matches previously published CSV rows.
"""
import csv
import hashlib
from pathlib import Path

from nyd_yeast_multivariate import (
    BASE, OUT, FILES, PEPTIDES,
    check_source, scaled, reordered_cluster,
)

TARGETS = (FILES[2], FILES[5])


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))


def main():
    if not OUT.is_dir():
        raise FileNotFoundError(f"Expected previous result folder: {OUT}")
    for path in TARGETS:
        if not (OUT/path).is_file():
            raise FileNotFoundError("Existing figure to replace not found: "+path)

    untouched={p.name:sha(p) for p in OUT.iterdir()
               if p.is_file() and p.name not in TARGETS}
    old_images={p:sha(OUT/p) for p in TARGETS}
    peptide_csv=read_csv(OUT/"peptide_sample_cluster_order.csv")
    protein_csv=read_csv(OUT/"protein_sample_cluster_order.csv")
    samples,logdata,_=check_source()
    peptidez,_,_=scaled(logdata[:,:2])
    proteinz,_,_=scaled(logdata[:,2:3])

    # Physically DELETE former PNGs, then render replacements at same paths.
    for filename in TARGETS:
        (OUT/filename).unlink()
    new_peptide=reordered_cluster(
        peptidez,samples,PEPTIDES,
        "Two yeast peptide features: sample hierarchical clustering (average / Euclidean)",
        FILES[2],
    )
    new_protein=reordered_cluster(
        proteinz,samples,["P00924"],
        "Yeast protein P00924: single-feature sample clustering (Euclidean)",
        FILES[5],
    )

    # CSV stores leaf ranks as strings; compare normalized records, not Python types.
    norm=lambda rows:[{**r,"cluster_leaf_order":str(r["cluster_leaf_order"])} for r in rows]
    if norm(new_peptide)!=peptide_csv or norm(new_protein)!=protein_csv:
        raise AssertionError("Sample cluster order changed unexpectedly in layout-only redraw")
    for filename in TARGETS:
        path=OUT/filename
        if not path.is_file() or path.stat().st_size<20000:
            raise AssertionError("Missing/undersized redrawn PNG: "+str(path))
        if old_images[filename]==sha(path):
            raise AssertionError("Redrawn PNG has same digest as prior overlapping figure")

    for filename,dig in untouched.items():
        if sha(OUT/filename)!=dig:
            raise AssertionError("An unrelated prior result was unexpectedly modified: "+filename)
    if len(list(OUT.glob("*.png")))!=6:
        raise AssertionError("Must retain exactly original six visualizations")
    print("RENDER PASS: old hierarchical PNGs physically removed and regenerated.")
    print("AUDIT PASS: numeric matrices, PCs, correlations, clustering CSVs, other four PNGs unchanged.")
    print("UPDATED",", ".join(TARGETS))


if __name__=="__main__":
    main()
