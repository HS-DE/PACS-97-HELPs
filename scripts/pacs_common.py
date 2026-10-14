#!/usr/bin/env python3
"""Common HELP reader, protein metadata and CV primitives for PACS pipelines.

No normalization or experimental correction is performed in this module.
"""
import math
import zipfile
from xml.etree import ElementTree as ET

import numpy as np

from filter_samples import readbook, tab
from yeast_protein_qc import COLORS, GROUPS, ROOT

PG = ROOT / "pg_matrix.tsv"
HELP = ROOT / "LH_ratio_tables_2026.9.22—ALL-HELP.xlsx"
ANN = ["Protein.Group", "Protein.Ids", "Protein.Names",
       "Genes", "First.Protein.Description"]
CV_MIN_DETECT_FRAC = 0.70
CV_MIN_N = 3
CV_MIN_MEAN = 1e-8


def load_ratio(samples):
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    with zipfile.ZipFile(HELP) as z:
        book = ET.fromstring(z.read("xl/workbook.xml"))
        names = [node.get("name") for node in book.find(ns + "sheets")]
    if names != ["L", "H", "ratio_LH"]:
        raise ValueError("Expected HELP sheets L, H, ratio_LH; saw " + str(names))
    shared, sheets = readbook(HELP)
    if len(sheets) != 3:
        raise ValueError("HELP sheet count differs from three")
    rows = tab(sheets[2], shared)
    if len(rows) != 98:
        raise ValueError("Expected 97 peptide rows; got " + str(len(rows)-1))
    header = rows[0]
    ids = {name: col for col, name in header.items() if col != "A"}
    if len(ids) != 75 or set(ids) != set(samples):
        raise ValueError("HELP ratio columns differ from metadata IDs")
    seq = [r.get("A", "") for r in rows[1:]]
    if len(set(seq)) != 97 or any(not x for x in seq):
        raise ValueError("HELP peptide sequences not unique or missing")
    medians, n_valid = {}, {}
    for sample, col in ids.items():
        values = []
        for row in rows[1:]:
            try:
                x = float(row.get(col, ""))
                if math.isfinite(x) and x > 0:
                    values.append(x)
            except (TypeError, ValueError):
                pass
        if not values:
            raise ValueError("No valid HELP ratios for " + sample)
        medians[sample] = float(np.median(values))
        n_valid[sample] = len(values)
    return medians, n_valid


def finite_positive(a):
    return np.isfinite(a) & (a > 0)


def cv_by_feature(matrix):
    """CV for each protein within one sample group, as in supplied R defaults."""
    good = finite_positive(matrix)
    n = good.sum(axis=1)
    valid = np.where(good, matrix, 0.0)
    total = valid.sum(axis=1)
    avg = np.divide(total, n, out=np.full(len(n), np.nan), where=n > 0)
    squares = (valid * valid).sum(axis=1)
    sample_var = np.divide(
        np.maximum(squares - np.divide(total * total, n,
                                      out=np.zeros_like(total), where=n > 0), 0.0),
        n - 1,
        out=np.full(len(n), np.nan), where=n > 1)
    cv = np.divide(100 * np.sqrt(sample_var), avg,
                   out=np.full(len(n), np.nan), where=avg >= CV_MIN_MEAN)
    needed = max(CV_MIN_N, math.ceil(matrix.shape[1] * CV_MIN_DETECT_FRAC))
    cv[(n < needed) | (~np.isfinite(cv))] = np.nan
    return cv, n, needed
