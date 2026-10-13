#!/usr/bin/env python3
"""Metadata-only filtration of reuploaded 2026-10-13 HELP and yeast XLSX.

Retain 75 metadata-listed sample columns (27 HC, 27 S, 6 LM Internal-QC,
15 pooled QC). Keep all L, H and ratio_LH peptide rows and both yeast peptides.
Never recompute, log-transform, impute or normalize any intensity or ratio.
Rewriting uses ZIP+XML so numeric cell text and retained formatting are preserved.
"""
import hashlib
import json
import re
import tempfile
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

from yeast_protein_qc import metadata, GROUPS

ROOT = Path(__file__).resolve().parents[1]
HELP = ROOT / "LH_ratio_tables_2026.10.13—ALL-HELP.xlsx"
YEAST = ROOT / "sample-QC-InternalQC-HELP_yeast.xlsx"
REPORT = ROOT / "HELP_YEAST_FILTER_REPORT.md"
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
ET.register_namespace("", NS[1:-1])
ET.register_namespace("r", "http://schemas.openxmlformats.org/officeDocument/2006/relationships")
ADDRESS = re.compile(r"^([A-Z]+)([0-9]+)$")


def column_number(s):
    i = 0
    for c in s:
        i = i * 26 + ord(c) - 64
    return i


def column_name(n):
    s = ""
    while n:
        n, k = divmod(n - 1, 26)
        s = chr(65 + k) + s
    return s


def hash_file(path):
    hasher = hashlib.sha256()
    with path.open("rb") as f:
        for piece in iter(lambda: f.read(1024*1024), b""):
            hasher.update(piece)
    return hasher.hexdigest()


def shared_strings(z):
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    root = ET.fromstring(z.read("xl/sharedStrings.xml"))
    return ["".join(t.text or "" for t in si.iter(NS+"t"))
            for si in root.findall(NS+"si")]


def cell_text(cell, strings):
    kind = cell.get("t")
    v = cell.find(NS+"v")
    if kind == "s":
        return strings[int(v.text)] if v is not None and v.text is not None else ""
    if kind == "inlineStr":
        return "".join(x.text or "" for x in cell.iter(NS+"t"))
    return v.text if v is not None and v.text is not None else ""


def sheet_values(root, strings):
    data = root.find(NS+"sheetData")
    if data is None:
        raise ValueError("Worksheet missing sheetData")
    lines = []
    for row in data.findall(NS+"row"):
        fields = {}
        rn = row.get("r")
        for cell in row.findall(NS+"c"):
            address = cell.get("r", "")
            m = ADDRESS.fullmatch(address)
            if not m or m.group(2) != rn:
                raise ValueError(f"Invalid XLSX cell coordinate {address}")
            if cell.find(NS+"f") is not None:
                raise ValueError("Spreadsheet formula found; cannot safely shift cell references")
            col = m.group(1)
            if col in fields:
                raise ValueError(f"Duplicate cell in {address}")
            fields[col] = cell_text(cell, strings)
        lines.append(fields)
    return lines


def workbook(path):
    with zipfile.ZipFile(path) as z:
        bad = z.testzip()
        if bad:
            raise ValueError(f"Corrupt XLSX archive: {bad}")
        x = ET.fromstring(z.read("xl/workbook.xml"))
        sheets = [e.get("name") for e in x.find(NS+"sheets")]
        names = [f"xl/worksheets/sheet{i}.xml" for i in range(1, len(sheets)+1)]
        if any(name not in z.namelist() for name in names):
            raise ValueError("Unexpected XLSX worksheet file paths")
        strings = shared_strings(z)
        roots = [ET.fromstring(z.read(name)) for name in names]
        records = [sheet_values(root, strings) for root in roots]
    return sheets, names, roots, records


def header(records):
    cols = sorted(((column_number(col), str(sid).strip())
                   for col, sid in records[0].items() if col != "A"),
                  key=lambda x: x[0])
    ids = [sid for _, sid in cols]
    if not ids or any(not sid for sid in ids) or len(ids) != len(set(ids)):
        raise ValueError("XLSX header has empty or duplicate sample IDs")
    return cols


def snapshot(rows, cols, accepted):
    return (
        [row.get("A", "") for row in rows[1:]],
        {sid: [r.get(column_name(col), "") for r in rows[1:]]
         for col, sid in cols if sid in accepted},
    )


def filter_plan(path, labels, sheets_expected, n_peptides):
    names, zip_paths, roots, records = workbook(path)
    if names != sheets_expected:
        raise ValueError(f"{path.name}: expected sheets {sheets_expected}, saw {names}")
    if any(len(rows) != n_peptides+1 for rows in records):
        raise ValueError(f"{path.name}: incorrect peptide row count {[len(x)-1 for x in records]}")
    columns = [header(x) for x in records]
    sid_lists = [[sid for _, sid in cols] for cols in columns]
    if sid_lists[1:] != sid_lists[:-1]:
        raise ValueError(f"{path.name}: sample order differs between sheets")
    src_ids = sid_lists[0]
    if not set(labels).issubset(src_ids):
        raise ValueError(f"{path.name}: metadata samples missing: {sorted(set(labels)-set(src_ids))}")
    peptide_lists = [[r.get("A","") for r in records_i[1:]] for records_i in records]
    if peptide_lists[1:] != peptide_lists[:-1]:
        raise ValueError(f"{path.name}: peptide rows differ between sheets")
    peptides = peptide_lists[0]
    if len(set(peptides)) != n_peptides or any(not x for x in peptides):
        raise ValueError(f"{path.name}: missing/duplicate peptide identifiers")
    if path == YEAST and set(peptides) != {"AADALLLK", "VNQIGTLSESIK"}:
        raise ValueError(f"{path.name}: unexpected yeast peptides {peptides}")
    removed = [sid for sid in src_ids if sid not in labels]
    kept = [sid for sid in src_ids if sid in labels]
    if len(kept) != 75 or len(set(kept)) != 75:
        raise ValueError("Expected 75 distinct retained sample IDs")
    original = [snapshot(rows, cols, labels) for rows, cols in zip(records, columns)]
    xml_updates = {}
    if removed:
        for key, root, cols in zip(zip_paths, roots, columns):
            mapping = {1: 1}
            current = 2
            for pos, sid in cols:
                if sid in labels:
                    mapping[pos] = current
                    current += 1
            if current != 77:
                raise AssertionError("Mapping does not have exactly 75 sample columns")
            for row in root.find(NS+"sheetData").findall(NS+"row"):
                rn = row.get("r")
                row.set("spans", "1:76")
                for cell in list(row.findall(NS+"c")):
                    match = ADDRESS.fullmatch(cell.get("r", ""))
                    if match is None:
                        raise ValueError("Unexpected cell address")
                    old_col = column_number(match.group(1))
                    if old_col in mapping:
                        cell.set("r", column_name(mapping[old_col]) + rn)
                    else:
                        row.remove(cell)
            dimension = root.find(NS+"dimension")
            if dimension is not None:
                dimension.set("ref", f"A1:BX{n_peptides+1}")
            xml_updates[key] = ET.tostring(root, encoding="utf-8", xml_declaration=True)
    return {
        "path": path, "sheet_names": names, "peptides": peptides,
        "old_ids": src_ids, "kept_ids": kept, "removed": removed,
        "original": original, "updates": xml_updates,
        "sha_before": hash_file(path), "bytes_before": path.stat().st_size,
    }


def create_and_check(plan, labels):
    if not plan["updates"]:
        return None, plan["sha_before"], plan["bytes_before"]
    with tempfile.NamedTemporaryFile(dir=ROOT, prefix="filtered_help_",
                                     suffix=".xlsx", delete=False) as fh:
        tmp = Path(fh.name)
    try:
        with zipfile.ZipFile(plan["path"]) as source, zipfile.ZipFile(tmp, "w") as target:
            for info in source.infolist():
                target.writestr(info, plan["updates"].get(info.filename,
                                                        source.read(info.filename)))
        sheets, _, roots, records = workbook(tmp)
        if sheets != plan["sheet_names"]:
            raise AssertionError("Sheet names changed")
        for i, rows in enumerate(records):
            cols = header(rows)
            if [sid for _, sid in cols] != plan["kept_ids"]:
                raise AssertionError(f"Sample order/number changed on sheet {i+1}")
            if snapshot(rows, cols, labels) != plan["original"][i]:
                raise AssertionError(f"Peptide/sample values changed on sheet {i+1}")
            if len(rows) != len(plan["peptides"])+1:
                raise AssertionError("Peptide row count changed")
        return tmp, hash_file(tmp), tmp.stat().st_size
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def main():
    labels = metadata()
    counts = {g: sum(x == g for x in labels.values()) for g in GROUPS}
    if counts != {"Internal-QC":6, "pooled QC":15, "HC":27, "S":27}:
        raise ValueError(f"Unexpected sample metadata: {counts}")
    with zipfile.ZipFile(YEAST) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        yeast_sheets = [x.get("name") for x in wb.find(NS+"sheets")]
    if len(yeast_sheets) != 1:
        raise ValueError("Yeast workbook should contain exactly one worksheet")

    plans = [
        filter_plan(HELP, labels, ["L","H","ratio_LH"], 97),
        filter_plan(YEAST, labels, yeast_sheets, 2),
    ]
    prepared = []
    try:
        for p in plans:
            prepared.append((p, *create_and_check(p, labels)))
        # Both have passed full validation; now replace the original uploads.
        for p, temp, digest, size in prepared:
            if temp is not None:
                temp.replace(p["path"])
        report = [
            "# Metadata filtering of 2026-10-13 HELP and yeast workbooks", "",
            "Using all_sample_metadata.xlsx and all_QC_metadata.xlsx, retain 75 samples:",
            "27 HC + 27 S + 6 Internal-QC (LM) + 15 pooled QC (QC1).",
            "Blank/unlisted sample columns are removed.",
            "",
            "Only sample columns are removed; no log2 transformation, normalization,",
            "imputation, ratio recomputation or peptide filtering is performed.",
            "Every retained sample/peptide value was compared to the original before",
            "modification; peptide sequences and three HELP sheets are preserved.",
            "Two yeast peptides, AADALLLK and VNQIGTLSESIK, remain independent.",
            "Original versions remain recoverable from the repository's Git history.",
            "",
        ]
        results = {}
        for p, temp, digest, size in prepared:
            results[p["path"].name] = {
                "before_samples": len(p["old_ids"]),
                "after_samples": 75, "removed_ids": p["removed"],
                "peptides": len(p["peptides"]), "sheets": p["sheet_names"],
                "changed": temp is not None,
                "sha256_before": p["sha_before"], "sha256_after": digest,
            }
            report += [
                "## " + p["path"].name, "",
                "- Worksheets: " + ", ".join(p["sheet_names"]),
                "- Peptide rows: " + str(len(p["peptides"])) + " (unchanged)",
                "- Sample columns: " + str(len(p["old_ids"])) + " -> 75",
                "- Removed columns: " + str(len(p["removed"])),
                "- Status: " + ("filtered" if temp is not None else "already filtered"),
                "- Previous SHA256: " + p["sha_before"],
                "- Current SHA256: " + digest,
                "- Previous size: " + str(p["bytes_before"]) + " bytes",
                "- Current size: " + str(size) + " bytes", "",
                "### Removed sample IDs", "",
            ]
            report.extend("- " + x for x in p["removed"])
            if not p["removed"]:
                report.append("- None")
            report += [""]
        report += [
            "## Important note", "",
            "Older PACS trial outputs were NOT recomputed with these updated files.",
            "Previous scripts still pointing to the deleted 2026.9.22 HELP filename",
            "must be updated before the next normalization analysis.", "",
        ]
        REPORT.write_text("\n".join(report), encoding="utf-8")
        print("PASS " + json.dumps({"groups": counts, "files": results},
                                   ensure_ascii=False))
    finally:
        for p, temp, digest, size in prepared:
            if temp is not None:
                temp.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
