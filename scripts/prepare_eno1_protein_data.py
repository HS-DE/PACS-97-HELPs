#!/usr/bin/env python3
"""Data-only cleaning: separate ENO1_YEAST and ENO2_YEAST from pg_matrix.tsv,
and delete metadata-unlisted sample columns from uploaded EN01.xlsx.

No analysis, normalization, conversion or reordering of remaining data.
Retains both yeast-enolase rows in the standalone Excel workbook.
Writes a concise provenance report; all inputs recoverable in Git history.
"""
import hashlib
import json
import os
import re
import tempfile
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

from yeast_protein_qc import metadata, GROUPS

ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "pg_matrix.tsv"
YEAST_PROTEIN = ROOT / "EN01.xlsx"
REPORT = ROOT / "ENO1_DATA_CLEANUP_REPORT.md"
ANNOTATIONS = ["Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
               "First.Protein.Description"]
EXPECTED_GROUPS = {"Internal-QC": 6, "pooled QC": 15, "HC": 27, "S": 27}
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
ADDR = re.compile(r"^([A-Z]+)([0-9]+)$")
ET.register_namespace("", NS[1:-1])
ET.register_namespace("r", "http://schemas.openxmlformats.org/officeDocument/2006/relationships")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda: f.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def to_number(label):
    n = 0
    for letter in label:
        n = 26 * n + ord(letter) - 64
    return n


def to_letters(n):
    text = ""
    while n:
        n, k = divmod(n - 1, 26)
        text = chr(65 + k) + text
    return text


def shared_strings(archive):
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return ["".join(node.text or "" for node in si.iter(NS + "t"))
            for si in root.findall(NS + "si")]


def read_value(cell, shared):
    val = cell.find(NS + "v")
    if cell.get("t") == "s":
        return shared[int(val.text)] if val is not None and val.text else ""
    if cell.get("t") == "inlineStr":
        return "".join(n.text or "" for n in cell.iter(NS + "t"))
    return val.text if val is not None and val.text is not None else ""


def sheet_rows(sheet, shared):
    data = sheet.find(NS + "sheetData")
    if data is None:
        raise ValueError("EN01 worksheet missing sheetData")
    result = []
    for row in data.findall(NS + "row"):
        vals = {}
        for cell in row.findall(NS + "c"):
            ref = cell.get("r", "")
            match = ADDR.fullmatch(ref)
            if match is None:
                raise ValueError(f"Unexpected EN01 cell ref: {ref}")
            if cell.find(NS + "f") is not None:
                raise ValueError("EN01 includes formulas; refusing unsafe column removal")
            col, rn = match.groups()
            if int(rn) != int(row.get("r")):
                raise ValueError("EN01 sheet inconsistent cell row refs")
            n = to_number(col)
            if n in vals:
                raise ValueError(f"Repeated cell at {ref}")
            vals[n] = read_value(cell, shared)
        result.append(vals)
    return result


def sheet_path(archive):
    wb = ET.fromstring(archive.read("xl/workbook.xml"))
    sheets = wb.find(NS + "sheets")
    if sheets is None or len(sheets) != 1:
        raise ValueError("Expected EN01.xlsx with exactly one worksheet")
    sheet_name = sheets[0].get("name")
    # The uploaded workbook uses standard worksheet entry; validate before editing.
    path = "xl/worksheets/sheet1.xml"
    if path not in archive.namelist():
        raise ValueError("Unexpected EN01.xlsx worksheet path")
    return sheet_name, path


def plan_matrix(source, accepted):
    if not source.is_file():
        raise FileNotFoundError(source)
    original_sha = digest(source)
    count = 0
    removed_lines = []
    old_header = None
    fd, name = tempfile.mkstemp(prefix="pacs_eno1_pg_", suffix=".tsv", dir=ROOT)
    os.close(fd)
    temp = Path(name)
    try:
        with source.open("r", encoding="utf-8-sig", newline="") as inp, \
                temp.open("w", encoding="utf-8", newline="") as out:
            hdr = inp.readline()
            if not hdr:
                raise ValueError("Empty pg_matrix.tsv")
            old_header = hdr.rstrip("\r\n").split("\t")
            if old_header[:5] != ANNOTATIONS:
                raise ValueError("Unexpected protein annotation columns: " + str(old_header[:5]))
            if len(old_header) != 80 or set(old_header[5:]) != set(accepted):
                raise ValueError("Protein matrix not matched to expected 75 metadata samples")
            name_idx = old_header.index("Protein.Names")
            out.write(hdr.lstrip("\ufeff"))
            for line_number, line in enumerate(inp, 2):
                cells = line.rstrip("\r\n").split("\t")
                if len(cells) != len(old_header):
                    raise ValueError(f"Invalid TSV width at line {line_number}")
                if cells[name_idx] in {"ENO1_YEAST", "ENO2_YEAST"}:
                    removed_lines.append((line_number, cells[0], cells[name_idx]))
                else:
                    out.write(line)
                count += 1
        if len(removed_lines) != 2 or {r[2] for r in removed_lines} != {"ENO1_YEAST", "ENO2_YEAST"}:
            raise ValueError(f"Expected one ENO1_YEAST and one ENO2_YEAST row: {removed_lines}")
        if count != 6677:
            raise ValueError(f"Unexpected source protein count {count}; expected 6677")
        with temp.open("r", encoding="utf-8", newline="") as f:
            retained_header = f.readline().rstrip("\r\n").split("\t")
            verified = 0
            for line in f:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != len(old_header) or fields[name_idx] in {"ENO1_YEAST", "ENO2_YEAST"}:
                    raise AssertionError("TSV cleaning output not valid")
                verified += 1
        if retained_header != old_header or verified != count - 2:
            raise AssertionError("Filtered TSV headers/row count changed unexpectedly")
        return {
            "temp": temp,
            "source_sha": original_sha,
            "filtered_sha": digest(temp),
            "old_rows": count,
            "new_rows": verified,
            "removed": removed_lines,
        }
    except BaseException:
        temp.unlink(missing_ok=True)
        raise


def plan_excel(source, accepted):
    if not source.is_file():
        raise FileNotFoundError(source)
    before = digest(source)
    with zipfile.ZipFile(source) as z:
        if z.testzip() is not None:
            raise ValueError("Damaged EN01.xlsx archive")
        title, key = sheet_path(z)
        strings = shared_strings(z)
        sheet = ET.fromstring(z.read(key))
        original = sheet_rows(sheet, strings)
        if not original:
            raise ValueError("Empty EN01 worksheet")
        headings = original[0]
        # EN01.xlsx is a small standalone table: sample IDs start in column B,
        # column A stores each yeast-enolase protein name; A1 may be blank.
        extra_columns = sorted(k for k in headings if k >= 2)
        extra_names = [str(headings[k]).strip() for k in extra_columns]
        if len(extra_names) != len(set(extra_names)) or not all(extra_names):
            raise ValueError("Duplicate or missing EN01 sample names")
        present = set(extra_names)
        missing = sorted(set(accepted) - present)
        if missing:
            raise ValueError(f"Missing metadata samples from EN01: {missing}")
        names_in_rows = [str(row.get(1, "")).strip() for row in original[1:]]
        if len(original) != 3 or set(names_in_rows) != {"ENO1_YEAST", "ENO2_YEAST"}:
            raise ValueError(
                f"EN01 should contain the two yeast-enolase rows; "
                f"observed {len(original)-1} rows with first-column IDs {names_in_rows}"
            )
        included = [n for n in extra_columns if headings[n].strip() in accepted]
        # Keep first protein identifier column and all 75 study/QC sample columns.
        keep_numbers = [1] + included
        if len(keep_numbers) != 76 or len(set(keep_numbers)) != 76:
            raise ValueError("Expected EN01 one protein-name column + 75 sample columns")
        map_index = {old: j for j, old in enumerate(keep_numbers, 1)}
        kept_names = [headings[k] for k in included]
        removed_names = [headings[k] for k in extra_columns if k not in included]
        # Compare the raw XLSX value/text strings for every retained cell.
        original_projection = [
            [r.get(index, "") for index in keep_numbers] for r in original
        ]
        if removed_names:
            if sheet.find(NS + "mergeCells") is not None:
                raise ValueError("EN01 contains merged cells; safe remapping not implemented")
            if sheet.find(NS + "autoFilter") is not None:
                raise ValueError("EN01 contains autoFilter; safe remapping not implemented")
            if sheet.find(NS + "tableParts") is not None:
                raise ValueError("EN01 contains tableParts; safe remapping not implemented")
            for row in sheet.find(NS + "sheetData").findall(NS + "row"):
                rn = row.get("r")
                row.set("spans", "1:76")
                for cell in list(row.findall(NS + "c")):
                    old = to_number(ADDR.fullmatch(cell.get("r")).group(1))
                    if old not in map_index:
                        row.remove(cell)
                    else:
                        cell.set("r", to_letters(map_index[old]) + rn)
            dimension = sheet.find(NS + "dimension")
            if dimension is not None:
                dimension.set("ref", "A1:BX3")  # 76 columns = BX
        # Validate changed workbook, including every preserved cell value and sequence.
        fd, name = tempfile.mkstemp(prefix="pacs_eno1_excel_", suffix=".xlsx", dir=ROOT)
        os.close(fd)
        temp = Path(name)
        if removed_names:
            new_sheet_data = ET.tostring(sheet, encoding="utf-8", xml_declaration=True)
            with zipfile.ZipFile(source) as inp, zipfile.ZipFile(temp, "w") as out:
                for info in inp.infolist():
                    out.writestr(info, new_sheet_data if info.filename == key else inp.read(info.filename))
        else:
            # Preserve the original XLSX byte-for-byte if already filtered.
            temp.write_bytes(source.read_bytes())
    try:
        with zipfile.ZipFile(temp) as zipped:
            if zipped.testzip() is not None:
                raise AssertionError("Output XLSX archive damaged")
            _, newkey = sheet_path(zipped)
            filtered = sheet_rows(ET.fromstring(zipped.read(newkey)), shared_strings(zipped))
            if len(filtered) != 3 or [r.get(1, "") for r in filtered[1:]] != names_in_rows:
                raise AssertionError("Filtered EN01 protein names or rows changed")
            header_ids = [filtered[0].get(i, "") for i in range(2, 77)]
            if header_ids != kept_names or set(header_ids) != set(accepted):
                raise AssertionError("Filtered EN01 sample names/order changed")
            if [[r.get(k, "") for k in range(1, 77)] for r in filtered] != original_projection:
                raise AssertionError("EN01 retained values changed during screening")
        return {
            "temp": temp,
            "old_samples": len(extra_names),
            "new_samples": 75,
            "removed_samples": removed_names,
            "sha_before": before,
            "sha_after": digest(temp),
            "sheet_name": title,
        }
    except BaseException:
        temp.unlink(missing_ok=True)
        raise


def main():
    labels = metadata()
    counts = {g: sum(t == g for t in labels.values()) for g in GROUPS}
    if counts != EXPECTED_GROUPS or len(labels) != 75:
        raise ValueError(f"Metadata group counts unexpected: {counts}")
    m = plan_matrix(MATRIX, labels)
    try:
        x = plan_excel(YEAST_PROTEIN, labels)
    except BaseException:
        m["temp"].unlink(missing_ok=True)
        raise
    try:
        # Both input sources have been independently validated first.
        m["temp"].replace(MATRIX)
        x["temp"].replace(YEAST_PROTEIN)
        report = [
            "# ENO1_YEAST separation and EN01.xlsx metadata filtering", "",
            "Data preparation only; no correction, statistical analysis,",
            "normalization, protein intensity conversion or new trial results.", "",
            "## Metadata selection", "",
            "- Study: HC 27, S 27",
            "- Internal-QC (LM): 6",
            "- pooled QC (QC1): 15",
            "- Total approved sample IDs: 75; Blank and other IDs removed", "",
            "## pg_matrix.tsv", "",
            f"- Protein rows: {m['old_rows']} → {m['new_rows']}.",
            "- Removed exactly two yeast-enolase rows, ENO1_YEAST and ENO2_YEAST.",
            "- All other protein rows, sample columns, and original TSV values preserved.",
            "- Removed Protein.Group identifiers and original line numbers:",
        ]
        report += [f"  - line {line}, Protein.Group: {group}, Protein.Names: {name}" for line, group, name in m["removed"]]
        report += [
            f"- Input SHA-256: {m['source_sha']}",
            f"- Output SHA-256: {m['filtered_sha']}", "",
            "## EN01.xlsx", "",
            f"- Sheet: {x['sheet_name']}",
            "- Protein rows: 2; first-column IDs ENO1_YEAST and ENO2_YEAST; retained.",
            f"- Sample columns: {x['old_samples']} → 75.",
            f"- Removed sample columns: {len(x['removed_samples'])}.",
            "- Preserved the protein identifier column and numeric value strings",
            "  of all 75 retained samples for both yeast-enolase protein rows.",
            f"- Input SHA-256: {x['sha_before']}",
            f"- Output SHA-256: {x['sha_after']}", "",
            "### Removed sample IDs", "",
        ]
        report += ["- " + str(sid) for sid in x["removed_samples"]] or ["- None"]
        report += [
            "", "## Other inputs and existing results", "",
            "- sample-QC-InternalQC-HELP_yeast.xlsx was not modified.",
            "- Other metadata and HELP files were not modified.",
            "- results/, results_v1/ and results_v2/ were not modified.",
            "- Prior result summaries describe the earlier protein matrix.",
            "- Original source versions remain available in Git history.", "",
        ]
        REPORT.write_text("\n".join(report), encoding="utf-8")
        print("PASS " + json.dumps({
            "pg_original_rows": m["old_rows"], "pg_filtered_rows": m["new_rows"],
            "eno1_rows_removed_from_pg": len(m["removed"]),
            "excel_before_samples": x["old_samples"], "excel_retained_samples": x["new_samples"],
            "excel_removed_sample_ids": x["removed_samples"],
            "groups": counts, "report": REPORT.name,
        }, ensure_ascii=False))
    finally:
        m["temp"].unlink(missing_ok=True)
        x["temp"].unlink(missing_ok=True)


if __name__ == "__main__":
    main()
