#!/usr/bin/env python3
"""Remove NYD sampleinfo rows absent from columns 6+ of NYD_pg_matrix.tsv.
Library ID is worksheet D column (文库编号); preserve 18 metadata columns.
"""
import csv, hashlib, json, re, tempfile, zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]/"NYD_data"
PG=ROOT/"NYD_pg_matrix.tsv"
INFO=ROOT/"NYD_sampleinfo.xlsx"
REPORT=ROOT/"NYD_sampleinfo_filter_report.md"
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
ET.register_namespace("",NS[1:-1])
ET.register_namespace("r","http://schemas.openxmlformats.org/officeDocument/2006/relationships")
REF=re.compile(r"^([A-Z]+)([0-9]+)$")

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def source_ids():
    with PG.open(encoding="utf-8-sig",newline="") as f:
        h=next(csv.reader(f,delimiter="\t"))
    if h[:5]!=["Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description"]:
        raise ValueError("Unexpected NYD protein matrix annotations")
    ids=h[5:]
    if len(ids)!=16 or len(set(ids))!=16 or not all(ids):
        raise ValueError("Expected 16 unique nonempty sample headers")
    return ids

def strings(z):
    if "xl/sharedStrings.xml" not in z.namelist(): return []
    r=ET.fromstring(z.read("xl/sharedStrings.xml"))
    return ["".join(t.text or "" for t in si.iter(NS+"t"))
            for si in r.findall(NS+"si")]

def val(c,ss):
    v=c.find(NS+"v")
    if c.get("t")=="s": return ss[int(v.text)] if v is not None and v.text else ""
    if c.get("t")=="inlineStr": return "".join(t.text or "" for t in c.iter(NS+"t"))
    return v.text if v is not None and v.text is not None else ""

def parsed(z):
    wb=ET.fromstring(z.read("xl/workbook.xml"))
    ssheets=wb.find(NS+"sheets")
    if ssheets is None or len(ssheets)!=1: raise ValueError("Expected 1 worksheet")
    key="xl/worksheets/sheet1.xml"
    if key not in z.namelist(): raise ValueError("Unexpected worksheet file")
    root=ET.fromstring(z.read(key))
    data=root.find(NS+"sheetData")
    if data is None: raise ValueError("Missing sheetData")
    st=strings(z)
    records=[]
    for row in data.findall(NS+"row"):
        values={}
        rn=int(row.get("r"))
        for c in row.findall(NS+"c"):
            m=REF.fullmatch(c.get("r",""))
            if not m or int(m.group(2))!=rn or c.find(NS+"f") is not None:
                raise ValueError("Unexpected cell reference or formula")
            values[m.group(1)]=val(c,st)
        records.append((rn,values))
    return key,root,data,st,records

def main():
    ids=source_ids()
    allowed=set(ids)
    before=sha(INFO)
    with zipfile.ZipFile(INFO) as z:
        if z.testzip(): raise ValueError("Corrupted workbook ZIP")
        key,root,data,ss,records=parsed(z)
    if len(records)!=41 or records[0][1].get("D")!="文库编号" or len(records[0][1])!=18:
        raise ValueError("Unexpected NYD sampleinfo layout")
    raw_ids=[r.get("D","").strip() for _,r in records[1:]]
    if not all(raw_ids) or len(set(raw_ids))!=40:
        raise ValueError("Blank or duplicate library IDs")
    kept=[x for x in raw_ids if x in allowed]
    removed=[x for x in raw_ids if x not in allowed]
    if len(kept)!=16 or set(kept)!=allowed:
        raise ValueError("Metadata sample set does not cover exactly 16 matrix IDs")
    risky=["mergeCells","autoFilter","tableParts","conditionalFormatting",
           "dataValidations","hyperlinks","drawing","legacyDrawing","rowBreaks"]
    if any(root.find(NS+x) is not None for x in risky):
        raise ValueError("Workbook contains row-dependent objects requiring separate handling: "+repr([x for x in risky if root.find(NS+x) is not None]))
    expected=[r for _,r in records if r.get("D") in allowed or r.get("D")=="文库编号"]
    nextrow=0
    for row in list(data.findall(NS+"row")):
        d=[c for c in row.findall(NS+"c") if REF.fullmatch(c.get("r","")).group(1)=="D"]
        identifier=val(d[0],ss).strip() if d else ""
        if row.get("r")!="1" and identifier not in allowed:
            data.remove(row)
            continue
        nextrow+=1
        row.set("r",str(nextrow))
        for c in row.findall(NS+"c"):
            m=REF.fullmatch(c.get("r",""))
            c.set("r",m.group(1)+str(nextrow))
    if nextrow!=17: raise AssertionError("Unexpected filtered row count")
    dim=root.find(NS+"dimension")
    if dim is not None: dim.set("ref","A1:R17")
    xml=ET.tostring(root,encoding="utf-8",xml_declaration=True)
    with tempfile.NamedTemporaryFile(dir=ROOT,prefix="clean_nyd_",suffix=".xlsx",delete=False) as f:
        tmp=Path(f.name)
    try:
        with zipfile.ZipFile(INFO) as zin,zipfile.ZipFile(tmp,"w") as out:
            for part in zin.infolist():
                out.writestr(part,xml if part.filename==key else zin.read(part.filename))
        with zipfile.ZipFile(tmp) as z:
            if z.testzip(): raise AssertionError("Corrupted output workbook")
            _,_,_,_,new=parsed(z)
            after_rows=[r for _,r in new]
            if after_rows!=expected or len(after_rows)!=17:
                raise AssertionError("Preserved metadata values/order differ from original")
            if any(len(r)!=18 for r in after_rows):
                raise AssertionError("Column widths lost")
            if [r["D"] for r in after_rows[1:]]!=kept:
                raise AssertionError("Incorrect filtered library IDs")
        result_sha=sha(tmp)
        tmp.replace(INFO)
        lines=[
            "# NYD sampleinfo metadata filtering","",
            "Based on NYD_pg_matrix.tsv sample columns 6 onward.",
            "Matching column in sampleinfo: D (文库编号), not E (样本编号).","",
            "- Matrix sample IDs: 16",
            "- Original metadata data rows: 40",
            "- Retained metadata data rows: 16",
            "- Removed metadata data rows: 24",
            "- The 18 metadata columns and retained cell data/order were preserved.",
            "- NYD_pg_matrix.tsv and NYD_pr_matrix.tsv were not changed.","",
            "## Retained library identifiers","",
        ]
        lines += ["- "+x for x in kept]
        lines += ["","## Removed library identifiers",""]
        lines += ["- "+x for x in removed]
        lines += ["","## Input/output SHA-256","",
                  "- Original: "+before,"- Filtered: "+result_sha,""]
        REPORT.write_text("\n".join(lines),encoding="utf-8")
        print("PASS "+json.dumps({"before":40,"retained":16,"removed":24,
              "removed_sample_ids":removed,"sha_before":before,"sha_after":result_sha},
              ensure_ascii=False),flush=True)
    finally:
        tmp.unlink(missing_ok=True)

if __name__=="__main__":
    main()
