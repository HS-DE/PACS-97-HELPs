#!/usr/bin/env python3
"""Inspect NYD yeast source spreadsheet layouts without modifying workbook files."""
import zipfile, re
from pathlib import Path
from xml.etree import ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]/"NYD_data"/"40samples"
N="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RNS="{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
PNS="{http://schemas.openxmlformats.org/package/2006/relationships}"
def examine(path):
    with zipfile.ZipFile(path) as z:
        ss=[]
        if "xl/sharedStrings.xml" in z.namelist():
            ss=["".join(x.text or "" for x in si.iter(N+"t")) for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(N+"si")]
        book=ET.fromstring(z.read("xl/workbook.xml"))
        sheets=book.find(N+"sheets").findall(N+"sheet")
        rels={x.get("Id"):x.get("Target") for x in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
        print("\nWORKBOOK",path.name,"SIZE",path.stat().st_size,"SHEETS",[(x.get("name"),x.get(RNS+"id")) for x in sheets],flush=True)
        for sh in sheets:
            target=rels[sh.get(RNS+"id")]
            sp=target.lstrip("/") if target.startswith("/") else "xl/"+target.lstrip("/")
            root=ET.fromstring(z.read(sp))
            data=root.find(N+"sheetData")
            rows=[]
            for xmlrow in data.findall(N+"row"):
                cells=[]
                for c in xmlrow.findall(N+"c"):
                    v=c.find(N+"v")
                    if v is not None and v.text is not None:
                        value=ss[int(v.text)] if c.get("t")=="s" else v.text
                    elif c.get("t")=="inlineStr":
                        value="".join(x.text or "" for x in c.iter(N+"t"))
                    else:value=""
                    if value!="" and value is not None: cells.append((c.get("r"),value))
                rows.append((xmlrow.get("r"),cells))
            print("SHEET",sh.get("name"),"ACTUAL ROWS",len(rows),"MAX_ROW",max(int(r[0]) for r in rows),"FIRST NONEMPTY CONTENT",flush=True)
            for r,cells in rows:
                if not cells:continue
                print("ROW",r,"COUNT",len(cells),"FIRST 9",cells[:9],"LAST 3",cells[-3:],flush=True)
                if len(rows)>150 and int(r)>110:
                    print("... (remaining rows truncated)");break
if __name__=="__main__":
    for f in ("NYD_sampleinfo2.xlsx","YEAST.xlsx"):examine(ROOT/f)
