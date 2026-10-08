#!/usr/bin/env python3
"""Metadata-only cleanup: keep all study + non-Blank QC sample_id columns."""
import re, zipfile, tempfile, json
from pathlib import Path
from xml.etree import ElementTree as E
ROOT=Path(__file__).resolve().parents[1]
N="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
T=lambda x:"{"+N+"}"+x
E.register_namespace("",N)
E.register_namespace("r","http://schemas.openxmlformats.org/officeDocument/2006/relationships")
def ci(s):
 n=0
 for c in s:n=n*26+ord(c)-64
 return n
def cn(n):
 s=""
 while n:n,v=divmod(n-1,26);s=chr(65+v)+s
 return s
def readbook(path):
 with zipfile.ZipFile(path) as z:
  strings=[]
  if "xl/sharedStrings.xml" in z.namelist():
   strings=["".join(t.text or "" for t in si.iter(T("t"))) for si in E.fromstring(z.read("xl/sharedStrings.xml")).findall(T("si"))]
  sheets=[E.fromstring(z.read("xl/worksheets/sheet"+str(i)+".xml")) for i in range(1,4) if "xl/worksheets/sheet"+str(i)+".xml" in z.namelist()]
 return strings,sheets
def tab(root,ss):
 records=[]
 for row in root.find(T("sheetData")).findall(T("row")):
  r={}
  for c in row.findall(T("c")):
   col=re.match(r"([A-Z]+)[0-9]+",c.get("r")).group(1);v=c.find(T("v"))
   r[col]=(ss[int(v.text)] if c.get("t")=="s" and v is not None else v.text if v is not None and v.text is not None else "")
  records.append(r)
 return records
def records(path):
 ss,sheets=readbook(path);rows=tab(sheets[0],ss);hdr={v:k for k,v in rows[0].items()}
 return [{name:row.get(c,"") for name,c in hdr.items()} for row in rows[1:]]
def rewrite(path,changes):
 with zipfile.ZipFile(path) as z:files=[(info,z.read(info.filename)) for info in z.infolist()]
 fd=tempfile.NamedTemporaryFile(dir=ROOT,suffix=".xlsx",delete=False);dst=Path(fd.name);fd.close()
 with zipfile.ZipFile(dst,"w") as z:
  for info,data in files:
   if info.filename in changes:data=E.tostring(changes[info.filename],encoding="utf-8",xml_declaration=True)
   z.writestr(info,data)
 with zipfile.ZipFile(dst) as z:assert z.testzip() is None
 dst.replace(path)
def main():
 study=records(ROOT/"all_sample_metadata.xlsx")
 qc=records(ROOT/"all_QC_metadata.xlsx")
 retainedqc=[q for q in qc if q["group1"].strip().lower()!="blank" and q["group2"].strip().lower()!="blank"]
 accepted=set([s["sample_id"] for s in study]+[s["sample_id"] for s in retainedqc])
 assert len(study)==54 and len(qc)==27 and len(retainedqc)==21 and len(accepted)==75
 assert sum(q["group1"]=="Internal-QC" for q in retainedqc)==6 and sum(q["group1"]=="QC1" for q in retainedqc)==15
 helpfile=ROOT/"LH_ratio_tables_2026.9.22—ALL-HELP.xlsx"
 ss,sheets=readbook(helpfile);changes={};removed=[];old_ids=None;peptides=None
 for i,root in enumerate(sheets,1):
  r=tab(root,ss);assert len(r)==98
  seq=[r0["A"] for r0 in r[1:]]
  if peptides is None:peptides=seq
  assert peptides==seq and len(set(seq))==97
  ids=[name for c,name in sorted(r[0].items(),key=lambda x:ci(x[0])) if c!="A"]
  if old_ids is None:old_ids=ids
  assert ids==old_ids and len(ids) in (75,109)
  keep={1:1};k=2
  for c,name in sorted(r[0].items(),key=lambda x:ci(x[0])):
   if c!="A" and name in accepted:keep[ci(c)]=k;k+=1
  assert len(keep)==76 and {x for x in ids if x in accepted}==accepted
  data=root.find(T("sheetData"))
  for row in data.findall(T("row")):
   rn=row.get("r");row.set("spans","1:76")
   for cell in list(row.findall(T("c"))):
    oc=ci(re.match(r"([A-Z]+)",cell.get("r")).group(1))
    if oc not in keep:row.remove(cell)
    else:cell.set("r",cn(keep[oc])+rn)
  dim=root.find(T("dimension"))
  if dim is not None:dim.set("ref","A1:BX98")
  changes["xl/worksheets/sheet"+str(i)+".xml"]=root
 removed=[x for x in old_ids if x not in accepted]
 if removed:rewrite(helpfile,changes)
 ssq,qroots=readbook(ROOT/"all_QC_metadata.xlsx")
 qroot=qroots[0];data=qroot.find(T("sheetData"));newno=0
 for row in list(data.findall(T("row"))):
  oldn=int(row.get("r"))
  if oldn>1 and qc[oldn-2] not in retainedqc:data.remove(row);continue
  newno+=1;row.set("r",str(newno))
  for cell in row.findall(T("c")):
   col=re.match(r"([A-Z]+)",cell.get("r")).group(1);cell.set("r",col+str(newno))
 qdim=qroot.find(T("dimension"))
 if qdim is not None:qdim.set("ref","A1:F22")
 if newno!=22:raise ValueError("Expected exactly 21 nonblank QC rows")
 if len(qc)>len(retainedqc):rewrite(ROOT/"all_QC_metadata.xlsx",{"xl/worksheets/sheet1.xml":qroot})
 pg=ROOT/"pg_matrix.tsv"
 with pg.open("r",encoding="utf8",newline="") as f:
  head=f.readline().rstrip("\r\n").split("\t")
  assert head[:5]==["Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description"]
  named=[x for x in head[5:] if x];extra=[x for x in named if x not in accepted]
  blank=sum(x=="" for x in head[5:])
  assert len(named) in (87,75)
  ix=[i for i,x in enumerate(head) if i<5 or x in accepted]
  assert len(ix)==80 and set(head[i] for i in ix[5:])==accepted
  t=tempfile.NamedTemporaryFile(mode="w",dir=ROOT,suffix=".tsv",encoding="utf8",delete=False,newline="")
  temp=Path(t.name);t.write("\t".join(head[i] for i in ix)+"\n");n=0
  for line in f:
   if not line.strip():continue
   cells=line.rstrip("\r\n").split("\t")
   assert len(cells)>ix[-1]
   t.write("\t".join(cells[i] for i in ix)+"\n");n+=1
  t.close()
 assert n==6613
 if len(head)!=80:temp.replace(pg)
 else:temp.unlink()
 hr=readbook(helpfile);tables=[tab(s,hr[0]) for s in hr[1]]
 cols=[[v for c,v in sorted(r[0].items(),key=lambda x:ci(x[0])) if c!="A"] for r in tables]
 assert len(cols)==3 and cols[0]==cols[1]==cols[2] and len(cols[0])==75 and set(cols[0])==accepted
 assert len(records(ROOT/"all_QC_metadata.xlsx"))==21
 with pg.open("r",encoding="utf8") as f:header=f.readline().rstrip("\r\n").split("\t")
 assert len(header)==80 and set(header[5:])==accepted
 report=["# 97-HELP metadata-based sample filtering","",
 "Retain all 54 study samples (27 HC + 27 S), 6 Internal-QC, and 15 QC1 pooled plasma; exclude 6 Blank QC records.",
 "All three HELP sheets retain 97 peptides x 75 sample columns.",
 "pg_matrix.tsv retains 6613 protein rows, 5 annotation columns + 75 sample columns.",
 "Original HELP sample columns: 109; removed 34. Original pg_matrix: 87 named sample columns plus 22 unnamed trailing columns; removed 12 named and 22 unnamed.",
 "This operation only deletes unused columns and Blank metadata rows; it does NOT alter retained intensity or ratio values.","",
 "## Removed HELP IDs",""]+[("- "+x) for x in removed]+["","## Removed pg_matrix IDs",""]+[("- "+x) for x in extra]
 (ROOT/"FILTERING_REPORT.md").write_text("\n".join(report)+"\n",encoding="utf8")
 readme=ROOT/"README.md";content=readme.read_text(encoding="utf8")
 if "## Filtered 97-HELP dataset" not in content:
  content+="\n## Filtered 97-HELP dataset\n\n54 study samples, 6 Internal-QC, and 15 QC1 pooled plasma samples are retained (75 total). Blank and metadata-unlisted samples have been excluded from L, H, ratio_LH, pg_matrix.tsv, and QC metadata. The original 97 peptides and retained intensity/ratio values are unchanged. See FILTERING_REPORT.md.\n"
  readme.write_text(content,encoding="utf8")
 print(json.dumps({"result":"PASS","samples":75,"help_rows":97,"protein_rows":n,"qc_rows":21,"removed_help":len(removed),"removed_pg":len(extra),"removed_blank_pg":blank}))
if __name__=="__main__":main()
