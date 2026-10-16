#!/usr/bin/env python3
"""Deliver original and two previously computed NYD 39-sample -fi PG matrices.

Outputs 3 same-order, same-scale log2 TSV matrices, a customer Excel workbook
with Readme plus 3 corresponding sheets, a sample-order key, and a ZIP.
No correction factors in the delivered matrices. Original results unchanged.
Uses Python's standard Open XML ZIP serialization for the workbook.
"""
from __future__ import annotations
import csv
import hashlib
import io
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path
from xml.sax.saxutils import escape
from zipfile import ZipFile,ZIP_DEFLATED

BASE=Path(__file__).resolve().parents[1]/"NYD_data"/"40samples"
SOURCE=BASE/"STGSWSTLK_Heavy_subtract_factor_39samples"
CLIENT=BASE/"client_visualizations_39samples"
OUT=BASE/"client_delivery_39samples_protein_matrices"
ANNOT=["Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description"]
EXCLUDED="MAP26092326E01"
SPEC=[
 ("Raw_log2","01_Raw_protein_matrix_log2.tsv",None),
 ("Overall_normalized","02_Overall_reference_normalized_protein_matrix_log2.tsv","global_subtract_fi"),
 ("Batch_normalized","03_Within_batch_normalized_protein_matrix_log2.tsv","batch_subtract_fi")
]
BOOK="NYD_39samples_Raw_and_Normalized_Protein_Matrices.xlsx"
BUNDLE="NYD_39samples_Protein_Matrices_Package.zip"
NS="http://schemas.openxmlformats.org/spreadsheetml/2006/main"

def read_tsv(path):
    with path.open(encoding="utf-8-sig",newline="") as file:
        data=csv.reader(file,delimiter="\t")
        return next(data),list(data)

def number(value):
    try:
        result=float(value)
        return result if math.isfinite(result) else None
    except (ValueError,TypeError):
        return None

def raw_log(value):
    value=number(value)
    return math.log2(value) if value is not None and value>0 else None

def hashfile(path):
    h=hashlib.sha256()
    with path.open("rb") as source:
        for part in iter(lambda:source.read(1024*1024),b""):
            h.update(part)
    return h.hexdigest()

def excel_column(idx):
    letters=""
    while idx:
        idx,v=divmod(idx-1,26)
        letters=chr(65+v)+letters
    return letters

def element(value,row,column,style):
    if value is None or value=="":
        return ""
    place=f"{excel_column(column)}{row}"
    attrs=f' r="{place}" s="{style}"'
    if isinstance(value,(int,float)) and not isinstance(value,bool):
        return f'<c{attrs}><v>{repr(value)}</v></c>'
    value=escape(str(value),{"\"":"&quot;","'":"&apos;"})
    return f'<c{attrs} t="inlineStr"><is><t xml:space="preserve">{value}</t></is></c>'

def worksheet(header,rows,is_info):
    file=io.StringIO()
    file.write('<?xml version="1.0" encoding="UTF-8"?>')
    file.write(f'<worksheet xmlns="{NS}">')
    if is_info:
        file.write('<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>')
        file.write('<cols><col min="1" max="1" width="27" customWidth="1"/>'
                   '<col min="2" max="2" width="88" customWidth="1"/></cols>')
    else:
        file.write('<sheetViews><sheetView workbookViewId="0"><pane xSplit="5" ySplit="1" topLeftCell="F2" activePane="bottomRight" state="frozen"/></sheetView></sheetViews>')
        file.write('<cols><col min="1" max="1" width="18" customWidth="1"/>'
                   '<col min="2" max="2" width="25" customWidth="1"/>'
                   '<col min="3" max="3" width="25" customWidth="1"/>'
                   '<col min="4" max="4" width="17" customWidth="1"/>'
                   '<col min="5" max="5" width="45" customWidth="1"/>'
                   '<col min="6" max="44" width="18" customWidth="1"/></cols>')
    file.write('<sheetData><row r="1" ht="25" customHeight="1">')
    for col,value in enumerate(header,1):
        file.write(element(value,1,col,1))
    file.write('</row>')
    for index,row in enumerate(rows,2):
        file.write(f'<row r="{index}">')
        for col,value in enumerate(row,1):
            file.write(element(value,index,col,3 if (not is_info and col>5) else 2))
        file.write('</row>')
    file.write('</sheetData>')
    if not is_info:
        file.write(f'<autoFilter ref="A1:{excel_column(len(header))}{len(rows)+1}"/>')
    file.write('</worksheet>')
    return file.getvalue().encode("utf-8")

def make_xlsx(path,sheets):
    count=len(sheets)
    doc={
      "[Content_Types].xml":('<?xml version="1.0" encoding="UTF-8"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        +''.join(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
                 for i in range(1,count+1))+
        '</Types>').encode(),
      "_rels/.rels":('<?xml version="1.0" encoding="UTF-8"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
        '</Relationships>').encode(),
      "xl/workbook.xml":('<?xml version="1.0" encoding="UTF-8"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<bookViews><workbookView/></bookViews><sheets>'+
        ''.join(f'<sheet name="{escape(s[0])}" sheetId="{i}" r:id="rId{i}"/>'
                for i,s in enumerate(sheets,1))+
        '</sheets></workbook>').encode(),
      "xl/_rels/workbook.xml.rels":('<?xml version="1.0" encoding="UTF-8"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'+
        ''.join(f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>'
                for i in range(1,count+1))+
        '<Relationship Id="rId9" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
        '</Relationships>').encode(),
      "xl/styles.xml":('<?xml version="1.0" encoding="UTF-8"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<numFmts count="1"><numFmt numFmtId="164" formatCode="0.000000"/></numFmts>'
        '<fonts count="3"><font><sz val="10"/><name val="Aptos"/></font>'
        '<font><b/><color rgb="FFFFFFFF"/><sz val="10"/><name val="Aptos"/></font>'
        '<font><sz val="10"/><color rgb="FF334155"/><name val="Aptos"/></font></fonts>'
        '<fills count="3"><fill><patternFill patternType="none"/></fill>'
        '<fill><patternFill patternType="gray125"/></fill>'
        '<fill><patternFill patternType="solid"><fgColor rgb="FF17324F"/></patternFill></fill></fills>'
        '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
        '<cellStyleXfs count="1"><xf fontId="0" fillId="0" borderId="0" numFmtId="0"/></cellStyleXfs>'
        '<cellXfs count="4">'
        '<xf fontId="0" fillId="0" borderId="0" numFmtId="0" xfId="0"/>'
        '<xf fontId="1" fillId="2" borderId="0" numFmtId="0" xfId="0" applyFont="1" applyFill="1"/>'
        '<xf fontId="2" fillId="0" borderId="0" numFmtId="0" xfId="0" applyFont="1"/>'
        '<xf fontId="0" fillId="0" borderId="0" numFmtId="164" xfId="0" applyNumberFormat="1"/>'
        '</cellXfs></styleSheet>').encode(),
    }
    for i,(name,header,data,summary) in enumerate(sheets,1):
        doc[f"xl/worksheets/sheet{i}.xml"]=worksheet(header,data,summary)
    with ZipFile(path,"w",compression=ZIP_DEFLATED,compresslevel=6) as z:
        for key,val in doc.items():
            z.writestr(key,val)
    with ZipFile(path,"r") as z:
        if z.testzip() is not None:
            raise AssertionError("XLSX archive integrity error")
        for name in z.namelist():
            ET.fromstring(z.read(name))
    return {x[0]:len(x[2]) for x in sheets}

def main():
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError("Client matrix package already exists; do not overwrite")
    with (CLIENT/"sample_label_key.csv").open(encoding="utf-8-sig",newline="") as f:
        samplekey=list(csv.DictReader(f))
    ids=[x["sample_id"] for x in samplekey]
    counts={}
    for row in samplekey:counts[row["batch"]]=counts.get(row["batch"],0)+1
    if len(ids)!=39 or len(set(ids))!=39 or EXCLUDED in ids or counts!={"T001":16,"T002":8,"T003":7,"T004":8}:
        raise AssertionError("Expected exactly the same 39-sample client cohort")
    input_header,input_rows=read_tsv(BASE/"NYD_pg_matrix2.tsv")
    if len(input_rows)!=4264 or len(input_header)!=45 or input_header[:5]!=ANNOT:
        raise AssertionError("Unexpected original NYD PG input")
    if set(input_header[5:])!=set(ids)|{EXCLUDED}:
        raise AssertionError("Unexpected original sample IDs")
    original_indices=[input_header.index(id) for id in ids]
    raw_rows=[
        source[:5]+[raw_log(source[j]) for j in original_indices]
        for source in input_rows
    ]
    header=ANNOT+ids
    normalized={}
    positives=sum(v is not None for row in raw_rows for v in row[5:])
    for name,_filename,branch in SPEC[1:]:
        path=SOURCE/branch/"corrected_NYD_pg_matrix2_log2.tsv"
        cols,rows=read_tsv(path)
        if len(cols)!=44 or len(rows)!=4264 or cols[:5]!=ANNOT or set(cols[5:])!=set(ids):
            raise AssertionError(name+": wrong corrected dimensions")
        colids=[cols.index(id) for id in ids]
        with (SOURCE/branch/"sample_factors_subtract_fi.csv").open(encoding="utf-8-sig",newline="") as f:
            f_by_sample={r["sample_id"]:float(r["fi_original_log2"]) for r in csv.DictReader(f)}
        if set(f_by_sample)!=set(ids):raise AssertionError("Incorrect factor cohort")
        batch=[]
        nvalid=0
        for i,(source,processed) in enumerate(zip(input_rows,rows)):
            if source[:5]!=processed[:5]:
                raise AssertionError("Protein-group annotation mismatch")
            converted=source[:5]
            for j,sid in enumerate(ids):
                before=raw_rows[i][5+j]
                raw_after=processed[colids[j]]
                after=number(raw_after)
                if before is None:
                    if raw_after!="":raise AssertionError("Nonpositive input was imputed")
                    converted.append(None)
                else:
                    expected=before-f_by_sample[sid]
                    if after is None or not math.isclose(after,expected,abs_tol=1e-10,rel_tol=1e-12):
                        raise AssertionError(f"Existing -fi correction mismatch {name}, row {i}, sample {sid}")
                    converted.append(after)
                    nvalid+=1
            batch.append(converted)
        if nvalid!=positives:raise AssertionError("Corrected detection rate changed")
        normalized[name]=batch
    OUT.mkdir(parents=True,exist_ok=False)
    files=[]
    for name,path,branch in SPEC:
        vals=raw_rows if branch is None else normalized[name]
        output=OUT/path
        with output.open("w",encoding="utf-8",newline="") as f:
            writer=csv.writer(f,delimiter="\t")
            writer.writerow(header)
            writer.writerows(vals)
        files.append(output)
    orderfile=OUT/"sample_order_and_batch.csv"
    with orderfile.open("w",encoding="utf-8-sig",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=("sample_id","batch","sample_label"))
        writer.writeheader()
        writer.writerows({k:r[k] for k in ("sample_id","batch","sample_label")} for r in samplekey)
    info=[
        ["Category","Description"],
        ["Study","NYD plasma protein group matrices"],
        ["Cohort","39 samples retained after sample quality control"],
        ["Batch counts","T001 16; T002 8; T003 7; T004 8"],
        ["Rows","4,264 protein groups"],
        ["Columns","5 protein annotations followed by 39 samples"],
        ["Units","log2 intensities, same scale in all three matrix sheets"],
        ["Raw_log2","Before normalization"],
        ["Overall_normalized","After overall-reference normalization"],
        ["Batch_normalized","After within-batch normalization"],
        ["Blanks","Missing or nonpositive original observations; no imputation"],
        ["Precision","Excel displays 6 decimal places; underlying values retained"],
        ["Matching","All matrices use identical protein and sample ordering"],
    ]
    sheets=[("Readme",info[0],info[1:],True)]
    for sheetname,path,branch in SPEC:
        sheets.append((sheetname,header,raw_rows if branch is None else normalized[sheetname],False))
    workbook=OUT/"NYD_39samples_Raw_and_Normalized_Protein_Matrices.xlsx"
    sheetcounts=make_xlsx(workbook,sheets)
    readme=OUT/"README.md"
    readme.write_text("\n".join([
        "# NYD 39-sample protein matrices · before and after normalization","",
        "These correspond exactly to the overall-reference and within-batch",
        "processed results used in the customer presentation figures.","",
        "All three matrices use 4,264 matching protein-group rows and",
        "39 identically ordered sample columns, plus five protein annotations.",
        "Every measurement is on the **log2 intensity scale**, including",
        "the original unnormalized data. Missing values are blank.","",
        "## Deliverables","",
        "- Excel workbook: "+workbook.name,
        "- Raw, before normalization: "+SPEC[0][1],
        "- After overall-reference normalization: "+SPEC[1][1],
        "- After within-batch normalization: "+SPEC[2][1],
        "- Sample and batch key: sample_order_and_batch.csv",
        "- Single-file download: NYD_39samples_Protein_Matrices_Package.zip","",
        "The same sample-quality-controlled cohort is used throughout.",
        "All numeric entries in the normalized sheets were checked against",
        "the completed underlying analyses. This package does not contain",
        "individual correction metadata or special target annotations.",""
    ])+"\n",encoding="utf-8")
    zip_path=OUT/"NYD_39samples_Protein_Matrices_Package.zip"
    with ZipFile(zip_path,"w",ZIP_DEFLATED,compresslevel=7) as z:
        for p in [readme,orderfile,workbook]+files:
            z.write(p,p.name)
    with ZipFile(zip_path,"r") as z:
        assert z.testzip() is None and len(z.namelist())==6
    manifest={
        "sample_count":len(ids),
        "protein_group_count":len(raw_rows),
        "sample_batch_counts":counts,
        "excluded_sample_id":EXCLUDED,
        "raw_valid_cells":positives,
        "measurement_scale":"log2 for all 3 matrices",
        "identical_sample_order":True,
        "identical_protein_annotation_order":True,
        "matching_missingness":True,
        "normalized_values_verified_against_existing_trial":True,
        "excel_sheet_row_counts":sheetcounts,
        "xlsx":workbook.name,
        "zip":zip_path.name,
        "matrix_files":[x[1] for x in SPEC],
        "input_PG_sha256":hashfile(BASE/"NYD_pg_matrix2.tsv")
    }
    (OUT/"verification_report.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    print("RESULT",json.dumps(manifest),flush=True)
    print("OUTPUT_SIZES",json.dumps({p.name:p.stat().st_size for p in [workbook,zip_path]+files}),flush=True)
    print("PASS 4264 x 39 matrices: Raw log2, Overall normalized, Batch normalized; Excel and ZIP",flush=True)

if __name__=="__main__":
    main()
