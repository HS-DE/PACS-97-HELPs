from collections import Counter
from yeast_protein_qc import ROOT,excel_rows,metadata
for f in ("all_QC_metadata.xlsx","all_sample_metadata.xlsx"):
 rows=excel_rows(ROOT/f)
 print("\nMETADATA",f,"n=",len(rows),"keys",list(rows[0]) if rows else [])
 for r in rows[:10]: print("ROW",r)
 print("GROUPS",Counter((x.get("group1"),x.get("group2")) for x in rows))
 print("PLATES",Counter((x.get("plate"),x.get("Plate"),x.get("plate_id"),x.get("Plate.ID")) for x in rows))
labels=metadata()
print("SAMPLE N",len(labels))
