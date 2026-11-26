"""Generate all intermediate model checkpoints from the three original inputs."""
import run_help_gmf as r
import sys
from pathlib import Path
import numpy as np,pandas as pd
p=r.p
root=Path(sys.argv[1]);inp=Path(sys.argv[2]);root.mkdir(parents=True,exist_ok=True)
try:m=pd.read_csv(inp/'metadata.csv',encoding='utf-8-sig')
except UnicodeDecodeError:m=pd.read_csv(inp/'metadata.csv',encoding='gb18030')
m.columns=m.columns.str.lower().str.strip();m=m[m.group.isin(['QC','LM-QC','sample'])][['sample_id','group','condition','batch','plate']].reset_index(drop=True);m.to_csv(root/'metadata_used.csv',index=False)
q=pd.read_csv(inp/'Plate1-16_QC_Hela_merged.pg_matrix.tsv',sep='\t');q.columns=[p.clean_sample_id(c) for c in q];ids=m.sample_id.tolist();pg=q['Protein.Group'].tolist();assert q.columns.is_unique and q['Protein.Group'].is_unique and m.sample_id.is_unique
y=p.positive_log(q[ids].to_numpy(float)).T;groups=m.group.to_numpy();plates=m.plate.to_numpy();report=p.read_help_report(inp/'Accurate_Quantification_Report_v4_苏创.csv')
for ref in ['LM-QC','QC']:
 folder=root/ref;folder.mkdir(exist_ok=True);tr=groups==ref;meta=m.copy();meta['group']=np.where(tr,'qc','sample');meta['condition']=np.where(tr,'QC',groups);meta['batch']=meta.batch.astype(str)
 heavy=p.extract_help_report(report,meta);h=p.positive_log(heavy.reindex(columns=ids).to_numpy(float)).T;tech,pc,t,co,keep,plate=r.help_model(y,h,plates,tr)
 np.savez_compressed(folder/'model.npz',technical_scores=t,coefficients=co,plate_effects=plate,technical_correction=tech,plate_correction=pc,protein_ids=pg,sample_ids=ids)
 if ref=='LM-QC':r.matrix(folder/'Gaussian_GMF_observed.tsv.gz',y-r.gmf(y,groups,plates),ids,pg)
