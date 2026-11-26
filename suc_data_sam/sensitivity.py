"""Exploratory conservative HELP response; fixed candidates, not independent model selection."""
import run_help_gmf as r
import numpy as np,pandas as pd,sys
from pathlib import Path
p=r.p;med=r.med;F=r.F;Z=r.Z
root=Path(sys.argv[1]);inp=Path(sys.argv[2]);m=pd.read_csv(root/'metadata_used.csv');ids=m.sample_id.tolist();groups=m.group.to_numpy();plates=m.plate.to_numpy()
q=pd.read_csv(inp/'Plate1-16_QC_Hela_merged.pg_matrix.tsv',sep='\t');q.columns=[p.clean_sample_id(c) for c in q];y=p.positive_log(q[ids].to_numpy(float)).T;report=p.read_help_report(inp/'Accurate_Quantification_Report_v4_苏创.csv')
rows=[];audit=[];diags={}
for ref in ['LM-QC','QC']:
 tr=groups==ref;meta=m.copy();meta['group']=np.where(tr,'qc','sample');meta['condition']=np.where(tr,'QC',groups);meta['batch']=meta.batch.astype(str)
 heavy=p.extract_help_report(report,meta);h=p.positive_log(heavy.reindex(columns=ids).to_numpy(float)).T;old,info=p.correct_data(meta,y.T,h.T,heavy.index.tolist());diags[ref]={'HELP_extraction':heavy.attrs['HELP_extraction'],'PACS_groups':{str(k):{z:v.get(z) for z in ['HELP','RMSE','Pearson','x','status','selection_reason']} for k,v in info['groups'].items()}}
 g=Z(med(h-med(h[tr],axis=0),axis=1));model=np.load(root/ref/'model.npz');t=model['technical_scores'];tc=model['technical_correction']
 for i in range(len(ids)):audit.append([ref,ids[i],groups[i],plates[i],int(F(h[i]).sum()),float(g[i]),float(2**(-g[i])),float(np.max(abs(t[i])))])
 for validation in ['cross_type','same_type_OOF']:
  eval_group=('QC' if ref=='LM-QC' else 'LM-QC') if validation=='cross_type' else ref
  evalmask=groups==eval_group
  pred={f'gated_{a}':y.copy() for a in [.25,.5,1.]};pred.update({'Raw':y.copy(),'HELP_GMF':y.copy()})
  if validation=='cross_type':folds=[(tr,evalmask)]
  elif ref=='LM-QC':folds=[(tr&(plates!=b),tr&(plates==b)) for b in np.unique(plates)]
  else:folds=[(tr&np.array([i.endswith(s) for i in ids]),tr&~np.array([i.endswith(s) for i in ids])) for s in ['-37','-73']]
  for tt,te in folds:
   tech,_,t,co,_,_=r.help_model(y,h,plates,tt);gate=np.max(abs(t),axis=1)<=2.5;pred['HELP_GMF'][te]=(y-tech)[te]
   for a in [.25,.5,1.]:pred[f'gated_{a}'][te]=(y-a*tech*gate[:,None])[te]
  for name,z in pred.items():
   cv,_=p.linear_cv(z[evalmask].T,3);rows.append([ref,validation,eval_group,name,int(F(cv).sum()),float(med(cv))])
pd.DataFrame(rows,columns=['reference','validation','evaluation_group','method','n_proteins','median_CV_percent']).to_csv(root/'gated_sensitivity.csv',index=False)
pd.DataFrame(audit,columns=['reference','sample_id','group','plate','HELP_detected','median_HELP_log2_delta','median_HELP_multiplier','max_abs_standardized_factor']).to_csv(root/'help_audit.csv',index=False)
import json
(root/'diagnostics.json').write_text(json.dumps(diags,indent=2,default=lambda o:o.item() if hasattr(o,'item') else str(o)),encoding='utf-8')
print(pd.read_csv(root/'gated_sensitivity.csv').to_string(index=False))
