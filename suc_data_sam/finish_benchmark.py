"""Checkpoint-aware benchmark: use completed observed outputs and model snapshots."""
import run_help_gmf as r
import numpy as np,pandas as pd,sys
from pathlib import Path
import matplotlib.pyplot as plt
p=r.p;med=r.med;F=r.F;Z=r.Z
root=Path(sys.argv[1]);inp=Path(sys.argv[2]);ref=sys.argv[3];folder=root/ref
m=pd.read_csv(root/'metadata_used.csv');ids=m.sample_id.tolist();groups=m.group.to_numpy();plates=m.plate.to_numpy();tr=groups==ref
q=pd.read_csv(inp/'Plate1-16_QC_Hela_merged.pg_matrix.tsv',sep='\t');q.columns=[p.clean_sample_id(c) for c in q];pg=q['Protein.Group'].tolist();y=p.positive_log(q[ids].to_numpy(float)).T
meta=m.copy();meta['group']=np.where(tr,'qc','sample');meta['condition']=np.where(tr,'QC',groups);meta['batch']=meta.batch.astype(str)
report=p.read_help_report(inp/'Accurate_Quantification_Report_v4_苏创.csv');heavy=p.extract_help_report(report,meta);h=p.positive_log(heavy.reindex(columns=ids).to_numpy(float)).T;seq=heavy.index.tolist();old,info=p.correct_data(meta,y.T,h.T,seq);old=old.T
model=np.load(folder/'model.npz');tech=model['technical_correction'];pc=model['plate_correction'];assert model['sample_ids'].tolist()==ids and model['protein_ids'].tolist()==pg
raw_gmf=p.positive_log(pd.read_csv(root/'LM-QC'/'Gaussian_GMF_observed.tsv.gz',sep='\t',index_col=0).to_numpy(float)).T;gc=Z(y-raw_gmf)
hm=Z(med(h-med(h[tr],axis=0),axis=1))[:,None];methods={'Raw':np.zeros_like(y),'Current_PACS':np.broadcast_to(Z(med(y-old,axis=1))[:,None],y.shape),'HELP_median':hm,'Soma_like':r.soma(y,plates,tr),'Gaussian_GMF':gc,'HELP_GMF':tech,'HELP_GMF_plate':tech+pc}
rows=[];cvrows=[];shifts=[];lopo=[];fold=[];imp=[]
for method,d in methods.items():
 z=y-d;assert np.array_equal(F(z),F(y))
 if method!='Raw':r.matrix(folder/(method+'_observed.tsv.gz'),z,ids,pg)
 for i in range(len(y)):shifts.append([ref,method,ids[i],groups[i],plates[i],float(med((y-z)[i])),float(med(abs((y-z)[i])))])
 for g0 in ['LM-QC','QC','sample']:
  ix=groups==g0;cv,n=p.linear_cv(z[ix].T,3);ok=F(cv);base,_=p.linear_cv(y[ix].T,3);rows.append([ref,method,g0,int(ix.sum()),int(ok.sum()),float(med(cv)),float((cv[ok]<20).mean()*100),float((cv[ok]<10).mean()*100),float((cv[ok]<base[ok]).mean()*100)]);cvrows.extend([ref,method,g0,pg[j],int(n[j]),float(cv[j])] for j in np.flatnonzero(ok))
summary=pd.DataFrame(rows,columns=['reference','method','evaluation_group','n_samples','n_proteins','median_CV_percent','pct_CV_lt20','pct_CV_lt10','pct_improved_vs_raw']);summary.to_csv(folder/'summary.csv',index=False,float_format='%.6g');pd.DataFrame(cvrows,columns=['reference','method','evaluation_group','Protein.Group','n_observed','CV_percent']).to_csv(folder/'protein_cv.csv.gz',index=False,compression='gzip');pd.DataFrame(shifts,columns=['reference','method','sample_id','group','plate','median_log2_removed','median_abs_log2_removed']).to_csv(folder/'sample_corrections.csv',index=False)
ax=summary.pivot(index='method',columns='evaluation_group',values='median_CV_percent').reindex(list(methods)).plot.bar(figsize=(11,5));ax.set_ylabel('Median protein CV (%)');ax.set_title(ref+' training reference');ax.tick_params(axis='x',rotation=25);plt.tight_layout();plt.savefig(folder/'cv_comparison.png',dpi=160);plt.close()
lpred=y.copy();lmed=y.copy()
for b in np.unique(plates):
 tt=tr&(plates!=b);te=plates==b;tc,*_=r.help_model(y,h,plates,tt);lpred[te]=y[te]-tc[te];lmed[te]=y[te]-Z(med(h[te]-med(h[tt],axis=0),axis=1))[:,None]
for name,z in [('Raw',y),('HELP_median',lmed),('HELP_GMF',lpred)]:
 for g0 in ['LM-QC','QC']:
  cv,_=p.linear_cv(z[groups==g0].T,3);lopo.append([ref,g0,name,int(F(cv).sum()),float(med(cv))])
pd.DataFrame(lopo,columns=['reference','evaluation_group','method','n_proteins','median_CV_percent']).to_csv(folder/'leave_one_plate_out.csv',index=False)
if ref=='QC':
 oo={k:y.copy() for k in ['Raw','HELP_median','Soma_like','HELP_GMF','HELP_GMF_plate']};pair=[]
 for suffix in ['-37','-73']:
  tt=tr&np.array([i.endswith(suffix) for i in ids]);te=tr&~tt;tc,pl,*_=r.help_model(y,h,plates,tt);hh=Z(med(h-med(h[tt],axis=0),axis=1))[:,None];cc={'Raw':np.zeros_like(y),'HELP_median':hh,'Soma_like':r.soma(y,plates,tt),'HELP_GMF':tc,'HELP_GMF_plate':tc+pl}
  for name,d in cc.items():oo[name][te]=(y-d)[te]
 for name,z in oo.items():
  cv,_=p.linear_cv(z[tr].T,3);pair.append([name,int(F(cv).sum()),float(med(cv))])
 pd.DataFrame(pair,columns=['method','n_proteins','median_CV_percent']).to_csv(folder/'same_type_QC_OOF.csv',index=False)
rng=np.random.default_rng(20261126);u=rng.random(y.shape);eligible=F(y)&(groups=='sample')[:,None]&(F(y).sum(0)>=8)[None];val=eligible&(u<.05);test=eligible&(u>=.05)&(u<.1);masked=y.copy();masked[val|test]=np.nan
for method in ['Raw','HELP_GMF','HELP_GMF_plate']:
 d=methods[method];zm=masked-d;scores=[]
 for k in [0,2,5]:
  yh=r.pred(zm,groups,k)+d;scores.append((float(np.sqrt(np.mean((yh[val]-y[val])**2))),k,yh))
 vr,k,yh=min(scores,key=lambda a:a[0]);e=yh[test]-y[test];imp.append([ref,method,k,vr,float(np.sqrt(np.mean(e**2))),float(np.mean(abs(e))),int(test.sum())])
 if method=='HELP_GMF_plate':
  z=y-d;elig=(~F(z))&(F(z).mean(0)>=.5)[None];yh=r.pred(z,groups,k);filled=np.where(elig,yh,z);assert np.array_equal(filled[F(z)],z[F(z)]);r.matrix(folder/'HELP_GMF_plate_imputed.tsv.gz',filled,ids,pg);pd.DataFrame(elig.T,index=pg,columns=ids).rename_axis('Protein.Group').to_csv(folder/'imputed_mask.tsv.gz',sep='\t',compression='gzip')
pd.DataFrame(imp,columns=['reference','method','rank','validation_RMSE_log2','test_RMSE_log2','test_MAE_log2','n_test']).to_csv(folder/'imputation_benchmark.csv',index=False)
yy=y.copy();sel=(groups=='sample')&(plates%2==0);yy[sel]+=np.log2(1.5)
for method,d in methods.items():
 if method=='Gaussian_GMF':continue
 z=y-d;zz=p.correct_data(meta,yy.T,h.T,seq)[0].T if method=='Current_PACS' else yy-d;change=zz-z;fold.append([ref,method,float(med(change[sel].ravel())),float(med(abs(change[sel]-np.log2(1.5)).ravel())),float(med(abs(change[~sel]).ravel()))])
pd.DataFrame(fold,columns=['reference','method','median_recovered_shift','median_abs_shift_error','median_abs_unperturbed_change']).to_csv(folder/'synthetic_shift.csv',index=False)
print(summary.to_string(index=False));print('complete',ref)
