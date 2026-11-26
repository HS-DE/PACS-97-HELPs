import run_help_gmf as r
import numpy as np,pandas as pd,sys,json,hashlib
from pathlib import Path
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
p=r.p;med=r.med;F=r.F;Z=r.Z
inp=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
pgfile=inp/'pg_matrix(2).tsv'; metafile=inp/'metadata(1).xlsx'; hfile=inp/'Accurate_Quantification_Report_v4_10ppm_no_msconvert_remove_8ScanPoints_transitionAgreement.xlsx'
protein,m,raw,heavy,seq,ids=p.load_inputs(pgfile,hfile,metafile,plate_col='Plate');y=raw.T;h=heavy.T;pg=protein['Protein.Group'].tolist();tr=m.group.eq('qc').to_numpy();plates=m.plate.to_numpy();groups=m.condition.to_numpy();report=p.read_help_report(hfile);m.to_csv(out/'metadata_used.csv',index=False)
old,info=p.correct_data(m,raw,heavy,seq);tech,pc,t,co,keep,plate=r.help_model(y,h,plates,tr);hm=Z(med(h-med(h[tr],axis=0),axis=1))[:,None]
methods={'Raw':y,'Current_PACS':old.T,'HELP_median':y-hm,'Soma_like':y-r.soma(y,plates,tr),'Gaussian_GMF':y-r.gmf(y,groups,plates),'HELP_GMF':y-tech,'HELP_GMF_plate':y-tech-pc}
for a in [.25,.5,1.]:methods[f'Gated_{a:g}']=y-a*tech*(np.max(abs(t),axis=1)<=2.5)[:,None]
cols=['method','scope','n_samples','n_proteins','raw_median_CV','median_CV','pct_CV_lt20','pct_CV_lt10','pct_improved'];rows=[];per=[];shifts=[]
def summarize(name,z,mask,scope,rows):
 cv,n=p.linear_cv(z[mask].T,3);base,_=p.linear_cv(y[mask].T,3);ok=F(cv);assert np.array_equal(ok,F(base));rows.append([name,scope,int(mask.sum()),int(ok.sum()),float(med(base)),float(med(cv)),float((cv[ok]<20).mean()*100),float((cv[ok]<10).mean()*100),float((cv[ok]<base[ok]).mean()*100)]);return cv,n
for name,z in methods.items():
 assert np.array_equal(F(y),F(z));r.matrix(out/(name+'_observed.tsv.gz'),z,ids,pg)
 for scope,mask in [('QC',tr),('HC',groups=='HC'),('S',groups=='S')]+[(f'QC_plate_{b}',tr&(plates==b)) for b in np.unique(plates)]:
  cv,n=summarize(name,z,mask,scope,rows);per.extend([name,scope,pg[j],int(n[j]),float(cv[j])] for j in np.flatnonzero(F(cv)))
 for i in range(len(ids)):shifts.append([name,ids[i],float(med(y[i]-z[i])),float(med(abs(y[i]-z[i])))])
sumdf=pd.DataFrame(rows,columns=cols);sumdf.to_csv(out/'summary.csv',index=False,float_format='%.8g');pd.DataFrame(per,columns=['method','scope','Protein.Group','n_observed','CV_percent']).to_csv(out/'protein_cv.csv.gz',index=False,compression='gzip');pd.DataFrame(shifts,columns=['method','sample_id','median_log2_removed','median_abs_log2_removed']).to_csv(out/'sample_corrections.csv',index=False)
valrows=[];assignments=[]
for mode in ['within_plate_5fold','leave_one_plate_out']:
 keys=['Raw','HELP_median','HELP_GMF','Gated_0.25','Gated_0.5','Gated_1']+(['Soma_like','HELP_GMF_plate'] if mode.startswith('within') else [])
 pred={k:y.copy() for k in keys}
 if mode.startswith('within'):
  ix=[sorted(np.flatnonzero(tr&(plates==b)),key=lambda i:ids[i]) for b in np.unique(plates)];assert all(len(a)==5 for a in ix);folds=[np.isin(np.arange(len(y)),[a[k] for a in ix]) for k in range(5)]
 else:folds=[tr&(plates==b) for b in np.unique(plates)]
 for k,te in enumerate(folds):
  tt=tr&~te;assignments.extend([mode,k+1,ids[i],plates[i]] for i in np.flatnonzero(te));mm=m.copy();mm.loc[te,'group']='sample';hh=p.positive_log(p.extract_help_report(report,mm).reindex(columns=ids).to_numpy(float)).T
  tc,pl,ts,*_=r.help_model(y,hh,plates,tt);ds={'Raw':np.zeros_like(y),'HELP_median':Z(med(hh-med(hh[tt],axis=0),axis=1))[:,None],'HELP_GMF':tc,'Soma_like':r.soma(y,plates,tt),'HELP_GMF_plate':tc+pl};ds.update({f'Gated_{a:g}':a*tc*(np.max(abs(ts),axis=1)<=2.5)[:,None] for a in [.25,.5,1.]})
  for name in keys:pred[name][te]=(y-ds[name])[te]
 for name,z in pred.items():summarize(name,z,tr,mode,valrows);r.matrix(out/(mode+'_'+name+'.tsv.gz'),z[tr],np.array(ids)[tr],pg)
valdf=pd.DataFrame(valrows,columns=cols);valdf.to_csv(out/'validation.csv',index=False,float_format='%.8g');pd.DataFrame(assignments,columns=['validation','fold','sample_id','plate']).to_csv(out/'folds.csv',index=False)
rng=np.random.default_rng(20261126);u=rng.random(y.shape);eligible=F(y)&~tr[:,None]&(F(y).sum(0)>=8)[None];val=eligible&(u<.05);test=eligible&(u>=.05)&(u<.1);masked=y.copy();masked[val|test]=np.nan;imp=[]
for name,d in [('Raw',np.zeros_like(y)),('HELP_GMF',tech),('HELP_GMF_plate',tech+pc)]:
 scores=[]
 for k in [0,2,5]:
  yh=r.pred(masked-d,groups,k)+d;scores.append((float(np.sqrt(np.mean((yh[val]-y[val])**2))),k,yh))
 vr,k,yh=min(scores,key=lambda a:a[0]);e=yh[test]-y[test];imp.append([name,k,vr,float(np.sqrt(np.mean(e**2))),float(np.mean(abs(e))),int(test.sum())])
 if name=='HELP_GMF_plate':
  z=y-d;fill=~F(z)&(F(z).mean(0)>=.5)[None];filled=np.where(fill,r.pred(z,groups,k),z);assert np.array_equal(filled[F(z)],z[F(z)]);r.matrix(out/'HELP_GMF_plate_imputed.tsv.gz',filled,ids,pg);pd.DataFrame(fill.T,index=pg,columns=ids).rename_axis('Protein.Group').to_csv(out/'imputed_mask.tsv.gz',sep='\t',compression='gzip')
pd.DataFrame(imp,columns=['method','rank','validation_RMSE_log2','test_RMSE_log2','test_MAE_log2','n_test']).to_csv(out/'imputation_benchmark.csv',index=False,float_format='%.8g')
shift=[];sel=groups=='S';yy=y.copy();yy[sel]+=np.log2(1.5)
for name,z in methods.items():
 if name=='Current_PACS':zz=p.correct_data(m,yy.T,heavy,seq)[0].T
 elif name=='Gaussian_GMF':zz=yy-r.gmf(yy,groups,plates)
 else:zz=z.copy();zz[sel]+=np.log2(1.5)
 change=zz-z;shift.append([name,float(med(change[sel])),float(med(abs(change[sel]-np.log2(1.5)))),float(med(abs(change[~sel])))])
pd.DataFrame(shift,columns=['method','median_recovered_shift','median_abs_shift_error','median_abs_unperturbed_change']).to_csv(out/'synthetic_shift.csv',index=False,float_format='%.8g')
np.savez_compressed(out/'model.npz',technical_scores=t,coefficients=co,plate_effects=plate,technical_correction=tech,plate_correction=pc,protein_ids=pg,sample_ids=ids)
cvh,_=p.linear_cv(h[tr].T,2);pd.DataFrame({'HELP':seq,'QC_CV':cvh,'selected_GMF':keep,'n_detected':F(h).sum(0)}).to_csv(out/'HELP_quality.csv',index=False);audit=pd.DataFrame({'sample_id':ids,'condition':groups,'plate':plates,'n_PG':F(y).sum(1),'n_HELP':F(h).sum(1),'HELP_median_multiplier':np.exp2(-hm[:,0]),'T1':t[:,0],'T2':t[:,1],'T3':t[:,2]});audit.to_csv(out/'HELP_audit.csv',index=False)
diag={'n_input_PG':len(pg),'n_retained':len(ids),'n_observed_PG':int(F(y).any(0).sum()),'n_complete_PG':int(F(y).all(0).sum()),'missing_fraction':float((~F(y)).mean()),'n_HELP':len(seq),'n_complete_HELP':int(F(h).all(0).sum()),'n_GMF_HELP':int(keep.sum()),'min_QC_HELP_CV':float(np.nanmin(cvh)),'PACS_groups':{str(k):{j:v.get(j) for j in ['HELP','RMSE','Pearson','x','status','selection_reason']} for k,v in info['groups'].items()},'sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [pgfile,metafile,hfile]}}
(out/'diagnostics.json').write_text(json.dumps(diag,indent=2,default=str))
for table,fn in [(sumdf.query("scope=='QC'"),'cv_comparison'),(valdf.query("scope=='within_plate_5fold'"),'heldout_cv')]:
 fig,ax=plt.subplots(figsize=(10,5));ax.bar(table.method,table.median_CV);ax.axhline(table.raw_median_CV.iloc[0],ls='--',color='gray');ax.set_ylabel('Median protein CV (%)');ax.tick_params(axis='x',rotation=35);fig.tight_layout();fig.savefig(out/(fn+'.png'),dpi=160);plt.close(fig)
complete=F(y).all(0);model=PCA(n_components=2).fit(y[:,complete]);fig,axes=plt.subplots(2,4,figsize=(15,8),sharex=True,sharey=True);coords=[]
for ax,(name,z) in zip(axes.flat,list(methods.items())[:7]):
 xy=model.transform(z[:,complete])
 for g in np.unique(groups):
  sel=groups==g;ax.scatter(xy[sel,0],xy[sel,1],label=g,s=25)
 ax.set_title(name);ax.set_xlabel('PC1');ax.set_ylabel('PC2');coords.extend([name,ids[i],groups[i],plates[i],*xy[i]] for i in range(len(ids)))
axes.flat[-1].axis('off');axes.flat[0].legend();fig.tight_layout();fig.savefig(out/'pca_comparison.png',dpi=160);plt.close(fig);pd.DataFrame(coords,columns=['method','sample_id','condition','plate','PC1','PC2']).to_csv(out/'pca_scores.csv',index=False)
print(json.dumps(diag,indent=2,default=str));print(sumdf.query("scope=='QC'").to_string(index=False));print(valdf.to_string(index=False));print(pd.DataFrame(imp).to_string(index=False));print(audit.sort_values('HELP_median_multiplier',ascending=False).head().to_string(index=False))
