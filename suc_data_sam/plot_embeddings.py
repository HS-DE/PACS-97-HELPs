"""Comparable PCA of observed complete proteins, and HELP factor scores."""
import run_help_gmf as r
import sys
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
p=r.p;root=Path(sys.argv[1]);inp=Path(sys.argv[2]);m=pd.read_csv(root/'metadata_used.csv');q=pd.read_csv(inp/'Plate1-16_QC_Hela_merged.pg_matrix.tsv',sep='\t');q.columns=[p.clean_sample_id(c) for c in q];ids=m.sample_id.tolist();y=p.positive_log(q[ids].to_numpy(float)).T
complete=np.isfinite(y).all(0);model=PCA(n_components=2).fit(y[:,complete]);rows=[];methods=['Raw','Current_PACS','Gaussian_GMF','HELP_GMF','HELP_GMF_plate','Soma_like'];colors={'LM-QC':'#2875a8','QC':'#de8743','sample':'#659b6a'}
for ref in ['LM-QC','QC']:
 folder=root/ref;fig,axes=plt.subplots(2,3,figsize=(14,9),sharex=True,sharey=True);allcoords=[]
 for name in methods:
  z=y if name=='Raw' else p.positive_log(pd.read_csv(folder/(name+'_observed.tsv.gz'),sep='\t',index_col=0).to_numpy(float)).T;assert np.array_equal(np.isfinite(z),np.isfinite(y));xy=model.transform(z[:,complete]);allcoords.append(xy)
  for i,(a,b) in enumerate(xy):rows.append([ref,name,ids[i],m.group[i],int(m.plate[i]),a,b,int(complete.sum())])
 for ax,name,xy in zip(axes.flat,methods,allcoords):
  for g,c in colors.items():ix=m.group.eq(g);ax.scatter(xy[ix,0],xy[ix,1],s=18,c=c,label=g,alpha=.8)
  ax.set_title(name);ax.set_xlabel('PC1');ax.set_ylabel('PC2');ax.grid(alpha=.2)
 axes.flat[0].legend();fig.suptitle(ref+' reference | common raw PCA axes | '+str(complete.sum())+' complete PGs');fig.tight_layout();fig.savefig(folder/'pca_comparison.png',dpi=150);plt.close(fig)
 t=np.load(folder/'model.npz')['technical_scores'];pd.DataFrame({'sample_id':ids,'group':m.group,'plate':m.plate,'HELP_global_std':t[:,0],'HELP_shape_PC1_std':t[:,1],'HELP_shape_PC2_std':t[:,2]}).to_csv(folder/'help_factor_scores.csv',index=False)
pd.DataFrame(rows,columns=['reference','method','sample_id','group','plate','PC1','PC2','n_complete_PGs']).to_csv(root/'pca_scores.csv',index=False)
