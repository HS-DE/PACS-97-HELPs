"""QC-anchored Gaussian factorization pilot; not original omicsGMF."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
import argparse, warnings, json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.linalg import svd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import PACS_norm as p
warnings.filterwarnings('ignore',message='All-NaN slice encountered')
warnings.filterwarnings('ignore',message='Degrees of freedom')
med=np.nanmedian
F=np.isfinite
Z=lambda a:np.nan_to_num(a,nan=0.)
def ridge(y,x,pen):
    a=np.einsum('ni,nj,np->pij',x,x,F(y),optimize=True)+np.diag(pen)[None]
    return np.linalg.solve(a,(x.T@Z(y)).T[...,None]).squeeze(-1).T

def lr(y,k,it=4):
    if k==0:return np.zeros_like(y)
    mask=F(y);z=Z(y)
    for _ in range(it):
        _,u=np.linalg.eigh(z@z.T);u=u[:,-k:];pred=u@(u.T@z);z=np.where(mask,y,pred)
    return pred

def help_model(y,h,plates,tr):
    cv,_=p.linear_cv(h[tr].T,2);keep=(F(h[tr]).mean(0)>=.8)&F(cv)
    dh=h[:,keep]-med(h[tr][:,keep],axis=0);g=Z(med(dh,axis=1));shape=Z(dh-g[:,None]);_,_,v=svd(shape[tr],full_matrices=False,check_finite=False)
    t=np.column_stack([g,shape@v[:2].T]);t=(t-t[tr].mean(0))/np.maximum(t[tr].std(0),.05)
    x=np.column_stack([np.ones(len(y)),t]);co=ridge(y[tr],x[tr],[1e-6,2,2,2]);tech=t@co[1:];support=F(y[tr]).sum(0);tech[:,support<4]=0
    r=y-tech;ref=Z(med(r[tr],axis=0));levels=np.unique(plates);local=np.full((len(levels),y.shape[1]),np.nan);n=np.zeros_like(local)
    for i,b in enumerate(levels):
        q=tr&(plates==b)
        if q.any():local[i]=med(r[q],axis=0)-ref;n[i]=F(r[q]).sum(0)
    gplate=Z(med(local,axis=1));res=local-gplate[:,None];shared=lr(res,2);plate=.5*gplate[:,None]+.5*shared+n/(n+4)*Z(res-shared);plate[:,support<4]=0;plate[~np.isin(levels,plates[tr])]=0
    return tech,plate[np.searchsorted(levels,plates)],t,co,keep,plate

def soma(y,plates,tr):
    out=np.zeros_like(y);ref=med(y[tr],axis=0)
    for b in np.unique(plates):
        q=tr&(plates==b)
        if q.any():out[plates==b]=Z(med(y[q],axis=0)-ref)
    return out

def design(g):return pd.get_dummies(g,dtype=float).to_numpy()
def gmf(y,groups,plates):
    g=design(groups);b=design(plates);x=np.column_stack([g,b]);latent=np.zeros_like(y)
    for _ in range(3):
        co=ridge(y-latent,x,[1e-6]*g.shape[1]+[2]*b.shape[1]);latent=lr(y-x@co,3,2)
    return b@co[g.shape[1]:]
def pred(y,groups,k):
    g=design(groups);mu=g@ridge(y,g,[1e-6]*g.shape[1]);return mu+lr(y-mu,k)
def matrix(path,y,ids,pg):
    pd.DataFrame(np.exp2(y.T),index=pg,columns=ids).rename_axis('Protein.Group').to_csv(path,sep='\t',compression='gzip',float_format='%.7g')
