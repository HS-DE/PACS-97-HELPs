#!/usr/bin/env python3
"""Targeted QC sensitivity analyses for SomaScan-like PACS calibration."""
import argparse,csv,math
from pathlib import Path
from collections import Counter
import numpy as np
from run_soma_like import load_metadata,load_pg,load_help,sample_help_offset,calibrate,get_cv,safe_row_csv,nanmed

def eval_cv(raw,corrected,indices,min_n=3):
    b,n=get_cv(raw,indices,min_n=min_n)
    a,_=get_cv(corrected,indices,min_n=min_n)
    mask=np.isfinite(a)&np.isfinite(b)
    if not mask.any():return (0,np.nan,np.nan,np.nan)
    return int(mask.sum()),float(np.median(b[mask])),float(np.median(a[mask])),float(np.mean(a[mask]<b[mask]-1e-8))

def main(args):
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    meta=load_metadata(args.metadata)
    _,_,ids,pg=load_pg(args.pg,meta)
    p=np.array([meta[s]['Plate'] for s in ids]);groups=[meta[s]['group'] for s in ids]
    keys,H,bad=load_help(args.help_report,ids)
    off,nhelp,sel=sample_help_offset(H,p)
    flagged=[ids[i] for i,v in enumerate(off) if abs(v)>math.log2(2.5)]
    print('flagged',flagged)
    inputs={
        'none':np.zeros_like(off),
        'HELP_full':off,
        'HELP_gate_2.5x':np.where(abs(off)>math.log2(2.5),0,off),
        'HELP_clip_2.5x':np.clip(off,-math.log2(2.5),math.log2(2.5)),
    }
    raw=np.log2(pg)
    records=[]
    for method,ofs in inputs.items():
        x1=raw+ofs[:,None]
        for cal in ('LM-QC','QC'):
            fit=calibrate(x1,meta,ids,cal)
            for eval_group in ('LM-QC','QC'):
                ix=[i for i,g in enumerate(groups) if g==eval_group]
                for min_n in (3,8,16,24):
                    if len(ix)<min_n:continue
                    for stage,arr in [('stage1',x1),('plate',fit['stage3_plate_scaling']),('protein',fit['stage4_per_protein'])]:
                        nn,b,a,imp=eval_cv(raw,arr,ix,min_n)
                        records.append([method,cal,eval_group,'full_fit',stage,min_n,nn,b,a,imp])
        # QC OOF by plate, two folds each using other QC replicate per plate
        qc=[i for i,g in enumerate(groups) if g=='QC'];qc=sorted(qc,key=lambda i:(p[i],ids[i]))
        oof3=np.full_like(raw,np.nan);oof4=np.full_like(raw,np.nan)
        for fold in (0,1):
            train=[];test=[]
            for pi in sorted(set(p)):
                pair=[i for i in qc if p[i]==pi]
                assert len(pair)==2
                train.append(pair[1-fold]);test.append(pair[fold])
            fit=calibrate(x1,meta,ids,'QC',train_indices=set(train))
            oof3[test]=fit['stage3_plate_scaling'][test]
            oof4[test]=fit['stage4_per_protein'][test]
            for min_n in (3,8,12):
                for stage,arr in [('plate',oof3),('protein',oof4)]:
                    nn,b,a,imp=eval_cv(raw,arr,test,min_n)
                    records.append([method,'QC','QC',f'OOF_fold_{fold+1}',stage,min_n,nn,b,a,imp])
        for min_n in (3,8,16,24):
            for stage,arr in [('plate',oof3),('protein',oof4)]:
                nn,b,a,imp=eval_cv(raw,arr,qc,min_n)
                records.append([method,'QC','QC','OOF_combined',stage,min_n,nn,b,a,imp])
        # QC independent LM Cal, exclude anomalous HELP samples from evaluation only
        f=calibrate(x1,meta,ids,'LM-QC')
        qc_valid=[i for i,g in enumerate(groups) if g=='QC' and ids[i] not in flagged]
        for min_n in (3,8,16,24):
            for stage,arr in [('stage1',x1),('plate',f['stage3_plate_scaling']),('protein',f['stage4_per_protein'])]:
                nn,b,a,imp=eval_cv(raw,arr,qc_valid,min_n)
                records.append([method,'LM-QC','QC','independent_without_HELP_outliers',stage,min_n,nn,b,a,imp])
    safe_row_csv(out/'sensitivity_CV_summary.tsv',[
        'sample_norm','calibrator','evaluation_group','validation','stage','min_detections','proteins',
        'median_cv_raw','median_cv_corrected','improved_fraction'],records)
    # With two QC replicates on each plate, CV within plate should be unaffected by
    # pure per-plate multiplicative correction. Only a sample-specific factor may change it.
    raw_cv=[]
    for pi in sorted(set(p)):
        ix=[i for i,g in enumerate(groups) if g=='QC' and p[i]==pi]
        cv,n=get_cv(raw,ix,2);l1,n1=get_cv(raw+off[:,None],ix,2)
        mask=np.isfinite(cv)&np.isfinite(l1)
        raw_cv.append([pi,len(ix),int(mask.sum()),float(nanmed(cv[mask])),float(nanmed(l1[mask]))])
    safe_row_csv(out/'QC_intraplate_pair_CV.tsv',
            ['plate','QC_replicates','PG_detected_both','median_pair_cv_raw','median_pair_cv_HELP'],raw_cv)
    print('selected rows')
    for r in records:
        if r[0] in ('none','HELP_full','HELP_gate_2.5x') and r[1]=='QC' and r[3]=='OOF_combined' and r[4]=='protein' and r[5]==3:
            print(r)
    print('DONE',out)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--metadata',type=Path,required=True);p.add_argument('--pg',type=Path,required=True)
    p.add_argument('--help-report',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    main(p.parse_args())