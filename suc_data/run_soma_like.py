#!/usr/bin/env python3
"""SomaScan-inspired HELP + calibrator normalization for the 16-plate PACS dataset.

Not the proprietary SomaScan pipeline. Steps 1,2,3,4 mimic sample scaling,
calibrator-only median normalization, plate scaling, and per-protein calibration.
Step 5 (ANML) is not run without its necessary external population reference.
All missing PG values remain missing. No response/clinical normalization.
Requires NumPy; XLSX reading uses only the Python standard library.
"""
import argparse
import csv
import json
import re
import zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
T = '{%s}' % NS


def canonical(value):
    m = re.search(r'P(\d+)-0*(\d+)(?:\.raw)?$', str(value), re.I)
    if not m:
        raise ValueError('Cannot parse sample name: '+str(value))
    return f'P{int(m[1])}-{int(m[2]):02d}'


def xlsx_rows(path, sheet_number=1):
    """Streaming basic XLSX worksheet values, without loading the entire workbook."""
    with zipfile.ZipFile(path) as z:
        if 'xl/sharedStrings.xml' in z.namelist():
            root = ET.fromstring(z.read('xl/sharedStrings.xml'))
            strings = [''.join(e.text or '' for e in si.iter(T+'t')) for si in root.findall(T+'si')]
        else:
            strings = []
        with z.open(f'xl/worksheets/sheet{sheet_number}.xml') as stream:
            for _, elem in ET.iterparse(stream, events=('end',)):
                if elem.tag != T+'row':
                    continue
                cells = {}
                for cell in elem.findall(T+'c'):
                    address = cell.attrib['r']
                    letters = re.match(r'[A-Z]+',address).group()
                    ci = 0
                    for letter in letters:
                        ci = ci*26 + ord(letter)-64
                    ci -= 1
                    node = cell.find(T+'v')
                    if node is None:
                        inline = cell.find(T+'is')
                        if inline is not None:
                            cells[ci] = ''.join(e.text or '' for e in inline.iter(T+'t'))
                        continue
                    val = node.text or ''
                    if cell.attrib.get('t') == 's':
                        val = strings[int(val)]
                    cells[ci] = val
                n = max(cells.keys(), default=-1)+1
                yield [cells.get(i,'') for i in range(n)]
                elem.clear()


def load_metadata(path):
    rows = xlsx_rows(path)
    head = next(rows)
    out = {}
    for row in rows:
        fields = dict(zip(head,row))
        if not fields.get('sample_id'):
            continue
        sid = canonical(fields['sample_id'])
        if sid in out:raise ValueError('Duplicate metadata ID '+sid)
        fields['sample_id'] = sid
        fields['Plate'] = int(float(fields['Plate']))
        out[sid] = fields
    return out


def load_pg(path, meta, include_groups=('LM-QC','QC','sample')):
    with path.open(newline='',encoding='utf-8-sig') as f:
        reader=csv.reader(f,delimiter='\t')
        head=next(reader)
        ids=[canonical(s) for s in head[5:]]
        assert len(ids)==len(set(ids)), 'Duplicate PG sample identifiers'
        bad=set(ids)-set(meta)
        if bad:raise ValueError('PG samples missing from metadata '+str(bad))
        want=[i for i,s in enumerate(ids) if meta[s]['group'] in include_groups]
        sample_ids=[ids[i] for i in want]
        rownames=[];annot=[];values=[]
        for j,row in enumerate(reader):
            if len(row)<5:continue
            annot.append(row[:5])
            rownames.append(f'PG{j+1:05d}')
            values.append([float(row[5+i]) if row[5+i] and float(row[5+i])>0 else np.nan for i in want])
    pg = np.asarray(values,dtype='float64').T
    if len(set(sample_ids)) != len(sample_ids): raise ValueError('PG duplicate sample')
    return head[:5], annot, sample_ids, pg


def load_help(path, sample_ids):
    rows = xlsx_rows(path,1)
    header=next(rows)
    col={h:i for i,h in enumerate(header)}
    needed=('Sample','Peptide_Key','Accurate_AUC','Peak_Scan_Points','N_Common_Transitions','Peak_Detection_QC_Flag')
    if not set(needed).issubset(col):raise ValueError('Missing HELP columns: '+str(set(needed)-set(col)))
    sidset=set(sample_ids)
    items={};bad=Counter();keys=set()
    for row in rows:
        def get(s):return row[col[s]] if col[s]<len(row) else ''
        sample=get('Sample')
        if not sample:continue
        sid=canonical(sample)
        if sid not in sidset:continue
        peptide=get('Peptide_Key');keys.add(peptide)
        try:
            intensity=float(get('Accurate_AUC'))
            nscan=float(get('Peak_Scan_Points'))
            ncommon=float(get('N_Common_Transitions'))
        except (ValueError,TypeError):
            bad['not_numeric']+=1;continue
        if get('Peak_Detection_QC_Flag')!='PASS' or nscan<8 or ncommon<2 or intensity<=0:
            bad['fail_qc']+=1;continue
        if (sid,peptide) in items:raise ValueError('Duplicate HELP '+str((sid,peptide)))
        items[(sid,peptide)] = intensity
    peptides=sorted(keys)
    out=np.full((len(sample_ids),len(peptides)),np.nan)
    iids={s:i for i,s in enumerate(sample_ids)};pids={k:i for i,k in enumerate(peptides)}
    for (sid,key),v in items.items():out[iids[sid],pids[key]]=v
    return peptides,out,bad


def nanmed(x,axis=None):
    with np.errstate(all='ignore'):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',category=RuntimeWarning)
            return np.nanmedian(x,axis=axis)


def get_cv(log, idx, min_n=3):
    """CV in *linear* intensities, detected values only, ddof=1."""
    a=np.exp2(log[np.asarray(idx)])
    n=np.sum(np.isfinite(a),axis=0)
    mn=np.nansum(a,axis=0)/np.maximum(n,1)
    v=np.nansum((a-mn)**2,axis=0)/np.maximum(n-1,1)
    cv=np.full(a.shape[1],np.nan)
    ok=(n>=min_n)&(mn>0)
    cv[ok]=100*np.sqrt(v[ok])/mn[ok]
    return cv,n


def median_cv(log,idx):
    cv,_=get_cv(log,idx)
    return float(nanmed(cv)) if np.isfinite(cv).any() else None


def qcs_summary(log_dict,group,metadata,ids,min_n=3):
    idx=[i for i,s in enumerate(ids) if metadata[s]['group']==group]
    cv_raw,n=get_cv(log_dict['Raw'],idx,min_n)
    out=[]
    for stage,log in log_dict.items():
        cv,_=get_cv(log,idx,min_n)
        mask=np.isfinite(cv)&np.isfinite(cv_raw)
        if not mask.any():continue
        b=cv_raw[mask];a=cv[mask]
        out.append({'evaluation_group':group,'stage':stage,'proteins':int(mask.sum()),
                    'median_cv_raw':float(np.median(b)),'median_cv_corrected':float(np.median(a)),
                    'cv_delta_percentage_points':float(np.median(a)-np.median(b)),
                    'improved_proteins':int((a<b-1e-9).sum()),
                    'improved_fraction':float((a<b-1e-9).mean())})
    return out


def safe_row_csv(path,header,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(header)
        for row in rows:w.writerow(row)


def write_matrix(path,head,annot,ids,log):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(head+ids)
        for j,an in enumerate(annot):
            w.writerow(an+[f'{2**v:.8g}' if np.isfinite(v) else '' for v in log[:,j]])


def sample_help_offset(help_a,plate,pep_min_frac=.75):
    # Only qualified heavy internal standards (fixed spike quantities) are used.
    ok=np.isfinite(help_a).mean(axis=0)>=pep_min_frac
    h=np.log2(help_a[:,ok]);offs=np.zeros(len(plate));npep=np.zeros(len(plate),int)
    for p in sorted(set(plate)):
        ins=np.where(plate==p)[0]
        # Plate-local HELP reference avoids conflating a systematic plate shift
        # with individual-sample readout effects.
        ref=nanmed(h[ins],axis=0)
        diff=h[ins]-ref
        for i,k in enumerate(ins):
            finite=np.isfinite(diff[i]);npep[k]=finite.sum()
            offs[k]=-float(nanmed(diff[i])) if finite.sum()>=5 else 0.
        offs[ins]-=np.median(offs[ins])
    return offs,npep,ok


def calibrate(log1,meta,ids,cal_group,train_indices=None,min_plates=8):
    """SomaScan-inspired steps 2,3,4. x: sample x protein log2."""
    plate=np.array([meta[s]['Plate'] for s in ids]);plates=sorted(set(plate))
    idx=[i for i,s in enumerate(ids) if meta[s]['group']==cal_group]
    if train_indices is not None:idx=[i for i in idx if i in train_indices]
    per_plate={p:[i for i in idx if plate[i]==p] for p in plates}
    if any(len(z)==0 for z in per_plate.values()):raise ValueError('Need calibrators on every plate')
    log2=log1.copy()
    c_off=np.zeros(len(ids))
    for p,ii in per_plate.items():
        if len(ii)<2:continue
        ref=nanmed(log1[ii],axis=0)
        for i in ii:
            v=log1[i]-ref
            off=-float(nanmed(v))
            log2[i]+=off
            c_off[i]=off
    # Only calibrators receive the intraplate median-normalization step.
    per_protein_plate=np.array([nanmed(log2[per_plate[p]],axis=0) for p in plates])
    global_ref=nanmed(per_protein_plate,axis=0)
    feature_plate_counts=np.sum(np.isfinite(per_protein_plate),axis=0)
    enough=feature_plate_counts>=min_plates
    plate_off=[]
    for j,p in enumerate(plates):
        d=global_ref-per_protein_plate[j]
        usable=np.isfinite(d)&enough
        plate_off.append(float(nanmed(d[usable])) if usable.sum()>=50 else 0.)
    plate_off=np.array(plate_off)
    log3=log2.copy()
    for j,p in enumerate(plates):log3[plate==p]+=plate_off[j]
    # Per-protein residual plate factors, evaluated after plate-scale correction.
    ref3=per_protein_plate+plate_off[:,None]
    global_ref3=nanmed(ref3,axis=0)
    shifts=np.full_like(ref3,np.nan)
    for j in range(len(plates)):
        shifts[j]=global_ref3-ref3[j]
    shifts[:,~enough]=np.nan
    log4=log3.copy()
    for j,p in enumerate(plates):
        shift=np.where(np.isfinite(shifts[j]),shifts[j],0.)
        log4[plate==p]+=shift
    return {'stage3_plate_scaling':log3,'stage4_per_protein':log4,
            'plate_off':dict(zip(plates,plate_off)),
            'protein_off':shifts, 'calibrator_off':c_off,
            'calibrator_index':idx, 'features_reference_eligible':enough,
            'n_per_plate':{str(p):len(q) for p,q in per_plate.items()}}


def produce(args):
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    meta=load_metadata(args.metadata)
    h,ann,ids,pg=load_pg(args.pg,meta)
    plate=np.array([meta[s]['Plate'] for s in ids]);groups=np.array([meta[s]['group'] for s in ids])
    peptides,help_a,bad=load_help(args.help_report,ids)
    log0=np.log2(pg)
    sample_off,npep,selected=sample_help_offset(help_a,plate)
    log1=log0+sample_off[:,None]
    print('DATA',len(ann),'PGs',len(ids),'eligible samples; HELP',len(peptides), 'selected',int(selected.sum()),flush=True)
    if (npep<5).any():print('WARNING low HELP coverage',[(ids[i],int(z)) for i,z in enumerate(npep) if z<5])
    summary=[];branch_data={}
    for cal_group in ('LM-QC','QC'):
        print('CALIBRATOR',cal_group,flush=True)
        fit=calibrate(log1,meta,ids,cal_group,min_plates=args.min_plates)
        stage_map={'Raw':log0,'stage1_HELP_sample':log1,
                   'stage3_plate_scaling':fit['stage3_plate_scaling'],
                   'stage4_per_protein':fit['stage4_per_protein']}
        branch_data[cal_group]=stage_map
        for eval_group in ('LM-QC','QC'):
            rs=qcs_summary(stage_map,eval_group,meta,ids)
            for r in rs:
                r['calibrator']=cal_group
                r['validation_type']='same_reference_in_sample' if cal_group==eval_group else 'cross_QC_type_independent'
            summary+=rs
        label=cal_group.replace('-','_')
        for name in ('stage3_plate_scaling','stage4_per_protein'):
            write_matrix(out/f'pg_matrix_{label}_{name}.tsv',h,ann,ids,stage_map[name])
        safe_row_csv(out/f'plate_factors_{label}.tsv', ['plate','log2_shift','linear_factor','calibrator_n'],
                     [(p,f'{v:.9f}',f'{2**v:.9f}',fit['n_per_plate'][str(p)]) for p,v in fit['plate_off'].items()])
        safe_row_csv(out/f'per_protein_plate_factors_{label}.tsv',h[:1]+[f'Plate{p}_log2_residual' for p in sorted(set(plate))],
                     [(ann[j][0],)+tuple('' if not np.isfinite(v) else f'{v:.8f}' for v in fit['protein_off'][:,j]) for j in range(len(ann))])
        safe_row_csv(out/f'calibrator_intraplate_factors_{label}.tsv',['sample_id','plate','group','log2_shift','linear_factor'],
                     [(ids[i],plate[i],groups[i],f'{fit["calibrator_off"][i]:.9f}',f'{2**fit["calibrator_off"][i]:.9f}') for i in fit['calibrator_index']])
    # Fair QC 2-fold OOF: one QC replicate/plate determines factors; the
    # other QC replicate/plate is held out. Both folds cover all 32 QCs.
    print('OUT OF FOLD QC',flush=True)
    qcs=sorted([i for i,g in enumerate(groups) if g=='QC'],key=lambda i:(plate[i],ids[i]))
    fold_map={};pred3=np.full_like(log0,np.nan);pred4=np.full_like(log0,np.nan)
    for fold in range(2):
        train=[];test=[]
        for p in sorted(set(plate)):
            p_qcs=[i for i in qcs if plate[i]==p]
            assert len(p_qcs)==2, f'plate {p} expected 2 QC replicates'
            train.append(p_qcs[1-fold]);test.append(p_qcs[fold])
        x=calibrate(log1,meta,ids,'QC',train_indices=set(train),min_plates=args.min_plates)
        pred3[test]=x['stage3_plate_scaling'][test]
        pred4[test]=x['stage4_per_protein'][test]
        for i in test:fold_map[ids[i]]=fold+1
    m={'Raw':log0,'stage1_HELP_sample':log1,'OOF_stage3_plate_scaling':pred3,'OOF_stage4_per_protein':pred4}
    oo=qcs_summary(m,'QC',meta,ids)
    for r in oo:
        r['calibrator']='QC';r['validation_type']='QC_within_plate_2fold_OOF'
    summary+=oo
    safe_row_csv(out/'QC_CV_summary.tsv', ['calibrator','validation_type','evaluation_group','stage','proteins','median_cv_raw','median_cv_corrected','cv_delta_percentage_points','improved_proteins','improved_fraction'],
                 [[r.get(k,'') for k in ('calibrator','validation_type','evaluation_group','stage','proteins','median_cv_raw','median_cv_corrected','cv_delta_percentage_points','improved_proteins','improved_fraction')] for r in summary])
    # QC CV per protein, with paired detection masking and original intensities.
    per_pg=[]
    for cal_group,stages in branch_data.items():
        for g in ('LM-QC','QC'):
            ii=[i for i,s in enumerate(ids) if meta[s]['group']==g]
            cv0,n=get_cv(log0,ii)
            c1,_=get_cv(log1,ii)
            c3,_=get_cv(stages['stage3_plate_scaling'],ii)
            c4,_=get_cv(stages['stage4_per_protein'],ii)
            for j,a in enumerate(ann):
                per_pg.append([cal_group,g,a[0],a[3],int(n[j])]+[f'{v[j]:.8g}' if np.isfinite(v[j]) else '' for v in (cv0,c1,c3,c4)])
    ii=[i for i,g in enumerate(groups) if g=='QC'];cv0,n=get_cv(log0,ii);c1,_=get_cv(log1,ii);c3,_=get_cv(pred3,ii);c4,_=get_cv(pred4,ii)
    for j,a in enumerate(ann):per_pg.append(['QC_OOF','QC',a[0],a[3],int(n[j])]+[f'{v[j]:.8g}' if np.isfinite(v[j]) else '' for v in (cv0,c1,c3,c4)])
    safe_row_csv(out/'QC_CV_per_protein.tsv',['calibrator','evaluation_group','Protein.Group','Genes','detected_n','raw_cv','help_cv','plate_cv','per_protein_cv'],per_pg)
    safe_row_csv(out/'sample_HELP_factors.tsv',['sample_id','plate','group','help_peptides_detected','log2_shift','linear_factor'],
                 [(s,plate[i],groups[i],npep[i],f'{v:.9f}',f'{2**v:.9f}') for i,(s,v) in enumerate(zip(ids,sample_off))])
    safe_row_csv(out/'HELP_coverage.tsv',['peptide_key','detected_n','detected_fraction','used_for_sample_normalization'],
                 [(pep,int(np.isfinite(help_a[:,j]).sum()),f'{np.isfinite(help_a[:,j]).mean():.6f}',int(selected[j])) for j,pep in enumerate(peptides)])
    safe_row_csv(out/'sample_manifest.tsv',['sample_id','plate','group','fold_if_QC'],
                 [(s,plate[i],groups[i],fold_map.get(s,'')) for i,s in enumerate(ids)])
    write_matrix(out/'pg_matrix_stage1_HELP_sample.tsv',h,ann,ids,log1)
    write_matrix(out/'pg_matrix_QC_OOF_final.tsv',h,ann,ids,pred4)
    info={'samples':dict(Counter(groups.tolist())),'plates':len(set(plate)),'proteins':len(ann),
          'help_peptides_all':len(peptides),'help_peptides_used':int(selected.sum()),
          'PG_missing_fraction':float(np.isnan(pg).mean()),
          'HELP_missing_fraction':float(np.isnan(help_a).mean()),
          'help_exclusion_reasons':dict(bad),'min_plates':args.min_plates,
          'stage5_ANML':'NOT RUN: requires healthy-population external reference and eligible matrix',
          'step2_LM_QC':'1 calibrator per plate; calibration-only step 2 mathematically identity',
          'validation_caveat':'Cross QC type tests assume materials comparable; QC OOF is like-for-like withheld replicate.'}
    (out/'run_manifest.json').write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf8')
    print('SUMMARY',flush=True)
    for row in summary:
        if row['stage'] in ('stage4_per_protein','OOF_stage4_per_protein'):
            print(row,flush=True)
    print('DONE',out,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--metadata',type=Path,required=True)
    ap.add_argument('--pg',type=Path,required=True)
    ap.add_argument('--help-report',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--min-plates',type=int,default=8)
    produce(ap.parse_args())