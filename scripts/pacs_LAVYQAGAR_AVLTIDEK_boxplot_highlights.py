#!/usr/bin/env python3
"""PACS original 75-sample raw HELP Heavy / matched PG boxplots with highlights."""
from __future__ import annotations
import csv
import hashlib
import json
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from filter_samples import readbook, tab
from yeast_protein_qc import GROUPS, COLORS
import pacs_two_reference_plate_trial as pacs

OUT=pacs.ROOT/"results"/"pacs_APOE_SERPINA1_two_HELP_boxplot_highlights"
TARGETS={
    "LAVYQAGAR":{"accession":"P02649","gene":"APOE","color":"#C54D47","marker":"D"},
    "AVLTIDEK":{"accession":"P01009","gene":"SERPINA1","color":"#087D94","marker":"o"},
}
PALETTE=dict(zip(GROUPS,COLORS))

def read_heavy(samples):
    strings,sheets=readbook(pacs.HELP)
    if len(sheets)!=3:raise AssertionError("Expected L/H/ratio_LH")
    rows=tab(sheets[1],strings)
    hdr={str(v).strip():k for k,v in rows[0].items() if k!="A"}
    if len(rows)!=98 or set(hdr)!=set(samples):
        raise ValueError("HELP H original sheet does not match 75 samples and 97 peptides")
    seq=[str(row.get("A","")).strip() for row in rows[1:]]
    if len(seq)!=len(set(seq)) or len(seq)!=97:raise AssertionError("HELP sequences are not unique")
    H=np.full((97,75),np.nan)
    for k,row in enumerate(rows[1:]):
        for j,sid in enumerate(samples):
            try:
                v=float(row.get(hdr[sid],"nan"))
                if np.isfinite(v) and v>0:H[k,j]=np.log2(v)
            except (ValueError,TypeError):pass
    for seq_id in TARGETS:
        if seq_id not in seq or not np.isfinite(H[seq.index(seq_id)]).all():
            raise ValueError(f"Heavy missing in at least one sample: {seq_id}")
    return seq,H

def pick_pg(pg,log_pg,seq_id,acc,gene):
    def tokenmatch(s):return acc in re.split(r"[;,|\s]+",str(s).strip())
    matches=[i for i in range(len(pg))
             if tokenmatch(pg.iloc[i]["Protein.Group"]) or tokenmatch(pg.iloc[i]["Protein.Ids"])]
    if len(matches)!=1:
        raise ValueError(f"Ambiguous PG source for {seq_id}/{acc}: {len(matches)} rows, indices {matches}")
    i=matches[0]
    row=pg.iloc[i]
    genenames=re.split(r"[;,|\s]+",str(row["Genes"]).strip())
    if gene not in genenames:
        raise ValueError(f"PG accession {acc} genes {row['Genes']} do not include {gene}")
    return i,{"HELP_sequence":seq_id,"UniProt_accession":acc,"Gene":gene,
              "source_pg_row_number_including_header":i+2,
              **{key:str(row[key]) for key in pacs.ANN},
              "n_positive_pg":int(np.isfinite(log_pg[i]).sum())}

def save_csv(path,records):
    with path.open("w",newline="",encoding="utf-8-sig") as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

def group_boundaries(ax,order,samples,info):
    for k in range(1,75):
        prev=info[samples[order[k-1]]];cur=info[samples[order[k]]]
        if cur["group"]!=prev["group"]:
            ax.axvline(k+.5,color="#546879",lw=1.25,alpha=.65)
        elif cur["plate"]!=prev["plate"]:
            ax.axvline(k+.5,color="#8795A5",lw=.85,ls=":",alpha=.65)
    for group in GROUPS:
        spots=[k+1 for k,j in enumerate(order) if info[samples[j]]["group"]==group]
        ax.text((spots[0]+spots[-1])/2,1.015,f"{group} (n={len(spots)})",
                fontsize=10.5,weight="bold",color=PALETTE[group],
                transform=ax.get_xaxis_transform(),ha="center",va="bottom",
                clip_on=False)
    ax.grid(axis="y",alpha=.15)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

def boxplot(source,markers,kind,order,samples,info,title,out):
    fig,ax=plt.subplots(figsize=(27,9.4))
    fig.subplots_adjust(left=.055,right=.99,top=.82,bottom=.33)
    pos=np.arange(1,76)
    vectors=[source[:,j][np.isfinite(source[:,j])] for j in order]
    if min(map(len,vectors))<2:raise ValueError(f"Empty raw boxplot {kind}")
    bp=ax.boxplot(vectors,positions=pos,widths=.65,patch_artist=True,
                  showfliers=False,manage_ticks=False,
                  medianprops={"color":"#203549","linewidth":1.5},
                  boxprops={"edgecolor":"#8795A4","linewidth":.7})
    for patch,j in zip(bp["boxes"],order):
        patch.set_facecolor(PALETTE[info[samples[j]]["group"]])
        patch.set_alpha(.48)
    group_boundaries(ax,order,samples,info)
    target_handles=[]
    for k,(seq,values) in enumerate(markers.items()):
        prop=TARGETS[seq]
        offsets=(-.15 if k==0 else .15) if len(markers)==2 else 0.
        y=values[order]
        good=np.isfinite(y)
        ax.scatter(pos[good]+offsets,y[good],
                   s=52,c=prop["color"],marker=prop["marker"],
                   edgecolor="white",lw=.55,zorder=7)
        label=(seq+" Heavy" if kind=="Heavy" else
               f"{prop['accession']} {prop['gene']} (HELP {seq})")
        target_handles.append(Line2D([0],[0],marker=prop["marker"],linestyle="None",
           color=prop["color"],markeredgecolor="white",markersize=9,label=label))
    ax.set_xlim(.35,75.65)
    ax.set_ylabel("Original log2 "+("HELP Heavy intensity" if kind=="Heavy"
                                   else "protein-group intensity"),fontsize=12)
    ax.set_xticks(pos,[samples[j] for j in order],rotation=90,fontsize=6.3)
    ax.set_title(title,loc="left",fontsize=14,fontweight="bold",pad=40)
    # Explicitly reserve figure bottom; no legends are inside plotting axes.
    fig.legend(handles=[Patch(color=PALETTE[g],alpha=.55,label=g) for g in GROUPS],
               loc="lower center",bbox_to_anchor=(.5,.16),ncol=4,
               frameon=False,fontsize=10)
    fig.legend(handles=target_handles,loc="lower center",bbox_to_anchor=(.5,.115),
               ncol=len(target_handles),frameon=False,fontsize=10)
    fig.savefig(out,dpi=185,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    if out.stat().st_size<23000:raise AssertionError(f"Bad plot size {out}")

def side_by_side(seqs,H,logPG,refs,order,samples,info,out):
    fig,axes=plt.subplots(2,1,figsize=(27,11.5),sharex=True)
    fig.subplots_adjust(left=.06,right=.985,top=.89,bottom=.22,hspace=.32)
    pos=np.arange(1,76)
    legend=[]
    for ax,seq in zip(axes,TARGETS):
        prop=TARGETS[seq]
        rawH=H[seqs.index(seq),order]
        rawPG=logPG[refs[seq],order]
        ax.plot(pos,rawH-np.nanmedian(rawH),c=prop["color"],marker=".",
                markersize=5,lw=1.35,label=seq+" Heavy (centered)")
        ax.plot(pos,rawPG-np.nanmedian(rawPG),c="#294356",marker="o",
                markersize=2.9,lw=1.3,label=prop["gene"]+" PG (centered)")
        ax.axhline(0,c="#A5B2BF",lw=.8)
        ax.set_ylabel("log2 deviation from own median")
        ax.set_title(f"{seq} (Heavy) vs {prop['gene']} / {prop['accession']} (PG)",
                     loc="left",pad=32,fontweight="bold",fontsize=12)
        group_boundaries(ax,order,samples,info)
        legend.extend(ax.get_legend_handles_labels()[0])
    axes[-1].set_xticks(pos,[samples[j] for j in order],rotation=90,fontsize=6)
    fig.legend(handles=legend,loc="lower center",bbox_to_anchor=(.5,.045),
               ncol=4,frameon=False)
    fig.suptitle("Matched raw HELP Heavy and endogenous PG intensity trends",
                 fontsize=15,fontweight="bold")
    fig.savefig(out,dpi=180,bbox_inches="tight",facecolor="white")
    plt.close(fig)

def main():
    if OUT.exists():raise FileExistsError(f"Independent result exists, no overwrite: {OUT}")
    pg,samples,info,raw,_,_=pacs.load()
    if raw.shape!=(6675,75):raise AssertionError(f"PG shape {raw.shape}")
    seqs,H=read_heavy(samples)
    ordered=np.array(sorted(range(75),key=lambda j:(
        GROUPS.index(info[samples[j]]["group"]),
        int(info[samples[j]]["plate"]),samples[j])),dtype=int)
    refs={};matches=[]
    for seq,prop in TARGETS.items():
        i,record=pick_pg(pg,raw,seq,prop["accession"],prop["gene"])
        refs[seq]=i;matches.append(record)
    full=int(np.isfinite(H).all(axis=1).sum())
    if full!=96:raise AssertionError(f"Unexpected all-75 Heavy coverage: {full}")
    OUT.mkdir(parents=True)
    save_csv(OUT/"target_PG_mapping.csv",matches)
    audit=[]
    for k,j in enumerate(ordered,1):
        sid=samples[j]
        row={"order":k,"sample_id":sid,"group":info[sid]["group"],
             "plate":info[sid]["plate"],"n_detected_PG":int(np.isfinite(raw[:,j]).sum())}
        for seq,prop in TARGETS.items():
            row[seq+"_Heavy_log2"]=float(H[seqs.index(seq),j])
            z=raw[refs[seq],j]
            row[prop["gene"]+"_PG_log2"]=float(z) if np.isfinite(z) else ""
        audit.append(row)
    save_csv(OUT/"75_samples_target_intensities.csv",audit)
    heavy={seq:H[seqs.index(seq)] for seq in TARGETS}
    proteins={seq:raw[refs[seq]] for seq in TARGETS}
    images=[
      ("01_97HELP_Heavy_both_highlights.png",H,heavy,"Heavy",
       "Original 97 HELP Heavy distributions: LAVYQAGAR and AVLTIDEK highlighted"),
      ("02_6675PG_APOE_SERPINA1_both_highlights.png",raw,proteins,"PG",
       "Original 6675 PG distributions: APOE/P02649 and SERPINA1/P01009 highlighted"),
    ]
    for seq,prop in TARGETS.items():
        images.extend([
         (f"03_Heavy_{seq}_highlight.png",H,{seq:heavy[seq]},"Heavy",
          f"Original 97 HELP Heavy distributions: highlight {seq}"),
         (f"04_PG_{prop['gene']}_{prop['accession']}_highlight.png",raw,
          {seq:proteins[seq]},"PG",
          f"Original 6675 PG distributions: highlight {prop['gene']} / {prop['accession']}"),
        ])
    for filename,matrix,highlights,kind,title in images:
        boxplot(matrix,highlights,kind,ordered,samples,info,title,OUT/filename)
    side_by_side(seqs,H,raw,refs,ordered,samples,info,
                 OUT/"05_two_HELP_and_PG_centered_trends.png")
    manifest={
      "target_map":TARGETS,"actual_PG_matching_rows":matches,
      "n_samples":75,"n_raw_PG":6675,"n_raw_HELP":97,
      "n_75_complete_Heavy_HELP":full,
      "original_sources_sha256":{
          x.name:hashlib.sha256(x.read_bytes()).hexdigest()
          for x in (pacs.PG,pacs.HELP,pacs.META_QC,pacs.META_STUDY)
      },"sample_order":"group -> plate -> sample ID",
      "H_sheet":"original H sheet, not L or ratio",
      "PG":"original pg_matrix.tsv, not corrected",
      "log2":True,"no_normalization":True,"no_imputation":True,
      "legends":"outside axes in reserved figure margin",
      "images":[r[0] for r in images]+["05_two_HELP_and_PG_centered_trends.png"]
    }
    (OUT/"manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n")
    readme=[
      "# PACS 75-sample HELP Heavy and matched protein highlights","",
      "Mappings supplied by project owner: LAVYQAGAR -> P02649 / APOE;",
      "AVLTIDEK -> P01009 / SERPINA1.","",
      "Original source H-sheet peptide values are plotted against 97 HELP Heavy",
      "per sample; original matched protein-group values are plotted against",
      "6675 PG per sample. Positive intensities transformed to log2.",
      "No normalization; missing values not imputed.",
      "All 75 samples use group -> plate -> sample ID ordering.",
      "Target protein row identity is verified via exact accession tokens and Gene.",
      "Legends are beneath axes in dedicated blank figure space.","",
      "## Matched protein rows","",
      "| HELP | UniProt | Gene | PG source row incl header | Positive PG samples |",
      "| --- | --- | --- | ---: | ---: |",
      *[f"| {r['HELP_sequence']} | {r['UniProt_accession']} | {r['Gene']} | "
        f"{r['source_pg_row_number_including_header']} | {r['n_positive_pg']} |"
        for r in matches],
      "","## Highlight charts","",
    ]
    for filename,*_ in images:
        readme.extend(["![highlight]("+filename+")",""])
    readme.extend([
       "![relative trends](05_two_HELP_and_PG_centered_trends.png)","",
       "## Tables","",
       "- target_PG_mapping.csv: exact source protein rows and annotations.",
       "- 75_samples_target_intensities.csv: ordered sample-level highlight values.",
       "- manifest.json: source hashes, mapping, plotting details."
    ])
    (OUT/"README.md").write_text("\n".join(readme)+"\n")
    print("MATCHED_PG",json.dumps(matches,ensure_ascii=False))
    print("PASS 75 samples, verified exact protein mappings, 7 figures, no modified input files")

if __name__=="__main__":
    main()
