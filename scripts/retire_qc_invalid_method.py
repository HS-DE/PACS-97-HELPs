#!/usr/bin/env python3
"""One-time safe retirement of a QC-invalid normalization experiment.

Delete only the rejected experiment outputs, entrypoint and workflow; move
shared numerical routines into pacs_common.py so the retained methods still
run. Keep all original raw inputs and the other experiment results.
"""
from pathlib import Path
import json
import re
import shutil

ROOT=Path(__file__).resolve().parents[1]

def read(path):
    return (ROOT/path).read_text(encoding="utf-8")

def write(path,value):
    target=ROOT/path
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(value,encoding="utf-8")

def replace(path,old,new):
    content=read(path)
    if old not in content:
        raise ValueError(f"Expected patch anchor not found in {path}: {old[:100]!r}")
    write(path,content.replace(old,new))

def replace_section(path,begin,end,replacement):
    content=read(path)
    a=content.index(begin)
    b=content.index(end,a)
    if b<=a:raise ValueError(f"Invalid section in {path}")
    write(path,content[:a]+replacement+content[b:])

def disable_push(path):
    content=read(path)
    new,n=re.subn(r"(?ms)^on:\n  push:\n.*?^  workflow_dispatch:\n",
                  "on:\n  workflow_dispatch:\n",content,count=1)
    if n!=1:raise ValueError(f"Push stanza unexpected in {path}")
    write(path,new)

def do_cleanup():
    source=read("scripts/pacs_trial1.py")
    a,b=source.index("def load_ratio("),source.index("def write_factors(")
    c,d=source.index("def finite_positive("),source.index("def plot_factors(")
    common='''#!/usr/bin/env python3
"""Shared HELP input reader, annotation constants and protein CV utilities.

This module has no normalization entrypoint; used by retained PACS methods.
"""
import math
import zipfile
from xml.etree import ElementTree as ET
import numpy as np
from filter_samples import readbook, tab
from yeast_protein_qc import COLORS, GROUPS, ROOT

PG = ROOT / "pg_matrix.tsv"
HELP = ROOT / "LH_ratio_tables_2026.9.22—ALL-HELP.xlsx"
ANN = ["Protein.Group", "Protein.Ids", "Protein.Names",
       "Genes", "First.Protein.Description"]
CV_MIN_DETECT_FRAC = 0.70
CV_MIN_N = 3
CV_MIN_MEAN = 1e-8

'''+source[a:b].strip()+"\n\n\n"+source[c:d].strip()+"\n"
    write("scripts/pacs_common.py",common)
    replace("scripts/pacs_trial1b.py","from pacs_trial1 import (","from pacs_common import (")
    for path in ("scripts/pacs_trial2.py","scripts/pacs_trial2_stepwise.py",
                 "scripts/pacs_eno1_protein_trial.py"):
        replace(path,"from pacs_trial1 import ","from pacs_common import ")

    b="scripts/pacs_trial1b.py"
    replace(b,
        "imputation is documented. Original files and Trial 1A are left untouched.",
        "imputation is documented. Original source files remain unchanged.")
    replace(b,"# same validated 97-HELP reader as Trial 1A",
              "# validated HELP reader; selected-peptide runs can override it")
    replace(b,
        '        "This exploratory test uses the **same division direction** as Trial 1A,",\n'
        '        "but does NOT correct any Internal-QC (LM) sample. Trial 1A and source files are unchanged.",',
        '        "This exploratory correction uses a pooled-QC-referenced HELP factor,",\n'
        '        "without correcting any Internal-QC (LM) sample. Source files are unchanged.",')
    replace(b,
        "- No original pg_matrix.tsv, HELP workbook, metadata or Trial 1A result was edited.",
        "- Original pg_matrix.tsv, HELP workbook and metadata were not edited.")
    replace("scripts/pacs_trial2.py",
        "Trial 1A / 1B and all original input files remain untouched.",
        "All original input files and previously published QC-safe results remain untouched.")

    p="scripts/pacs_help5_reruns.py"
    replace(p,"Reruns: Trial 1A, Trial 1B, Trial 2 stepwise yeast-peptide a/b, and",
              "Reruns: QC-safe HELP-only division, Trial 2 yeast-peptide a/b, and")
    replace(p,"import pacs_trial1 as t1","import pacs_common as common")
    replace(p,'METHODS = ("trial1a", "trial1b", "trial2_peptide", "trial3_eno1")',
              'METHODS = ("trial1b", "trial2_peptide", "trial3_eno1")')
    replace(p,"t1.ANN","common.ANN")
    replace_section(p,"def run_1a():","def run_1b():","")
    replace(p,
        '    for name in ("trial1a", "trial1b"):\n'
        '        fs = pd.read_csv(OUT / name / "normalization_factors.csv")\n'
        '        if fs.n_valid_HELP_ratios.min() < 1 or fs.n_valid_HELP_ratios.max() > 5:\n'
        '            raise AssertionError(f"Incorrect five HELP count in {name}")\n'
        '        for _, r in fs.iterrows():\n'
        '            if abs(float(r["sample_mR"]) - help_cache["raw_medians"][r["sample_id"]]) > 1e-9:\n'
        '                raise AssertionError(f"Incorrect five HELP raw median in {name}")',
        '    fs = pd.read_csv(OUT / "trial1b" / "normalization_factors.csv")\n'
        '    if fs.n_valid_HELP_ratios.min() < 1 or fs.n_valid_HELP_ratios.max() > 5:\n'
        '        raise AssertionError("Unexpected five-HELP count in HELP-only correction")\n'
        '    for _, r in fs.iterrows():\n'
        '        if abs(float(r["sample_mR"]) - help_cache["raw_medians"][r["sample_id"]]) > 1e-9:\n'
        '            raise AssertionError("Unexpected five-HELP raw median in HELP-only correction")')
    replace(p,"four types of experiments","three types of experiments")
    replace(p,"Four historical algorithms re-executed","Three retained algorithms re-executed")
    replace(p,
        '        "| [trial1a](trial1a/README.md) | median(raw L/H), pooled-QC ref; divide protein by sample/ref factor; includes multiplied-direction sensitivity | All groups corrected, including QC |",\n'
        '        "| [trial1b](trial1b/README.md) | Same as 1A, two PCA views | Internal-QC unchanged, pooled QC corrected |",',
        '        "| [trial1b](trial1b/README.md) | HELP-only median(raw L/H), pooled-QC reference; two PCA views | Internal-QC unchanged, pooled QC corrected |",')
    replace(p,
        '        "Do not equate raw-ratio medians used in Trial 1A/1B with median(log2 ratio)",\n'
        '        "used in Trial2/3; the old procedures explicitly differed in this respect.",',
        '        "HELP-only correction uses median(raw L/H), whereas two-step Trial2/3",\n'
        '        "methods use median(log2 L/H); these are distinct calculations.",')
    replace(p,'"description": "Trial1A/1B: raw L/H median. Trial2/3: log2 L/H median.",',
              '"description": "QC-safe HELP-only: raw L/H median. Trial2/3: log2 L/H median.",')
    replace(p,"    run_1a()\n","")
    replace(p,"1A 1B Trial2 peptide and Trial3 ENO1",
              "HELP-only Trial2 peptide and Trial3 ENO1")

    p="results/README.md"
    replace(p,"four types of experiments","three types of experiments")
    replace(p,"Four historical algorithms re-executed","Three retained algorithms re-executed")
    replace(p,
        "| [trial1a](trial1a/README.md) | median(raw L/H), pooled-QC ref; divide protein by sample/ref factor; includes multiplied-direction sensitivity | All groups corrected, including QC |\n"
        "| [trial1b](trial1b/README.md) | Same as 1A, two PCA views | Internal-QC unchanged, pooled QC corrected |",
        "| [trial1b](trial1b/README.md) | HELP-only median(raw L/H), pooled-QC reference; two PCA views | Internal-QC unchanged, pooled QC corrected |")
    replace(p,
        "Do not equate raw-ratio medians used in Trial 1A/1B with median(log2 ratio)\n"
        "used in Trial2/3; the old procedures explicitly differed in this respect.",
        "HELP-only correction uses median(raw L/H), whereas two-step Trial2/3\n"
        "methods use median(log2 L/H); these are distinct calculations.")

    manifest=json.loads(read("results/run_manifest.json"))
    manifest["methods"]=[x for x in manifest["methods"] if x!="trial1a"]
    manifest["description"]="QC-safe HELP-only: raw L/H median. Trial2/3: log2 L/H median."
    write("results/run_manifest.json",json.dumps(manifest,ensure_ascii=False,indent=2)+"\n")
    for p in ("results/trial1b/README.md","results_v1/pacs_trial1b/README.md"):
        replace(p,
            "This exploratory test uses the **same division direction** as Trial 1A,\n"
            "but does NOT correct any Internal-QC (LM) sample. Trial 1A and source files are unchanged.",
            "This exploratory HELP-only correction uses a pooled-QC-referenced L/H ratio\n"
            "to scale pooled QC, HC and S. Internal-QC (LM) remains exactly unchanged.")
        replace(p,"- No original pg_matrix.tsv, HELP workbook, metadata or Trial 1A result was edited.",
                  "- Original pg_matrix.tsv, HELP workbook and metadata were not edited.")
    replace("results_v1/pacs_trial2/README.md",
        "Trial 1A / 1B and all original input files remain untouched.",
        "All original input files and other retained results remain untouched.")

    workflows=(
        ".github/workflows/pacs-help5-reruns.yml",
        ".github/workflows/pacs-trial1b.yml",
        ".github/workflows/pacs-trial2.yml",
        ".github/workflows/pacs-trial2-stepwise.yml",
        ".github/workflows/pacs-eno1-protein-trial.yml",
    )
    for p in workflows:disable_push(p)
    p=".github/workflows/pacs-help5-reruns.yml"
    replace(p,"rerun-four-historical-methods:","rerun-three-retained-methods:")
    replace(p,'for trial in ("trial1a","trial1b"):', 'for trial in ("trial1b",):')
    replace(p,
        '              if trial=="trial1b":\n'
        '                  for sid in help.loc[help.group=="Internal-QC","sample_id"]:\n'
        '                      assert original[sid].equals(corrected[sid]),f"Internal-QC unexpectedly modified {sid}"',
        '              for sid in help.loc[help.group=="Internal-QC","sample_id"]:\n'
        '                  assert original[sid].equals(corrected[sid]),f"Internal-QC unexpectedly modified {sid}"')
    replace(p,"PASS: all four methods, exactly five HELP peptides",
              "PASS: three retained methods, exactly five HELP peptides")
    replace(p,'          assert (base/"run_manifest.json").exists()\n',
              '          assert (base/"run_manifest.json").exists()\n'
              '          assert not (base/"trial1a").exists()\n')
    # Stop future regeneration; preserve only workflows whose QC logic is valid.
    for path in ("results/trial1a","results_v1/pacs_trial1"):
        target=ROOT/path
        if not target.is_dir():raise ValueError(f"Expected obsolete folder: {target}")
        shutil.rmtree(target)
    for path in ("scripts/pacs_trial1.py",".github/workflows/pacs-trial1.yml"):
        target=ROOT/path
        if not target.is_file():raise ValueError(f"Expected obsolete file: {target}")
        target.unlink()
    for path in ("results/trial1b","results/trial2_peptide","results/trial3_eno1",
                 "results_v1/pacs_trial1b","results_v1/pacs_trial2"):
        if not (ROOT/path).is_dir():
            raise AssertionError(f"Retained result directory unexpectedly missing: {path}")
    for path in ("results/README.md","scripts/pacs_help5_reruns.py",
                 "scripts/pacs_trial1b.py","results/trial1b/README.md"):
        content=read(path)
        if re.search(r"trial\s*1a|trial1a|pacs_trial1(?!b)",content,re.I):
            raise AssertionError(f"Obsolete references remain in {path}")
    print("CLEANED: invalid method outputs/entrypoint/old workflow removed; valid methods retained.")

if __name__=="__main__":
    do_cleanup()
