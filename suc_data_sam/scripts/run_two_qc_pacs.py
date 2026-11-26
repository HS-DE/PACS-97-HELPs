#!/usr/bin/env python3
"""Run the supplied PACS workflow twice: LM-QC or QC is the only pooled QC."""
from __future__ import annotations
import argparse
import importlib.util
import tempfile
from pathlib import Path

import pandas as pd


def read_metadata(path: Path) -> pd.DataFrame:
    try:
        meta = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    except UnicodeDecodeError:
        # The supplied CSV had a stray non-UTF8 byte in an unnamed trailing column.
        meta = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="latin1")
    meta.columns = [str(c).strip().casefold() for c in meta.columns]
    meta = meta.loc[:, [c for c in meta.columns if not c.startswith("unnamed:")]].copy()
    required = ["sample_id", "group", "condition", "batch"]
    missing = [c for c in required if c not in meta.columns]
    if missing:
        raise ValueError(f"metadata missing required columns: {missing}")
    if "plate" not in meta.columns:
        meta["plate"] = meta["batch"]
    meta["sample_id"] = meta.sample_id.astype(str).str.strip()
    if meta.sample_id.eq("").any() or not meta.sample_id.is_unique:
        raise ValueError("sample_id must be non-empty and unique")
    return meta


def map_qc_strategy(meta: pd.DataFrame, selected_qc: str) -> pd.DataFrame:
    """Only selected_qc becomes pooled QC; the other QC types are excluded."""
    result = meta.copy()
    mapped = []
    for value in result.group.astype(str).str.strip():
        key = value.casefold()
        if key == selected_qc.casefold():
            mapped.append("qc")
        elif key == "sample":
            mapped.append("sample")
        elif key == "blank":
            mapped.append("blank")
        else:
            mapped.append("internal-calibrator")
    result["group"] = mapped
    return result[["sample_id", "group", "condition", "batch", "plate"]]


def load_pacs_module(script_dir: Path):
    path = script_dir / "PACS_norm.py"
    spec = importlib.util.spec_from_file_location("PACS_norm", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load PACS script: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-table", type=Path, required=True)
    parser.add_argument("--help-table", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    pacs = load_pacs_module(script_dir)
    meta = read_metadata(args.metadata)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    for selected in ("LM-QC", "QC"):
        mapped = map_qc_strategy(meta, selected)
        out = args.out_dir / selected
        with tempfile.TemporaryDirectory(prefix=f"metadata_{selected}_") as temp:
            meta_path = Path(temp) / "metadata.csv"
            mapped.to_csv(meta_path, index=False, encoding="utf-8-sig")
            pacs.run(
                args.protein_table,
                args.help_table,
                meta_path,
                out,
                plate_col="plate",
                min_samples=2,
            )
        print(f"Completed {selected}: {out}")


if __name__ == "__main__":
    main()
