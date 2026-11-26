#!/usr/bin/env python3
"""Exploratory HELP latent-factor diagnostics; not a validated normalization."""
from __future__ import annotations
import argparse
import importlib.util
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler


def load_pacs_module(script_dir: Path):
    spec = importlib.util.spec_from_file_location("PACS_norm", script_dir / "PACS_norm.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_metadata(path: Path) -> pd.DataFrame:
    try:
        meta = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    except UnicodeDecodeError:
        meta = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="latin1")
    meta.columns = [str(c).strip().casefold() for c in meta.columns]
    meta = meta.loc[:, [c for c in meta.columns if not c.startswith("unnamed:")]].copy()
    if "plate" not in meta.columns:
        meta["plate"] = meta["batch"]
    return meta


def mapped_metadata(meta: pd.DataFrame, selected_qc: str) -> pd.DataFrame:
    out = meta.copy()
    groups = []
    for value in out.group.astype(str).str.strip():
        key = value.casefold()
        if key == selected_qc.casefold():
            groups.append("qc")
        elif key == "sample":
            groups.append("sample")
        elif key == "blank":
            groups.append("blank")
        else:
            groups.append("internal-calibrator")
    out["group"] = groups
    return out[["sample_id", "group", "condition", "batch", "plate"]]


def cv_by_protein(log_values: np.ndarray) -> np.ndarray:
    observed = np.isfinite(log_values)
    linear = np.exp2(log_values)
    n = observed.sum(axis=1)
    mean = np.divide(
        np.where(observed, linear, 0).sum(axis=1), n,
        out=np.full(len(n), np.nan), where=n > 0
    )
    ss = np.where(observed, linear - mean[:, None], 0) ** 2
    var = np.divide(
        ss.sum(axis=1), n - 1,
        out=np.full(len(n), np.nan), where=n > 1
    )
    cv = 100 * np.sqrt(var) / mean
    cv[n < 3] = np.nan
    return cv


def grouped_global_factor_cv(H: np.ndarray, y: np.ndarray, groups: np.ndarray, k: int):
    splitter = GroupKFold(n_splits=min(8, len(np.unique(groups))))
    pred = np.full(len(y), np.nan)
    baseline = np.full(len(y), np.nan)
    for train, test in splitter.split(H, y, groups):
        scaler = StandardScaler()
        x_train = scaler.fit_transform(H[train])
        x_test = scaler.transform(H[test])
        pca = PCA(n_components=k, svd_solver="full")
        scores_train = pca.fit_transform(x_train)
        scores_test = pca.transform(x_test)
        model = LinearRegression().fit(scores_train, y[train])
        pred[test] = model.predict(scores_test)
        baseline[test] = np.mean(y[train])
    return pred, baseline


def proteinwise_oof(H: np.ndarray, Y: np.ndarray, groups: np.ndarray, k: int = 1,
                    min_train_obs: int = 6):
    """QC-only per-protein regressions, evaluated with plate-grouped OOF predictions."""
    splitter = GroupKFold(n_splits=min(8, len(np.unique(groups))))
    out = np.full_like(Y, np.nan)
    for train, test in splitter.split(H, groups=groups):
        scaler = StandardScaler()
        z_train = scaler.fit_transform(H[train])
        z_test = scaler.transform(H[test])
        pca = PCA(n_components=k, svd_solver="full")
        s_train = pca.fit_transform(z_train)
        s_test = pca.transform(z_test)
        x_train = np.column_stack([np.ones(len(train)), s_train])
        x_test = np.column_stack([np.ones(len(test)), s_test])
        for protein in range(Y.shape[0]):
            values = Y[protein, train]
            valid = np.isfinite(values)
            if valid.sum() < max(min_train_obs, k + 2):
                continue
            beta = np.linalg.lstsq(x_train[valid], values[valid], rcond=None)[0]
            predicted = x_test @ beta
            reference = np.median(values[valid])
            out[protein, test] = Y[protein, test] - (predicted - reference)
    out[~np.isfinite(Y)] = np.nan
    return out


def run_one(pacs, protein, help_report, metadata, selected_qc, out_dir):
    mapped = mapped_metadata(metadata, selected_qc)
    with tempfile.TemporaryDirectory(prefix=f"metadata_{selected_qc}_") as temp:
        meta_path = Path(temp) / "metadata.csv"
        mapped.to_csv(meta_path, index=False, encoding="utf-8-sig")
        _, meta, raw, heavy, sequences, _ = pacs.load_inputs(
            protein, help_report, meta_path, plate_col="plate"
        )

    meta = meta.reset_index(drop=True)
    qc = meta.group.astype(str).str.casefold().eq("qc").to_numpy()
    qc_indices = np.flatnonzero(qc)
    plate_groups = meta.loc[qc, "plate"].astype(str).to_numpy()

    # Complete HELP features only; no missing-value imputation.
    H_all = heavy.T
    complete = np.isfinite(H_all).all(axis=0)
    H = H_all[:, complete]
    help_names = np.asarray(sequences)[complete]
    y_all = np.nanmedian(raw, axis=0)
    y_qc = y_all[qc] - np.median(y_all[qc])
    H_qc = H[qc]

    out_dir.mkdir(parents=True, exist_ok=True)
    score_rows = []
    prediction_rows = []
    for k in (1, 2, 3):
        pred, baseline = grouped_global_factor_cv(H_qc, y_qc, plate_groups, k)
        score_rows.append({
            "strategy": selected_qc, "n_QC": len(qc_indices),
            "n_HELP_features": H.shape[1], "n_components": k,
            "grouped_CV_RMSE": np.sqrt(mean_squared_error(y_qc, pred)),
            "baseline_RMSE": np.sqrt(mean_squared_error(y_qc, baseline)),
            "grouped_CV_R2": r2_score(y_qc, pred),
            "RMSE_improvement_pct": 100 * (
                1 - np.sqrt(mean_squared_error(y_qc, pred)) /
                np.sqrt(mean_squared_error(y_qc, baseline))
            )
        })
        if k == 1:
            for i, idx in enumerate(qc_indices):
                prediction_rows.append({
                    "sample_id": meta.iloc[idx].sample_id,
                    "plate": meta.iloc[idx].plate,
                    "observed_PG_global_shift_log2": y_qc[i],
                    "predicted_PG_global_shift_log2_OOF": pred[i],
                    "baseline_prediction_log2": baseline[i],
                })

    # PCA fit on all selected QC samples; scores are diagnostic, not a correction.
    scaler = StandardScaler().fit(H_qc)
    Z_qc = scaler.transform(H_qc)
    pca = PCA(n_components=min(5, len(qc_indices) - 1, H.shape[1]), svd_solver="full")
    pca.fit(Z_qc)
    scores = pca.transform(scaler.transform(H))
    score_table = meta.copy()
    for j in range(scores.shape[1]):
        score_table[f"PC{j+1}"] = scores[:, j]
    score_table["is_selected_qc"] = qc
    help_global = np.median(H, axis=1)
    score_table["HELP_global_shift_log2"] = help_global - np.median(help_global[qc])
    score_table["PG_global_shift_log2"] = y_all - np.median(y_all[qc])
    score_table.to_csv(out_dir / "HELP_latent_scores.tsv", sep="\t", index=False)
    pd.DataFrame(score_rows).to_csv(out_dir / "HELP_global_factor_grouped_CV.tsv", sep="\t", index=False)
    pd.DataFrame(prediction_rows).to_csv(out_dir / "HELP_global_factor_OOF_predictions.tsv", sep="\t", index=False)

    loadings = pd.DataFrame({"HELP": help_names})
    for j, component in enumerate(pca.components_):
        loadings[f"PC{j+1}_loading"] = component
    loadings.to_csv(out_dir / "HELP_PCA_loadings.tsv", sep="\t", index=False)

    # Global-shift correlation is diagnostic only.
    r, p = pearsonr(score_table.loc[qc, "HELP_global_shift_log2"],
                    score_table.loc[qc, "PG_global_shift_log2"])
    pd.DataFrame([{
        "strategy": selected_qc, "n_QC": len(qc_indices),
        "HELP_PG_global_shift_Pearson_r": r, "p_value": p
    }]).to_csv(out_dir / "HELP_global_shift_correlation.tsv", sep="\t", index=False)

    # Protein-wise factor regression is evaluated OOF; no final matrix is promoted.
    proteinwise = proteinwise_oof(H_qc, raw[:, qc], plate_groups, k=1)
    cv_raw = cv_by_protein(raw[:, qc])
    cv_factor = cv_by_protein(proteinwise)
    valid = np.isfinite(cv_raw) & np.isfinite(cv_factor)
    pd.DataFrame([{
        "strategy": selected_qc, "n_components": 1,
        "n_proteins": int(valid.sum()),
        "raw_median_CV_pct": float(np.median(cv_raw[valid])),
        "HELP_factor_proteinwise_OOF_median_CV_pct": float(np.median(cv_factor[valid])),
        "proteins_improved_n": int((cv_factor[valid] < cv_raw[valid]).sum()),
    }]).to_csv(out_dir / "HELP_factor_proteinwise_OOF_validation.tsv", sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-table", type=Path, required=True)
    parser.add_argument("--help-table", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    pacs = load_pacs_module(script_dir)
    metadata = read_metadata(args.metadata)
    for selected in ("LM-QC", "QC"):
        run_one(
            pacs, args.protein_table, args.help_table, metadata, selected,
            args.out_dir / selected
        )


if __name__ == "__main__":
    main()
