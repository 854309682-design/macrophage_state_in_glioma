"""
Step 2 — Build AnnData objects for the three scRNA-seq discovery/validation cohorts.

Follows the pipeline in analysis_strategy_1.md §4 Step 1. Produces per-cohort AnnData
(genes x cells) saved to objects/, printing dims + a quick sanity view BEFORE saving
(inspect-then-save).

Cohorts / formats (all VERIFIED this session):
  - GSE103224 (Yuan 2018) : headerless gene x cell UMI .txt.gz
        col1=Ensembl(versioned), col2=symbol, col3+=UMI. No cell barcode -> name cells
        "<sample>_<i>".
  - GSE131928 (Neftel 2019): processed TPM, gene x cell, header GENE + barcodes
        (Smart-seq2 + 10x), extracted from GSE131928_RAW.tar.
  - GSE163120 (Pombo Antunes): 10x-style CSV gene x cell (feature col1 = gene id/symbol,
        header = barcodes), + an annot csv (barcode -> cluster/sample). Human GBM matrices
        extracted from GSE163120_RAW.tar.

Run:  env/.venv/bin/python scripts/00_download/step02_build_anndata.py
"""
import gzip
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
DATASETS = os.environ.get("GLIOMA_DATA") or os.path.join(os.path.dirname(ROOT), "Datasets")
OBJ = os.path.join(ROOT, "objects")
os.makedirs(OBJ, exist_ok=True)

sc.settings.verbosity = 2


def _save(adata: ad.AnnData, name: str):
    # outer-join concat leaves missing genes as NaN -> replace with 0 (absent = 0 counts)
    if np.isnan(adata.X).any():
        adata.X = np.nan_to_num(adata.X)
    out = os.path.join(OBJ, f"step02_{name}.h5ad")
    adata.write_h5ad(out)
    print(f"  saved -> {out}")


def read_gse103224():
    """Headerless gene x cell UMI matrix (col1 Ensembl, col2 symbol, col3+ counts)."""
    path = os.path.join(DATASETS, "GSE103224")
    files = sorted(
        f for f in os.listdir(path) if f.endswith(".filtered.matrix.txt.gz")
    )
    var = None
    mats = []
    obs = []
    for f in files:
        sample = f.split(".")[0].split("_")[1]
        full = os.path.join(path, f)
        df = pd.read_csv(full, sep="\t", header=None, compression="gzip")
        if var is None:
            var = pd.DataFrame({"ensembl": df[0].astype(str),
                                "symbol": df[1].astype(str)})
            var.index = var["ensembl"]
            var = var[["symbol"]]  # keep ensembl as index, symbol as column
        mat = df.iloc[:, 2:].values.T.astype(np.float32)  # cells x genes
        mats.append(mat)
        obs += [f"{sample}_{i}" for i in range(mat.shape[0])]
    X = np.vstack(mats)
    adata = ad.AnnData(X=X, var=var)
    adata.obs_names = obs
    adata.obs["sample"] = [o.split("_")[0] for o in obs]
    adata.obs["dataset"] = "GSE103224"
    print(f"GSE103224 combined: {adata.shape}")
    _save(adata, "gse103224")
    return adata


def read_gse131928():
    """Processed TPM, gene x cell (header GENE + barcodes), Smart-seq2 + 10x."""
    path = os.path.join(DATASETS, "GSE131928")
    adatas = []
    for f in sorted(os.listdir(path)):
        if not f.endswith(".tsv.gz"):
            continue
        df = pd.read_csv(os.path.join(path, f), sep="\t", index_col=0, compression="gzip")
        a = ad.AnnData(X=df.values.T.astype(np.float32), var=pd.DataFrame(index=df.index))
        a.obs_names = [str(c) for c in df.columns]
        a.obs["dataset"] = "GSE131928"
        a.obs["platform"] = "Smartseq2" if "Smartseq2" in f else "10X"
        print(f"  GSE131928 {os.path.basename(f)}: {a.shape}")
        adatas.append(a)
    adata = ad.concat(adatas, join="outer")
    adata.var_names_make_unique()
    print(f"GSE131928 combined: {adata.shape}")
    _save(adata, "gse131928")
    return adata


def read_gse163120():
    """Human GBM matrices (10x-style CSV) + annot (barcode -> cluster/sample)."""
    path = os.path.join(DATASETS, "GSE163120")
    adatas = []
    for f in sorted(os.listdir(path)):
        if ".filtered.gene.bc.matrix.csv.gz" not in f or "Mouse" in f:
            continue  # human RNA matrices only (skip mouse) here
        # annot file for this GSM
        gsm = f.split("_")[0]
        annot = next(
            (x for x in os.listdir(path) if x.startswith(gsm) and "annot" in x), None
        )
        df = pd.read_csv(os.path.join(path, f), index_col=0, compression="gzip")
        a = ad.AnnData(X=df.values.T.astype(np.float32),
                       var=pd.DataFrame(index=[str(i) for i in df.index]))
        a.obs_names = [str(c) for c in df.columns]
        a.obs["dataset"] = "GSE163120"
        label = "recurrent" if ".R" in f else ("newly-diagnosed" if ".ND" in f else gsm)
        a.obs["disease_phase"] = label
        if annot:
            ann = pd.read_csv(os.path.join(path, annot), compression="gzip")
            ann["cell"] = ann["cell"].astype(str)
            ann = ann.set_index("cell").reindex(a.obs_names)
            a.obs["cluster"] = ann["cluster"].values
            a.obs["sample"] = ann["sample"].values
            a.obs["tam_label"] = ann["ident"].values
        print(f"  GSE163120 {gsm} ({label}): {a.shape}")
        adatas.append(a)
    adata = ad.concat(adatas, join="outer")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()  # same 10x barcodes recur across R/ND matrices
    print(f"GSE163120 combined: {adata.shape}")
    _save(adata, "gse163120")
    return adata


if __name__ == "__main__":
    print("=== GSE103224 ===")
    read_gse103224()
    print("=== GSE131928 ===")
    read_gse131928()
    print("=== GSE163120 (human) ===")
    read_gse163120()
    print("DONE")
