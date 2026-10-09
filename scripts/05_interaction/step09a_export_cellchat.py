"""
Step 5a — Export a labeled, downsampled matrix for CellChat (R).

From step04_annotated.h5ad (log-normalized, lineage labels), downsample per lineage
and write:
  results/cellchat/expr.mtx     (genes x cells, MatrixMarket)
  results/cellchat/genes.txt    (gene symbols)
  results/cellchat/cells.txt    (cell ids)
  results/cellchat/labels.csv   (cell -> lineage)

Run: env/.venv/bin/python scripts/05_interaction/step09a_export_cellchat.py
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io as sio

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
OBJ = os.path.join(ROOT, "objects")
OUT = os.path.join(ROOT, "results/cellchat")
os.makedirs(OUT, exist_ok=True)
PER = 1500

if __name__ == "__main__":
    a = ad.read_h5ad(os.path.join(OBJ, "step04_annotated.h5ad"))
    print(f"loaded annotated: {a.shape}", flush=True)
    keep_lineages = ["malignant", "monocyte_mac", "microglia", "T_cell", "endothelial", "NK", "neutrophil"]
    idx = []
    rng = np.random.default_rng(42)
    for lin in keep_lineages:
        cells = np.where(a.obs["lineage"].values == lin)[0]
        if len(cells) > PER:
            cells = rng.choice(cells, PER, replace=False)
        idx.append(cells)
        print(f"  {lin}: {len(cells)}", flush=True)
    idx = np.concatenate(idx)
    sub = a[idx].copy()
    X = sub.X.toarray().T if sp.issparse(sub.X) else np.asarray(sub.X).T  # genes x cells
    sio.mmwrite(os.path.join(OUT, "expr.mtx"), sp.csr_matrix(X))
    pd.Series(sub.var_names).to_csv(os.path.join(OUT, "genes.txt"), index=False, header=False)
    pd.Series(sub.obs_names).to_csv(os.path.join(OUT, "cells.txt"), index=False, header=False)
    sub.obs["lineage"].to_csv(os.path.join(OUT, "labels.csv"), index=True, header=["lineage"])
    print(f"exported {X.shape} (genes x cells) to {OUT}", flush=True)
    print(f"lineage counts: {sub.obs['lineage'].value_counts().to_dict()}", flush=True)
    print("DONE", flush=True)
