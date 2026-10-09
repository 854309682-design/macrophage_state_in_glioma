"""
Step 3a — scVI integration (correct gene-space strategy).

Earlier attempt intersected per-cohort top-3000 HVG lists -> only 54 common genes
(too few; technologies barely overlap). Correct approach:
  1. Convert each cohort to gene symbols.
  2. Build the COMMON gene universe (genes present in all 3 cohorts, by symbol).
  3. Subset each cohort to the common universe, concatenate.
  4. Compute HVGs on the concatenated common-gene object (top N_HVG).
  5. scVI (batch=cohort), UMAP + Leiden. Save objects/step04_integrated.h5ad.

Run: env/.venv/bin/python scripts/02_integration/step04_integrate.py
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
OBJ = os.path.join(ROOT, "objects")
N_HVG = 3000

sc.settings.verbosity = 2
sc.settings.n_jobs = 8

COHORTS = ["gse103224", "gse131928", "gse163120"]


def load_symbols(name):
    a = ad.read_h5ad(os.path.join(OBJ, f"step03_{name}.h5ad"))
    print(f"  {name}: {a.shape}", flush=True)
    if "symbol" in a.var.columns:
        sym = a.var["symbol"].astype(str).values
        a.var = a.var.drop(columns=["symbol"])
        a.var_names = sym
    a.var_names_make_unique()
    a.X = np.nan_to_num(a.X)
    return a


if __name__ == "__main__":
    adatas = {c: load_symbols(c) for c in COHORTS}
    # common gene universe (present in all cohorts, by symbol)
    common = sorted(set.intersection(*[set(a.var_names) for a in adatas.values()]))
    print(f"common gene universe: {len(common)}", flush=True)
    if len(common) < 2000:
        print("WARNING: small common universe; using it as-is", flush=True)

    concats = []
    for c in COHORTS:
        sub = adatas[c][:, adatas[c].var_names.isin(common)].copy()
        sub.obs["cohort"] = c
        concats.append(sub)
    adata = ad.concat(concats, join="inner")
    adata.var_names_make_unique()
    adata.X = np.nan_to_num(adata.X)
    adata.layers["counts"] = adata.X.copy()
    print(f"concatenated (common universe): {adata.shape}", flush=True)
    adata.write_h5ad(os.path.join(OBJ, "step04_common.h5ad"))  # full common-gene set for annotation

    # HVG on the common-gene concat
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat",
                                layer="counts", subset=False)
    adata = adata[:, adata.var["highly_variable"]].copy()
    adata.layers["counts"] = adata.X.copy()
    print(f"after HVG (on common universe): {adata.shape}", flush=True)

    scvi.model.SCVI.setup_anndata(adata, batch_key="cohort", layer="counts")
    model = scvi.model.SCVI(adata, n_latent=30, n_layers=2)
    model.train(max_epochs=30, early_stopping=True, plan_kwargs={"lr": 1e-3})
    adata.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.8, flavor="igraph", n_iterations=2)

    adata.write_h5ad(os.path.join(OBJ, "step04_integrated.h5ad"))
    print(f"saved -> objects/step04_integrated.h5ad ({adata.shape})", flush=True)
    print("DONE", flush=True)
