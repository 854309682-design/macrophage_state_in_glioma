"""
Step 1b — Build the two extra discovery cohorts for H1 replication.

GSE89567 (Venteicher 2017) and GSE70630 (Tirosh 2016) are processed
gene(symbol) x cell TSV matrices with a header (genes as index). Both are built to
objects/step02_<cohort>.h5ad (cells x genes), matching the other cohorts.

Run: env/.venv/bin/python scripts/00_download/step02b_build_extra_cohorts.py
"""
import os

import anndata as ad
import numpy as np
import pandas as pd

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
DATASETS = os.environ.get("GLIOMA_DATA") or os.path.join(os.path.dirname(ROOT), "Datasets")
OBJ = os.path.join(ROOT, "objects")

FILES = {
    "gse89567": "GSE89567_IDH_A_processed_data.txt.gz",
    "gse70630": "GSE70630_OG_processed_data_v2.txt.gz",
}


if __name__ == "__main__":
    for name, fn in FILES.items():
        df = pd.read_csv(os.path.join(DATASETS, fn), sep="\t", index_col=0, compression="gzip")
        print(f"  {name}: raw gene x cell {df.shape}", flush=True)
        X = df.values.T.astype(np.float32)  # cells x genes
        X = np.nan_to_num(X)
        a = ad.AnnData(X=X, var=pd.DataFrame(
            index=[str(g).strip("'\"") for g in df.index]))
        a.obs_names = [str(c) for c in df.columns]
        a.obs["cohort"] = name
        a.obs["dataset"] = name
        a.var_names_make_unique()
        a.obs_names_make_unique()
        out = os.path.join(OBJ, f"step02_{name}.h5ad")
        a.write_h5ad(out)
        print(f"  saved {out} ({a.shape})", flush=True)
    print("DONE", flush=True)
