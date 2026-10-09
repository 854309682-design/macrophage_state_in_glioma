"""
Step 4a — TAM-IS state discovery (first pass).

Loads the integrated object (has X_scVI/X_umap) + the full-gene myeloid object
(step04_myeloid.h5ad), subsets to myeloid, re-clusters on the scVI latent, scores
immunosuppressive TAM-IS vs microglia-like signatures (full genes), and flags the
myeloid subclusters that look like TAM-IS (high TAM-IS signature, low microglia).

Signatures (strategy §4 Step 4 + literature):
  TAM-IS (Mo-Mac / immunosuppressive): SPP1, TREM2, APOE, LGALS3, CD163, C1QA, C1QB,
                                          LYZ, FTL, FTH1, IL1B, CCL3
  microglia-like (MG): P2RY12, TMEM119, CX3CR1, SALL1, C3, CSF1R

Outputs: objects/step05_myeloid_reclustered.h5ad + results/tables/step05_tamis_candidates.tsv

Run: env/.venv/bin/python scripts/04_myeloid_states/step05_tamis_discovery.py
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
OBJ = os.path.join(ROOT, "objects")

sc.settings.verbosity = 1

TAMIS_GENES = ["SPP1", "TREM2", "APOE", "LGALS3", "CD163", "C1QA", "C1QB",
               "LYZ", "FTL", "FTH1", "IL1B", "CCL3", "TGFB1", "GPNMB"]
MICRO_GENES = ["P2RY12", "TMEM119", "CX3CR1", "SALL1", "C3", "CSF1R"]


if __name__ == "__main__":
    intg = ad.read_h5ad(os.path.join(OBJ, "step04_integrated.h5ad"))
    myel_full = ad.read_h5ad(os.path.join(OBJ, "step04_myeloid.h5ad"))
    print(f"integrated {intg.shape}; myeloid-full {myel_full.shape}", flush=True)

    # subset integrated to myeloid cells (match by obs_names)
    myo = intg[intg.obs_names.isin(set(myel_full.obs_names))].copy()
    print(f"myeloid in integrated: {myo.shape}", flush=True)

    # re-cluster myeloid on scVI latent
    sc.pp.neighbors(myo, use_rep="X_scVI", n_neighbors=15)
    sc.tl.umap(myo)
    sc.tl.leiden(myo, resolution=1.0, flavor="igraph", n_iterations=2)
    print(f"myeloid subclusters: {myo.obs['leiden'].nunique()}", flush=True)

    # full-gene myeloid object aligned to myo cell order
    f = myel_full[myel_full.obs_names.isin(set(myo.obs_names))].copy()
    f = f[list(myo.obs_names)].copy()
    print(f"aligned full-gene myeloid: {f.shape}", flush=True)

    def score(a, genes):
        g = [x for x in genes if x in a.var_names]
        return g

    tg = score(f, TAMIS_GENES)
    mg = score(f, MICRO_GENES)
    print(f"scoring TAM-IS ({len(tg)}) microglia ({len(mg)})", flush=True)
    sc.tl.score_genes(f, tg, score_name="tamis_score", use_raw=False)
    sc.tl.score_genes(f, mg, score_name="micro_score", use_raw=False)
    myo.obs["tamis_score"] = f.obs["tamis_score"].values
    myo.obs["micro_score"] = f.obs["micro_score"].values

    # per-cluster mean
    cl = myo.obs.groupby("leiden")[["tamis_score", "micro_score"]].mean()
    cl["tamis_minus_micro"] = cl["tamis_score"] - cl["micro_score"]
    cl["n_cells"] = myo.obs["leiden"].value_counts().reindex(cl.index)
    cl = cl.sort_values("tamis_minus_micro", ascending=False)
    cl["is_tamis"] = (cl["tamis_minus_micro"] > 0) & (cl["n_cells"] >= 100)
    print("\nmyeloid subcluster TAM-IS score table:", flush=True)
    print(cl.round(3).to_string(), flush=True)

    myo.write_h5ad(os.path.join(OBJ, "step05_myeloid_reclustered.h5ad"))
    os.makedirs(os.path.join(ROOT, "results/tables"), exist_ok=True)
    cl.to_csv(os.path.join(ROOT, "results/tables/step05_tamis_candidates.tsv"), sep="\t")
    print("saved step05_myeloid_reclustered.h5ad + step05_tamis_candidates.tsv", flush=True)
    print("DONE", flush=True)
