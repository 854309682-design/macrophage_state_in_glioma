"""
Step 5f — Export the stratified myeloid subsample as a loom for pySCENIC (GRN step).

Stratified subsample of the myeloid compartment (<= PER_CAP cells per state x cohort,
<= TOTAL_CAP cells), log-normalized X, written as a SCOPE/pySCENIC-compatible loom
(rows = genes, columns = cells; row_attrs Gene, col_attrs CellID).

Writes:
  results/scenic/step09f_myeloid.loom     genes x cells (float32)
  results/scenic/step09f_cellmeta.csv     CellID, state, cohort, leiden
  results/scenic/step09f_subsample.tsv    cell counts per state x cohort

Run: env/.venv/bin/python scripts/05_interaction/step09f_export_scenic.py
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
import loompy as lp

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
OBJ = os.path.join(ROOT, "objects")
OUT = os.path.join(ROOT, "results/scenic")
os.makedirs(OUT, exist_ok=True)

PER_CAP = 2000
TOTAL_CAP = 15000
SEED = 42


def dense(sub):
    x = sub.X
    return np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)


if __name__ == "__main__":
    state = ad.read_h5ad(os.path.join(OBJ, "step05_myeloid_reclustered.h5ad"), backed="r")
    lab = state.obs["state"].copy()
    state.file.close()

    mye = ad.read_h5ad(os.path.join(OBJ, "step04_myeloid.h5ad"), backed="r")
    obs = mye.obs[["cohort", "leiden"]].copy()
    obs["state"] = lab.reindex(obs.index).values
    obs = obs.dropna(subset=["state"])
    print(f"myeloid cells with state: {obs.shape[0]}; states: {obs['state'].value_counts().to_dict()}", flush=True)

    rng = np.random.default_rng(SEED)
    picked = []
    for (st, co), grp in obs.groupby(["state", "cohort"], observed=True):
        ii = grp.index.to_numpy()
        if len(ii) > PER_CAP:
            ii = rng.choice(ii, PER_CAP, replace=False)
        picked.extend(ii.tolist())
        print(f"  {st:14s} {co:10s} -> {len(ii)}", flush=True)
    if len(picked) > TOTAL_CAP:
        picked = rng.choice(np.array(picked), TOTAL_CAP, replace=False).tolist()
    picked = list(dict.fromkeys(picked))
    print(f"subsample: {len(picked)} cells", flush=True)

    sub = mye[mye.obs_names.isin(set(picked))].to_memory()[picked].copy()
    mye.file.close()
    x = dense(sub).astype(np.float32)
    print(f"matrix cells x genes: {x.shape}", flush=True)

    lp.create(
        os.path.join(OUT, "step09f_myeloid.loom"),
        x.T,
        row_attrs={"Gene": np.array(sub.var_names, dtype=str)},
        col_attrs={"CellID": np.array(sub.obs_names, dtype=str)},
    )

    meta = obs.loc[picked, ["state", "cohort", "leiden"]].copy()
    meta.index.name = "CellID"
    meta.to_csv(os.path.join(OUT, "step09f_cellmeta.csv"))

    obs.groupby(["state", "cohort"], observed=True).size().rename("n_cells").reset_index() \
        .to_csv(os.path.join(OUT, "step09f_subsample.tsv"), sep="\t", index=False)

    print("state x cohort in subsample:", flush=True)
    print(meta.groupby(["state", "cohort"], observed=True).size().to_string(), flush=True)
    print("DONE", flush=True)
