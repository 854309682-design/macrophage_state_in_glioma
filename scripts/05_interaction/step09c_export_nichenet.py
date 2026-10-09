"""
Step 5c — Export sender/receiver pseudobulk for NicheNet (Step 5, ligand-target layer).

Sender   = malignant cells in the hypoxic niche (top-quartile 15-gene hypoxia score, per cohort).
Receiver = TAM-IS myeloid cells (state == "TAM-IS", Step 4 labels).

Writes to results/nichenet/:
  step09c_sender_expr.csv       genes x samples ("all" pooled + one column per cohort)
  step09c_receiver_expr.csv     genes x samples ("all" pooled + one column per cohort)
  step09c_sender_detect.csv     genes x samples, fraction of cells with expression > 0
  step09c_receiver_detect.csv   genes x samples, fraction of cells with expression > 0
  step09c_geneset_tamis.txt     TAM-IS gene set (curated signature + top Step 4 markers)
  step09c_cells_summary.tsv     cell counts per group

Run: env/.venv/bin/python scripts/05_interaction/step09c_export_nichenet.py
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
TBL = os.path.join(ROOT, "results/tables")
OUT = os.path.join(ROOT, "results/nichenet")
os.makedirs(OUT, exist_ok=True)

sc.settings.verbosity = 0

# 15-gene hypoxia signature and TAM-IS signature — identical to H1 (step06) and H2 (step06b).
HYPOXIA = ["VEGFA", "SLC2A1", "LDHA", "PGK1", "PDK1", "HK2", "BNIP3", "NDRG1", "ENO1",
           "PDK3", "TPI1", "P4HA1", "ADM", "MRPS17", "ANKZF1"]
TAMIS = ["SPP1", "TREM2", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "LGALS3", "CD163",
         "CTSL", "CCL3", "CCL4", "GPNMB", "LPL", "HLA-DRA", "HLA-DRB1", "S100A9",
         "HMOX1", "FTL", "FTH1"]
# housekeeping / high-abundance genes excluded from the marker-derived gene set
HK_BLOCK = {"B2M", "MALAT1", "GAPDH", "ACTB", "TMSB4X", "EEF1A1", "TPT1", "VIM",
            "ANXA2", "LGALS1", "S100A11", "HSP90AA1", "HSPA1A", "BASP1", "FTL", "FTH1"}


def dense(sub):
    x = sub.X
    return np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x)


def pseudobulk(a, groups):
    """Pooled ("all") and per-level mean expression + detection fraction (genes as rows)."""
    x = dense(a)
    levels = ["all"] + sorted(a.obs[groups].unique())
    mean_rows, det_rows, counts = [], [], []
    for g in levels:
        m = np.ones(a.n_obs, dtype=bool) if g == "all" else (a.obs[groups] == g).values
        mean_rows.append(x[m].mean(axis=0))
        det_rows.append((x[m] > 0).mean(axis=0))
        counts.append(int(m.sum()))
    means = pd.DataFrame(np.vstack(mean_rows).T, index=a.var_names, columns=levels)
    detect = pd.DataFrame(np.vstack(det_rows).T, index=a.var_names, columns=levels)
    return means, detect, pd.Series(counts, index=levels)


def tamis_geneset(universe):
    mk = pd.read_csv(os.path.join(TBL, "step05_tamis_markers.tsv"), sep="\t", index_col=0)
    keep = mk[(mk["logfoldchanges"] >= 1) & (mk["pvals_adj"] < 0.01)]
    keep = keep[~keep.index.str.startswith(("RPL", "RPS", "MT-"))]
    keep = keep[~keep.index.isin(HK_BLOCK)]
    top = keep.sort_values("scores", ascending=False).head(150).index.tolist()
    return sorted((set(TAMIS) | set(top)) & set(universe))


if __name__ == "__main__":
    # --- receiver: TAM-IS myeloid cells ---
    state = ad.read_h5ad(os.path.join(OBJ, "step05_myeloid_reclustered.h5ad"), backed="r")
    tamis_cells = set(state.obs_names[state.obs["state"] == "TAM-IS"])
    state.file.close()
    print(f"TAM-IS cells (labels): {len(tamis_cells)}", flush=True)

    mye = ad.read_h5ad(os.path.join(OBJ, "step04_myeloid.h5ad"), backed="r")
    rec = mye[mye.obs_names.isin(tamis_cells)].to_memory()
    mye.file.close()
    print(f"receiver (full genes): {rec.shape}", flush=True)
    rec_means, rec_detect, rec_counts = pseudobulk(rec, "cohort")
    universe = list(rec.var_names)
    del rec

    # --- sender: malignant cells, hypoxic-high (top quartile per cohort) ---
    ann = ad.read_h5ad(os.path.join(OBJ, "step04_annotated.h5ad"), backed="r")
    sen = ann[(ann.obs["lineage"] == "malignant").values].to_memory()
    ann.file.close()
    print(f"sender (all malignant): {sen.shape}", flush=True)

    hyp = [g for g in HYPOXIA if g in sen.var_names]
    sc.tl.score_genes(sen, hyp, score_name="hypoxia", use_raw=False, random_state=0)
    thr = sen.obs.groupby("cohort", observed=True)["hypoxia"].transform(lambda s: s.quantile(0.75))
    sen = sen[(sen.obs["hypoxia"] >= thr).values].copy()
    print(f"sender (hypoxic-high): {sen.shape}; hypoxia genes used: {len(hyp)}", flush=True)
    sen_means, sen_detect, sen_counts = pseudobulk(sen, "cohort")
    del sen

    # --- align gene space and write ---
    genes = [g for g in universe if g in sen_means.index]
    sen_means, sen_detect = sen_means.loc[genes], sen_detect.loc[genes]
    rec_means, rec_detect = rec_means.loc[genes], rec_detect.loc[genes]

    sen_means.to_csv(os.path.join(OUT, "step09c_sender_expr.csv"))
    rec_means.to_csv(os.path.join(OUT, "step09c_receiver_expr.csv"))
    sen_detect.to_csv(os.path.join(OUT, "step09c_sender_detect.csv"))
    rec_detect.to_csv(os.path.join(OUT, "step09c_receiver_detect.csv"))

    geneset = tamis_geneset(genes)
    pd.Series(geneset).to_csv(os.path.join(OUT, "step09c_geneset_tamis.txt"),
                              index=False, header=False)

    summary = pd.DataFrame({
        "role": ["sender"] * sen_means.shape[1] + ["receiver"] * rec_means.shape[1],
        "group": list(sen_means.columns) + list(rec_means.columns),
        "n_cells": list(sen_counts.reindex(sen_means.columns).values)
        + list(rec_counts.reindex(rec_means.columns).values),
    })
    summary.to_csv(os.path.join(OUT, "step09c_cells_summary.tsv"), sep="\t", index=False)

    print(f"expressed genes: {len(genes)}; TAM-IS gene set: {len(geneset)}", flush=True)
    print(summary.to_string(index=False), flush=True)
    print("DONE", flush=True)
