"""
Step 3b — inferCNV-based malignant calling (infercnvpy), to validate the marker-based labels.

The pipeline's malignant label was marker-based (SOX2/OLIG2/EGFR). Here we infer CNVs with
infercnvpy using the non-malignant compartment (T/NK/endothelial/myeloid) as the diploid
reference, score each cell (cnv_score), and ask how well that score agrees with the marker-based
malignant label (AUROC), plus which marker lineages score highest.

A per-cohort subsample (<= PER_COHORT cells) is used for tractability.

Outputs:
  results/tables/step03b_infercnv_summary.tsv   per cohort x lineage mean cnv_score
  results/tables/step03b_infercnv_auroc.tsv     AUROC(cnv_score -> marker-malignant)
  results/figures/supplementary/step03b_infercnv_heatmap.{pdf,jpeg}

Run: env/.venv/bin/python scripts/03_infercnv/step03b_infercnv.py
"""
import os

import anndata as ad
import numpy as np
import pandas as pd

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
OBJ = os.path.join(ROOT, "objects")
TBL = os.path.join(ROOT, "results/tables")
FIG = os.path.join(ROOT, "results/figures/supplementary")
GTF = os.path.join(DATASETS, "raw/gencode/gencode.v44.basic.annotation.gtf.gz")
POS_CACHE = os.path.join(TBL, "step03b_gene_positions.tsv")
os.makedirs(FIG, exist_ok=True)

PER_COHORT = 3000
SEED = 42
REF_CATS = ["T_cell", "NK", "endothelial", "monocyte_mac", "microglia", "neutrophil"]


def save(fig, name):
    for ext in ("pdf", "jpeg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    print(f"  saved {name}", flush=True)


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import infercnvpy as cnv
    import scanpy as sc
    from sklearn.metrics import roc_auc_score

    a = ad.read_h5ad(os.path.join(OBJ, "step04_annotated.h5ad"))
    rng = np.random.default_rng(SEED)
    keep = []
    for c, grp in a.obs.groupby("cohort", observed=True):
        idx = grp.index.to_numpy()
        if len(idx) > PER_COHORT:
            idx = rng.choice(idx, PER_COHORT, replace=False)
        keep.extend(idx.tolist())
    a = a[keep].copy()
    del a.obsm  # drop inherited embeddings to save RAM
    print(f"subsample: {a.shape}; cohorts: {a.obs['cohort'].value_counts().to_dict()}", flush=True)

    # ---- gene positions (cached) ----
    if os.path.exists(POS_CACHE):
        pos = pd.read_csv(POS_CACHE, sep="\t", index_col=0)
        a.var = a.var.join(pos[["chromosome", "start", "end"]], how="left")
    else:
        cnv.io.genomic_position_from_gtf(GTF, a, inplace=True)
        a.var[["chromosome", "start", "end"]].to_csv(POS_CACHE, sep="\t")
    a = a[:, a.var["chromosome"].notna()].copy()
    a.var["chromosome"] = a.var["chromosome"].astype(str)
    print(f"genes with position: {a.shape[1]}", flush=True)

    # ---- infer CNV ----
    cnv.tl.infercnv(a, reference_key="lineage", reference_cat=REF_CATS,
                    window_size=100, step=10, n_jobs=8)
    cnv.tl.pca(a)
    try:
        cnv.pp.neighbors(a)
    except Exception:
        sc.pp.neighbors(a, use_rep="X_cnv")
    cnv.tl.leiden(a, key_added="cnv_leiden", resolution=1.0, flavor="igraph", n_iterations=2)
    cnv.tl.cnv_score(a)
    print("cnv_score computed", flush=True)

    # ---- agreement with marker-based malignant label ----
    y = (a.obs["lineage"] == "malignant").astype(int).values
    auc = roc_auc_score(y, a.obs["cnv_score"].values)
    print(f"AUROC(cnv_score vs marker-malignant) overall = {auc:.3f}", flush=True)
    auroc_rows = [{"group": "all", "n": int(a.shape[0]), "n_malignant": int(y.sum()),
                   "auroc": round(float(auc), 3)}]
    for c in sorted(a.obs["cohort"].unique()):
        m = (a.obs["cohort"] == c).values
        if len(set(y[m])) == 2:
            auroc_rows.append({"group": c, "n": int(m.sum()), "n_malignant": int(y[m].sum()),
                               "auroc": round(float(roc_auc_score(y[m], a.obs["cnv_score"].values[m])), 3)})

    summ = (a.obs.groupby(["cohort", "lineage"], observed=True)["cnv_score"]
            .agg(["count", "mean", "median"]).reset_index())
    summ.to_csv(os.path.join(TBL, "step03b_infercnv_summary.tsv"), sep="\t", index=False)
    pd.DataFrame(auroc_rows).to_csv(os.path.join(TBL, "step03b_infercnv_auroc.tsv"), sep="\t", index=False)

    # ---- figure: infercnvpy's chromosome-level heatmap (returns dict of Axes) ----
    try:
        res = cnv.pl.chromosome_heatmap(a, groupby="lineage", show=False, figsize=(13, 4))
        if isinstance(res, dict):
            fig = next(iter(res.values())).figure
        else:
            fig = getattr(res, "figure", res)
        save(fig, "step03b_infercnv_heatmap")
    except Exception as e:  # noqa: BLE001
        print(f"  [skip figure] {type(e).__name__}: {e}", flush=True)

    print("\n=== AUROC(cnv_score -> marker-malignant) ===", flush=True)
    print(pd.DataFrame(auroc_rows).to_string(index=False), flush=True)
    print("\n=== mean cnv_score by lineage ===", flush=True)
    print(summ.sort_values("mean", ascending=False).to_string(index=False), flush=True)
    print("DONE", flush=True)
