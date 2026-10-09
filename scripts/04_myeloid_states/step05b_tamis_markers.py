"""
Step 4b (marker) — Define the TAM-IS marker signature.

Loads the re-clustered myeloid object (scores present). Defines:
  TAM-IS cells  = subclusters with high TAM-IS signature (tamis_minus_micro > 1.0)
  Microglia     = subclusters with negative tamis_minus_micro (< 0)
Then runs Wilcoxon rank_genes_groups (TAM-IS vs Microglia) to derive TAM-IS markers.
Writes objects/step05_myeloid_reclustered.h5ad (with tamis_label) + a marker table.
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

sc.settings.verbosity = 1

TAMIS_THRESH = 1.0


if __name__ == "__main__":
    a = ad.read_h5ad(os.path.join(OBJ, "step05_myeloid_reclustered.h5ad"))
    cl = a.obs.groupby("leiden")[["tamis_score", "micro_score"]].mean()
    cl["diff"] = cl["tamis_score"] - cl["micro_score"]
    tn = set(cl[cl["diff"] > TAMIS_THRESH].index)  # TAM-IS candidates
    mn = set(cl[cl["diff"] < 0].index)             # microglia-like
    print(f"TAM-IS clusters: {sorted(tn)} | microglia clusters: {sorted(mn)}", flush=True)

    a.obs["state"] = "myeloid_other"
    a.obs.loc[a.obs["leiden"].astype(str).isin(tn), "state"] = "TAM-IS"
    a.obs.loc[a.obs["leiden"].astype(str).isin(mn), "state"] = "microglia"
    print(a.obs["state"].value_counts().to_string(), flush=True)

    # Wilcoxon markers: TAM-IS vs microglia
    sub = a[a.obs["state"].isin(["TAM-IS", "microglia"])].copy()
    sc.tl.rank_genes_groups(sub, groupby="state", groups=["TAM-IS"],
                            reference="microglia", method="wilcoxon", use_raw=False)
    res = sc.get.rank_genes_groups_df(sub, group="TAM-IS")
    res = res[res["logfoldchanges"] > 0.5]
    res = res[res["pvals_adj"] < 0.05]
    # exclude the scoring genes themselves for a clean "marker" definition
    res = res.sort_values("scores", ascending=False)
    print(f"TAM-IS markers (log2FC>0.5, adjP<0.05): {len(res)}", flush=True)
    print(res.head(30).to_string(), flush=True)

    res.to_csv(os.path.join(TBL, "step05_tamis_markers.tsv"), sep="\t", index=False)
    a.write_h5ad(os.path.join(OBJ, "step05_myeloid_reclustered.h5ad"))
    print("saved step05_tamis_markers.tsv" + " + updated reclustered object", flush=True)
    print("DONE", flush=True)
