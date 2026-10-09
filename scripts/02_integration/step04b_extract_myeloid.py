"""
Step 3c — Annotate the integrated object and extract the myeloid compartment.

Loads objects/step04_integrated.h5ad (scVI latent, UMAP, Leiden), scores established
lineage marker sets, annotates each Leiden cluster by its dominant lineage, then
subsets to myeloid (microglia-like + monocyte-derived macrophages) and saves
objects/step04_myeloid.h5ad. Returns the full annotated object too.

Markers (glioma TME; strategy §4 Step 3-4 + literature):
  - microglia(MG-like): P2RY12, TMEM119, CX3CR1, SALL1, TREM2
  - monocyte-macrophage(Mo-Mac/MDM): CD14, LYZ, SPP1, APOE, LGALS3, CD68, C1QA/C1QB
  - T cell: CD3D, CD3E, CD8A, IL7R, TRAC
  - NK: NKG7, GNLY, PRF1
  - endothelial: PECAM1, VWF, CLDN5
  - malignant-ish (GBM): GFAP, SOX2, OLIG1, OLIG2, SOX4
  - neutrophil: FCGR3B, CSF3R, S100A8, S100A9

Run: env/.venv/bin/python scripts/02_integration/step04b_extract_myeloid.py
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

MARKERS = {
    "microglia": ["P2RY12", "TMEM119", "CX3CR1", "SALL1", "TREM2"],
    "monocyte_mac": ["CD14", "LYZ", "SPP1", "APOE", "LGALS3", "CD68", "C1QA", "C1QB"],
    "T_cell": ["CD3D", "CD3E", "CD8A", "IL7R", "TRAC"],
    "NK": ["NKG7", "GNLY", "PRF1"],
    "endothelial": ["PECAM1", "VWF", "CLDN5"],
    "malignant": ["GFAP", "SOX2", "OLIG1", "OLIG2"],
    "neutrophil": ["FCGR3B", "CSF3R", "S100A8", "S100A9"],
    "prolif": ["MKI67", "TOP2A"],
}


def score_all(a):
    scores = {}
    for lab, genes in MARKERS.items():
        genes = [g for g in genes if g in a.var_names]
        if not genes:
            print(f"  [warn] no genes for {lab}", flush=True)
            a.obs[f"score_{lab}"] = 0.0
            continue
        sc.tl.score_genes(a, genes, score_name=f"score_{lab}", use_raw=False)
        scores[lab] = genes
    return scores


if __name__ == "__main__":
    # full common-gene object for scoring (markers present); integrated for leiden/cohort
    a = ad.read_h5ad(os.path.join(OBJ, "step04_common.h5ad"))
    print(f"loaded common: {a.shape}", flush=True)
    integ = ad.read_h5ad(os.path.join(OBJ, "step04_integrated.h5ad"))
    a.obs["leiden"] = integ.obs["leiden"].values
    a.obs["cohort"] = integ.obs["cohort"].values
    print(f"common {a.shape} vs integrated {integ.shape} (rows align)", flush=True)
    score_all(a)

    # for each Leiden cluster -> mean lineage score
    cluster_scores = a.obs.groupby("leiden")[
        [f"score_{k}" for k in MARKERS]].mean()
    # annotate cluster by argmax over the four "solid" lineages
    lineages = ["microglia", "monocyte_mac", "T_cell", "NK", "endothelial",
                "malignant", "neutrophil"]
    cluster_scores["annotation"] = cluster_scores[[f"score_{k}" for k in lineages]].idxmax(
        axis=1).str.replace("score_", "")
    print("cluster -> lineage:", flush=True)
    print(cluster_scores["annotation"].to_string(), flush=True)

    cn_map = cluster_scores["annotation"].to_dict()
    a.obs["lineage"] = a.obs["leiden"].map(cn_map)

    # myeloid = microglia + monocyte_mac
    myel = a[a.obs["lineage"].isin(["microglia", "monocyte_mac"])].copy()
    print(f"myeloid cells: {myel.n_obs} ({myel.n_obs / a.n_obs * 100:.1f}%)", flush=True)

    # save annotated + myeloid subsets
    a.write_h5ad(os.path.join(OBJ, "step04_annotated.h5ad"))
    myel.write_h5ad(os.path.join(OBJ, "step04_myeloid.h5ad"))
    print("saved step04_annotated.h5ad + step04_myeloid.h5ad", flush=True)

    print("lineage distribution:", flush=True)
    print(a.obs["lineage"].value_counts().to_string(), flush=True)
    print("DONE", flush=True)
