"""
H1 rigor — scArches reference mapping.

Reference = myeloid cells from the discovery cohorts (GSE103224 + GSE131928).
Query     = myeloid cells from the independent cohort (GSE163120, Pombo Antunes).
Train scVI on the reference, map the query via scArches surgery (load_query_data),
then kNN-transfer the reference TAM-IS label to query cells in the shared latent and
quantify TAM-IS recovery in the held-out cohort.

Outputs: results/tables/step11_arches_h1.tsv, objects/step11_query_mapped.h5ad

Run: env/.venv/bin/python scripts/02_integration/step11_arches_h1.py
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from sklearn.neighbors import NearestNeighbors

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
OBJ = os.path.join(ROOT, "objects")
TBL = os.path.join(ROOT, "results/tables")

TAMIS = ["SPP1", "TREM2", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "LGALS3", "CD163",
         "CTSL", "CCL3", "CCL4", "GPNMB", "LPL", "HLA-DRA", "HLA-DRB1", "S100A9",
         "HMOX1", "FTL", "FTH1"]
MICRO = ["P2RY12", "TMEM119", "CX3CR1", "SALL1", "CRYBB1", "CST3", "PROS1", "MERTK"]

REF_COHORTS = ["gse103224", "gse131928"]
QUERY_COHORT = "gse163120"

sc.settings.verbosity = 1


def score(a, genes, name):
    g = [x for x in genes if x in a.var_names]
    sc.tl.score_genes(a, g, score_name=name, use_raw=False)
    return len(g)


if __name__ == "__main__":
    integ = ad.read_h5ad(os.path.join(OBJ, "step04_integrated.h5ad"))
    myo_names = set(ad.read_h5ad(os.path.join(OBJ, "step05_myeloid_reclustered.h5ad"), backed="r").obs_names)
    myo = integ[integ.obs_names.isin(myo_names)].copy()
    print(f"myeloid: {myo.shape}", flush=True)

    ref = myo[myo.obs["cohort"].isin(REF_COHORTS)].copy()
    qry = myo[myo.obs["cohort"] == QUERY_COHORT].copy()
    print(f"reference {ref.shape} | query {qry.shape}", flush=True)

    # reference TAM-IS label
    ng = score(ref, TAMIS, "tamis")
    score(ref, MICRO, "micro")
    thr = ref.obs["tamis"].quantile(0.66)
    ref.obs["is_tamis"] = (ref.obs["tamis"] > thr).astype(int)
    print(f"reference TAM-IS genes found: {ng}; TAM-IS prevalence: {ref.obs['is_tamis'].mean():.2f}", flush=True)

    scvi.model.SCVI.setup_anndata(ref, batch_key="cohort", layer="counts")
    model = scvi.model.SCVI(ref, n_latent=30, n_layers=2)
    model.train(max_epochs=30, early_stopping=True, plan_kwargs={"lr": 1e-3})
    ref.obsm["X_scvi"] = model.get_latent_representation()

    qmodel = scvi.model.SCVI.load_query_data(qry, model)
    qmodel.train(max_epochs=30, plan_kwargs={"lr": 1e-3})
    qry.obsm["X_scvi"] = qmodel.get_latent_representation()

    # kNN label transfer in the reference latent
    nn = NearestNeighbors(n_neighbors=30).fit(ref.obsm["X_scvi"])
    _, idx = nn.kneighbors(qry.obsm["X_scvi"])
    ref_label = ref.obs["is_tamis"].values
    qry.obs["tamis_frac"] = ref_label[idx].mean(axis=1)
    qry.obs["mapped_tamis"] = (qry.obs["tamis_frac"] > 0.5).astype(int)
    qry.write_h5ad(os.path.join(OBJ, "step11_query_mapped.h5ad"))

    frac = qry.obs["mapped_tamis"].mean()
    print(f"\nH1 scArches: query {QUERY_COHORT} n={qry.n_obs}; "
          f"cells mapping to TAM-IS (kNN majority) = {frac*100:.1f}%", flush=True)
    res = pd.DataFrame([{
        "reference": "gse103224,gse131928", "query": QUERY_COHORT,
        "n_query": qry.n_obs, "ref_tamis_prevalence": round(ref.obs["is_tamis"].mean(), 3),
        "query_mapped_tamis_frac": round(frac, 3),
        "mean_tamis_frac": round(float(qry.obs["tamis_frac"].mean()), 3),
    }])
    res.to_csv(os.path.join(TBL, "step11_arches_h1.tsv"), sep="\t", index=False)
    print(res.to_string(index=False), flush=True)
    print("DONE", flush=True)
