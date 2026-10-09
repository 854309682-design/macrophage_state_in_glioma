"""
Step 4c — H1 replication (signature-projection reproducibility).

For each of the 5 scRNA cohorts, score the published TAM-IS vs microglia signatures
(derived from Cheng 2021 Cell / Friebel 2020 Immunity / Ochocka 2021 Nat Commun +
our markers), identify the myeloid compartment by a myeloid-marker score, and report
per-cohort whether a TAM-IS-high (Mo-Mac) population is present.

H1 is supported if TAM-IS-high Mo-Mac are detected in every cohort.

Signatures:
  TAM-IS (SPP1+ TAM / TREM2+ LAM + complement/MHC-II/inflammatory):
    SPP1, TREM2, APOE, APOC1, C1QA, C1QB, C1QC, LGALS3, CD163, CTSL, CCL3, CCL4,
    GPNMB, LPL, HLA-DRA, HLA-DRB1, S100A9, HMOX1, FTL, FTH1
  microglia (HomMG): P2RY12, TMEM119, CX3CR1, SALL1, CRYBB1, CST3, PROS1, MERTK
  myeloid: LYZ, CD14, ITGAM, AIF1, CSF1R, C1QA, TYROBP, FCER1G

Run: env/.venv/bin/python scripts/04_myeloid_states/step06_h1_replication.py
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

sc.settings.verbosity = 0

COHORTS = ["gse103224", "gse131928", "gse163120", "gse89567", "gse70630"]
TAMIS = ["SPP1", "TREM2", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "LGALS3",
         "CD163", "CTSL", "CCL3", "CCL4", "GPNMB", "LPL", "HLA-DRA", "HLA-DRB1",
         "S100A9", "HMOX1", "FTL", "FTH1"]
MICRO = ["P2RY12", "TMEM119", "CX3CR1", "SALL1", "CRYBB1", "CST3", "PROS1", "MERTK"]
MYELOID = ["LYZ", "CD14", "ITGAM", "AIF1", "CSF1R", "C1QA", "TYROBP", "FCER1G"]


def load(name):
    for pref in ("step03", "step02"):
        p = os.path.join(OBJ, f"{pref}_{name}.h5ad")
        if os.path.exists(p):
            a = ad.read_h5ad(p)
            if "symbol" in a.var.columns:  # GSE103224 kept ensembl index + symbol col
                sym = a.var["symbol"].astype(str).values
                a.var = a.var.drop(columns=["symbol"])
                a.var_names = sym
            a.var_names_make_unique()
            return a, pref
    raise FileNotFoundError(name)


def score(a, genes, name):
    g = [x for x in genes if x in a.var_names]
    sc.tl.score_genes(a, g, score_name=name, use_raw=False)
    return g


rows = []
for c in COHORTS:
    a, pref = load(c)
    a.X = np.nan_to_num(a.X)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    gt = score(a, TAMIS, "tamis")
    gm = score(a, MICRO, "micro")
    gy = score(a, MYELOID, "myeloid")
    mycol = a.obs["myeloid"] > 0.0
    my = a.obs[mycol]
    tamis_hi = (my["tamis"] > 0.0) & (my["micro"] < 0.0)
    print(f"{c} [{pref}]: {a.n_obs} cells, myeloid={int(mycol.sum())} "
          f"({mycol.mean()*100:.1f}%), TAM-IS-high my={int(tamis_hi.sum())} "
          f"({tamis_hi.mean()*100:.1f}% of myeloid)", flush=True)
    rows.append({
        "cohort": c, "source_obj": pref, "n_cells": a.n_obs,
        "n_myeloid": int(mycol.sum()), "pct_myeloid": round(mycol.mean() * 100, 1),
        "n_tamis_high_myeloid": int(tamis_hi.sum()),
        "pct_tamis_high_of_myeloid": round(tamis_hi.mean() * 100, 1),
        "tamis_genes_found": len(gt), "micro_genes_found": len(gm),
        "myeloid_genes_found": len(gy),
    })

df = pd.DataFrame(rows)
os.makedirs(TBL, exist_ok=True)
df.to_csv(os.path.join(TBL, "step06_h1_replication.tsv"), sep="\t", index=False)
print("\n=== H1 replication summary ===", flush=True)
print(df.to_string(index=False), flush=True)
n_rep = (df["pct_tamis_high_of_myeloid"] > 2).sum()
print(f"\ncohorts with TAM-IS-high myeloid >2%: {n_rep}/5", flush=True)
print("DONE", flush=True)
