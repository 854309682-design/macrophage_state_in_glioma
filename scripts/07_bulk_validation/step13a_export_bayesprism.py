"""
H3 rigor (b) — export raw-count reference + CGGA bulk for BayesPrism.

Rebuilds raw counts for the annotated cells from the step02 (pre-normalization) objects,
attaches lineage labels, downsamples per lineage, and exports genes x cells counts.
Also exports CGGA_325/693 bulk (genes x samples) on the same gene space.

Outputs: results/bayesprism/{ref.mtx, ref_genes.txt, ref_cells.txt, ref_labels.txt,
         bulk_325.tsv, bulk_693.tsv, common_genes.txt}

Run: env/.venv/bin/python scripts/07_bulk_validation/step13a_export_bayesprism.py
"""
import os
import zipfile

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io as sio

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
DATASETS = os.environ.get("GLIOMA_DATA") or os.path.join(os.path.dirname(ROOT), "Datasets")
OBJ = os.path.join(ROOT, "objects")
CGGA = os.path.join(DATASETS, "CGGA")
OUT = os.path.join(ROOT, "results/bayesprism")
os.makedirs(OUT, exist_ok=True)
PER = 1000
COHORTS = ["gse103224", "gse131928", "gse163120"]


def load_raw(name):
    a = ad.read_h5ad(os.path.join(OBJ, f"step02_{name}.h5ad"))
    if "symbol" in a.var.columns:
        sym = a.var["symbol"].astype(str).values
        a.var = a.var.drop(columns=["symbol"])
        a.var_names = sym
    a.var_names_make_unique()
    return a


def read_cgga_rsem(tag):
    with zipfile.ZipFile(os.path.join(CGGA, f"CGGA.mRNAseq_{tag}.RSEM-genes.zip")) as zf:
        n = [x for x in zf.namelist() if x.endswith(".txt") and "__MACOSX" not in x][0]
        with zf.open(n) as fh:
            df = pd.read_csv(fh, sep="\t", index_col=0)
    df.index = [str(g).strip("'\"") for g in df.index]
    return df


if __name__ == "__main__":
    lin = ad.read_h5ad(os.path.join(OBJ, "step04_annotated.h5ad"), backed="r").obs["lineage"]
    print(f"lineage labels: {len(lin)}", flush=True)

    adatas = []
    for c in COHORTS:
        a = load_raw(c)
        adatas.append(a)
    ref = ad.concat(adatas, join="outer")  # same names as integration (no suffix)
    ref.var_names_make_unique()
    ref = ref[ref.obs_names.isin(lin.index)].copy()
    ref.obs["lineage"] = lin.reindex(ref.obs_names).values
    print(f"raw reference: {ref.shape}", flush=True)

    # downsample per lineage
    rng = np.random.default_rng(42)
    idx = []
    for L in ["malignant", "monocyte_mac", "microglia", "T_cell", "endothelial", "NK", "neutrophil"]:
        cells = np.where(ref.obs["lineage"].values == L)[0]
        if len(cells) > PER:
            cells = rng.choice(cells, PER, replace=False)
        idx.append(cells)
    idx = np.concatenate(idx)
    sub = ref[idx].copy()
    print(f"downsampled reference: {sub.shape}", flush=True)

    # bulk (CGGA) on shared genes
    bulk = {}
    for tag in ["325", "693"]:
        b = read_cgga_rsem(tag)
        bulk[tag] = b
    common = set(sub.var_names)
    for tag in ["325", "693"]:
        common &= set(bulk[tag].index)
    common = sorted(common)
    print(f"common genes: {len(common)}", flush=True)
    pd.Series(common).to_csv(os.path.join(OUT, "common_genes.txt"), index=False, header=False)

    sub = sub[:, common].copy()
    X = sub.X
    X = sp.csr_matrix(X).T  # genes x cells
    sio.mmwrite(os.path.join(OUT, "ref.mtx"), X)
    np.savetxt(os.path.join(OUT, "ref_genes.txt"), sub.var_names, fmt="%s")
    np.savetxt(os.path.join(OUT, "ref_cells.txt"), sub.obs_names, fmt="%s")
    np.savetxt(os.path.join(OUT, "ref_labels.txt"), sub.obs["lineage"].values, fmt="%s")
    for tag in ["325", "693"]:
        bulk[tag].loc[common].to_csv(os.path.join(OUT, f"bulk_{tag}.tsv"), sep="\t")
    print(f"exported reference {X.shape} + bulk to {OUT}", flush=True)
    print("DONE", flush=True)
