"""
H3 rigor — deconvolution (NNLS) of CGGA bulk with the scRNA reference.

Builds a cell-type signature matrix (mean linear counts per lineage) from the annotated
scRNA object, then non-negative least squares (NNLS) deconvolves each CGGA bulk sample
into lineage fractions. Re-tests H3 with the deconvolved monocyte_mac (TAM-IS) fraction.

Outputs: results/tables/step12_cgga_deconv.tsv (fractions), step12_cgga_deconv_cox.tsv

Run: env/.venv/bin/python scripts/07_bulk_validation/step12_deconv.py
"""
import os
import zipfile

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import nnls
import statsmodels.api as sm

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
DATASETS = os.environ.get("GLIOMA_DATA") or os.path.join(os.path.dirname(ROOT), "Datasets")
OBJ = os.path.join(ROOT, "objects")
CGGA = os.path.join(DATASETS, "CGGA")
TBL = os.path.join(ROOT, "results/tables")

LINEAGES = ["malignant", "monocyte_mac", "microglia", "T_cell", "endothelial", "NK", "neutrophil"]


def read_cgga_rsem(tag):
    with zipfile.ZipFile(os.path.join(CGGA, f"CGGA.mRNAseq_{tag}.RSEM-genes.zip")) as zf:
        n = [x for x in zf.namelist() if x.endswith(".txt") and "__MACOSX" not in x][0]
        with zf.open(n) as fh:
            df = pd.read_csv(fh, sep="\t", index_col=0)
    df.index = [str(g).strip("'\"") for g in df.index]
    return df


def read_cgga_clin(tag):
    with zipfile.ZipFile(os.path.join(CGGA, f"CGGA.mRNAseq_{tag}_clinical.zip")) as zf:
        n = [x for x in zf.namelist() if x.endswith(".txt") and "__MACOSX" not in x][0]
        with zf.open(n) as fh:
            return pd.read_csv(fh, sep="\t").set_index("CGGA_ID")


if __name__ == "__main__":
    # ---- reference signature matrix (genes x lineages), linear counts ----
    a = ad.read_h5ad(os.path.join(OBJ, "step04_annotated.h5ad"))
    C = a.layers["counts"]
    C = C.toarray() if sp.issparse(C) else np.asarray(C)
    genes = np.array(a.var_names)
    sig = {}
    for lin in LINEAGES:
        m = (a.obs["lineage"].values == lin)
        sig[lin] = C[m].mean(axis=0)
    S = pd.DataFrame(sig, index=genes)
    S = S[~S.index.duplicated()]
    S = S.div(S.sum(axis=0), axis=1) * 1e6          # CPM per lineage
    print(f"signature matrix (CPM): {S.shape}", flush=True)

    rows, cox_rows = [], []
    for tag in ["325", "693"]:
        bulk = read_cgga_rsem(tag)                      # genes x samples (RSEM, already linear)
        bulk = bulk.div(bulk.sum(axis=0), axis=1) * 1e6  # CPM per sample
        common = S.index.intersection(bulk.index)
        Sg = np.nan_to_num(S.loc[common].values.astype(float))
        Bg = np.nan_to_num(bulk.loc[common].values.astype(float))
        fractions = []
        for j in range(Bg.shape[1]):
            b = np.clip(Bg[:, j], 0, None)
            w, _ = nnls(Sg, b)
            fractions.append(w / (w.sum() + 1e-9))
        F = pd.DataFrame(fractions, index=bulk.columns, columns=LINEAGES)
        F["tag"] = tag
        rows.append(F)
        print(f"CGGA_{tag}: deconvolved {F.shape}", flush=True)

        cl = read_cgga_clin(tag)
        d = F.join(cl, how="inner")
        d["OS"] = pd.to_numeric(d["OS"], errors="coerce")
        d["event"] = pd.to_numeric(d["Censor (alive=0; dead=1)"], errors="coerce")
        d = d[d["OS"].notna() & d["event"].notna() & (d["OS"] > 0)].copy()
        d["mac_frac"] = d["monocyte_mac"]
        for model, covs in [("mac_only", ["mac_frac"]),
                            ("mac_plus_micro", ["mac_frac", "microglia"])]:
            c = d.dropna(subset=covs)
            if c["mac_frac"].std() < 1e-6:
                print(f"  {tag} {model}: mac_frac has no variance; skipped", flush=True)
                continue
            X = sm.add_constant(c[covs].astype(float))
            try:
                # ties='efron': the default (breslow) can give a degenerate Hessian.
                res = sm.PHReg(c["OS"], X, status=c["event"].astype(int), ties="efron").fit()
            except Exception as e:
                print(f"  {tag} {model}: fit failed ({e})", flush=True)
                continue
            for i, v in enumerate(X.columns):
                cox_rows.append({"cohort": f"CGGA_{tag}", "model": model, "covariate": v,
                                 "HR": round(float(np.exp(res.params[i])), 3),
                                 "p": float(res.pvalues[i]), "n": len(c)})

    pd.concat(rows).to_csv(os.path.join(TBL, "step12_cgga_deconv.tsv"), sep="\t")
    cox = pd.DataFrame(cox_rows)
    cox.to_csv(os.path.join(TBL, "step12_cgga_deconv_cox.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 200)
    print("\n=== H3 with deconvolved TAM-IS (monocyte_mac) fraction ===", flush=True)
    print(cox.to_string(index=False), flush=True)
    print("DONE", flush=True)
