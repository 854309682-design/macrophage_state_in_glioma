"""
H3 rigor (b) — survival with BayesPrism-deconvolved fractions.

Reads results/bayesprism/fraction_<tag>.csv (samples x cell types), joins CGGA clinical,
and runs Cox on the deconvolved monocyte_mac (TAM-IS) fraction (and mac + microglia).
Resolves the signature-vs-NNLS contradiction with the proper method.

Outputs: results/tables/step13_bayesprism_cox.tsv

Run: env/.venv/bin/python scripts/07_bulk_validation/step13c_bayesprism_cox.py
"""
import os
import zipfile

import numpy as np
import pandas as pd
import statsmodels.api as sm

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
DATASETS = os.environ.get("GLIOMA_DATA") or os.path.join(os.path.dirname(ROOT), "Datasets")
BP = os.path.join(ROOT, "results/bayesprism")
CGGA = os.path.join(DATASETS, "CGGA")
TBL = os.path.join(ROOT, "results/tables")


def read_clin(tag):
    with zipfile.ZipFile(os.path.join(CGGA, f"CGGA.mRNAseq_{tag}_clinical.zip")) as zf:
        n = [x for x in zf.namelist() if x.endswith(".txt") and "__MACOSX" not in x][0]
        with zf.open(n) as fh:
            return pd.read_csv(fh, sep="\t").set_index("CGGA_ID")


rows = []
for tag in ["325", "693"]:
    f = os.path.join(BP, f"fraction_{tag}.csv")
    if not os.path.exists(f):
        print(f"missing {f}; skip", flush=True)
        continue
    theta = pd.read_csv(f, index_col=0)          # samples x cell types (BayesPrism 2.2.3 get.fraction)
    F = theta                                     # already sample-major — do NOT transpose
    cl = read_clin(tag)
    d = F.join(cl, how="inner")
    d["OS"] = pd.to_numeric(d["OS"], errors="coerce")
    d["event"] = pd.to_numeric(d["Censor (alive=0; dead=1)"], errors="coerce")
    d = d[d["OS"].notna() & d["event"].notna() & (d["OS"] > 0)].copy()
    d["mac_frac"] = d["monocyte_mac"]
    print(f"\n=== CGGA_{tag}: n={len(d)}, events={int(d['event'].sum())} ===", flush=True)
    print("mean fractions:", {c: round(float(d[c].mean()), 3) for c in F.columns}, flush=True)
    for model, covs in [("mac_only", ["mac_frac"]),
                        ("mac_plus_micro", ["mac_frac", "microglia"])]:
        c = d.dropna(subset=covs)
        if c["mac_frac"].std() < 1e-9:
            print(f"  {model}: no variance", flush=True); continue
        X = sm.add_constant(c[covs])
        try:
            res = sm.PHReg(c["OS"], X, status=c["event"].astype(int), ties="efron").fit()
        except Exception as e:
            print(f"  {model}: fit failed {e}", flush=True); continue
        for i, v in enumerate(X.columns):
            rows.append({"cohort": f"CGGA_{tag}", "model": model, "covariate": v,
                         "HR": round(float(np.exp(res.params[i])), 3),
                         "p": float(res.pvalues[i]), "n": len(c)})

out = pd.DataFrame(rows)
out.to_csv(os.path.join(TBL, "step13_bayesprism_cox.tsv"), sep="\t", index=False)
pd.set_option("display.width", 200)
print("\n=== BayesPrism H3 (deconvolved TAM-IS fraction) ===", flush=True)
print(out.to_string(index=False), flush=True)
print("DONE", flush=True)
