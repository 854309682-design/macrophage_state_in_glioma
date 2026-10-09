"""
Step 17 — Deconvolution method consistency (§8.5): BayesPrism vs NNLS vs MuSiC.

`analysis_strategy_1.md` §8.5 requires >=2 deconvolution methods to agree (Spearman > 0.7) on
the TAM-IS (monocyte_mac) fraction, else report the inconsistency and adjudicate with a third.
Here we compare the CGGA-325/693 `monocyte_mac` and `microglia` fractions across three methods:
  - BayesPrism : results/bayesprism/fraction_<tag>.csv            (samples x cell types)
  - NNLS       : results/tables/step12_cgga_deconv.tsv            (step12, scipy nnls)
  - MuSiC      : results/bayesprism/fraction_music_<tag>.csv      (step14, pseudo-samples)
and re-test the H3 association of the monocyte_mac fraction per method (Cox, ties='efron').

Outputs: results/tables/step17_deconv_consistency.tsv

Run: env/.venv/bin/python scripts/07_bulk_validation/step17_deconv_consistency.py
"""
import os
import zipfile

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

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


if __name__ == "__main__":
    nnls_all = pd.read_csv(os.path.join(TBL, "step12_cgga_deconv.tsv"), sep="\t", index_col=0)
    corr_rows, cox_rows = [], []

    for tag in ["325", "693"]:
        methods = {}
        bp = pd.read_csv(os.path.join(BP, f"fraction_{tag}.csv"), index_col=0)   # samples x types
        methods["BayesPrism"] = bp
        nn = nnls_all[nnls_all["tag"].astype(str) == str(tag)].drop(columns=["tag"])
        methods["NNLS"] = nn
        mf = os.path.join(BP, f"fraction_music_{tag}.csv")
        if os.path.exists(mf):
            mu = pd.read_csv(mf, index_col=0)
            # MuSiC columns may be prefixed; keep the monocyte_mac / microglia ones
            methods["MuSiC"] = mu

        # pairwise Spearman on the monocyte_mac fraction
        names = list(methods)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = names[i], names[j]
                ca, cb = methods[a], methods[b]
                if "monocyte_mac" not in ca.columns or "monocyte_mac" not in cb.columns:
                    continue
                common = ca.index.intersection(cb.index)
                rho, p = spearmanr(ca.loc[common, "monocyte_mac"], cb.loc[common, "monocyte_mac"])
                corr_rows.append({"cohort": f"CGGA_{tag}", "method_a": a, "method_b": b,
                                  "cell_type": "monocyte_mac", "n": len(common),
                                  "spearman_rho": round(float(rho), 3), "p": float(p)})

        # Cox per method (ties='efron')
        cl = read_clin(tag)
        for name, frac in methods.items():
            if "monocyte_mac" not in frac.columns:
                continue
            d = frac.join(cl, how="inner")
            d["OS"] = pd.to_numeric(d["OS"], errors="coerce")
            d["event"] = pd.to_numeric(d["Censor (alive=0; dead=1)"], errors="coerce")
            d = d[d["OS"].notna() & d["event"].notna() & (d["OS"] > 0)].copy()
            d["mac_frac"] = d["monocyte_mac"]
            c = d.dropna(subset=["mac_frac"])
            if len(c) < 10 or c["mac_frac"].std() < 1e-9 or c["event"].sum() < 5:
                continue
            X = sm.add_constant(c[["mac_frac"]].astype(float))
            res = sm.PHReg(c["OS"], X, status=c["event"].astype(int), ties="efron").fit()
            sd = float(c["mac_frac"].std())
            cox_rows.append({"cohort": f"CGGA_{tag}", "method": name, "n": len(c),
                             "events": int(c["event"].sum()),
                             "mac_frac_mean": round(float(c["mac_frac"].mean()), 4),
                             "mac_frac_sd": round(sd, 4),
                             "HR_per_unit": round(float(np.exp(res.params[1])), 3),
                             "HR_per_SD": round(float(np.exp(res.params[1] * sd)), 3),
                             "p": float(res.pvalues[1])})

    corr = pd.DataFrame(corr_rows)
    cox = pd.DataFrame(cox_rows)
    out = os.path.join(TBL, "step17_deconv_consistency.tsv")
    with open(out, "w") as fh:
        fh.write("# pairwise Spearman of the monocyte_mac fraction\n")
        corr.to_csv(fh, sep="\t", index=False)
        fh.write("\n# H3 association of the monocyte_mac fraction, per method (Cox, ties=efron)\n")
        cox.to_csv(fh, sep="\t", index=False)
    pd.set_option("display.width", 200)
    print("=== method agreement (monocyte_mac, Spearman) ===", flush=True)
    print(corr.to_string(index=False), flush=True)
    print("\n=== H3 per method (mac fraction) ===", flush=True)
    print(cox.to_string(index=False), flush=True)
    print("\nDONE", flush=True)
