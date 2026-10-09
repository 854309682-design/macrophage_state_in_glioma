"""
Step 16 — Cross-cohort random-effects meta-analysis of the TAM-IS hazard ratio (H3).

Pools the multivariable-adjusted TAM-IS HR across TCGA and CGGA-325/693, overall and within
GBM/G4, LGG/G2-G3 and IDH-mut / IDH-wt strata, using a DerSimonian–Laird random-effects model
on ln(HR).

Sources:
  - TCGA: results/tables/step07_tcga_idh_mgmt_stratified.tsv (step07f, cBioPortal covariates;
    SE recovered from HR + two-sided p).
  - CGGA: re-fit here from the CGGA RSEM + clinical (local) so every stratum has a proper SE.

METHOD NOTE (important): all Cox fits here use ties='efron'. This statsmodels build defaults to
ties='breslow', which yields a NEGATIVE-DEFINITE Hessian (bse=0, p=0.0) for several of these
subsets — the older step07c/d/e tables therefore contain spurious p=0.0 rows (e.g. TCGA IDH-mut,
CGGA G4/G2-G3). Those are recomputed correctly here.

Outputs: results/tables/step16_meta_analysis.tsv

Run: env/.venv/bin/python scripts/07_bulk_validation/step16_meta_analysis.py
"""
import os
import zipfile

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
DATASETS = os.environ.get("GLIOMA_DATA") or os.path.join(os.path.dirname(ROOT), "Datasets")
CGGA = os.path.join(DATASETS, "CGGA")
TBL = os.path.join(ROOT, "results/tables")

TAMIS = ["SPP1", "TREM2", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "LGALS3", "CD163",
         "CTSL", "CCL3", "CCL4", "GPNMB", "LPL", "HLA-DRA", "HLA-DRB1", "S100A9",
         "HMOX1", "FTL", "FTH1"]


def se_from_hr_p(hr, p):
    z = norm.isf(max(min(p, 1 - 1e-16), 1e-300) / 2.0)
    return abs(np.log(hr)) / z


def read_cgga(tag):
    def rd(z, kind):
        with zipfile.ZipFile(z) as zf:
            n = [x for x in zf.namelist() if x.endswith(".txt") and "__MACOSX" not in x][0]
            with zf.open(n) as fh:
                return pd.read_csv(fh, sep="\t", index_col=0 if kind == "expr" else None)
    e = rd(os.path.join(CGGA, f"CGGA.mRNAseq_{tag}.RSEM-genes.zip"), "expr")
    e.index = [str(g).strip("'\"") for g in e.index]
    cl = rd(os.path.join(CGGA, f"CGGA.mRNAseq_{tag}_clinical.zip"), "clin").set_index("CGGA_ID")
    return e, cl


def cgga_scores(tag):
    e, cl = read_cgga(tag)
    g = [x for x in TAMIS if x in e.index]
    z = e.loc[g]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, np.nan), axis=0)
    s = pd.DataFrame({"tamis": z.mean(0)}).join(cl)
    s["OS"] = pd.to_numeric(s["OS"], errors="coerce")
    s["event"] = pd.to_numeric(s["Censor (alive=0; dead=1)"], errors="coerce")
    s = s[s["OS"].notna() & s["event"].notna() & (s["OS"] > 0)].copy()
    s["age"] = pd.to_numeric(s["Age"], errors="coerce")
    s["grade_num"] = s["Grade"].astype(str).str.extract(r"(IV|III|II)")[0].map({"II": 2, "III": 3, "IV": 4})
    s["idh"] = (s["IDH_mutation_status"].astype(str).str.lower().str.contains("mut")).astype(float)
    s["mgmt"] = (s["MGMTp_methylation_status"].astype(str).str.lower().str.contains("methyl")).astype(float)
    return s


def fit(s, covs, tag, stratum, studies):
    c = s.dropna(subset=covs)
    if len(c) < 30 or c["event"].sum() < 10:
        return
    X = sm.add_constant(c[covs].astype(float))
    res = sm.PHReg(c["OS"], X, status=c["event"].astype(int), ties="efron").fit()
    i = list(X.columns).index("tamis")
    studies.append({"study": tag, "stratum": stratum, "n": len(c), "events": int(c["event"].sum()),
                    "HR": float(np.exp(res.params[i])), "se_logHR": float(res.bse[i]),
                    "p": float(res.pvalues[i])})


def dl_meta(rows, label):
    y = np.log(np.array([r["HR"] for r in rows], float))
    se = np.array([r["se_logHR"] for r in rows], float)
    k = len(rows)
    w = 1.0 / se ** 2
    yfix = (w * y).sum() / w.sum()
    Q = float((w * (y - yfix) ** 2).sum())
    C = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    ws = 1.0 / (se ** 2 + tau2)
    mu = (ws * y).sum() / ws.sum()
    se_mu = 1.0 / np.sqrt(ws.sum())
    p = 2 * (1 - norm.cdf(abs(mu / se_mu)))
    i2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    return {"pool": label, "k": k, "pooled_HR": round(float(np.exp(mu)), 3),
            "lo": round(float(np.exp(mu - 1.96 * se_mu)), 3),
            "hi": round(float(np.exp(mu + 1.96 * se_mu)), 3), "p": float(p),
            "Q": round(Q, 2), "I2_pct": round(i2, 1), "tau2": round(tau2, 4)}


if __name__ == "__main__":
    studies = []

    # ---- TCGA (unbiased, step07f) ----
    t = pd.read_csv(os.path.join(TBL, "step07_tcga_idh_mgmt_stratified.tsv"), sep="\t")
    tcga_map = [("TCGA_all_age_grade_idh", "TCGA_all", "all"),
                ("TCGA_GBM_G4_age_idh", "TCGA_G4", "G4"),
                ("TCGA_LGG_G2G3_age_idh", "TCGA_G2G3", "G2G3"),
                ("TCGA_IDHmut", "TCGA_IDHmut", "IDHmut"),
                ("TCGA_IDHwt", "TCGA_IDHwt", "IDHwt")]
    for ds, tag, stratum in tcga_map:
        r = t[(t["dataset"] == ds) & (t["covariate"] == "tamis")].iloc[0]
        studies.append({"study": tag, "stratum": stratum, "n": int(r["n"]), "events": int(r["events"]),
                        "HR": float(r["HR"]), "se_logHR": se_from_hr_p(float(r["HR"]), float(r["p"])),
                        "p": float(r["p"])})

    # ---- CGGA (re-fit with efron ties) ----
    for tag in ["325", "693"]:
        s = cgga_scores(tag)
        fit(s, ["tamis", "age", "grade_num", "idh", "mgmt"], f"CGGA_{tag}_all", "all", studies)
        fit(s[s["grade_num"] == 4], ["tamis"], f"CGGA_{tag}_G4", "G4", studies)
        fit(s[s["grade_num"].isin([2, 3])], ["tamis"], f"CGGA_{tag}_G2G3", "G2G3", studies)
        fit(s[s["idh"] == 1], ["tamis"], f"CGGA_{tag}_IDHmut", "IDHmut", studies)
        fit(s[s["idh"] == 0], ["tamis"], f"CGGA_{tag}_IDHwt", "IDHwt", studies)

    df = pd.DataFrame(studies)

    pools = [dl_meta([r for r in studies if r["stratum"] == "all"], "PRIMARY all cohorts (adjusted)")]
    for grp, lab in [("G4", "G4/GBM strata"), ("G2G3", "G2-G3/LGG strata"),
                     ("IDHmut", "IDH-mutant strata"), ("IDHwt", "IDH-wildtype strata")]:
        rows = [r for r in studies if r["stratum"] == grp]
        if len(rows) >= 2:
            pools.append(dl_meta(rows, lab))
    pools = pd.DataFrame(pools)

    df.to_csv(os.path.join(TBL, "step16_meta_analysis.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 200)
    print("=== per-study (ties=efron) ===", flush=True)
    print(df[["study", "stratum", "n", "events", "HR", "se_logHR", "p"]].to_string(index=False), flush=True)
    print("\n=== random-effects pools (DerSimonian-Laird) ===", flush=True)
    print(pools.to_string(index=False), flush=True)
    print("\nDONE", flush=True)
