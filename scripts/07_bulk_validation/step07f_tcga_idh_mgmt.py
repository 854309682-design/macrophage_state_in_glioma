"""
Step 7f (P2) — H3 with complete, uniform TCGA clinical covariates (incl. IDH/MGMT).

Motivation / bug fixed: step07b/d/e took TCGA survival from the GDC `cases` endpoint,
which (a) does NOT expose IDH or MGMT, and (b) returns NO `days_to_last_follow_up` for
GBM (3/293 non-null) -> every LIVING GBM case got OS=NaN and was dropped, silently
turning the GBM Cox into a dead-only analysis (the "spurious all-Dead cohort" trap).
Verified: GDC GBM subset showed 205/207 events.

Fix: take ALL covariates from cBioPortal study `lgggbm_tcga_pub` (Brennan/Ceccarelli
TCGA LGG+GBM marker paper), which provides uniform, complete patient-level
OS_MONTHS + OS_STATUS (+ AGE) and sample-level GRADE, IDH_STATUS, MGMT_PROMOTER_STATUS.
Sample attributes are aggregated to patient, then the H3 Cox is refit adjusted for
age + grade + IDH (+ MGMT), plus G4/GBM vs G2-G3 and IDH-stratified models — mirroring
step07e's CGGA analysis so TCGA and CGGA are compared on equal footing.

Outputs:
  Datasets/TCGA/TCGA_lgggbm_clinical.tsv                 (cached patient-level covariates)
  results/tables/step07_tcga_idh_mgmt_stratified.tsv

Run: env/.venv/bin/python scripts/07_bulk_validation/step07f_tcga_idh_mgmt.py
"""
import json
import os
import urllib.request

import numpy as np
import pandas as pd
import statsmodels.api as sm

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
DATASETS = os.environ.get("GLIOMA_DATA") or os.path.join(os.path.dirname(ROOT), "Datasets")
TCGA = os.path.join(DATASETS, "TCGA")
TBL = os.path.join(ROOT, "results/tables")
CACHE = os.path.join(TCGA, "TCGA_lgggbm_clinical.tsv")
STUDY = "lgggbm_tcga_pub"

TAMIS = ["SPP1", "TREM2", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "LGALS3", "CD163",
         "CTSL", "CCL3", "CCL4", "GPNMB", "LPL", "HLA-DRA", "HLA-DRB1", "S100A9",
         "HMOX1", "FTL", "FTH1"]
MYELOID = ["PTPRC", "ITGAM", "CSF1R", "AIF1", "TYROBP", "FCER1G", "LYZ", "CD68", "CD14"]

# attribute -> cBioPortal clinicalDataType
ATTRS = {"os_status": ("OS_STATUS", "PATIENT"), "os_months": ("OS_MONTHS", "PATIENT"),
         "age": ("AGE", "PATIENT"), "grade": ("GRADE", "SAMPLE"),
         "idh": ("IDH_STATUS", "SAMPLE"), "mgmt": ("MGMT_PROMOTER_STATUS", "SAMPLE")}

rows = []


def cbio(attribute, level):
    url = ("https://www.cbioportal.org/api/studies/%s/clinical-data"
           "?clinicalDataType=%s&attributeId=%s&projection=SUMMARY" % (STUDY, level, attribute))
    with urllib.request.urlopen(url, timeout=180) as r:
        d = json.load(r)
    out = {}
    for x in d:
        p = x.get("patientId")
        if p and p not in out:          # first non-null per patient (sample attrs)
            out[p] = x.get("value")
    return out


def covariates():
    if os.path.exists(CACHE):
        print("using cached covariates", CACHE, flush=True)
        return pd.read_csv(CACHE, sep="\t")
    cols = {k: cbio(a, lv) for k, (a, lv) in ATTRS.items()}
    pats = sorted(set().union(*[set(v) for v in cols.values()]))
    df = pd.DataFrame({"case": pats})
    for k, v in cols.items():
        df[k] = df["case"].map(v)
    df.to_csv(CACHE, sep="\t", index=False)
    print("wrote covariates:", CACHE, df.shape, flush=True)
    return df


def cox(df, covs, tag):
    c = df.dropna(subset=covs)
    if len(c) < 30 or c["event"].sum() < 10:
        rows.append({"dataset": tag, "model": "|".join(covs), "n": len(c),
                     "events": int(c["event"].sum()), "covariate": "tamis",
                     "HR": np.nan, "p": np.nan})
        return
    X = sm.add_constant(c[covs].astype(float))
    try:
        # ties='efron' explicitly: the default ('breslow' here) can give a degenerate Hessian.
        res = sm.PHReg(c["OS"], X, status=c["event"].astype(int), ties="efron").fit()
    except Exception as e:  # e.g. LinAlgError from collinear covariates -> record and continue
        print(f"  [{tag}] fit failed: {type(e).__name__}", flush=True)
        rows.append({"dataset": tag, "model": "|".join(covs), "n": len(c),
                     "events": int(c["event"].sum()), "covariate": "tamis",
                     "HR": np.nan, "p": np.nan})
        return
    for i, v in enumerate(X.columns):
        rows.append({"dataset": tag, "model": "|".join(covs), "n": len(c),
                     "events": int(c["event"].sum()), "covariate": v,
                     "HR": round(float(np.exp(res.params[i])), 3), "p": float(res.pvalues[i])})


if __name__ == "__main__":
    m = pd.read_csv(os.path.join(TBL, "step07_tcga_signature_matrix.tsv"), sep="\t")
    tam = m.groupby("case")[[g for g in TAMIS if g in m.columns]].mean()
    mye = m.groupby("case")[[g for g in MYELOID if g in m.columns]].mean()
    tz = (tam - tam.mean()) / tam.std()
    mz = (mye - mye.mean()) / mye.std()
    t = pd.DataFrame({"case": tz.index, "tamis": tz.mean(axis=1).values,
                      "myeloid": mz.mean(axis=1).values})
    t = t.merge(covariates(), on="case", how="inner")

    t["OS"] = pd.to_numeric(t["os_months"], errors="coerce")
    t["event"] = t["os_status"].astype(str).str.contains("DECEASED").astype(int)
    t = t[t["OS"].notna() & (t["OS"] > 0)].copy()
    t["age_y"] = pd.to_numeric(t["age"], errors="coerce")
    t["grade_num"] = t["grade"].astype(str).str.extract(r"G(\d)").astype(float)
    t["idh_mut"] = np.where(t["idh"].isin(["Mutant", "WT"]), t["idh"].eq("Mutant"), np.nan)
    t["mgmt_meth"] = np.where(t["mgmt"].isin(["Methylated", "Unmethylated"]),
                              t["mgmt"].eq("Methylated"), np.nan)
    t["gbm"] = (t["grade_num"] == 4)

    print("grade:", t["grade"].value_counts(dropna=False).to_dict(), flush=True)
    print("age range: %.0f-%.0f" % (t["age_y"].min(), t["age_y"].max()), flush=True)
    print("n with OS: %d; events: %d; grade: %d; idh: %d; mgmt: %d" % (
        len(t), int(t["event"].sum()), t["grade_num"].notna().sum(),
        t["idh_mut"].notna().sum(), t["mgmt_meth"].notna().sum()), flush=True)
    print("G4 events/alive: %d/%d" % (
        int(t.loc[t["gbm"], "event"].sum()), int((1 - t.loc[t["gbm"], "event"]).sum())), flush=True)

    cox(t, ["tamis", "age_y", "grade_num"], "TCGA_all_age_grade")
    cox(t, ["tamis", "age_y", "grade_num", "idh_mut"], "TCGA_all_age_grade_idh")
    cox(t, ["tamis", "age_y", "grade_num", "idh_mut", "mgmt_meth"], "TCGA_all_age_grade_idh_mgmt")
    cox(t, ["tamis", "myeloid", "age_y", "grade_num", "idh_mut"], "TCGA_all_state_over_count_idh")
    cox(t[t["gbm"]], ["tamis", "age_y", "idh_mut"], "TCGA_GBM_G4_age_idh")
    cox(t[~t["gbm"]], ["tamis", "age_y", "idh_mut"], "TCGA_LGG_G2G3_age_idh")
    cox(t[t["idh_mut"] == True], ["tamis"], "TCGA_IDHmut")
    cox(t[t["idh_mut"] == False], ["tamis"], "TCGA_IDHwt")
    cox(t[t["idh_mut"] == True], ["tamis", "myeloid"], "TCGA_IDHmut_state_over_count")
    cox(t[t["idh_mut"] == False], ["tamis", "myeloid"], "TCGA_IDHwt_state_over_count")

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(TBL, "step07_tcga_idh_mgmt_stratified.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_rows", 300)
    print("\n=== tamis rows ===", flush=True)
    print(out[out["covariate"] == "tamis"].to_string(index=False), flush=True)
    print("\n=== all rows ===", flush=True)
    print(out.to_string(index=False), flush=True)
    print("DONE", flush=True)
