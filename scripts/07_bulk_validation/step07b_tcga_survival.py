"""
Step 7b — TCGA-GBM/LGG survival analysis (H3).

Loads the TCGA signature matrix (TPM of TAM-IS/microglia/immune genes), aggregates to
case level, fetches clinical (OS, age, grade, IDH) from GDC, scores TAM-IS vs microglia
(mean z), then:
  - Kaplan-Meier (TAM-IS high vs low by median) + log-rank,
  - multivariable Cox (statsmodels PHReg) adjusting for age and grade (and IDH if present).
Outputs: results/tables/step07_tcga_survival_cox.tsv, results/figures/main/step07_tcga_km.{pdf,jpeg}

Run: env/.venv/bin/python scripts/07_bulk_validation/step07b_tcga_survival.py
"""
import json
import os
import urllib.request

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.duration.survfunc import SurvfuncRight, survdiff

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
TBL = os.path.join(ROOT, "results/tables")
FIG = os.path.join(ROOT, "results/figures/main")
os.makedirs(FIG, exist_ok=True)

TAMIS = ["SPP1", "TREM2", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "LGALS3", "CD163",
         "CTSL", "CCL3", "CCL4", "GPNMB", "LPL", "HLA-DRA", "HLA-DRB1", "S100A9",
         "HMOX1", "FTL", "FTH1"]
MICRO = ["P2RY12", "TMEM119", "CX3CR1", "SALL1", "CRYBB1", "CST3", "PROS1", "MERTK"]


def clinical(cases):
    """GDC cases endpoint (verified field nesting):
    vital_status + days_to_death under demographic; days_to_last_follow_up + age + grade under diagnoses."""
    payload = {"filters": {"op": "in",
                           "content": {"field": "project.project_id",
                                       "value": ["TCGA-GBM", "TCGA-LGG"]}},
               "fields": "submitter_id,demographic.vital_status,demographic.days_to_death,"
                         "diagnoses.days_to_last_follow_up,diagnoses.age_at_diagnosis,"
                         "diagnoses.tumor_grade",
               "format": "json", "size": 5000}
    recs = []
    try:
        req = urllib.request.Request("https://api.gdc.cancer.gov/cases",
                                     data=json.dumps(payload).encode(),
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=120) as r:
            d = json.load(r)
        for h in d["data"]["hits"]:
            diag = (h.get("diagnoses") or [{}])[0]
            recs.append({
                "case": h["submitter_id"],
                "vital": (h.get("demographic") or {}).get("vital_status"),
                "dtd": (h.get("demographic") or {}).get("days_to_death"),
                "dtf": diag.get("days_to_last_follow_up"),
                "age": diag.get("age_at_diagnosis"),
                "grade": diag.get("tumor_grade"),
            })
    except Exception as e:
        print("clinical err", e, flush=True)
    return pd.DataFrame(recs).drop_duplicates("case")


def km_plot(df, x, fn):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 5))
    for grp, sub in df.groupby(x):
        sf = SurvfuncRight(sub["os"], sub["event"])
        ax.step(sf.surv_times, sf.surv_prob, where="post", label=f"{grp} (n={len(sub)})")
    ax.set_xlabel("Overall survival (days)")
    ax.set_ylabel("Survival probability")
    ax.set_title("TCGA-GBM/LGG: TAM-IS high vs low")
    ax.legend()
    fig.tight_layout()
    fig.savefig(fn + ".pdf")
    fig.savefig(fn + ".jpeg", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    m = pd.read_csv(os.path.join(TBL, "step07_tcga_signature_matrix.tsv"), sep="\t")
    # aggregate to case level (mean TPM across aliquots)
    genes = [g for g in TAMIS + MICRO if g in m.columns]
    g_mean = m.groupby("case")[genes].mean()
    z = (g_mean - g_mean.mean()) / g_mean.std()
    z.columns = genes
    score = pd.DataFrame({
        "tamis": z[[g for g in TAMIS if g in z.columns]].mean(axis=1),
        "micro": z[[g for g in MICRO if g in z.columns]].mean(axis=1),
    })
    score.index.name = "case"
    score = score.reset_index()

    cl = clinical(set(m["case"].dropna()))
    df = score.merge(cl, on="case", how="inner")
    df["os"] = df["dtd"].fillna(df["dtf"])
    df["event"] = (df["vital"] == "Dead").astype(int)
    df = df[df["os"].notna() & (df["os"] > 0)].copy()
    df["age_y"] = pd.to_numeric(df["age"], errors="coerce") / 365.25
    df["grade_num"] = df["grade"].astype(str).str.extract(r"G(\d)").astype(float)
    df["tamis_high"] = (df["tamis"] > df["tamis"].median()).astype(int)
    print(f"n cases with OS: {len(df)}; events: {df['event'].sum()}", flush=True)

    # log-rank
    grp = df["tamis_high"].values
    print("vital_status:", df["vital"].value_counts().to_dict(), flush=True)
    chi2, p = survdiff(df["os"].values, df["event"].values, grp)
    print(f"log-rank (TAM-IS high vs low): chi2={chi2:.2f}, p={p:.3g}", flush=True)

    # Cox multivariable (use available covariates; grade may be missing in GDC)
    print("non-null: age_y=%d grade_num=%d" % (df["age_y"].notna().sum(),
                                               df["grade_num"].notna().sum()), flush=True)
    cov_cols = ["tamis"] + [c for c in ["age_y", "grade_num"] if df[c].notna().sum() > 50]
    cox = df.dropna(subset=cov_cols)
    print(f"Cox n={len(cox)} covariates={cov_cols}", flush=True)
    X = sm.add_constant(cox[cov_cols])
    res = sm.PHReg(cox["os"], X, status=cox["event"]).fit()
    print(res.summary(), flush=True)
    coef = pd.DataFrame({"covariate": X.columns, "coef": res.params,
                         "HR": np.exp(res.params), "p": res.pvalues})
    coef.to_csv(os.path.join(TBL, "step07_tcga_survival_cox.tsv"), sep="\t", index=False)
    print(coef.to_string(index=False), flush=True)

    km_plot(df, "tamis_high", os.path.join(FIG, "step07_tcga_km"))
    print("DONE", flush=True)
