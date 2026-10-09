"""
Step 7c — CGGA-325 / CGGA-693 survival analysis (H3 validation) + state-over-count test.

For each cohort:
  - read RSEM expression (genes x samples) + clinical (OS, Censor, Age, Grade, IDH, MGMT),
  - score TAM-IS (immunosuppressive state), microglia, and a broad myeloid-content score,
  - Kaplan-Meier (TAM-IS high vs low) + log-rank,
  - multivariable Cox: TAM-IS adjusted for age, grade, IDH, MGMT,
  - STATE-OVER-COUNT: Cox with TAM-IS + myeloid_content -> does the state remain prognostic
    after accounting for total myeloid abundance?

Outputs: results/tables/step07_cgga_survival.tsv, results/tables/step07_cgga_cox.tsv,
         results/figures/main/step07_cgga_{325,693}_km.{pdf,jpeg}

Run: env/.venv/bin/python scripts/07_bulk_validation/step07c_cgga_survival.py
"""
import os
import zipfile

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.duration.survfunc import SurvfuncRight, survdiff

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
DATASETS = os.environ.get("GLIOMA_DATA") or os.path.join(os.path.dirname(ROOT), "Datasets")
CGGA = os.path.join(DATASETS, "CGGA")
TBL = os.path.join(ROOT, "results/tables")
FIG = os.path.join(ROOT, "results/figures/main")
os.makedirs(FIG, exist_ok=True)

GENES = {
    "TAMIS": ["SPP1", "TREM2", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "LGALS3",
              "CD163", "CTSL", "CCL3", "CCL4", "GPNMB", "LPL", "HLA-DRA", "HLA-DRB1",
              "S100A9", "HMOX1", "FTL", "FTH1"],
    "MICRO": ["P2RY12", "TMEM119", "CX3CR1", "SALL1", "CRYBB1", "CST3", "PROS1", "MERTK"],
    "MYELOID_CONTENT": ["PTPRC", "ITGAM", "CSF1R", "AIF1", "TYROBP", "FCER1G", "LYZ",
                        "CD68", "CD14"],
}
ALL = sorted({g for v in GENES.values() for g in v})

COHORTS = {"CGGA_325": "325", "CGGA_693": "693"}


def read_zip_matrix(tag):
    z = os.path.join(CGGA, f"CGGA.mRNAseq_{tag}.RSEM-genes.zip")
    with zipfile.ZipFile(z) as zf:
        name = [n for n in zf.namelist() if n.endswith(".txt") and "__MACOSX" not in n][0]
        with zf.open(name) as fh:
            df = pd.read_csv(fh, sep="\t", index_col=0)
    df.index = [str(g).strip("'\"") for g in df.index]
    return df


def read_zip_clinical(tag):
    z = os.path.join(CGGA, f"CGGA.mRNAseq_{tag}_clinical.zip")
    with zipfile.ZipFile(z) as zf:
        name = [n for n in zf.namelist() if n.endswith(".txt") and "__MACOSX" not in n][0]
        with zf.open(name) as fh:
            df = pd.read_csv(fh, sep="\t")
    return df


def km_plot(df, fn, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 5))
    for grp, sub in df.groupby("tamis_high"):
        sf = SurvfuncRight(sub["OS"], sub["event"])
        ax.step(sf.surv_times, sf.surv_prob, where="post", label=f"{'high' if grp else 'low'} (n={len(sub)})")
    ax.set_xlabel("Overall survival (days)")
    ax.set_ylabel("Survival probability")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(fn + ".pdf")
    fig.savefig(fn + ".jpeg", dpi=300)
    plt.close(fig)


summary = []
cox_rows = []
for label, tag in COHORTS.items():
    expr = read_zip_matrix(tag)          # genes x samples
    cl = read_zip_clinical(tag)
    cl = cl.set_index("CGGA_ID")
    genes = [g for g in ALL if g in expr.index]
    sub = expr.loc[genes]
    # z-score per gene across samples
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0)
    score = pd.DataFrame({
        "tamis": z.loc[[g for g in GENES["TAMIS"] if g in z.index]].mean(axis=0),
        "micro": z.loc[[g for g in GENES["MICRO"] if g in z.index]].mean(axis=0),
        "myeloid": z.loc[[g for g in GENES["MYELOID_CONTENT"] if g in z.index]].mean(axis=0),
    })
    df = score.join(cl, how="inner")
    df["OS"] = pd.to_numeric(df["OS"], errors="coerce")
    df["event"] = pd.to_numeric(df["Censor (alive=0; dead=1)"], errors="coerce")
    df = df[df["OS"].notna() & df["event"].notna() & (df["OS"] > 0)].copy()
    df["age"] = pd.to_numeric(df["Age"], errors="coerce")
    df["grade_num"] = df["Grade"].astype(str).str.extract(r"(IV|III|II)")
    df["grade_num"] = df["grade_num"].map({"II": 2, "III": 3, "IV": 4})
    df["idh"] = (df["IDH_mutation_status"].astype(str).str.lower().str.contains("mut")).astype(float)
    df["mgmt"] = (df["MGMTp_methylation_status"].astype(str).str.lower().str.contains("methylated")).astype(float)
    df["tamis_high"] = (df["tamis"] > df["tamis"].median()).astype(int)
    print(f"\n=== {label}: n={len(df)}, events={int(df['event'].sum())} ===", flush=True)

    chi2, p = survdiff(df["OS"].values, df["event"].values, df["tamis_high"].values)
    print(f"  log-rank TAM-IS high vs low: chi2={chi2:.2f}, p={p:.3g}", flush=True)
    summary.append({"cohort": label, "n": len(df), "events": int(df["event"].sum()),
                    "logrank_chi2": round(chi2, 2), "logrank_p": p})

    cov = ["tamis", "age", "grade_num", "idh", "mgmt"]
    c = df.dropna(subset=cov)
    X = sm.add_constant(c[cov])
    res = sm.PHReg(c["OS"], X, status=c["event"].astype(int), ties="efron").fit()
    for i, v in enumerate(X.columns):
        cox_rows.append({"cohort": label, "model": "multivariable", "covariate": v,
                         "HR": float(np.exp(res.params[i])), "p": float(res.pvalues[i])})
    # fair comparison with TCGA: age + grade only
    c1 = df.dropna(subset=["tamis", "age", "grade_num"])
    X1 = sm.add_constant(c1[["tamis", "age", "grade_num"]])
    res1 = sm.PHReg(c1["OS"], X1, status=c1["event"].astype(int), ties="efron").fit()
    for i, v in enumerate(X1.columns):
        cox_rows.append({"cohort": label, "model": "age_grade_only", "covariate": v,
                         "HR": float(np.exp(res1.params[i])), "p": float(res1.pvalues[i])})
    # state-over-count
    c2 = df.dropna(subset=["tamis", "myeloid"])
    X2 = sm.add_constant(c2[["tamis", "myeloid"]])
    res2 = sm.PHReg(c2["OS"], X2, status=c2["event"].astype(int), ties="efron").fit()
    for i, v in enumerate(X2.columns):
        cox_rows.append({"cohort": label, "model": "state_over_count", "covariate": v,
                         "HR": float(np.exp(res2.params[i])), "p": float(res2.pvalues[i])})
    print(pd.DataFrame(cox_rows).query("cohort==@label").to_string(index=False), flush=True)
    km_plot(df, os.path.join(FIG, f"step07_cgga_{tag}_km"), f"CGGA-{tag}: TAM-IS high vs low")

pd.DataFrame(summary).to_csv(os.path.join(TBL, "step07_cgga_survival.tsv"), sep="\t", index=False)
pd.DataFrame(cox_rows).to_csv(os.path.join(TBL, "step07_cgga_cox.tsv"), sep="\t", index=False)
print("\nDONE", flush=True)
