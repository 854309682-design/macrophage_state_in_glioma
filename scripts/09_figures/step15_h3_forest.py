"""
Step 15 — H3 forest (Fig 5), rebuilt from the UNBIASED TCGA numbers.

The original Fig 5 (step10_figures.fig5) was drawn from `step07_stratified_stateovercount.tsv`,
whose TCGA rows came from the GDC `cases` endpoint — which has no IDH/MGMT and, critically, no
`days_to_last_follow_up` for GBM, so the GBM Cox was dead-only and the "merged" HR was an
LGG-dominated subset (see handoff §5). This script redraws Fig 5 from:
  - results/tables/step07_tcga_idh_mgmt_stratified.tsv  (cBioPortal-based, unbiased; step07f)
  - results/tables/step07_cgga_cox.tsv                  (CGGA multivariable, IDH/MGMT-adjusted)
95% CIs are recovered from the reported HR + two-sided p (SE = |ln HR| / z, z = Phi^-1(1-p/2)),
so no re-fit is needed.

Outputs: results/figures/main/fig05_h3_forest.{pdf,jpeg}

Run: env/.venv/bin/python scripts/09_figures/step15_h3_forest.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm

plt.rcParams["axes.grid"] = False
plt.rcParams["figure.dpi"] = 100

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
TBL = os.path.join(ROOT, "results/tables")
FIG = os.path.join(ROOT, "results/figures/main")
os.makedirs(FIG, exist_ok=True)

# (display label, source, selector)  — selector picks one row from the table
ROWS = [
    ("TCGA all (age+grade)",              "tcga", "TCGA_all_age_grade"),
    ("TCGA all (+IDH)",                   "tcga", "TCGA_all_age_grade_idh"),
    ("TCGA all (+IDH+MGMT)",              "tcga", "TCGA_all_age_grade_idh_mgmt"),
    ("TCGA GBM/G4 (age+IDH)",             "tcga", "TCGA_GBM_G4_age_idh"),
    ("TCGA LGG/G2-G3 (age+IDH)",          "tcga", "TCGA_LGG_G2G3_age_idh"),
    ("TCGA state-over-count (+myeloid)",  "tcga", "TCGA_all_state_over_count_idh"),
    ("CGGA-325 (age+grade+IDH+MGMT)",     "cgga", "CGGA_325"),
    ("CGGA-693 (age+grade+IDH+MGMT)",     "cgga", "CGGA_693"),
]


def save(fig, name):
    for ext in ("pdf", "jpeg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name}")


def lookup(tcga, cgga, surv, source, key):
    if source == "tcga":
        r = tcga[(tcga["dataset"] == key) & (tcga["covariate"] == "tamis")].iloc[0]
        return float(r["HR"]), float(r["p"]), int(r["n"]), int(r["events"])
    r = cgga[(cgga["cohort"] == key) & (cgga["model"] == "multivariable")
             & (cgga["covariate"] == "tamis")].iloc[0]
    s = surv[surv["cohort"] == key].iloc[0]     # n / events live in the survival table
    return float(r["HR"]), float(r["p"]), int(s["n"]), int(s["events"])


def main():
    tcga = pd.read_csv(os.path.join(TBL, "step07_tcga_idh_mgmt_stratified.tsv"), sep="\t")
    cgga = pd.read_csv(os.path.join(TBL, "step07_cgga_cox.tsv"), sep="\t")
    surv = pd.read_csv(os.path.join(TBL, "step07_cgga_survival.tsv"), sep="\t")

    rows = []
    for label, source, key in ROWS:
        hr, p, n, ev = lookup(tcga, cgga, surv, source, key)
        if p <= 0 or p >= 1:            # cannot derive SE from a degenerate p
            p = min(max(p, 1e-300), 1 - 1e-16)
        z = norm.isf(p / 2.0)
        se = abs(np.log(hr)) / z if z > 0 else np.nan
        lo, hi = np.exp(np.log(hr) - 1.96 * se), np.exp(np.log(hr) + 1.96 * se)
        rows.append({"label": label, "HR": hr, "lo": lo, "hi": hi, "p": p, "n": n, "events": ev})
        print(f"  {label:36s} HR={hr:.3f} [{lo:.3f},{hi:.3f}] p={p:.3g} (n={n}, ev={ev})", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(ROOT, "results/tables", "step15_h3_forest.tsv"), sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    y = np.arange(len(df))[::-1]
    for yi, (_, r) in zip(y, df.iterrows()):
        sig = r["p"] < 0.05
        col = "#B2182B" if (r["HR"] > 1 and sig) else ("#2166AC" if (r["HR"] < 1 and sig) else "#999999")
        ax.errorbar(r["HR"], yi, xerr=[[r["HR"] - r["lo"]], [r["hi"] - r["HR"]]],
                    fmt="o", color=col, ecolor=col, capsize=3, ms=7, lw=1.5)
        ax.text(5.9, yi, f"{r['HR']:.2f} [{r['lo']:.2f}-{r['hi']:.2f}]  p={r['p']:.1e}"
                         f"  n={r['n']}", va="center", ha="right", fontsize=7.5,
                color="black" if sig else "#777777")
    ax.axvline(1.0, color="black", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(df["label"], fontsize=8.5)
    ax.set_xlabel("TAM-IS hazard ratio (95% CI, log scale)")
    ax.set_xlim(0.5, 6.2)
    ax.set_title("Fig 5 | TAM-IS prognostic value is context-dependent\n"
                 "(unbiased TCGA covariates; CGGA IDH/MGMT-adjusted)")
    fig.tight_layout()
    save(fig, "fig05_h3_forest")


if __name__ == "__main__":
    try:
        main()
        print("DONE", flush=True)
    except Exception as e:  # noqa: BLE001
        print(f"FAILED: {type(e).__name__}: {e}", flush=True)
        raise
