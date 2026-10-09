"""
Step 19 — Fig 3 (H1): per-cohort TAM-IS prevalence + cross-cohort transfer AUROC.

Replaces the old Fig 3 so it presents the *rigorous* H1 evidence (step18 classifier transfer), not
just prevalence and not the ambiguous scArches result (step11_arches_h1).

Panel A: % TAM-IS-high of myeloid per cohort (step06_h1_replication.tsv).
Panel B: TAM-IS vs myeloid_other classifier AUROC — within-cohort 5-fold CV and cross-cohort
         transfer (step18_h1_rigor.tsv), with the 0.5 chance line.

Outputs: results/figures/main/fig03_h1_replication.{pdf,jpeg}
Run: env/.venv/bin/python scripts/09_figures/step19_h1_fig.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams["axes.grid"] = False
plt.rcParams["figure.dpi"] = 100

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
TBL = os.path.join(ROOT, "results/tables")
FIG = os.path.join(ROOT, "results/figures/main")
os.makedirs(FIG, exist_ok=True)


def save(fig, name):
    for ext in ("pdf", "jpeg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name}", flush=True)


if __name__ == "__main__":
    h1 = pd.read_csv(os.path.join(TBL, "step06_h1_replication.tsv"), sep="\t")
    rig = pd.read_csv(os.path.join(TBL, "step18_h1_rigor.tsv"), sep="\t")

    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 4.8))

    # ---- Panel A: per-cohort prevalence ----
    a0.bar(h1["cohort"], h1["pct_tamis_high_of_myeloid"], color="steelblue")
    for i, v in enumerate(h1["pct_tamis_high_of_myeloid"]):
        a0.text(i, v + 0.6, f"{v:.1f}%", ha="center", fontsize=8.5)
    a0.set_ylabel("% TAM-IS-high of myeloid")
    a0.set_title("A  Per-cohort TAM-IS prevalence")
    a0.set_ylim(0, max(h1["pct_tamis_high_of_myeloid"]) * 1.25)
    a0.tick_params(axis="x", labelrotation=30)

    # ---- Panel B: classifier AUROC (step18) ----
    lab, auc, col = [], [], []
    order = {"within_cohort_CV": 0, "cross_cohort_transfer": 1}
    for _, r in rig.sort_values("kind", key=lambda s: s.map(order)).iterrows():
        if r["kind"] == "within_cohort_CV":
            lab.append(f"{r['train']} (5-fold CV)"); col.append("steelblue")
        else:
            lab.append(f"{r['train']} \u2192 {r['test']}"); col.append("darkorange")
        auc.append(r["auroc"])
    a1.barh(range(len(auc))[::-1], auc, color=col)
    for yi, v in zip(range(len(auc))[::-1], auc):
        a1.text(v + 0.01, yi, f"{v:.3f}", va="center", fontsize=8.5)
    a1.axvline(0.5, color="black", ls="--", lw=1)
    a1.set_yticks(range(len(auc))[::-1])
    a1.set_yticklabels(lab, fontsize=8.5)
    a1.set_xlim(0.4, 1.03)
    a1.set_xlabel("AUROC (TAM-IS vs myeloid_other)")
    a1.set_title("B  Reproducibility (classifier transfer)")

    fig.suptitle("Fig 3 | H1: TAM-IS recurs across cohorts (blue = within-cohort; orange = cross-cohort)",
                 y=1.02)
    fig.tight_layout()
    save(fig, "fig03_h1_replication")
    print("DONE", flush=True)
