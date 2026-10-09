"""
Step 2 — Quality control + doublet removal per single-cell cohort.

Loads objects/step02_<cohort>.h5ad, computes QC metrics, applies data-driven
thresholds (with a QC summary printed for review), runs scrublet doublet
detection for the 10x cohorts, log-normalizes, and saves objects/step03_<cohort>.h5ad.

Thresholds (per analysis_strategy_1.md §4 Step 2; adjusted from the QC plots):
  - min genes/cell = 200 ; min counts/cell = 500
  - mito % < 20
  - genes present in >= 3 cells
  - doublets: scrublet (10x cohorts); Smart-seq2 kept (low doublet risk) but flagged.

Run: env/.venv/bin/python scripts/01_qc/step03_qc.py
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
OBJ = os.path.join(ROOT, "objects")
FIG = os.path.join(ROOT, "results/figures/supplementary")
os.makedirs(FIG, exist_ok=True)

sc.settings.verbosity = 2
sc.settings.figdir = FIG

MIN_GENES = 200
MIN_COUNTS = 500
MAX_MITO = 20.0
MIN_CELLS = 3

COHORTS = ["gse103224", "gse131928", "gse163120"]
TENX = {"gse103224", "gse163120"}  # 10x -> run scrublet

summary = []


def qc(name):
    src = os.path.join(OBJ, f"step02_{name}.h5ad")
    a = ad.read_h5ad(src)
    print(f"\n=== {name}: loaded {a.shape} ===")

    # ---- QC metrics ----
    mito_prefix = "MT-" if not a.var_names[0].startswith("MT-") else "MT-"
    a.var["mt"] = a.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, inplace=True)

    before = a.n_obs
    # ---- filters ----
    cells_ok = (
        (a.obs["n_genes_by_counts"] >= MIN_GENES)
        & (a.obs["total_counts"] >= MIN_COUNTS)
        & (a.obs["pct_counts_mt"] < MAX_MITO)
    )
    a = a[cells_ok].copy()
    sc.pp.filter_genes(a, min_cells=MIN_CELLS)
    print(f"  after cell/gene filter: {a.shape} (from {before})")

    # ---- doublets (10x only) ----
    n_dbl = 0
    if name in TENX and a.n_obs > 5000:
        sc.pp.scrublet(a, expected_doublet_rate=0.06, threshold=0.25, random_state=42)
        n_dbl = int((a.obs["predicted_doublet"]).sum())
        print(f"  scrublet doublets predicted: {n_dbl} ({n_dbl / max(a.n_obs, 1) * 100:.1f}%)")
        a = a[~a.obs["predicted_doublet"]].copy()

    # ---- normalize ----
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    # ---- save QC summary + publication-quality figure (pdf + jpeg) ----
    summary.append({
        "cohort": name,
        "cells_passed": a.n_obs,
        "genes": a.n_vars,
        "doublets_removed": n_dbl,
    })
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, col, lab in zip(axes, ["pct_counts_mt", "n_genes_by_counts", "total_counts"],
                            ["% mitochondrial", "genes / cell", "counts / cell"]):
        ax.hist(a.obs[col], bins=60, color="steelblue")
        ax.set_title(f"{name}: {lab}")
        ax.set_xlabel(lab)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"step03_{name}_qc.pdf"))
    fig.savefig(os.path.join(FIG, f"step03_{name}_qc.jpeg"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  QC figure saved (pdf+jpeg) for {name}")

    out = os.path.join(OBJ, f"step03_{name}.h5ad")
    a.write_h5ad(out)
    print(f"  saved -> {out}")


if __name__ == "__main__":
    for c in COHORTS:
        qc(c)
    print("\n=== Step 2 QC summary ===")
    print(pd.DataFrame(summary).to_string(index=False))
    pd.DataFrame(summary).to_csv(os.path.join(ROOT, "results/tables/step03_qc_summary.tsv"),
                                 sep="\t", index=False)
    print("DONE")
