"""
Step 5j — AUCell regulon-activity heatmap across TAM-IS states (pySCENIC, tail).

Reads the AUCell matrix (step09h_aucell.csv, cells x regulons) + the cell metadata
(step09f_cellmeta.csv: CellID, state, cohort, leiden), collapses to per-state mean
AUC, z-scores each regulon across states, and plots the most state-variable
regulons as a heatmap. Also writes the underlying state-mean table.

Expected biology: TAM-IS-enriched regulons of the hypoxic/inflammatory programme
(HIF1A, STAT3, CEBPB/CEBPD, MAFB, NFKB-type) vs microglia/homeostatic regulons.

Outputs:
  results/scenic/step09j_state_regulon_auc.tsv        state x regulon mean AUC (raw)
  results/figures/main/fig07_tamis_regulon_auc.{pdf,jpeg}

Run: env/.venv/bin/python scripts/05_interaction/step09j_aucell_heatmap.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["axes.grid"] = False
plt.rcParams["figure.dpi"] = 100

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
SCENIC = os.path.join(ROOT, "results/scenic")
FIG = os.path.join(ROOT, "results/figures/main")
os.makedirs(FIG, exist_ok=True)

TOP_N = 40
STATE_ORDER = ["TAM-IS", "myeloid_other", "microglia"]


def save(fig, name):
    for ext in ("pdf", "jpeg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name}")


def main():
    auc = pd.read_csv(os.path.join(SCENIC, "step09h_aucell.csv"), index_col=0)
    meta = pd.read_csv(os.path.join(SCENIC, "step09f_cellmeta.csv"), index_col=0)
    print(f"AUC matrix: {auc.shape[0]} cells x {auc.shape[1]} regulons", flush=True)

    meta = meta.reindex(auc.index)
    keep = meta["state"].notna().values
    auc, meta = auc.loc[keep], meta.loc[keep]
    print(f"cells with state: {auc.shape[0]}", flush=True)

    states = [s for s in STATE_ORDER if s in set(meta["state"])]
    states += [s for s in meta["state"].unique() if s not in states]
    mean_auc = auc.groupby(meta["state"]).mean().reindex(states)
    mean_auc.index.name = "state"
    mean_auc.to_csv(os.path.join(SCENIC, "step09j_state_regulon_auc.tsv"), sep="\t")
    print(f"state-mean table: {mean_auc.shape}", flush=True)

    # pick the most state-variable regulons
    var = mean_auc.var(axis=0)
    top = var.sort_values(ascending=False).head(min(TOP_N, mean_auc.shape[1])).index
    z = mean_auc[top].apply(lambda r: (r - r.mean()) / (r.std() if r.std() else 1.0), axis=0)
    mat = z.T  # regulons x states

    fig, ax = plt.subplots(figsize=(1.1 * len(states) + 3.0, 0.28 * len(mat) + 2.0))
    im = ax.imshow(mat.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(states)))
    ax.set_xticklabels(states, rotation=30, ha="right")
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels([r.replace("(+)", "").strip() for r in mat.index], fontsize=8)
    for i in range(len(mat)):
        for j in range(len(states)):
            ax.text(j, i, f"{mat.values[i, j]:.1f}", ha="center", va="center",
                    fontsize=5, color="black")
    cb = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("regulon AUC (z across states)")
    ax.set_title("Fig 7 | TAM-IS regulon activity across myeloid states (pySCENIC AUCell)")
    ax.set_xlabel("state")
    save(fig, "fig07_tamis_regulon_auc")

    if "TAM-IS" in mean_auc.index:
        tam = mean_auc.loc["TAM-IS"].sort_values(ascending=False).head(15)
        print("top TAM-IS regulons:\n" + tam.to_string(), flush=True)


if __name__ == "__main__":
    try:
        main()
        print("DONE", flush=True)
    except Exception as e:  # noqa: BLE001
        print(f"FAILED: {type(e).__name__}: {e}", flush=True)
        raise
