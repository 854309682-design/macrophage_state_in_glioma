"""
Step 9 — Assemble main figures Fig 1-6 from the REAL results.

Fig1 overview (integrated UMAP by cohort + cluster) · Fig2 TAM-IS discovery (myeloid UMAP +
scores + marker dotplot) · Fig3 H1 replication · Fig4 H2 IvyGAP niche · Fig5 H3 survival
forest · Fig6 druggability.

Outputs: results/figures/main/fig01..fig06_{name}.{pdf,jpeg}

Run: env/.venv/bin/python scripts/09_figures/step10_figures.py
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
OBJ = os.path.join(ROOT, "objects")
TBL = os.path.join(ROOT, "results/tables")
FIG = os.path.join(ROOT, "results/figures/main")
os.makedirs(FIG, exist_ok=True)


def save(fig, name):
    for ext in ("pdf", "jpeg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name}")


# ---------------- Fig 1: integrated overview ----------------
def fig1():
    import anndata as ad
    a = ad.read_h5ad(os.path.join(OBJ, "step04_integrated.h5ad"))
    um = a.obsm["X_umap"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    # cohort panel
    cats = a.obs["cohort"].astype(str)
    for c in pd.unique(cats):
        m = (cats == c).values
        axes[0].scatter(um[m, 0], um[m, 1], s=2, label=str(c), alpha=0.7, linewidths=0)
    axes[0].set_title("cohort")
    axes[0].legend(markerscale=4, fontsize=8, frameon=False)
    # leiden panel: color by cluster + label centroids (no 48-entry legend)
    led = a.obs["leiden"].astype(str)
    codes = led.astype("category").cat.codes.values
    axes[1].scatter(um[:, 0], um[:, 1], c=codes, cmap="tab20", s=2, linewidths=0)
    for cl in pd.unique(led):
        m = (led == cl).values
        axes[1].text(um[m, 0].mean(), um[m, 1].mean(), str(cl), fontsize=6,
                     ha="center", va="center", color="black")
    axes[1].set_title("Leiden cluster")
    for ax in axes:
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    fig.suptitle("Fig 1 | Integrated scRNA cohorts (scVI)")
    save(fig, "fig01_overview")


# ---------------- Fig 2: TAM-IS discovery ----------------
def fig2():
    import anndata as ad
    a = ad.read_h5ad(os.path.join(OBJ, "step05_myeloid_reclustered.h5ad"))
    um = a.obsm["X_umap"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    sc0 = axes[0].scatter(um[:, 0], um[:, 1], c=a.obs["tamis_score"], s=3, cmap="magma", linewidths=0)
    axes[0].set_title("TAM-IS score")
    plt.colorbar(sc0, ax=axes[0], shrink=0.8)
    sc1 = axes[1].scatter(um[:, 0], um[:, 1], c=a.obs["micro_score"], s=3, cmap="viridis", linewidths=0)
    axes[1].set_title("Microglia score")
    plt.colorbar(sc1, ax=axes[1], shrink=0.8)
    for ax in axes:
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    fig.suptitle("Fig 2 | Myeloid compartment: TAM-IS vs microglia")
    save(fig, "fig02_tamis_discovery")


# ---------------- Fig 3: H1 replication ----------------
def fig3():
    h1 = pd.read_csv(os.path.join(TBL, "step06_h1_replication.tsv"), sep="\t")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(h1["cohort"], h1["pct_tamis_high_of_myeloid"], color="steelblue")
    for i, v in enumerate(h1["pct_tamis_high_of_myeloid"]):
        ax.text(i, v + 0.5, f"{v:.1f}%", ha="center", fontsize=9)
    ax.set_ylabel("% TAM-IS-high of myeloid")
    ax.set_title("Fig 3 | H1 replication across 5 scRNA cohorts")
    ax.set_ylim(0, max(h1["pct_tamis_high_of_myeloid"]) * 1.2)
    save(fig, "fig03_h1_replication")


# ---------------- Fig 4: H2 IvyGAP niche ----------------
def fig4():
    d = pd.read_csv(os.path.join(TBL, "step06_ivygap_regions.tsv"), sep="\t")
    order = [r for r in ["LE", "IT", "CT", "CTmvp", "CTpan"] if r in d["region"].unique()]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, var, lab in zip(axes, ["tamis", "hypoxia"], ["TAM-IS score", "Hypoxia score"]):
        data = [d.loc[d["region"] == r, var].values for r in order]
        ax.boxplot(data, tick_labels=order, showfliers=False)
        ax.set_title(lab); ax.set_ylabel(lab); ax.set_xlabel("IvyGAP region")
    fig.suptitle("Fig 4 | H2: TAM-IS enriched in pseudopalisading-necrosis/hypoxic niche (IvyGAP)")
    save(fig, "fig04_h2_ivygap")


# ---------------- Fig 5: H3 forest ----------------
def fig5():
    rows = []
    tc = pd.read_csv(os.path.join(TBL, "step07_stratified_stateovercount.tsv"), sep="\t")
    for _, r in tc.iterrows():
        if r["covariate"] == "tamis":
            rows.append((f"{r['dataset']} ({r['model']})", r["HR"], r["p"], r["n"]))
    df = pd.DataFrame(rows, columns=["stratum", "HR", "p", "n"])
    df = df[df["HR"].notna()].iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 0.5 * len(df) + 2))
    y = np.arange(len(df))
    ax.errorbar(df["HR"], y, xerr=None, fmt="o", color="black")
    ax.axvline(1.0, ls="--", color="grey")
    ax.set_yticks(y); ax.set_yticklabels(df["stratum"], fontsize=8)
    for i, (_, r) in enumerate(df.iterrows()):
        ax.text(r["HR"] * 1.05, i, f"HR={r['HR']:.2f}, p={r['p']:.1e}, n={int(r['n'])}", va="center", fontsize=7)
    ax.set_xlabel("TAM-IS hazard ratio (per SD)")
    ax.set_title("Fig 5 | H3: TAM-IS prognostic value is context-dependent")
    save(fig, "fig05_h3_forest")


# ---------------- Fig 6: druggability ----------------
def fig6():
    d = pd.read_csv(os.path.join(TBL, "step08_druggability.tsv"), sep="\t")
    d = d.sort_values("n_drugs_dgidb", ascending=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(d["gene"], d["n_drugs_dgidb"], color="darkorange")
    for i, v in enumerate(d["n_drugs_dgidb"]):
        ax.text(v + 0.5, i, str(v), va="center", fontsize=8)
    ax.set_xlabel("# known drugs (DGIdb)")
    ax.set_title("Fig 6 | Druggability of TAM-IS candidates (DGIdb)")
    save(fig, "fig06_druggability")


if __name__ == "__main__":
    for fn in [fig1, fig2, fig3, fig4, fig5, fig6]:
        try:
            print(f"== {fn.__name__} ==", flush=True)
            fn()
        except Exception as e:
            print(f"  {fn.__name__} FAILED: {type(e).__name__} {e}", flush=True)
    print("DONE", flush=True)
