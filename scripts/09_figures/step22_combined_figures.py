"""
Step 22 — compose the 27 main panels into 5 combined figures + 8 supplementary panels into 2.

Re-uses the panel draw functions from step21 (`pNN(ax)`) and places them into GridSpec layouts
with A/B/C... panel letters, producing true vector combined figures:

  figure1_atlas          | 5 panels — data & integrated myeloid atlas
  figure2_definition     | 5 panels — defining the TAM-IS state
  figure3_repro_spatial  | 6 panels — H1 reproducibility & H2 spatial niche
  figure4_clinical       | 5 panels — H3 clinical value
  figure5_mechanism      | 6 panels — deconvolution, mechanism & druggability
  figureS1_validation    | 4 panels — technical & malignant-label validation
  figureS2_stats         | 4 panels — supporting statistics

Outputs: results/figures/combined/*.pdf + *.jpeg
Run: env/.venv/bin/python scripts/09_figures/step22_combined_figures.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import step21_main_figures as s21

ROOT = s21.ROOT
TBL = s21.TBL
SCENIC = s21.SCENIC
OBJ = s21.OBJ
COMBINED = os.path.join(ROOT, "results/figures/combined")
os.makedirs(COMBINED, exist_ok=True)


def save(fig, name):
    for ext in ("pdf", "jpeg"):
        fig.savefig(os.path.join(COMBINED, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name}", flush=True)


# ---------------- supplementary panel draw functions ----------------
def s_qc_dist(ax):
    import anndata as ad
    parts, labels = [], []
    for c in ["gse103224", "gse131928", "gse163120"]:
        a = ad.read_h5ad(os.path.join(OBJ, f"step03_{c}.h5ad"), backed="r")
        parts.append(a.obs["n_genes_by_counts"].values)
        labels.append(c)
    ax.boxplot(parts, tick_labels=labels, showfliers=False)
    ax.set_ylabel("genes / cell")
    ax.set_title("Per-cohort QC: genes per cell")


def s_infercnv_img(ax):
    img = plt.imread(os.path.join(ROOT, "results/figures/supplementary/step03b_infercnv_heatmap.jpeg"))
    ax.imshow(img)
    ax.axis("off")
    ax.set_title("inferCNV heatmap (malignant vs reference)", fontsize=9)


def s_infercnv_auroc(ax):
    d = pd.read_csv(os.path.join(TBL, "step03b_infercnv_auroc.tsv"), sep="\t")
    ax.bar(d["group"], d["auroc"], color="slateblue")
    for i, v in enumerate(d["auroc"]):
        ax.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=8)
    ax.axhline(0.5, ls="--", color="grey")
    ax.set_ylim(0.4, 1.05)
    ax.set_ylabel("AUROC (cnv score \u2192 malignant)")
    ax.set_title("inferCNV AUROC")


def s_deconv_mean(ax):
    d = s21._deconv_blocks()[1]
    lab = d["cohort"].str.replace("CGGA_", "") + " / " + d["method"]
    ax.bar(np.arange(len(d)), d["mac_frac_mean"],
           color=["darkorange" if m == "NNLS" else "steelblue" for m in d["method"]])
    ax.set_xticks(np.arange(len(d)))
    ax.set_xticklabels(lab, rotation=30, fontsize=7)
    ax.set_ylabel("mean monocyte_mac fraction")
    ax.set_title("Deconvolution: mean macrophage fraction by method")


def s_cellchat_paths(ax):
    d = pd.read_csv(os.path.join(TBL, "step09_cellchat_summary.csv")).head(15).iloc[::-1]
    ax.barh(np.arange(len(d)), d["prob"], color="mediumpurple")
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels(d["pathway"], fontsize=8)
    ax.set_xlabel("total communication probability")
    ax.set_title("CellChat: top pathways")


def s_nichenet_lt(ax):
    d = pd.read_csv(os.path.join(TBL, "step09d_nichenet_ligand_target.tsv"), sep="\t")
    lig = pd.read_csv(os.path.join(TBL, "step09d_nichenet_ligand_activities.tsv"), sep="\t")
    lig = lig.sort_values("pearson", ascending=False).head(8)["test_ligand"].tolist()
    d = d[d["ligand"].isin(lig)]
    top_t = d.groupby("target")["weight"].max().sort_values(ascending=False).head(15).index
    piv = d.pivot_table(index="target", columns="ligand", values="weight").reindex(index=top_t, columns=lig)
    im = ax.imshow(piv.values, cmap="viridis", aspect="auto")
    plt.colorbar(im, ax=ax, shrink=0.8, label="regulatory weight")
    ax.set_xticks(np.arange(len(lig)))
    ax.set_xticklabels(lig, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(top_t)))
    ax.set_yticklabels(top_t, fontsize=6.5)
    ax.set_title("NicheNet: ligand \u2192 target weights")


def s_scenic_all(ax):
    d = pd.read_csv(os.path.join(SCENIC, "step09j_state_regulon_auc.tsv"), sep="\t").set_index("state")
    d = d.reindex([s for s in ["TAM-IS", "myeloid_other", "microglia"] if s in d.index])
    z = (d - d.mean(axis=1).values[:, None]) / d.std(axis=1).values[:, None]
    top = z.var(axis=0).sort_values(ascending=False).head(30).index
    piv = z[top]
    v = np.nanmax(np.abs(piv.values))
    im = ax.imshow(piv.values, cmap="RdBu_r", aspect="auto", vmin=-v, vmax=v)
    plt.colorbar(im, ax=ax, shrink=0.8, label="z(mean AUC)")
    ax.set_xticks(np.arange(len(top)))
    ax.set_xticklabels(top, rotation=90, fontsize=6.5)
    ax.set_yticks(np.arange(piv.shape[0]))
    ax.set_yticklabels(piv.index, fontsize=8)
    ax.set_title("pySCENIC regulon activity (all states)")


def s_cgga_forest(ax):
    d = pd.read_csv(os.path.join(TBL, "step07_cgga_cox.tsv"), sep="\t")
    d = d[(d["covariate"] == "tamis") & d["p"].notna() & (d["p"] > 0)].iloc[::-1]
    y = np.arange(len(d))
    ax.barh(y, d["HR"], color="steelblue")
    for yi, (_, r) in zip(y, d.iterrows()):
        ax.text(r["HR"] + 0.02, yi, f"p={r['p']:.3f}", va="center", fontsize=7)
    ax.axvline(1.0, ls="--", color="black", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(d["cohort"].str.replace("CGGA_", "CGGA-") + " / " + d["model"], fontsize=7)
    ax.set_xlabel("TAM-IS HR")
    ax.set_title("CGGA: TAM-IS HR by model (state-over-count)")


SUPP = {
    "supp_s1a_qc_dist": s_qc_dist,
    "supp_s1b_infercnv_img": s_infercnv_img,
    "supp_s1c_infercnv_auroc": s_infercnv_auroc,
    "supp_s1d_deconv_mean": s_deconv_mean,
    "supp_s2a_cellchat_paths": s_cellchat_paths,
    "supp_s2b_nichenet_lt": s_nichenet_lt,
    "supp_s2c_scenic_all": s_scenic_all,
    "supp_s2d_cgga_forest": s_cgga_forest,
}

DRAW = {name: fn for name, (_num, fn, _size) in s21.PANELS.items()}
DRAW.update(SUPP)

LAYOUTS = {
    "figure1_atlas": (
        "Figure 1 | Study design & integrated myeloid atlas", (16, 9), (2, 3),
        [("fig02_umap_cohort", 0, 0), ("fig03_umap_lineage", 0, 1), ("fig05_umap_tamis_score", 0, 2),
         ("fig01_qc_cells_passed", 1, 0), ("fig04_umap_state", 1, 1)],
    ),
    "figure2_definition": (
        "Figure 2 | Defining the TAM-IS state", (12, 14), (3, 2),
        [("fig06_umap_microglia_score", 0, 0), ("fig08_marker_volcano", 0, 1),
         ("fig07_marker_dotplot", 1, 0, 1, 2),
         ("fig09_prevalence", 2, 0), ("fig10_score_violin", 2, 1)],
    ),
    "figure3_repro_spatial": (
        "Figure 3 | Cross-cohort reproducibility (H1) & spatial niche (H2)", (16, 9), (2, 3),
        [("fig11_step18_auroc", 0, 0), ("fig12_step20_hvg_heatmap", 0, 1), ("fig13_state_composition", 0, 2),
         ("fig14_ivygap_tamis", 1, 0), ("fig15_ivygap_hypoxia", 1, 1), ("fig16_ivygap_colocalization", 1, 2)],
    ),
    "figure4_clinical": (
        "Figure 4 | Clinical value (H3)", (16, 9), (2, 3),
        [("fig17_km_tcga", 0, 0), ("fig18_km_cgga325", 0, 1), ("fig19_km_cgga693", 0, 2),
         ("fig20_forest_stratified", 1, 0), ("fig21_forest_meta", 1, 1)],
    ),
    "figure5_mechanism": (
        "Figure 5 | Deconvolution, mechanism & druggability", (16, 9), (2, 3),
        [("fig22_deconv_agreement", 0, 0), ("fig23_deconv_hr", 0, 1), ("fig24_cellchat", 0, 2),
         ("fig25_nichenet", 1, 0), ("fig26_scenic_regulon_heatmap", 1, 1), ("fig27_druggability", 1, 2)],
    ),
    "figureS1_validation": (
        "Figure S1 | Technical & malignant-label validation", (12, 10), (2, 2),
        [("supp_s1a_qc_dist", 0, 0), ("supp_s1b_infercnv_img", 0, 1),
         ("supp_s1c_infercnv_auroc", 1, 0), ("supp_s1d_deconv_mean", 1, 1)],
    ),
    "figureS2_stats": (
        "Figure S2 | Supporting statistics", (12, 10), (2, 2),
        [("supp_s2a_cellchat_paths", 0, 0), ("supp_s2b_nichenet_lt", 0, 1),
         ("supp_s2c_scenic_all", 1, 0), ("supp_s2d_cgga_forest", 1, 1)],
    ),
}


def compose(fname, title, figsize, grid, panels):
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(grid[0], grid[1], hspace=0.55, wspace=0.38)
    for i, item in enumerate(panels):
        name, r, c = item[0], item[1], item[2]
        rs = item[3] if len(item) > 3 else 1
        cs = item[4] if len(item) > 4 else 1
        ax = fig.add_subplot(gs[r:r + rs, c:c + cs])
        DRAW[name](ax)
        ax.set_title(ax.get_title(), fontsize=9)
        ax.text(-0.14, 1.08, "ABCDEFGH"[i], transform=ax.transAxes,
                fontsize=15, fontweight="bold", va="top")
    fig.suptitle(title, fontsize=12, y=0.995)
    save(fig, fname)


if __name__ == "__main__":
    import sys
    want = sys.argv[1:]
    for fname, (title, figsize, grid, panels) in LAYOUTS.items():
        if want and fname not in want:
            continue
        try:
            print(f"== {fname} ==", flush=True)
            compose(fname, title, figsize, grid, panels)
        except Exception as e:
            print(f"  {fname} FAILED: {type(e).__name__} {e}", flush=True)
    print("DONE", flush=True)
