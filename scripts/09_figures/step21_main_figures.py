"""
Step 21 — main-figure panels as importable draw functions + 27 standalone files.

Each panel is a `pNN(ax)` function that draws onto a supplied matplotlib Axes (no figure
creation, no saving) so the same code can be re-used by step22 to compose combined figures.
The standalone writer runs each panel into its own figure and saves
`results/figures/main/fig01..fig27_*.{pdf,jpeg}`.

Evidence chain (one plot per file):
  fig01 QC per cohort · fig02/03 integrated UMAP (cohort/lineage) · fig04-06 myeloid UMAP
  (state / TAM-IS score / microglia score) · fig07/08 marker dotplot & volcano · fig09-13 H1
  · fig14-16 H2 IvyGAP · fig17-19 KM · fig20/21 H3 forest & meta · fig22/23 deconvolution
  · fig24/25/26 CellChat / NicheNet / pySCENIC · fig27 druggability.

Run (all):        env/.venv/bin/python scripts/09_figures/step21_main_figures.py
Run (subset):     env/.venv/bin/python scripts/09_figures/step21_main_figures.py 03 18 19
"""
import io
import os
import sys
import zipfile

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
SCENIC = os.path.join(ROOT, "results/scenic")
FIG = os.path.join(ROOT, "results/figures/main")
DATASETS = os.environ.get("GLIOMA_DATA") or os.path.join(os.path.dirname(ROOT), "Datasets")
CGGA = os.path.join(DATASETS, "CGGA")
os.makedirs(FIG, exist_ok=True)

C_TAMIS, C_MICRO = "magma", "viridis"
STATE_ORDER = ["TAM-IS", "myeloid_other", "microglia"]

CANON = ["SPP1", "TREM2", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "LGALS3",
         "CD163", "GPNMB", "LPL", "P2RY12", "TMEM119", "CX3CR1", "SALL1", "CSF1R"]

KM_GENES = {
    "TAMIS": ["SPP1", "TREM2", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "LGALS3",
              "CD163", "CTSL", "CCL3", "CCL4", "GPNMB", "LPL", "HLA-DRA", "HLA-DRB1",
              "S100A9", "HMOX1", "FTL", "FTH1"],
}


def umap_scatter(ax, um, c, title, cmap=None, cats=None, legend=False):
    if cats is None:
        sc0 = ax.scatter(um[:, 0], um[:, 1], c=c, cmap=cmap, s=2, linewidths=0)
        plt.colorbar(sc0, ax=ax, shrink=0.8)
    else:
        cats = np.asarray(cats)
        for v in pd.unique(cats):
            m = cats == v
            ax.scatter(um[m, 0], um[m, 1], s=2, label=str(v), alpha=0.7, linewidths=0)
        if legend:
            ax.legend(markerscale=4, fontsize=8, frameon=False)
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")


def km_curve(ax, df, title, xlabel):
    from statsmodels.duration.survfunc import SurvfuncRight
    for grp, sub in df.groupby("tamis_high"):
        sf = SurvfuncRight(sub["OS"].values, sub["event"].values)
        ax.step(sf.surv_times, sf.surv_prob, where="post",
                label=f"{'high' if grp else 'low'} (n={len(sub)})")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Survival probability")
    ax.set_title(title)
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, fontsize=8)


def forest(ax, df, title):
    df = df.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(df))
    ax.errorbar(df["HR"], y, xerr=[df["HR"] - df["lo"], df["hi"] - df["HR"]],
                fmt="o", color="black", capsize=2)
    ax.axvline(1.0, ls="--", color="grey", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(df["label"], fontsize=8)
    for i, r in df.iterrows():
        ax.text(r["hi"] * 1.02, i, f"{r['HR']:.2f} ({r['lo']:.2f}-{r['hi']:.2f}), p={r['p']:.1e}",
                va="center", fontsize=6.5)
    ax.set_xlim(0, max(df["hi"]) * 1.6)
    ax.set_xlabel("Hazard ratio")
    ax.set_title(title)


def _deconv_blocks():
    txt = open(os.path.join(TBL, "step17_deconv_consistency.tsv")).read().strip()
    blocks = [b for b in txt.split("\n\n") if b.strip()]
    out = []
    for b in blocks:
        rows = [r for r in b.splitlines() if not r.startswith("#")]
        out.append(pd.read_csv(io.StringIO("\n".join(rows)), sep="\t"))
    return out


def _cgga_df(tag):
    def read_zip(pat, index_col):
        z = os.path.join(CGGA, pat)
        with zipfile.ZipFile(z) as zf:
            name = [n for n in zf.namelist() if n.endswith(".txt") and "__MACOSX" not in n][0]
            with zf.open(name) as fh:
                return pd.read_csv(fh, sep="\t", index_col=index_col)
    expr = read_zip(f"CGGA.mRNAseq_{tag}.RSEM-genes.zip", 0)
    expr.index = [str(g).strip("'\"") for g in expr.index]
    cl = read_zip(f"CGGA.mRNAseq_{tag}_clinical.zip", None).set_index("CGGA_ID")
    genes = [g for g in KM_GENES["TAMIS"] if g in expr.index]
    sub = expr.loc[genes]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0)
    df = pd.DataFrame({"tamis": z.mean(axis=0)}).join(cl, how="inner")
    df["OS"] = pd.to_numeric(df["OS"], errors="coerce")
    df["event"] = pd.to_numeric(df["Censor (alive=0; dead=1)"], errors="coerce")
    df = df[df["OS"].notna() & df["event"].notna() & (df["OS"] > 0)].copy()
    df["tamis_high"] = (df["tamis"] > df["tamis"].median()).astype(int)
    return df


# ---------------- panel draw functions (ax only, no saving) ----------------
def p01(ax):
    d = pd.read_csv(os.path.join(TBL, "step03_qc_summary.tsv"), sep="\t")
    ax.bar(d["cohort"], d["cells_passed"], color="steelblue")
    for i, r in d.iterrows():
        ax.text(i, r["cells_passed"] + 400,
                f"{int(r['cells_passed']):,}\n{int(r['genes']):,} genes", ha="center", fontsize=7)
    ax.set_ylabel("Cells passing QC")
    ax.set_title("Cells passing QC per cohort")
    ax.set_ylim(0, d["cells_passed"].max() * 1.28)


def p02(ax):
    import anndata as ad
    a = ad.read_h5ad(os.path.join(OBJ, "step04_integrated.h5ad"))
    umap_scatter(ax, a.obsm["X_umap"], None, "Integrated scRNA by cohort",
                 cats=a.obs["cohort"].astype(str), legend=True)


def p03(ax):
    import anndata as ad
    a = ad.read_h5ad(os.path.join(OBJ, "step04_integrated.h5ad"))
    ann = ad.read_h5ad(os.path.join(OBJ, "step04_annotated.h5ad"), backed="r")
    lin = ann.obs["lineage"].reindex(a.obs_names).astype(str).values
    umap_scatter(ax, a.obsm["X_umap"], None, "Integrated scRNA by lineage", cats=lin, legend=True)


def p04(ax):
    import anndata as ad
    a = ad.read_h5ad(os.path.join(OBJ, "step05_myeloid_reclustered.h5ad"))
    umap_scatter(ax, a.obsm["X_umap"], None, "Myeloid states",
                 cats=a.obs["state"].astype(str), legend=True)


def p05(ax):
    import anndata as ad
    a = ad.read_h5ad(os.path.join(OBJ, "step05_myeloid_reclustered.h5ad"))
    umap_scatter(ax, a.obsm["X_umap"], a.obs["tamis_score"].values, "TAM-IS score", cmap=C_TAMIS)


def p06(ax):
    import anndata as ad
    a = ad.read_h5ad(os.path.join(OBJ, "step05_myeloid_reclustered.h5ad"))
    umap_scatter(ax, a.obsm["X_umap"], a.obs["micro_score"].values, "Microglia score", cmap=C_MICRO)


def p07(ax):
    import anndata as ad
    st = ad.read_h5ad(os.path.join(OBJ, "step05_myeloid_reclustered.h5ad"), backed="r")
    state = st.obs["state"].astype(str)
    am = ad.read_h5ad(os.path.join(OBJ, "step04_myeloid.h5ad"), backed="r")
    genes = [g for g in CANON if g in am.var_names]
    X = np.asarray(am[:, genes].to_memory().X)
    state = state.reindex(am.obs_names).values
    mean, frac = {}, {}
    for s in STATE_ORDER:
        m = state == s
        mean[s] = X[m].mean(axis=0)
        frac[s] = (X[m] > 0).mean(axis=0)
    mean = np.vstack([mean[s] for s in STATE_ORDER])
    frac = np.vstack([frac[s] for s in STATE_ORDER])
    mean = (mean - mean.min(axis=0)) / (mean.max(axis=0) - mean.min(axis=0) + 1e-9)
    xg, yg = np.meshgrid(np.arange(len(genes)), np.arange(len(STATE_ORDER)))
    sc0 = ax.scatter(xg.ravel(), yg.ravel(), c=mean.ravel(), s=frac.ravel() * 200 + 10,
                     cmap="Reds", edgecolors="grey", linewidths=0.4)
    plt.colorbar(sc0, ax=ax, shrink=0.8, label="scaled mean expr")
    ax.set_xticks(np.arange(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=8)
    ax.set_yticks(np.arange(len(STATE_ORDER)))
    ax.set_yticklabels(STATE_ORDER)
    ax.set_title("Canonical marker expression by state")


def p08(ax):
    d = pd.read_csv(os.path.join(TBL, "step05_tamis_markers.tsv"), sep="\t")
    d["nlogp"] = -np.log10(d["pvals_adj"].clip(lower=1e-300))
    sig = d["pvals_adj"] < 0.05
    ax.scatter(d.loc[~sig, "logfoldchanges"], d.loc[~sig, "nlogp"], s=6, color="lightgrey")
    ax.scatter(d.loc[sig, "logfoldchanges"], d.loc[sig, "nlogp"], s=8, color="crimson")
    for _, r in d[sig].nlargest(12, "scores").iterrows():
        ax.annotate(r["names"], (r["logfoldchanges"], r["nlogp"]), fontsize=6)
    ax.axvline(0, ls="--", color="grey", lw=0.8)
    ax.set_xlabel("log2 fold change (TAM-IS vs microglia)")
    ax.set_ylabel("-log10 adj. p")
    ax.set_title("TAM-IS marker volcano")


def p09(ax):
    d = pd.read_csv(os.path.join(TBL, "step06_h1_replication.tsv"), sep="\t")
    ax.bar(d["cohort"], d["pct_tamis_high_of_myeloid"], color="steelblue")
    for i, v in enumerate(d["pct_tamis_high_of_myeloid"]):
        ax.text(i, v + 0.6, f"{v:.1f}%", ha="center", fontsize=8)
    ax.set_ylabel("% TAM-IS-high of myeloid")
    ax.set_title("TAM-IS prevalence per cohort")
    ax.set_ylim(0, d["pct_tamis_high_of_myeloid"].max() * 1.25)
    ax.tick_params(axis="x", labelrotation=30)


def p10(ax):
    import anndata as ad
    a = ad.read_h5ad(os.path.join(OBJ, "step05_myeloid_reclustered.h5ad"), backed="r")
    d = a.obs[["cohort", "tamis_score"]].copy()
    order = sorted(d["cohort"].unique())
    ax.violinplot([d.loc[d["cohort"] == c, "tamis_score"].values for c in order], showextrema=False)
    ax.set_xticks(np.arange(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=20, fontsize=8)
    ax.set_ylabel("TAM-IS score")
    ax.set_title("TAM-IS score distribution per cohort")


def p11(ax):
    d = pd.read_csv(os.path.join(TBL, "step18_h1_rigor.tsv"), sep="\t")
    lab, auc, col = [], [], []
    for _, r in d.iterrows():
        if r["kind"] == "within_cohort_CV":
            lab.append(f"{r['train']} (5-fold CV)"); col.append("steelblue")
        else:
            lab.append(f"{r['train']} \u2192 {r['test']}"); col.append("darkorange")
        auc.append(r["auroc"])
    ax.barh(np.arange(len(auc))[::-1], auc, color=col)
    for yi, v in zip(np.arange(len(auc))[::-1], auc):
        ax.text(v + 0.01, yi, f"{v:.3f}", va="center", fontsize=8)
    ax.axvline(0.5, color="black", ls="--", lw=1)
    ax.set_yticks(np.arange(len(auc))[::-1])
    ax.set_yticklabels(lab, fontsize=8)
    ax.set_xlim(0.4, 1.03)
    ax.set_xlabel("AUROC (TAM-IS vs myeloid_other)")
    ax.set_title("H1 reproducibility (step18, 3k-HVG)")


def p12(ax):
    d = pd.read_csv(os.path.join(TBL, "step20_h1_hvg_independent.tsv"), sep="\t")
    t = d[d["kind"] == "cross_cohort_transfer"].copy()
    t["pair"] = t["train"].str.replace("gse", "") + "\u2192" + t["test"].str.replace("gse", "")
    piv = t.pivot_table(index=["label_def", "feature_set"], columns="pair", values="auroc")
    im = ax.imshow(piv.values, cmap="RdYlBu_r", vmin=0.4, vmax=1.0, aspect="auto")
    plt.colorbar(im, ax=ax, shrink=0.8, label="transfer AUROC")
    ax.set_xticks(np.arange(piv.shape[1]))
    ax.set_xticklabels(piv.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(np.arange(piv.shape[0]))
    ax.set_yticklabels([f"{a} / {b}" for a, b in piv.index], fontsize=6.5)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6)
    ax.set_title("H1 HVG-independent re-test (step20)")


def p13(ax):
    d = pd.read_csv(os.path.join(TBL, "step20_label_summary.tsv"), sep="\t")
    d = d[d["label_def"] == "L1_cell_marker"].set_index("cohort")
    bottom = np.zeros(len(d))
    for s, col in zip(STATE_ORDER, ["crimson", "steelblue", "seagreen"]):
        ax.bar(d.index, d[s], bottom=bottom, label=s, color=col)
        bottom += d[s].values
    ax.set_ylabel("Cells")
    ax.set_title("Myeloid state composition per cohort")
    ax.legend(frameon=False, fontsize=8)
    ax.tick_params(axis="x", labelrotation=20)


def p14(ax):
    d = pd.read_csv(os.path.join(TBL, "step06_ivygap_regions.tsv"), sep="\t")
    order = [r for r in ["LE", "IT", "CT", "CTmvp", "CTpan"] if r in d["region"].unique()]
    ax.boxplot([d.loc[d["region"] == r, "tamis"].values for r in order],
               tick_labels=order, showfliers=False)
    ax.set_xlabel("IvyGAP region")
    ax.set_ylabel("TAM-IS score")
    ax.set_title("H2: TAM-IS by IvyGAP region")


def p15(ax):
    d = pd.read_csv(os.path.join(TBL, "step06_ivygap_regions.tsv"), sep="\t")
    order = [r for r in ["LE", "IT", "CT", "CTmvp", "CTpan"] if r in d["region"].unique()]
    ax.boxplot([d.loc[d["region"] == r, "hypoxia"].values for r in order],
               tick_labels=order, showfliers=False)
    ax.set_xlabel("IvyGAP region")
    ax.set_ylabel("Hypoxia score")
    ax.set_title("H2: hypoxia by IvyGAP region")


def p16(ax):
    from scipy.stats import spearmanr
    d = pd.read_csv(os.path.join(TBL, "step06_ivygap_regions.tsv"), sep="\t")
    rho, p = spearmanr(d["tamis"], d["hypoxia"])
    ax.scatter(d["hypoxia"], d["tamis"], s=10, color="steelblue", alpha=0.7)
    ax.set_xlabel("Hypoxia score")
    ax.set_ylabel("TAM-IS score")
    ax.set_title(f"H2: TAM-IS vs hypoxia (rho={rho:.2f}, p={p:.1e})")


def p17(ax):
    sm = pd.read_csv(os.path.join(TBL, "step07_tcga_signature_matrix.tsv"), sep="\t")
    cl = pd.read_csv(os.path.join(DATASETS, "TCGA/TCGA_lgggbm_clinical.tsv"), sep="\t")
    genes = [g for g in KM_GENES["TAMIS"] if g in sm.columns]
    z = sm[genes].apply(lambda c: (c - c.mean()) / c.std(), axis=0)
    sm = sm.assign(tamis=z.mean(axis=1)).groupby("case", as_index=False)["tamis"].mean()
    df = sm.merge(cl, on="case", how="inner")
    df["OS"] = pd.to_numeric(df["os_months"], errors="coerce")
    df["event"] = df["os_status"].astype(str).str.startswith("1").astype(int)
    df = df[df["OS"].notna() & (df["OS"] > 0)].copy()
    df["tamis_high"] = (df["tamis"] > df["tamis"].median()).astype(int)
    km_curve(ax, df, "TCGA: TAM-IS high vs low", "Overall survival (months)")


def p18(ax):
    df = _cgga_df("325")
    km_curve(ax, df, f"CGGA-325: TAM-IS high vs low (n={len(df)})", "Overall survival (days)")


def p19(ax):
    df = _cgga_df("693")
    km_curve(ax, df, f"CGGA-693: TAM-IS high vs low (n={len(df)})", "Overall survival (days)")


def p20(ax):
    d = pd.read_csv(os.path.join(TBL, "step15_h3_forest.tsv"), sep="\t")
    forest(ax, d, "H3: TAM-IS HR (stratified)")


def p21(ax):
    d = pd.read_csv(os.path.join(TBL, "step16_meta_analysis.tsv"), sep="\t")
    d = d.rename(columns={"study": "label"})
    d["lo"] = d["HR"] * np.exp(-1.96 * d["se_logHR"])
    d["hi"] = d["HR"] * np.exp(1.96 * d["se_logHR"])
    forest(ax, d[["label", "HR", "lo", "hi", "p"]], "H3: HR by study & stratum (meta input)")


def p22(ax):
    d = _deconv_blocks()[0]
    d["pair"] = d["method_a"] + " vs " + d["method_b"]
    x = np.arange(len(d))
    ax.bar(x, d["spearman_rho"], color="steelblue")
    for i, v in enumerate(d["spearman_rho"]):
        ax.text(i, v + 0.015, f"{v:.2f}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{a}\n{b}" for a, b in zip(d["cohort"], d["pair"])], fontsize=7)
    ax.set_ylabel("Spearman rho (monocyte_mac fraction)")
    ax.set_ylim(0, 1.05)
    ax.set_title("Deconvolution method agreement")


def p23(ax):
    d = _deconv_blocks()[1]
    d["pair"] = d["cohort"] + " / " + d["method"]
    y = np.arange(len(d))[::-1]
    col = ["darkorange" if m == "NNLS" else "steelblue" for m in d["method"]]
    ax.barh(y, d["HR_per_SD"], color=col)
    for yi, v in zip(y, d["HR_per_SD"]):
        ax.text(v + 0.01, yi, f"{v:.2f}", va="center", fontsize=8)
    ax.axvline(1.0, ls="--", color="black", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(d["pair"], fontsize=8)
    ax.set_xlabel("HR per SD of monocyte_mac fraction")
    ax.set_title("Macrophage-fraction HR (per SD), by method")


def p24(ax):
    d = pd.read_csv(os.path.join(TBL, "step09_cellchat_interactions.csv"))
    d = d[(d["source"] == "monocyte_mac") | (d["target"] == "monocyte_mac")]
    d = d.sort_values("prob", ascending=False).head(20).iloc[::-1]
    lab = [f"{r.ligand}\u2192{r.receptor} ({r.source}\u2192{r.target})" for r in d.itertuples()]
    ax.barh(np.arange(len(d)), d["prob"], color="mediumpurple")
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels(lab, fontsize=6.5)
    ax.set_xlabel("Communication probability")
    ax.set_title("CellChat: top L-R interactions with Mo-Mac")


def p25(ax):
    d = pd.read_csv(os.path.join(TBL, "step09d_nichenet_ligand_activities.tsv"), sep="\t")
    d = d.sort_values("pearson", ascending=False).head(20).iloc[::-1]
    ax.barh(np.arange(len(d)), d["pearson"], color="teal")
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels(d["test_ligand"], fontsize=8)
    ax.set_xlabel("Pearson ligand-activity score")
    ax.set_title("NicheNet: top ligands for the TAM-IS program")


def p26(ax):
    d = pd.read_csv(os.path.join(SCENIC, "step09j_state_regulon_auc.tsv"), sep="\t").set_index("state")
    d = d.reindex([s for s in STATE_ORDER if s in d.index])
    z = (d - d.mean(axis=1).values[:, None]) / d.std(axis=1).values[:, None]
    top = z.loc["TAM-IS"].sort_values(ascending=False).head(25).index
    piv = z[top]
    v = np.nanmax(np.abs(piv.values))
    im = ax.imshow(piv.values, cmap="RdBu_r", aspect="auto", vmin=-v, vmax=v)
    plt.colorbar(im, ax=ax, shrink=0.8, label="z(mean AUC)")
    ax.set_xticks(np.arange(len(top)))
    ax.set_xticklabels(top, rotation=90, fontsize=6.5)
    ax.set_yticks(np.arange(piv.shape[0]))
    ax.set_yticklabels(piv.index, fontsize=8)
    ax.set_title("TAM-IS top regulons (pySCENIC)")


def p27(ax):
    d = pd.read_csv(os.path.join(TBL, "step08_druggability.tsv"), sep="\t").sort_values("n_drugs_dgidb")
    ax.barh(d["gene"], d["n_drugs_dgidb"], color="darkorange")
    for i, v in enumerate(d["n_drugs_dgidb"]):
        ax.text(v + 0.3, i, str(v), va="center", fontsize=8)
    ax.set_xlabel("# known drugs (DGIdb)")
    ax.set_title("Druggability of TAM-IS candidate genes")


# name -> (number, draw function, standalone figsize)
PANELS = {
    "fig01_qc_cells_passed":        ("fig01", p01, (6, 4.5)),
    "fig02_umap_cohort":            ("fig02", p02, (6, 5.5)),
    "fig03_umap_lineage":           ("fig03", p03, (6, 5.5)),
    "fig04_umap_state":             ("fig04", p04, (6, 5.5)),
    "fig05_umap_tamis_score":       ("fig05", p05, (6, 5.5)),
    "fig06_umap_microglia_score":   ("fig06", p06, (6, 5.5)),
    "fig07_marker_dotplot":         ("fig07", p07, (10, 3.2)),
    "fig08_marker_volcano":         ("fig08", p08, (6.5, 5)),
    "fig09_prevalence":             ("fig09", p09, (6, 4.5)),
    "fig10_score_violin":           ("fig10", p10, (6, 4.5)),
    "fig11_step18_auroc":           ("fig11", p11, (6.5, 4.5)),
    "fig12_step20_hvg_heatmap":     ("fig12", p12, (10, 4.5)),
    "fig13_state_composition":      ("fig13", p13, (6.5, 4.5)),
    "fig14_ivygap_tamis":           ("fig14", p14, (6.5, 4.5)),
    "fig15_ivygap_hypoxia":         ("fig15", p15, (6.5, 4.5)),
    "fig16_ivygap_colocalization":  ("fig16", p16, (6, 5)),
    "fig17_km_tcga":                ("fig17", p17, (6, 5)),
    "fig18_km_cgga325":             ("fig18", p18, (6, 5)),
    "fig19_km_cgga693":             ("fig19", p19, (6, 5)),
    "fig20_forest_stratified":      ("fig20", p20, (8.5, 5)),
    "fig21_forest_meta":            ("fig21", p21, (8.5, 6)),
    "fig22_deconv_agreement":       ("fig22", p22, (6.5, 4.5)),
    "fig23_deconv_hr":              ("fig23", p23, (6.5, 4.5)),
    "fig24_cellchat":               ("fig24", p24, (7.5, 6)),
    "fig25_nichenet":               ("fig25", p25, (6.5, 6)),
    "fig26_scenic_regulon_heatmap": ("fig26", p26, (9, 2.2)),
    "fig27_druggability":           ("fig27", p27, (6.5, 5)),
}


def save(fig, name):
    for ext in ("pdf", "jpeg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name}", flush=True)


def make_single(name, number, fn, size):
    fig, ax = plt.subplots(figsize=size)
    fn(ax)
    ax.set_title(f"{number} | {ax.get_title()}")
    save(fig, name)


if __name__ == "__main__":
    want = [w.zfill(2) for w in sys.argv[1:]]
    for name, (number, fn, size) in PANELS.items():
        if want and number[-2:] not in want:
            continue
        try:
            print(f"== {number} ==", flush=True)
            make_single(name, number, fn, size)
        except Exception as e:
            print(f"  {number} FAILED: {type(e).__name__} {e}", flush=True)
    print("DONE", flush=True)
