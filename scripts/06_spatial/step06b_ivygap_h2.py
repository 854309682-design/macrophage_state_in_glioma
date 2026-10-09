"""
Step 6 — H2 spatial/niche validation with IvyGAP (GSE107559).

IvyGAP = region-annotated bulk RNA-seq (leading edge LE, infiltrating tumor IT,
cellular tumor CT, microvascular proliferation CTmvp, pseudopalisading necrosis CTpan).
Scores TAM-IS and a hypoxia signature per sample, maps samples to regions, and tests
whether TAM-IS is enriched in hypoxia/necrosis niches (Kruskal-Wallis + pairwise) and
co-localizes with hypoxia (Spearman).

Outputs: results/tables/step06_ivygap_regions.tsv, results/figures/main/step06_ivygap_h2.{pdf,jpeg}

Run: env/.venv/bin/python scripts/06_spatial/step06b_ivygap_h2.py
"""
import os

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
DATASETS = os.environ.get("GLIOMA_DATA") or os.path.join(os.path.dirname(ROOT), "Datasets")
IVY = os.path.join(DATASETS, "IvyGAP")
TBL = os.path.join(ROOT, "results/tables")
FIG = os.path.join(ROOT, "results/figures/main")
os.makedirs(FIG, exist_ok=True)

TAMIS = ["SPP1", "TREM2", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "LGALS3", "CD163",
         "CTSL", "CCL3", "CCL4", "GPNMB", "LPL", "HLA-DRA", "HLA-DRB1", "S100A9",
         "HMOX1", "FTL", "FTH1"]
MICRO = ["P2RY12", "TMEM119", "CX3CR1", "SALL1", "CRYBB1", "CST3", "PROS1", "MERTK"]
HYPOXIA = ["VEGFA", "SLC2A1", "LDHA", "PGK1", "PDK1", "HK2", "BNIP3", "NDRG1", "ENO1",
           "PDK3", "TPI1", "P4HA1", "ADM", "MRPS17", "ANKZF1"]

REGION_ORDER = ["LE", "IT", "CT", "CTmvp", "CTpan"]


def score_mean_z(expr, genes):
    genes = [g for g in genes if g in expr.index]
    sub = expr.loc[genes]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


if __name__ == "__main__":
    # FPKM matrix: gene_id (numeric) x samples; map gene_id -> symbol via rows-genes
    expr = pd.read_csv(os.path.join(IVY, "GSE107559_ivygap_fpkm_table.csv.gz"),
                       index_col=0, compression="gzip")
    rg = pd.read_csv(os.path.join(IVY, "GSE107559_ivygap_rows-genes.csv.gz"), compression="gzip")
    id2sym = dict(zip(rg["gene_id"].astype(str), rg["gene_symbol"].astype(str)))
    expr.index = [id2sym.get(str(g), str(g)) for g in expr.index]
    expr = expr[~expr.index.duplicated(keep="first")]
    expr.columns = [str(c) for c in expr.columns]
    print(f"FPKM: {expr.shape}", flush=True)

    # annotation: rna_well_id -> structure_abbreviation
    from openpyxl import load_workbook
    wb = load_workbook(os.path.join(IVY, "GSE107559_ivygap_columns-samples.xlsx"), read_only=True)
    ws = wb[wb.sheetnames[0]]
    rows = list(ws.iter_rows(values_only=True))
    hdr = [str(h) if h else "" for h in rows[0]]
    ann = pd.DataFrame(rows[1:], columns=hdr)
    ann["rna_well_id"] = ann["rna_well_id"].astype(str)

    # derive region code from study_structure_specimen_abbreviation (e.g. ASTR_CTmvp-...)
    def region(s):
        s = str(s)
        for code in ["CTpan", "CTmvp", "CT", "LE", "IT"]:
            if f"_{code}-" in s or f"_{code}_" in s:
                return code
        return np.nan
    ann["region"] = ann["study_structure_specimen_abbreviation"].map(region)
    reg_map = ann.set_index("rna_well_id")["region"].to_dict()
    print("region counts:", ann["region"].value_counts().to_dict(), flush=True)

    scores = pd.DataFrame({
        "tamis": score_mean_z(expr, TAMIS),
        "micro": score_mean_z(expr, MICRO),
        "hypoxia": score_mean_z(expr, HYPOXIA),
    })
    scores["region"] = [reg_map.get(c) for c in scores.index]
    scores = scores.dropna(subset=["region"])
    print(f"samples with region: {len(scores)}", flush=True)

    # Kruskal-Wallis across regions for TAM-IS and hypoxia
    for var in ["tamis", "hypoxia"]:
        groups = [scores.loc[scores["region"] == r, var] for r in REGION_ORDER if (scores["region"] == r).sum() > 2]
        H, p = stats.kruskal(*groups)
        print(f"Kruskal-Wallis {var} across regions: H={H:.2f}, p={p:.3g}", flush=True)

    # Spearman TAM-IS x hypoxia
    rho, p = stats.spearmanr(scores["tamis"], scores["hypoxia"])
    print(f"Spearman TAM-IS x hypoxia: rho={rho:.3f}, p={p:.3g}", flush=True)

    # region means + pairwise vs CT (cellular tumor = reference core)
    reg_mean = scores.groupby("region")[["tamis", "hypoxia", "micro"]].agg(["mean", "count"])
    print(reg_mean.to_string(), flush=True)
    scores.to_csv(os.path.join(TBL, "step06_ivygap_regions.tsv"), sep="\t")

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    order = [r for r in REGION_ORDER if r in scores["region"].unique()]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, var, lab in zip(axes, ["tamis", "hypoxia"], ["TAM-IS score", "Hypoxia score"]):
        data = [scores.loc[scores["region"] == r, var].values for r in order]
        ax.boxplot(data, tick_labels=order, showfliers=False)
        ax.set_title(f"IvyGAP: {lab} by region")
        ax.set_ylabel(lab)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "step06_ivygap_h2.pdf"))
    fig.savefig(os.path.join(FIG, "step06_ivygap_h2.jpeg"), dpi=300)
    plt.close(fig)
    print("saved step06_ivygap_h2 figure + step06_ivygap_regions.tsv", flush=True)
    print("DONE", flush=True)
