"""
Step 20 — H1 rigor with an HVG-independent TAM-IS definition.

step18 (step18_h1_rigor.py) tested whether the TAM-IS vs myeloid_other classifier
transfers across cohorts (within-cohort 5-fold CV AUROC ~0.95; cross-cohort transfer
~0.70).  But BOTH sides of that test were tied to the 3,000-HVG selection:

  * labels   = myeloid subclusters found on the scVI latent, which was built from the
               3,000 HVGs (step04_integrate -> step05_tamis_discovery -> step05b), and
  * features = the same 3,000-HVG expression matrix (step05_myeloid_reclustered.X).

This step removes that tie and re-runs the SAME protocol on the full 17,526-gene space
of objects/step04_myeloid.h5ad (the common-gene universe, log-normalised).

Labels (all HVG-free, derived on full genes):
  L1  cell-level marker rule  : TAM-IS = (tamis_score > 0) & (micro_score < 0);
                                microglia = (micro_score > 0) & (tamis_score < 0);
                                myeloid_other = the rest.  (matches step06 H1 prevalence)
      Rationale: the signature score carries a strong cohort offset on this object
      (mean tamis-micro: gse103224 0.01 vs gse131928/gse163120 ~1.16), so an absolute
      threshold gives gse103224 zero TAM-IS cells; the sign rule is cohort-robust.
  L2  full-gene re-clustering : PCA on all 17,526 genes (NO HVG ranking) -> neighbours
                                -> leiden 1.0, then per-cluster mean (tamis-micro) > 1.0
                                marks TAM-IS / < 0 marks microglia (the step05b rule).
  L3  calibration-free        : per-cohort tertiles of (tamis-micro) - TAM-IS = top tertile,
                                microglia = bottom tertile, myeloid_other = middle.  Removes the
                                cohort score-offset, to separate a genuine program difference from
                                gse103224's calibration artifact.

Features:
  F1  all 17,526 genes.
  F2  F1 minus the 20 signature genes                        (leakage control - mandatory,
      because the L1 label is a function of those genes; a ceiling AUROC on F1 is circular).
  F3  genes detected in >= 5 % of cells in EVERY cohort       (HVG-free detection filter).

Protocol identical to step18:
  binary task = TAM-IS vs myeloid_other (microglia dropped);
  StandardScaler + LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced");
  (a) within-cohort StratifiedKFold(5) AUROC, (b) cross-cohort train->test both ways.
  The step18 result is carried in as the baseline row set (label_def=step18_3kHVG).

Outputs:
  results/tables/step20_h1_hvg_independent.tsv   (AUROC grid, long format)
  results/tables/step20_label_summary.tsv        (per-cohort state counts per label def)

Run: env/.venv/bin/python scripts/04_myeloid_states/step20_h1_hvg_independent.py
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
OBJ = os.path.join(ROOT, "objects")
TBL = os.path.join(ROOT, "results/tables")

sc.settings.verbosity = 1

COHORTS = ["gse103224", "gse131928", "gse163120"]
TAMIS = ["SPP1", "TREM2", "APOE", "LGALS3", "CD163", "C1QA", "C1QB", "LYZ",
         "FTL", "FTH1", "IL1B", "CCL3", "TGFB1", "GPNMB"]
MICRO = ["P2RY12", "TMEM119", "CX3CR1", "SALL1", "C3", "CSF1R"]
SIG = sorted(set(TAMIS) | set(MICRO))

TAMIS_THRESH = 1.0        # cluster-mean (tamis - micro) threshold, as in step05b
DETECT_MIN = 0.05         # F3: gene must be detected in >= 5 % of cells in every cohort
MIN_CELLS = 50
MIN_POS = 20

LABEL_DEFS = ["L1_cell_marker", "L2_fullgene_cluster", "L3_quantile"]
FEATURE_SETS = ["F1_allgenes", "F2_no_signature", "F3_detection"]


def dense(x):
    return np.asarray(x.todense() if hasattr(x, "todense") else x, dtype=np.float32)


def model():
    return LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced", n_jobs=-1)


def eval_combo(label_def, feature_set, X, y, coh, rows):
    """Run within-cohort 5-fold CV + cross-cohort transfer for one (label, feature) combo."""
    # ---- (a) within-cohort 5-fold CV ----
    for c in COHORTS:
        m = coh == c
        if m.sum() < MIN_CELLS * 2 or min(int(y[m].sum()), int((~y[m].astype(bool)).sum())) < MIN_CELLS:
            print(f"  [skip CV] {label_def}/{feature_set} {c}: n={int(m.sum())}, pos={int(y[m].sum())}",
                  flush=True)
            continue
        aucs = []
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=42).split(X[m], y[m]):
            sc_ = StandardScaler().fit(X[m][tr])
            clf = model().fit(sc_.transform(X[m][tr]), y[m][tr])
            aucs.append(roc_auc_score(y[m][te], clf.predict_proba(sc_.transform(X[m][te]))[:, 1]))
        rows.append({"label_def": label_def, "feature_set": feature_set, "kind": "within_cohort_CV",
                     "train": c, "test": f"{c} (5-fold CV)", "n_test": int(m.sum()),
                     "n_pos": int(y[m].sum()), "auroc": round(float(np.mean(aucs)), 3),
                     "auroc_sd": round(float(np.std(aucs)), 3)})
        print(f"  CV {label_def}/{feature_set} {c}: AUROC {np.mean(aucs):.3f} +/- {np.std(aucs):.3f}",
              flush=True)

    # ---- (b) cross-cohort transfer ----
    for tr_c in COHORTS:
        mtr = coh == tr_c
        if mtr.sum() < MIN_CELLS * 2 or min(int(y[mtr].sum()), int((~y[mtr].astype(bool)).sum())) < MIN_CELLS:
            continue
        sc_ = StandardScaler().fit(X[mtr])
        clf = model().fit(sc_.transform(X[mtr]), y[mtr])
        for te_c in COHORTS:
            if te_c == tr_c:
                continue
            mte = coh == te_c
            if mte.sum() < MIN_CELLS or min(int(y[mte].sum()), int((~y[mte].astype(bool)).sum())) < MIN_POS:
                print(f"  [skip transfer] {label_def}/{feature_set} {tr_c}->{te_c}: insufficient",
                      flush=True)
                continue
            auc = roc_auc_score(y[mte], clf.predict_proba(sc_.transform(X[mte]))[:, 1])
            rows.append({"label_def": label_def, "feature_set": feature_set,
                         "kind": "cross_cohort_transfer", "train": tr_c, "test": te_c,
                         "n_test": int(mte.sum()), "n_pos": int(y[mte].sum()),
                         "auroc": round(float(auc), 3), "auroc_sd": np.nan})
            print(f"  transfer {label_def}/{feature_set} {tr_c} -> {te_c}: AUROC {auc:.3f}", flush=True)


if __name__ == "__main__":
    a = ad.read_h5ad(os.path.join(OBJ, "step04_myeloid.h5ad"))
    genes = np.asarray(a.var_names)
    coh = a.obs["cohort"].values.astype(str)
    print(f"myeloid object: {a.shape}; cohorts {dict(pd.Series(coh).value_counts())}", flush=True)

    # ---- full-gene signature scores (the same 14 / 6 gene sets as step05) ----
    sc.tl.score_genes(a, TAMIS, score_name="tamis_score", use_raw=False)
    sc.tl.score_genes(a, MICRO, score_name="micro_score", use_raw=False)
    ts = a.obs["tamis_score"].values
    ms = a.obs["micro_score"].values
    delta = ts - ms

    Xfull = dense(a.X)                     # 51,402 x 17,526, log-normalised
    print(f"X {Xfull.shape}; signature genes present: "
          f"{sum(g in set(genes) for g in SIG)}/{len(SIG)}", flush=True)

    # ---- label definition L1: cell-level sign rule ----
    tamis1 = (ts > 0) & (ms < 0)
    micro1 = (ms > 0) & (ts < 0)
    other1 = ~(tamis1 | micro1)
    lab1 = np.where(tamis1, "TAM-IS", np.where(micro1, "microglia", "myeloid_other"))

    # ---- label definition L2: full-gene (HVG-free) re-clustering ----
    cl = ad.AnnData(X=Xfull.copy(), obs=pd.DataFrame({"cohort": coh}, index=a.obs_names))
    sc.pp.pca(cl, n_comps=50, svd_solver="randomized", random_state=42)
    sc.pp.neighbors(cl, n_neighbors=15, use_rep="X_pca")
    sc.tl.leiden(cl, resolution=1.0, flavor="igraph", n_iterations=2, key_added="leiden_full")
    cl_mean = pd.Series(delta, index=a.obs_names).groupby(cl.obs["leiden_full"].values).mean()
    tamis_cl = set(cl_mean[cl_mean > TAMIS_THRESH].index)
    micro_cl = set(cl_mean[cl_mean < 0].index)
    leiden_full = cl.obs["leiden_full"].values.astype(str)
    lab2 = np.where(np.isin(leiden_full, list(tamis_cl)), "TAM-IS",
                    np.where(np.isin(leiden_full, list(micro_cl)), "microglia", "myeloid_other"))
    print(f"L2: {len(tamis_cl)} TAM-IS clusters, {len(micro_cl)} microglia clusters "
          f"(of {cl_mean.size})", flush=True)

    # ---- label definition L3: calibration-free within-cohort tertiles ----
    # (separates a genuine program difference from gse103224's score-offset artifact)
    q = pd.DataFrame({"cohort": coh, "delta": delta}).groupby("cohort", observed=True)[
        "delta"].quantile([1 / 3, 2 / 3]).unstack()
    lab3 = np.empty(len(delta), dtype=object)
    for c in COHORTS:
        m = coh == c
        lab3[m] = np.where(delta[m] >= q.loc[c, 2 / 3], "TAM-IS",
                           np.where(delta[m] < q.loc[c, 1 / 3], "microglia", "myeloid_other"))
    print("L3: per-cohort tertiles of (tamis - micro) [calibration-free]", flush=True)

    labels = {"L1_cell_marker": lab1, "L2_fullgene_cluster": lab2, "L3_quantile": lab3}

    # ---- label summary ----
    summ = []
    for ld, lab in labels.items():
        for c in COHORTS:
            m = coh == c
            for st in ["TAM-IS", "microglia", "myeloid_other"]:
                summ.append({"label_def": ld, "cohort": c, "state": st, "n": int((lab[m] == st).sum())})
    summ = pd.DataFrame(summ).pivot_table(index=["label_def", "cohort"], columns="state",
                                          values="n", fill_value=0).reset_index()
    summ.to_csv(os.path.join(TBL, "step20_label_summary.tsv"), sep="\t", index=False)
    print("\n=== label summary ===\n" + summ.to_string(index=False), flush=True)

    # ---- feature-space masks ----
    sig_idx = np.array([i for i, g in enumerate(genes) if g in set(SIG)])
    det = np.ones(Xfull.shape[1], dtype=bool)
    for c in COHORTS:
        m = coh == c
        det &= (Xfull[m] > 0).mean(axis=0) >= DETECT_MIN
    feat_masks = {
        "F1_allgenes": np.ones(Xfull.shape[1], dtype=bool),
        "F2_no_signature": np.ones(Xfull.shape[1], dtype=bool),
        "F3_detection": det,
    }
    feat_masks["F2_no_signature"][sig_idx] = False
    for k, v in feat_masks.items():
        print(f"  feature space {k}: {int(v.sum())} genes", flush=True)

    # ---- run the grid ----
    rows = []
    for ld in LABEL_DEFS:
        lab = labels[ld]
        keep = lab != "microglia"                    # binary task: TAM-IS vs myeloid_other
        y = (lab == "TAM-IS").astype(int)
        for fs in FEATURE_SETS:
            print(f"\n--- {ld} / {fs} ---", flush=True)
            eval_combo(ld, fs, Xfull[np.ix_(keep, feat_masks[fs])], y[keep], coh[keep], rows)

    out = pd.DataFrame(rows)
    # carry in the step18 3,000-HVG baseline for direct comparison
    try:
        base = pd.read_csv(os.path.join(TBL, "step18_h1_rigor.tsv"), sep="\t")
        base = base.assign(label_def="step18_3kHVG", feature_set="F0_3kHVG")
        out = pd.concat([base[out.columns], out], ignore_index=True)
    except FileNotFoundError:
        print("  [warn] step18_h1_rigor.tsv not found; baseline row set omitted", flush=True)

    out.to_csv(os.path.join(TBL, "step20_h1_hvg_independent.tsv"), sep="\t", index=False)
    print("\n=== H1: HVG-independent vs step18 baseline ===", flush=True)
    print(out.to_string(index=False), flush=True)
    print("DONE", flush=True)
