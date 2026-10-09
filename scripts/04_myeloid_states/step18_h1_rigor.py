"""
Step 18 — H1 rigor: is TAM-IS a reproducible, transferable state across cohorts?

The handoff flags H1 as "suggestive, not rigorous" (the scArches result assigned 96.8 % of query
cells to TAM-IS — too high to trust). Here we test reproducibility with a simple, honest
supervised criterion instead:

  * binary task = TAM-IS vs myeloid_other (the two main states; microglia excluded to avoid the
    resident-vs-monocyte confound),
  * features = the 3,000 HVG expression space,
  * (a) within-cohort 5-fold CV AUROC, and
  * (b) cross-cohort transfer: train on one cohort, test on another (both directions).

If the TAM-IS expression program is cohort-independent, a classifier trained on one cohort should
separate TAM-IS in an unseen cohort well above chance.

Outputs: results/tables/step18_h1_rigor.tsv
Run: env/.venv/bin/python scripts/04_myeloid_states/step18_h1_rigor.py
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
OBJ = os.path.join(ROOT, "objects")
TBL = os.path.join(ROOT, "results/tables")

COHORTS = ["gse131928", "gse163120", "gse103224"]
MIN_CELLS = 50


def dense(mat):
    return np.asarray(mat.todense()) if hasattr(mat, "todense") else np.asarray(mat)


def model():
    return LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced", n_jobs=-1)


if __name__ == "__main__":
    a = ad.read_h5ad(os.path.join(OBJ, "step05_myeloid_reclustered.h5ad"))
    keep = a.obs["state"].isin(["TAM-IS", "myeloid_other"]).values
    X = dense(a[keep].X).astype(np.float32)
    y = (a.obs["state"].values[keep] == "TAM-IS").astype(int)
    coh = a.obs["cohort"].values[keep]
    print(f"cells: {X.shape[0]} x {X.shape[1]} genes; TAM-IS={int(y.sum())}", flush=True)

    rows = []
    # ---- (a) within-cohort 5-fold CV AUROC ----
    for c in COHORTS:
        m = coh == c
        if m.sum() < MIN_CELLS * 2 or y[m].min(initial=0) == y[m].max():
            print(f"  [skip CV] {c}: n={int(m.sum())}, TAM-IS={int(y[m].sum())}", flush=True)
            continue
        if min(y[m].sum(), (~y[m].astype(bool)).sum()) < MIN_CELLS:
            print(f"  [skip CV] {c}: minority class too small", flush=True)
            continue
        aucs = []
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=42).split(X[m], y[m]):
            sc = StandardScaler().fit(X[m][tr])
            clf = model().fit(sc.transform(X[m][tr]), y[m][tr])
            aucs.append(roc_auc_score(y[m][te], clf.predict_proba(sc.transform(X[m][te]))[:, 1]))
        rows.append({"test": f"{c} (5-fold CV)", "train": c, "n_test": int(m.sum()),
                     "n_pos": int(y[m].sum()), "auroc": round(float(np.mean(aucs)), 3),
                     "auroc_sd": round(float(np.std(aucs)), 3), "kind": "within_cohort_CV"})
        print(f"  CV {c}: AUROC {np.mean(aucs):.3f} +/- {np.std(aucs):.3f}", flush=True)

    # ---- (b) cross-cohort transfer ----
    for tr_c in COHORTS:
        mtr = coh == tr_c
        if mtr.sum() < MIN_CELLS * 2 or min(y[mtr].sum(), (~y[mtr].astype(bool)).sum()) < MIN_CELLS:
            print(f"  [skip train] {tr_c}: insufficient", flush=True)
            continue
        sc = StandardScaler().fit(X[mtr])
        clf = model().fit(sc.transform(X[mtr]), y[mtr])
        for te_c in COHORTS:
            if te_c == tr_c:
                continue
            mte = coh == te_c
            if mte.sum() < MIN_CELLS or min(y[mte].sum(), (~y[mte].astype(bool)).sum()) < 20:
                print(f"  [skip test] {tr_c}->{te_c}: insufficient", flush=True)
                continue
            auc = roc_auc_score(y[mte], clf.predict_proba(sc.transform(X[mte]))[:, 1])
            rows.append({"test": te_c, "train": tr_c, "n_test": int(mte.sum()),
                         "n_pos": int(y[mte].sum()), "auroc": round(float(auc), 3),
                         "auroc_sd": np.nan, "kind": "cross_cohort_transfer"})
            print(f"  transfer {tr_c} -> {te_c}: AUROC {auc:.3f}", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(TBL, "step18_h1_rigor.tsv"), sep="\t", index=False)
    print("\n=== H1 rigor (AUROC) ===", flush=True)
    print(out.to_string(index=False), flush=True)
    print("DONE", flush=True)
