"""
Step 7d — H3 refinement: grade-stratified Cox + state-over-count (TCGA and CGGA).

For TCGA-GBM/LGG and CGGA_325/693:
  - score TAM-IS, microglia, and broad myeloid-content from expression,
  - STATE-OVER-COUNT: Cox(tamis + myeloid_content) -> does the state survive adjustment?
  - GRADE-STRATIFIED: Cox(tamis) within each grade (G2/G3/G4) -> is the discordance a
    grade-composition effect?
Outputs: results/tables/step07_stratified_stateovercount.tsv

Run: env/.venv/bin/python scripts/07_bulk_validation/step07d_stratified.py
"""
import json
import os
import urllib.request
import zipfile

import numpy as np
import pandas as pd
import statsmodels.api as sm

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
DATASETS = os.environ.get("GLIOMA_DATA") or os.path.join(os.path.dirname(ROOT), "Datasets")
TBL = os.path.join(ROOT, "results/tables")
CGGA = os.path.join(DATASETS, "CGGA")

TAMIS = ["SPP1", "TREM2", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "LGALS3", "CD163",
         "CTSL", "CCL3", "CCL4", "GPNMB", "LPL", "HLA-DRA", "HLA-DRB1", "S100A9",
         "HMOX1", "FTL", "FTH1"]
MYELOID = ["PTPRC", "ITGAM", "CSF1R", "AIF1", "TYROBP", "FCER1G", "LYZ", "CD68", "CD14"]

rows = []


def cox(df, covs, tag):
    c = df.dropna(subset=covs)
    if len(c) < 30 or c["event"].sum() < 10:
        return
    X = sm.add_constant(c[covs])
    res = sm.PHReg(c["OS"], X, status=c["event"]).fit()
    for i, v in enumerate(X.columns):
        rows.append({"dataset": tag, "model": "|".join(covs), "n": len(c),
                     "events": int(c["event"].sum()), "covariate": v,
                     "HR": round(float(np.exp(res.params[i])), 3),
                     "p": float(res.pvalues[i])})


def scores_from_z(expr, genes, sample_index):
    genes = [g for g in genes if g in expr.index]
    sub = expr.loc[genes]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0).reindex(sample_index)


# ---------- TCGA ----------
def gdc_clinical():
    payload = {"filters": {"op": "in", "content": {"field": "project.project_id",
               "value": ["TCGA-GBM", "TCGA-LGG"]}},
               "fields": "submitter_id,demographic.vital_status,demographic.days_to_death,"
                         "diagnoses.days_to_last_follow_up,diagnoses.age_at_diagnosis,"
                         "diagnoses.tumor_grade", "format": "json", "size": 5000}
    req = urllib.request.Request("https://api.gdc.cancer.gov/cases",
                                 data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json"})
    d = json.load(urllib.request.urlopen(req, timeout=120))
    out = []
    for h in d["data"]["hits"]:
        dg = (h.get("diagnoses") or [{}])[0]
        out.append({"case": h["submitter_id"],
                    "vital": (h.get("demographic") or {}).get("vital_status"),
                    "dtd": (h.get("demographic") or {}).get("days_to_death"),
                    "dtf": dg.get("days_to_last_follow_up"), "age": dg.get("age_at_diagnosis"),
                    "grade": dg.get("tumor_grade")})
    return pd.DataFrame(out).drop_duplicates("case")


m = pd.read_csv(os.path.join(TBL, "step07_tcga_signature_matrix.tsv"), sep="\t")
tam = m.groupby("case")[[g for g in TAMIS if g in m.columns]].mean()
mye = m.groupby("case")[[g for g in MYELOID if g in m.columns]].mean()
tz = (tam - tam.mean()) / tam.std()
mz = (mye - mye.mean()) / mye.std()
tcga = pd.DataFrame({"case": tz.index, "tamis": tz.mean(axis=1).values,
                     "myeloid": mz.mean(axis=1).values})
tcga = tcga.merge(gdc_clinical(), on="case", how="inner")
tcga["OS"] = tcga["dtd"].fillna(tcga["dtf"])
tcga["event"] = (tcga["vital"] == "Dead").astype(int)
tcga = tcga[tcga["OS"].notna() & (tcga["OS"] > 0)].copy()
tcga["age"] = pd.to_numeric(tcga["age"], errors="coerce") / 365.25
tcga["grade_num"] = tcga["grade"].astype(str).str.extract(r"G(\d)").astype(float)
print(f"TCGA n={len(tcga)}", flush=True)

cox(tcga, ["tamis"], "TCGA")
cox(tcga, ["tamis", "myeloid"], "TCGA_state_over_count")
cox(tcga, ["tamis", "age", "grade_num"], "TCGA_age_grade")
for g in [2, 3, 4]:
    cox(tcga[tcga["grade_num"] == g], ["tamis"], f"TCGA_gradeG{g}")

# ---------- CGGA ----------
def read_cgga(tag):
    def rd(z, kind):
        with zipfile.ZipFile(z) as zf:
            n = [x for x in zf.namelist() if x.endswith(".txt") and "__MACOSX" not in x][0]
            with zf.open(n) as fh:
                return pd.read_csv(fh, sep="\t", index_col=0 if kind == "expr" else None)
    expr = rd(os.path.join(CGGA, f"CGGA.mRNAseq_{tag}.RSEM-genes.zip"), "expr")
    expr.index = [str(g).strip("'\"") for g in expr.index]
    cl = rd(os.path.join(CGGA, f"CGGA.mRNAseq_{tag}_clinical.zip"), "clin").set_index("CGGA_ID")
    return expr, cl


for tag in ["325", "693"]:
    expr, cl = read_cgga(tag)
    s = pd.DataFrame({
        "tamis": scores_from_z(expr, TAMIS, cl.index),
        "myeloid": scores_from_z(expr, MYELOID, cl.index),
    }).join(cl)
    s["OS"] = pd.to_numeric(s["OS"], errors="coerce")
    s["event"] = pd.to_numeric(s["Censor (alive=0; dead=1)"], errors="coerce")
    s = s[s["OS"].notna() & s["event"].notna() & (s["OS"] > 0)].copy()
    s["age"] = pd.to_numeric(s["Age"], errors="coerce")
    s["grade_num"] = s["Grade"].astype(str).str.extract(r"(IV|III|II)")[0].map({"II": 2, "III": 3, "IV": 4})
    tagc = f"CGGA_{tag}"
    print(f"{tagc} n={len(s)}", flush=True)
    cox(s, ["tamis"], tagc)
    cox(s, ["tamis", "myeloid"], f"{tagc}_state_over_count")
    cox(s, ["tamis", "age", "grade_num"], f"{tagc}_age_grade")
    for g in [2, 3, 4]:
        cox(s[s["grade_num"] == g], ["tamis"], f"{tagc}_gradeG{g}")

out = pd.DataFrame(rows)
out.to_csv(os.path.join(TBL, "step07_stratified_stateovercount.tsv"), sep="\t", index=False)
pd.set_option("display.width", 200)
print("\n=== key rows (tamis HR / state-over-count / grade-stratified) ===", flush=True)
pd.set_option("display.max_rows", 200)
print(out[out["covariate"].isin(["tamis", "myeloid"])].to_string(index=False), flush=True)
print("DONE", flush=True)
