"""
Step 7e — H3 decisive sub-question: GBM-only & IDH-stratified.

TCGA: split by project (TCGA-GBM vs TCGA-LGG); Cox(tamis) + state-over-count within each.
CGGA: split by grade (G4=GBM vs G2/G3) and by IDH (mut vs wt); Cox(tamis) within each.
Tests whether the "state over count" holds in the paper's main setting (IDH-wt GBM).
Outputs: results/tables/step07_gbm_idh_stratified.tsv
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
        rows.append({"dataset": tag, "model": "|".join(covs), "n": len(c),
                     "events": int(c["event"].sum()), "covariate": "tamis",
                     "HR": np.nan, "p": np.nan}); return
    X = sm.add_constant(c[covs])
    res = sm.PHReg(c["OS"], X, status=c["event"]).fit()
    for i, v in enumerate(X.columns):
        if v in ("tamis",):
            rows.append({"dataset": tag, "model": "|".join(covs), "n": len(c),
                         "events": int(c["event"].sum()), "covariate": v,
                         "HR": round(float(np.exp(res.params[i])), 3), "p": float(res.pvalues[i])})


def gdc():
    payload = {"filters": {"op": "in", "content": {"field": "project.project_id",
               "value": ["TCGA-GBM", "TCGA-LGG"]}},
               "fields": "submitter_id,project.project_id,demographic.vital_status,"
                         "demographic.days_to_death,diagnoses.days_to_last_follow_up,"
                         "diagnoses.age_at_diagnosis,diagnoses.tumor_grade",
               "format": "json", "size": 5000}
    d = json.load(urllib.request.urlopen(urllib.request.Request(
        "https://api.gdc.cancer.gov/cases", data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"}), timeout=120))
    out = []
    for h in d["data"]["hits"]:
        dg = (h.get("diagnoses") or [{}])[0]
        out.append({"case": h["submitter_id"], "project": h["project"]["project_id"],
                    "vital": (h.get("demographic") or {}).get("vital_status"),
                    "dtd": (h.get("demographic") or {}).get("days_to_death"),
                    "dtf": dg.get("days_to_last_follow_up"), "age": dg.get("age_at_diagnosis")})
    return pd.DataFrame(out).drop_duplicates("case")


m = pd.read_csv(os.path.join(TBL, "step07_tcga_signature_matrix.tsv"), sep="\t")
tam = m.groupby("case")[[g for g in TAMIS if g in m.columns]].mean()
mye = m.groupby("case")[[g for g in MYELOID if g in m.columns]].mean()
tz = (tam - tam.mean()) / tam.std()
mz = (mye - mye.mean()) / mye.std()
t = pd.DataFrame({"case": tz.index, "tamis": tz.mean(axis=1).values,
                  "myeloid": mz.mean(axis=1).values}).merge(gdc(), on="case", how="inner")
t["OS"] = t["dtd"].fillna(t["dtf"]); t["event"] = (t["vital"] == "Dead").astype(int)
t = t[t["OS"].notna() & (t["OS"] > 0)]
print("TCGA by project:", t["project"].value_counts().to_dict(), flush=True)
for proj, tag in [("TCGA-GBM", "TCGA_GBM"), ("TCGA-LGG", "TCGA_LGG")]:
    cox(t[t["project"] == proj], ["tamis"], tag)
    cox(t[t["project"] == proj], ["tamis", "myeloid"], tag + "_state_over_count")


def read_cgga(tag):
    def rd(z, kind):
        with zipfile.ZipFile(z) as zf:
            n = [x for x in zf.namelist() if x.endswith(".txt") and "__MACOSX" not in x][0]
            with zf.open(n) as fh:
                return pd.read_csv(fh, sep="\t", index_col=0 if kind == "expr" else None)
    e = rd(os.path.join(CGGA, f"CGGA.mRNAseq_{tag}.RSEM-genes.zip"), "expr")
    e.index = [str(g).strip("'\"") for g in e.index]
    cl = rd(os.path.join(CGGA, f"CGGA.mRNAseq_{tag}_clinical.zip"), "clin").set_index("CGGA_ID")
    return e, cl


for tag in ["325", "693"]:
    e, cl = read_cgga(tag)
    g = [x for x in TAMIS if x in e.index]; mg = [x for x in MYELOID if x in e.index]
    zt = e.loc[g]; zt = zt.sub(zt.mean(axis=1), axis=0).div(zt.std(axis=1).replace(0, np.nan), axis=0)
    zm = e.loc[mg]; zm = zm.sub(zm.mean(axis=1), axis=0).div(zm.std(axis=1).replace(0, np.nan), axis=0)
    s = pd.DataFrame({"tamis": zt.mean(0), "myeloid": zm.mean(0)}).join(cl)
    s["OS"] = pd.to_numeric(s["OS"], errors="coerce")
    s["event"] = pd.to_numeric(s["Censor (alive=0; dead=1)"], errors="coerce")
    s = s[s["OS"].notna() & s["event"].notna() & (s["OS"] > 0)].copy()
    s["grade_num"] = s["Grade"].astype(str).str.extract(r"(IV|III|II)")[0].map({"II": 2, "III": 3, "IV": 4})
    s["idh_mut"] = s["IDH_mutation_status"].astype(str).str.lower().str.contains("mut")
    tc = f"CGGA_{tag}"
    cox(s[s["grade_num"] == 4], ["tamis"], tc + "_G4_GBM")
    cox(s[s["grade_num"].isin([2, 3])], ["tamis"], tc + "_G2G3")
    cox(s[s["idh_mut"]], ["tamis"], tc + "_IDHmut")
    cox(s[~s["idh_mut"]], ["tamis"], tc + "_IDHwt")
    cox(s[s["idh_mut"]], ["tamis", "myeloid"], tc + "_IDHmut_state_over_count")
    cox(s[~s["idh_mut"]], ["tamis", "myeloid"], tc + "_IDHwt_state_over_count")

out = pd.DataFrame(rows)
out.to_csv(os.path.join(TBL, "step07_gbm_idh_stratified.tsv"), sep="\t", index=False)
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 200)
print(out.to_string(index=False), flush=True)
print("DONE", flush=True)
