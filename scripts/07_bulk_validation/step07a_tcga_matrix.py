"""
Step 7a — Build the TCGA-GBM/LGG sample x signature matrix.

For each of the 925 STAR-Counts files (Datasets/TCGA/star_counts/<file_id>.tsv):
  - map file_id -> case_submitter_id via the GDC API,
  - extract the TAM-IS / microglia / immune signature genes (gene_name) and their
    TPM (tpm_unstranded),
  - assemble a samples x genes matrix.
Also fetch the per-case clinical (already saved) and merge.

Outputs: results/tables/step07_tcga_signature_matrix.tsv + step07_tcga_samples.tsv

Run: env/.venv/bin/python scripts/07_bulk_validation/step07a_tcga_matrix.py
"""
import json
import os
import urllib.request

import numpy as np
import pandas as pd

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
DATASETS = os.environ.get("GLIOMA_DATA") or os.path.join(os.path.dirname(ROOT), "Datasets")
TCGA = os.path.join(DATASETS, "TCGA")
IDS = "/tmp/tcga_ids.txt"
OUT = os.path.join(ROOT, "results/tables")

SIG = {
    "TAMIS": ["SPP1", "TREM2", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "LGALS3",
              "CD163", "CTSL", "CCL3", "CCL4", "GPNMB", "LPL", "HLA-DRA", "HLA-DRB1",
              "S100A9", "HMOX1", "FTL", "FTH1"],
    "MICRO": ["P2RY12", "TMEM119", "CX3CR1", "SALL1", "CRYBB1", "CST3", "PROS1", "MERTK"],
    "CD8T": ["CD8A", "CD8B", "GZMB", "PRF1", "IFNG"],
    "CHECKPOINT": ["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
    "MES": ["VIM", "CD44", "ANXA2", "LGALS1"],
    "MYELOID_CONTENT": ["PTPRC", "ITGAM", "CSF1R", "AIF1", "TYROBP", "FCER1G",
                        "LYZ", "CD68", "CD14"],
    "HN": ["ACTB", "GAPDH", "RPL13A", "B2M"],
}
ALLGENES = sorted({g for v in SIG.values() for g in v})


def gdc_map():
    from datetime import datetime
    payload = {"filters": {"op": "in",
                           "content": {"field": "file_id",
                                       "value": [l.split("\t")[0] for l in open(IDS) if l.strip()]}},
               "fields": "file_id,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type",
               "format": "json", "size": 2000}
    req = urllib.request.Request("https://api.gdc.cancer.gov/files",
                                 data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        d = json.load(r)
    rows = {}
    for h in d["data"]["hits"]:
        fid = h["file_id"]
        case = h["cases"][0]["submitter_id"] if h.get("cases") else None
        rows[fid] = case
    return rows


if __name__ == "__main__":
    mp = gdc_map()
    print(f"mapped {len(mp)} files -> cases", flush=True)
    recs = []
    for fn in sorted(os.listdir(os.path.join(TCGA, "star_counts"))):
        fid = fn.replace(".tsv", "")
        path = os.path.join(TCGA, "star_counts", fn)
        # STAR augmented gene counts: skip '#' comment lines, then header
        df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
        df = df.rename(columns={df.columns[0]: "gene_id"})
        gcol = "gene_name" if "gene_name" in df.columns else df.columns[1]
        tcol = "tpm_unstranded" if "tpm_unstranded" in df.columns else df.columns[-3]
        sub = df[df[gcol].isin(ALLGENES)][[gcol, tcol]].drop_duplicates(gcol)
        vals = dict(zip(sub[gcol], pd.to_numeric(sub[tcol], errors="coerce")))
        rec = {"file_id": fid, "case": mp.get(fid)}
        rec.update({g: vals.get(g, np.nan) for g in ALLGENES})
        recs.append(rec)
    mat = pd.DataFrame(recs)
    mat.to_csv(os.path.join(OUT, "step07_tcga_signature_matrix.tsv"), sep="\t", index=False)
    print(f"matrix: {mat.shape}; cases with data: {mat['case'].notna().sum()}", flush=True)
    print(mat.head(3).iloc[:, :6].to_string(), flush=True)
    print("DONE", flush=True)
