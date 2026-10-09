"""
Step 1 — Inspect the local GSE103224 raw matrix format.

This resolves the format question flagged in analysis_strategy_1.md §4 Step 1 and handoff.md §3:
for each of the 8 GSE103224 patients, PJDATA contains three files:
  - GSM*_PJ0XX.filtered.matrix.txt.gz   -> gene x cell, headerless, col1=Ensembl(versioned), col2=symbol, col3+=UMI values  [the one we parse]
  - GSM*_PJ0XX.filtered.matrix.csv      -> 10x-style CSV (larger)
  - GSM*_PJ0XX.filtered.matrix-2.csv    -> 10x-style CSV (larger)

The strategy doc and the .txt.gz agree (gene x cell, headerless UMI). We read the .txt.gz.

Output
  - logs/step01_gse103224_inspect.log
  - results/tables/step01_gse103224_matrix_dims.tsv   (per-patient genes x cells)
  - writes a manifest of per-patient .txt.gz paths to Datasets/GSE103224_file_manifest.tsv

Run:  env/.venv/bin/python scripts/00_download/step01_inspect_gse103224.py
"""
import gzip, glob, os, sys

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
DATASETS = os.environ.get("GLIOMA_DATA") or os.path.join(os.path.dirname(ROOT), "Datasets")
DATA_ROOT = os.path.join(DATASETS, "GSE103224")
OUT_DIR = ROOT

log_lines = []
def log(msg):
    log_lines.append(msg)
    print(msg)

def count_fields_first_line(path):
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n")
    return len(header.split("\t")), header[:60]

log("=== GSE103224 raw matrix inspection ===")
log(f"data root: {DATA_ROOT}")

manifests = sorted(glob.glob(os.path.join(DATA_ROOT, "*filtered.matrix.txt.gz")))
log(f"found {len(manifests)} .txt.gz matrices")

rows = []
for p in manifests:
    base = os.path.basename(p)
    # e.g. GSM2758471_PJ016.filtered.matrix.txt.gz -> PJ016
    sample = base.split(".")[0].split("_")[1] if "_" in base else base
    # counts
    with gzip.open(p, "rt") as fh:
        n_genes = sum(1 for _ in fh)
    n_fields, first = count_fields_first_line(p)
    n_cells = n_fields - 2  # minus Ensembl, minus symbol
    rows.append([sample, base, n_genes, n_cells, n_fields, first])
    log(f"{sample:8s} genes={n_genes:6d} cells={n_cells:6d} header_sample='{first}'")

# Write dims table
tbl_out = os.path.join(OUT_DIR, "results/tables/step01_gse103224_matrix_dims.tsv")
os.makedirs(os.path.dirname(tbl_out), exist_ok=True)
with open(tbl_out, "w") as fh:
    fh.write("sample\tfile\tn_genes\tn_cells\tn_fields\theader_sample\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")
log(f"wrote {tbl_out}")

# Write manifest
man_out = os.path.join(DATASETS, "GSE103224_file_manifest.tsv")
os.makedirs(os.path.dirname(man_out), exist_ok=True)
with open(man_out, "w") as fh:
    fh.write("sample\tfile\tpath\n")
    for r in rows:
        fh.write(f"{r[0]}\t{r[1]}\t{os.path.normpath(os.path.join(DATA_ROOT, r[1]))}\n")
log(f"wrote {man_out}")

# Save log
os.makedirs(os.path.join(OUT_DIR, "logs"), exist_ok=True)
with open(os.path.join(OUT_DIR, "logs/step01_gse103224_inspect.log"), "w") as fh:
    fh.write("\n".join(log_lines) + "\n")

n_cells_total = sum(r[3] for r in rows)
log(f"TOTAL cells across 8 samples: {n_cells_total}")
log("DONE")
