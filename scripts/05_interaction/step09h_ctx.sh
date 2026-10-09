#!/usr/bin/env bash
# Step 5h — pySCENIC cisTarget motif pruning (ctx). Records the exact parameters.
# Run from anywhere:  bash scripts/05_interaction/step09h_ctx.sh
set -euo pipefail

# --- paths: portable — override with GLIOMA_DATA (defaults: a Datasets/ dir beside this repo) ---
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATASETS="${GLIOMA_DATA:-$(cd "${REPO_ROOT}/.." && pwd)/Datasets}"   # shared, outside the repo (see DATA_MANIFEST.md)
RES="${REPO_ROOT}/Datasets/scenic_resources"               # fallback if DATASETS moved
if [ ! -d "${RES}" ]; then RES="${DATASETS}/scenic_resources"; fi

ADJ="${REPO_ROOT}/results/scenic/step09g_adj.tsv"
LOOM="${REPO_ROOT}/results/scenic/step09f_myeloid.loom"
OUT="${REPO_ROOT}/results/scenic/step09h_regulons.csv"

if [ ! -s "${ADJ}" ]; then
  echo "ERROR: ${ADJ} missing/empty — run step09g (GRNBoost2) first." >&2
  exit 1
fi

echo "[step09h] ctx: ${ADJ} -> ${OUT}"
"${REPO_ROOT}/env/.venv/bin/python" "${REPO_ROOT}/scripts/05_interaction/step09h_ctx.py" ctx \
  "${ADJ}" \
  "${RES}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather" \
  "${RES}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather" \
  --annotations_fname "${RES}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl" \
  --expression_mtx_fname "${LOOM}" \
  --output "${OUT}" \
  --num_workers 16 \
  --nes_threshold 3.0 \
  --rank_threshold 5000
echo "[step09h] DONE -> ${OUT}"
