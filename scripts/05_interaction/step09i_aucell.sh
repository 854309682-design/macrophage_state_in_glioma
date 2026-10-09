#!/usr/bin/env bash
# Step 5i — pySCENIC AUCell regulon activity. Records the exact parameters.
# Run from anywhere:  bash scripts/05_interaction/step09i_aucell.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

LOOM="${REPO_ROOT}/results/scenic/step09f_myeloid.loom"
REG="${REPO_ROOT}/results/scenic/step09h_regulons.csv"
OUT="${REPO_ROOT}/results/scenic/step09h_aucell.csv"

if [ ! -s "${REG}" ]; then
  echo "ERROR: ${REG} missing/empty — run step09h (ctx) first." >&2
  exit 1
fi

echo "[step09i] aucell: ${LOOM} + ${REG} -> ${OUT}"
"${REPO_ROOT}/env/.venv/bin/python" "${REPO_ROOT}/scripts/05_interaction/step09i_aucell.py" aucell \
  "${LOOM}" \
  "${REG}" \
  --output "${OUT}" \
  --num_workers 16 \
  --seed 42
echo "[step09i] DONE -> ${OUT}"
