#!/usr/bin/env bash
# Step 5e — Download the pySCENIC / cisTarget hg38 resources (run from a network that can reach aertslab).
#
# WHY: resources.aertslab.org is unreachable from the lab machine (TCP timeout), so run this with a
#      VPN on, or on another machine / a relay. Everything else (pip/conda/github/bioconductor) is fine.
#
# Files (v9 / mc9nr gene-based — the combo proven with pySCENIC 0.12.1):
#   hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather       ~1.2 GB
#   hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather ~1.2 GB
#   motifs-v9-nr.hgnc-m0.001-o0.0.tbl                                                    ~99 MB
#   allTFs_hg38.txt                                                                      ~30 KB
#
# Usage:  bash scripts/05_interaction/step09e_download_scenic_resources.sh
set -euo pipefail

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ="${GLIOMA_PROJ:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
DATASETS="${GLIOMA_DATA:-$(cd "${PROJ}/.." && pwd)/Datasets}"
DEST="${DATASETS}/scenic_resources"
BASE="https://resources.aertslab.org/cistarget"
DBDIR="$BASE/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based"
PY="${PROJ}/env/.venv/bin/python"

FEATHERS=(
  "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
  "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
)
OTHERS=(
  "$BASE/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
  "$BASE/tf_lists/allTFs_hg38.txt"
)

mkdir -p "$DEST"
cd "$DEST"

for f in "${FEATHERS[@]}"; do
  echo "==> $f"
  aria2c -x8 -s8 -c --file-allocation=none --retry-wait=5 --max-tries=0 \
         --summary-interval=15 -d "$DEST" -o "$f" "$DBDIR/$f"
done

for url in "${OTHERS[@]}"; do
  echo "==> $(basename "$url")"
  aria2c -x4 -s4 -c --file-allocation=none --retry-wait=5 --max-tries=0 \
         -d "$DEST" -o "$(basename "$url")" "$url"
done

echo
echo "==> files:"; ls -la "$DEST"

echo
echo "==> verifying sha1 of the two feather DBs against the official per-file manifests"
# aertslab serves one '<sha1>  <filename>' manifest per feather (no global manifest).
for f in "${FEATHERS[@]}"; do
  want=$(curl -sSL --max-time 30 "$DBDIR/$f.sha1sum.txt" | awk '{print $1}' | head -1)
  have=$(sha1sum "$f" | awk '{print $1}')
  if [ -n "$want" ] && [ "$want" = "$have" ]; then
    echo "OK sha1      $f"
  else
    echo "!! SHA1 MISMATCH for $f (want=$want have=$have) — re-run; aria2c resumes."
  fi
done

echo
echo "==> validating the feather DBs open in ctxcore (real integrity check):"
GLIOMA_SCENIC_DEST="$DEST" "$PY" - <<'PY'
import os

from ctxcore.rnkdb import FeatherRankingDatabase

d = os.environ["GLIOMA_SCENIC_DEST"]
for f in ["hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
          "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"]:
    db = FeatherRankingDatabase(d + "/" + f, name=f)
    print("OK", f, "| genes:", db.total_genes,
          "| motifs:", len(db.ct_db.all_motif_or_track_ids))
PY
echo "DONE — tell the assistant the files are in place to resume pySCENIC (Step 5, step09f/step09g)."
