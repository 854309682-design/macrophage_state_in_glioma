"""
Step 5h — pySCENIC cisTarget motif pruning (pySCENIC step 2), NumPy-2 shim.

Thin wrapper around the pySCENIC 0.12.1 CLI. `pyscenic ctx` (and every other
sub-command) crashes at import time on this env because pyscenic/transform.py
still references np.object (and pyscenic/rss.py + diptest.py reference np.float)
— aliases removed in NumPy 2.0. The shim below restores both aliases (1:1 with
the Python builtins, so behaviour-preserving) before importing the CLI, then
forwards the arguments verbatim.

Canonical invocation — prefer the sibling shell wrapper (records the exact
parameters in version control):

  bash scripts/05_interaction/step09h_ctx.sh

Equivalent direct call (run from the repo root):

  env/.venv/bin/python scripts/05_interaction/step09h_ctx.py ctx \
    results/scenic/step09g_adj.tsv \
    <DATASETS>/scenic_resources/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
    <DATASETS>/scenic_resources/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
    --annotations_fname <DATASETS>/scenic_resources/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname results/scenic/step09f_myeloid.loom \
    --output results/scenic/step09h_regulons.csv \
    --num_workers 16 --nes_threshold 3.0 --rank_threshold 5000

(Note: the adjacencies input requires --expression_mtx_fname so the modules can
be rebuilt inside ctx; the loom is genes x cells, so no -t/--transpose.)
"""
import sys

import numpy as np

for _alias, _target in (("object", object), ("float", float)):
    if not hasattr(np, _alias):
        setattr(np, _alias, _target)

from pyscenic.cli.pyscenic import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
