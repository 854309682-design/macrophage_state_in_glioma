"""
Step 5g — GRNBoost2 (pySCENIC step 1) with a NumPy-2 compatibility shim.

pySCENIC 0.12.1 (2022) still references np.object / np.float, aliases removed in
NumPy 2.0 (this env runs numpy 2.4). Both map 1:1 onto the Python builtins, so
restoring them before importing pySCENIC is behaviour-preserving.

All arguments are forwarded verbatim to pySCENIC's arboreto_with_multiprocessing
(the single-node multiprocessing entry point, recommended over the dask backend).

Usage:
  env/.venv/bin/python scripts/05_interaction/step09g_grn.py \
    results/scenic/step09f_myeloid.loom \
    $GLIOMA_DATA/scenic_resources/allTFs_hg38.txt \
    --method grnboost2 --output results/scenic/step09g_adj.tsv \
    --num_workers 16 --seed 42
"""
import sys

import numpy as np

for _alias, _target in (("object", object), ("float", float)):
    if not hasattr(np, _alias):
        setattr(np, _alias, _target)

from pyscenic.cli.arboreto_with_multiprocessing import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
