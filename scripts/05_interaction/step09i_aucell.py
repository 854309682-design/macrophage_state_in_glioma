"""
Step 5i — pySCENIC AUCell regulon activity (pySCENIC step 3), NumPy-2 shim.

Thin wrapper around the pySCENIC 0.12.1 CLI. `pyscenic aucell` crashes at import
time on this env because pyscenic still references np.object / np.float — aliases
removed in NumPy 2.0. The shim below restores both (1:1 with the Python builtins,
behaviour-preserving) before importing the CLI, then forwards the arguments
verbatim.

Canonical invocation — prefer the sibling shell wrapper:

  bash scripts/05_interaction/step09i_aucell.sh

Equivalent direct call (run from the repo root):

  env/.venv/bin/python scripts/05_interaction/step09i_aucell.py aucell \
    results/scenic/step09f_myeloid.loom \
    results/scenic/step09h_regulons.csv \
    --output results/scenic/step09h_aucell.csv \
    --num_workers 16 --seed 42

(Note: aucell's signature loader accepts the ctx CSV of enriched motifs —
load_signatures() dispatches .csv/.tsv to df2regulons(load_motifs(...)) — so the
ctx output feeds aucell directly; no GMT conversion needed.)
"""
import sys

import numpy as np

for _alias, _target in (("object", object), ("float", float)):
    if not hasattr(np, _alias):
        setattr(np, _alias, _target)

from pyscenic.cli.pyscenic import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
