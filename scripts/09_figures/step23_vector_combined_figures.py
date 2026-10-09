"""
Step 23 - true-vector (PDF) versions of the 7 combined manuscript figures.

Background
----------
`manuscript/Figures/*.pdf` were previously produced by *rasterising* the step-22
JPEG composites and stamping a text-title band on top with Pillow. The resulting
PDFs therefore contain a single 300-dpi bitmap and no vector content (no embedded
fonts, no selectable text). Journal submission wants vector line art, so this
step re-emits the same seven layouts as genuine vector PDFs.

The layouts, panel draw functions and A/B/C panel letters are imported verbatim
from step 22, so the geometry matches the reviewed figures; only the output
format (vector PDF) and the figure title (the manuscript figure-legend sentence,
drawn as vector text) change.

Outputs (vector PDF + 300-dpi JPEG):
  manuscript/Figures/figure1_atlas.{pdf,jpeg}
  manuscript/Figures/figure2_definition.{pdf,jpeg}
  manuscript/Figures/figure3_repro_spatial.{pdf,jpeg}
  manuscript/Figures/figure4_clinical.{pdf,jpeg}
  manuscript/Figures/figure5_mechanism.{pdf,jpeg}
  manuscript/Figures/figureS1_validation.{pdf,jpeg}
  manuscript/Figures/figureS2_stats.{pdf,jpeg}

Run (all):     env/.venv/bin/python scripts/09_figures/step23_vector_combined_figures.py
Run (subset):  env/.venv/bin/python scripts/09_figures/step23_vector_combined_figures.py figure1_atlas
"""
import os
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import step22_combined_figures as s22

OUT = os.path.join(s22.ROOT, "manuscript/Figures")
os.makedirs(OUT, exist_ok=True)

# Figure titles = the manuscript figure-legend sentences, verbatim from
# manuscript/manuscript.md > Figure Legends.
TITLES = {
    "figure1_atlas": "Figure 1. Single-cell atlas of human glioma and "
                     "identification of the TAM-IS macrophage state.",
    "figure2_definition": "Figure 2. Molecular definition of the TAM-IS state.",
    "figure3_repro_spatial": "Figure 3. Cross-cohort reproducibility and "
                             "hypoxic-niche localization of TAM-IS.",
    "figure4_clinical": "Figure 4. Prognostic association of the TAM-IS signature "
                        "in bulk glioma cohorts.",
    "figure5_mechanism": "Figure 5. Deconvolution robustness, cell\u2013cell "
                         "communication, regulatory programs, and druggability of TAM-IS.",
    "figureS1_validation": "Figure S1. Technical validation of cell annotation and deconvolution.",
    "figureS2_stats": "Figure S2. Extended communication, regulatory, and clinical statistics.",
}

# Titles fit on one line on the wide (12-16 in) canvases; only wrap if a title
# would otherwise be unreasonably long.
WRAP = 140

# Per-figure canvas + inter-panel spacing (2026-09-13, 2nd pass).
# The step-22 default (hspace=0.55, wspace=0.38) let panels whose tick labels /
# colorbars extend outside the axes bleed into the neighbouring cell (e.g.
# fig3 A/B, fig5 B/C, figS2 panel D vs panel C's colorbar). These wider canvases
# and larger gaps were verified overlap-free with a tight-bbox check
# (/tmp/check_layout2.py). Values: (width_in, height_in, hspace, wspace).
LAYOUT = {
    "figure1_atlas":         (19.0, 10.0, 0.75, 0.80),
    "figure2_definition":    (13.0, 15.0, 0.90, 0.50),
    "figure3_repro_spatial": (19.0, 10.0, 0.75, 0.80),
    "figure4_clinical":      (19.0, 10.0, 0.75, 0.80),
    "figure5_mechanism":     (19.5, 10.0, 0.75, 1.05),
    "figureS1_validation":   (14.0, 11.5, 0.90, 0.70),
    "figureS2_stats":        (14.0, 11.0, 0.80, 0.85),
}


def compose(fname, title, figsize, grid, panels, hspace, wspace):
    """Same layout as step22.compose; title is the manuscript legend sentence and
    both outputs are written to manuscript/Figures/ (vector PDF + 300-dpi JPEG)."""
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(grid[0], grid[1], hspace=hspace, wspace=wspace)
    for i, item in enumerate(panels):
        name, r, c = item[0], item[1], item[2]
        rs = item[3] if len(item) > 3 else 1
        cs = item[4] if len(item) > 4 else 1
        ax = fig.add_subplot(gs[r:r + rs, c:c + cs])
        s22.DRAW[name](ax)
        ax.set_title(ax.get_title(), fontsize=9)
        ax.text(-0.14, 1.08, "ABCDEFGH"[i], transform=ax.transAxes,
                fontsize=15, fontweight="bold", va="top")
    fig.suptitle("\n".join(textwrap.wrap(title, WRAP)), fontsize=11, y=1.0,
                 va="bottom")
    fig.savefig(os.path.join(OUT, f"{fname}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, f"{fname}.jpeg"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {fname} (vector PDF + JPEG)", flush=True)


if __name__ == "__main__":
    import sys
    want = sys.argv[1:]
    for fname, (title22, figsize0, grid, panels) in s22.LAYOUTS.items():
        if want and fname not in want:
            continue
        w, h, hspace, wspace = LAYOUT[fname]
        try:
            print(f"== {fname} ==", flush=True)
            compose(fname, TITLES[fname], (w, h), grid, panels, hspace, wspace)
        except Exception as e:
            print(f"  {fname} FAILED: {type(e).__name__} {e}", flush=True)
    print("DONE", flush=True)
