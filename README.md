# macrophage_state_in_glioma

Analysis code accompanying the manuscript

> **A reproducible, interpretable multi-omics integration framework for clinical decision-making in glioma: a hypoxic-niche macrophage state with context-dependent prognostic value**
>
> *(submitted to Frontiers in Big Data — Research Topic 78756, "Multi-Omics Integration for Clinical Decision-Making")*

Public glioma single-cell, spatially resolved and bulk transcriptomes are integrated to define
**TAM-IS** — an immunosuppressive tumour-associated macrophage **state** — and to test three
pre-registered hypotheses. The pipeline is entirely computational (pure dry-lab).

| | Hypothesis | Outcome reported in the manuscript |
|---|---|---|
| **H1** | the TAM-IS state is reproducible and transferable across independent cohorts | **partially supported** — cross-cohort transfer at moderate-to-good accuracy (AUROC 0.68–0.91; headline 0.712 / 0.674), but at chance level for one composition-shifted cohort |
| **H2** | TAM-IS is coupled to the hypoxic / peri-necrotic niche | **supported** — enrichment in hypoxic, perinecrotic tumour regions across all 13 regional sections and all 3 cohorts |
| **H3** | TAM-IS abundance is prognostic in bulk cohorts | **context-dependent** — associated with worse survival in IDH-mutant, grade 2–3 glioma (meta-analytic HR 1.91, p = 1.4 × 10⁻⁸), not in IDH-wildtype glioblastoma; the independence claim for IDH-wt GBM is withdrawn |

Scale: 110,462 quality-controlled cells across 3 discovery cohorts, plus independent single-cell,
spatial (IvyGAP) and bulk (TCGA, CGGA) validation layers.

---

## Repository layout

The numeric prefixes follow the order in which analyses were added during the study — they are
*not* a strict topological order. Each script is standalone and reads/writes versioned
intermediate objects (`objects/stepXX_*`) and result tables (`results/tables/stepXX_*`).

| Directory | Script | What it does |
|---|---|---|
| `00_download/` | `step01_inspect_gse103224.py` | inspect the local GSE103224 raw matrix format |
| | `step02_build_anndata.py` | build AnnData objects for the three discovery cohorts |
| | `step02b_build_extra_cohorts.py` | build the two extra cohorts used for H1 replication |
| `01_qc/` | `step03_qc.py` | per-cohort QC + doublet removal (Scrublet) |
| `02_integration/` | `step04_integrate.py` | scVI integration (raw counts, `n_latent=30`) |
| | `step04b_extract_myeloid.py` | annotate lineages and extract the myeloid compartment |
| | `step11_arches_h1.py` | scArches reference mapping — **exploratory, withdrawn; see Caveats** |
| `03_infercnv/` | `step03b_infercnv.py` | inferCNV-based malignant calling (infercnvpy), validating the marker-based labels |
| `04_myeloid_states/` | `step05_tamis_discovery.py` | TAM-IS state discovery (first pass) |
| | `step05b_tamis_markers.py` | define the fixed TAM-IS marker signature |
| | `step06_h1_replication.py` | H1 replication (signature-projection reproducibility) |
| | `step18_h1_rigor.py` | H1 rigor: leave-one-dataset-out / cross-cohort transfer |
| | `step20_h1_hvg_independent.py` | H1 rigor with an HVG-independent TAM-IS definition |
| `05_interaction/` | `step09a_export_cellchat.py` + `step09b_cellchat.R` | CellChat v2 cell–cell communication |
| | `step09c_export_nichenet.py` + `step09d_nichenet.R` | NicheNet ligand→target analysis |
| | `step09e_download_scenic_resources.sh` | fetch the pySCENIC / cisTarget hg38 resources |
| | `step09f_export_scenic.py`, `step09g_grn.py`, `step09h_ctx.py` (+`.sh`), `step09i_aucell.py` (+`.sh`), `step09j_aucell_heatmap.py` | pySCENIC GRN → cisTarget → AUCell regulon activity |
| `06_spatial/` | `step06b_ivygap_h2.py` | H2 spatial/niche validation with IvyGAP (GSE107559) |
| `07_bulk_validation/` | `step07a_tcga_matrix.py`, `step07b_tcga_survival.py` | TCGA-GBM/LGG signature matrix + survival |
| | `step07c_cgga_survival.py`, `step07d_stratified.py`, `step07e_gbm_idh.py` | CGGA + grade/IDH-stratified survival, state-over-count control |
| | `step07f_tcga_idh_mgmt.py` | TCGA OS with uniform cBioPortal covariates (IDH/MGMT) |
| | `step12_deconv.py`, `step13a_export_bayesprism.py`, `step13b_bayesprism.R`, `step13c_bayesprism_cox.py`, `step14_music.R`, `step17_deconv_consistency.py` | deconvolution (NNLS / BayesPrism / MuSiC) + method agreement |
| | `step16_meta_analysis.py` | cross-cohort random-effects meta-analysis of the TAM-IS HR |
| `08_causal_druggability/` | `step08_druggability.py` | druggability of TAM-IS genes (DepMap/CRISPR evidence, DGIdb, Open Targets) |
| `09_figures/` | `step10_figures.py`, `step15_h3_forest.py`, `step19_h1_fig.py` | individual figure panels |
| | `step21_main_figures.py` | all main panels as importable draw functions (27 standalone files) |
| | `step22_combined_figures.py` | compose the panels into 5 combined figures + 2 supplementary |
| | `step23_vector_combined_figures.py` | true-vector (PDF) + 300-dpi JPEG versions of the 7 combined figures |

`DATA_MANIFEST.md` catalogues every input and reference resource (accession, local filename,
purpose, verification status).

## Environment

### Python (all steps except the R steps)

Python **3.11** (`numpy` 2.4.6, `pandas` 2.3.3, `scipy` 1.17.1, `scikit-learn` 1.9.0,
`anndata` 0.12.19, `scanpy` 1.11.5, `scvi-tools` 1.4.2, `pySCENIC` 0.12.1, `cell2location` 0.1.5,
`infercnvpy`, `scrublet` 0.2.3, `matplotlib` 3.11.1, `seaborn` 0.13.2).
`torch` is the **CPU** build — the pipeline requires no GPU.

```bash
uv venv env/.venv --python 3.11
uv pip install --python env/.venv/bin/python -r env/requirements.txt
# or, without uv:
python3.11 -m venv env/.venv && env/.venv/bin/pip install -r env/requirements.txt
```

The scripts' docstrings and the `Run:` lines assume this interpreter at `env/.venv/bin/python`.

### R (the `*.R` steps only)

R **4.6.1** with: `CellChat`, `nichenetr`, `BayesPrism`, `MuSiC`, `SingleCellExperiment`,
`Matrix`, `dplyr`, `tidyr`, `tibble`, `ggplot2`.
(`inferCNV` is done in Python via `infercnvpy`, so no R inferCNV install is needed.)

## Configuration — where the data and outputs live

**No path editing is required.** Every script resolves two roots from its own location at run time:

| Variable | Meaning | Default | Override |
|---|---|---|---|
| `ROOT` | repository root (intermediate objects, result tables, figures) | two levels above `scripts/<step>/` | `GLIOMA_PROJ` |
| `DATASETS` | shared raw-data / reference-resource directory | `<ROOT>/../Datasets` | `GLIOMA_DATA` |

The expected default layout is therefore:

```
glioma/
├── Datasets/            # raw GEO/TCGA/CGGA/IvyGAP downloads + reference resources
└── macrophage_state_in_glioma/     # this repository
    ├── scripts/         # analysis scripts (this repo)
    ├── env/             # Python environment + requirements.txt
    ├── objects/         # versioned intermediate objects (h5ad/loom/rds) — created at run time
    ├── results/         # tables/ and figures/ — created at run time
    └── logs/            # per-step run logs — created at run time
```

If your data live elsewhere, export the overrides instead of editing code:

```bash
export GLIOMA_PROJ=/path/to/this/repo
export GLIOMA_DATA=/path/to/Datasets
```

## Data

All input data are public; nothing is redistributed in this repository. See `DATA_MANIFEST.md`
for file-level detail.

| Dataset | Type | Role in the study | Source |
|---|---|---|---|
| GSE103224 | 10x scRNA-seq (8 patients) | discovery | GEO |
| GSE131928 | Smart-seq2 scRNA-seq (28 patients) | discovery | GEO |
| GSE163120 | 10x scRNA-seq (4 patients) | discovery | GEO |
| GSE89567 | scRNA-seq | H1 validation (IDH-mutant) | GEO |
| GSE70630 | scRNA-seq | H1 validation | GEO |
| GSE107559 (IvyGAP) | regionally dissected spatial profiles | H2 spatial validation | GEO / Ivy GAP |
| TCGA-GBM / TCGA-LGG (`lgggbm_tcga_pub`) | bulk RNA-seq + clinical | H3 prognosis | GDC / cBioPortal |
| CGGA `mRNAseq_325`, `mRNAseq_693` | bulk RNA-seq + clinical | H3 validation | CGGA portal |
| NicheNet, pySCENIC/cisTarget hg38, DGIdb, Open Targets, DepMap | reference resources | mechanism / druggability | upstream providers |

## Running the pipeline

Run in the order below; each step consumes the previous steps' objects under `objects/` and
writes tables to `results/tables/`. Long steps print progress; logs are written to `logs/`.

```bash
P=env/.venv/bin/python

# --- discovery atlas ---
$P scripts/00_download/step01_inspect_gse103224.py
$P scripts/00_download/step02_build_anndata.py
$P scripts/00_download/step02b_build_extra_cohorts.py
$P scripts/01_qc/step03_qc.py
$P scripts/03_infercnv/step03b_infercnv.py
$P scripts/02_integration/step04_integrate.py
$P scripts/02_integration/step04b_extract_myeloid.py

# --- TAM-IS definition + H1 ---
$P scripts/04_myeloid_states/step05_tamis_discovery.py
$P scripts/04_myeloid_states/step05b_tamis_markers.py
$P scripts/04_myeloid_states/step06_h1_replication.py
$P scripts/04_myeloid_states/step18_h1_rigor.py
$P scripts/04_myeloid_states/step20_h1_hvg_independent.py

# --- mechanism: CellChat, NicheNet, pySCENIC ---
$P scripts/05_interaction/step09a_export_cellchat.py && Rscript scripts/05_interaction/step09b_cellchat.R
$P scripts/05_interaction/step09c_export_nichenet.py && Rscript scripts/05_interaction/step09d_nichenet.R
bash scripts/05_interaction/step09e_download_scenic_resources.sh   # one-off resource download
$P scripts/05_interaction/step09f_export_scenic.py
$P scripts/05_interaction/step09g_grn.py results/scenic/step09f_myeloid.loom \
     "$GLIOMA_DATA/scenic_resources/allTFs_hg38.txt" \
     --method grnboost2 --output results/scenic/step09g_adj.tsv --num_workers 16 --seed 42
bash scripts/05_interaction/step09h_ctx.sh
bash scripts/05_interaction/step09i_aucell.sh
$P scripts/05_interaction/step09j_aucell_heatmap.py

# --- H2 spatial ---
$P scripts/06_spatial/step06b_ivygap_h2.py

# --- H3 bulk ---
$P scripts/07_bulk_validation/step07a_tcga_matrix.py
$P scripts/07_bulk_validation/step07b_tcga_survival.py
$P scripts/07_bulk_validation/step07c_cgga_survival.py
$P scripts/07_bulk_validation/step07d_stratified.py
$P scripts/07_bulk_validation/step07e_gbm_idh.py
$P scripts/07_bulk_validation/step07f_tcga_idh_mgmt.py
$P scripts/07_bulk_validation/step12_deconv.py
$P scripts/07_bulk_validation/step13a_export_bayesprism.py && Rscript scripts/07_bulk_validation/step13b_bayesprism.R
$P scripts/07_bulk_validation/step13c_bayesprism_cox.py
Rscript scripts/07_bulk_validation/step14_music.R
$P scripts/07_bulk_validation/step16_meta_analysis.py
$P scripts/07_bulk_validation/step17_deconv_consistency.py

# --- druggability + figures ---
$P scripts/08_causal_druggability/step08_druggability.py
$P scripts/09_figures/step21_main_figures.py
$P scripts/09_figures/step22_combined_figures.py
$P scripts/09_figures/step23_vector_combined_figures.py
```

Notes:
* `step22` / `step23` import `step21_main_figures` and must be run from their own directory
  (or with `scripts/09_figures` on `PYTHONPATH`) so the import resolves.
* Random seed `42` throughout; the statistical unit is the patient; multiple testing uses BH-FDR.

## Figures

`step21_main_figures.py` holds one draw function per panel; `step22_combined_figures.py` places
them into the 5 main + 2 supplementary combined figures; `step23_vector_combined_figures.py`
writes the final **true-vector PDF** plus 300-dpi JPEG versions, each carrying the manuscript
figure-legend sentence as vector text. Figures are written to `results/figures/` and
`manuscript/Figures/`.

## Caveats

* `02_integration/step11_arches_h1.py` — scArches reference mapping was run and then
  **withdrawn**: its internal quality metric (`query_mapped_tamis_frac` ≈ 0.97) was implausible
  on inspection. The script is kept for transparency; its output is **not** reported in the
  manuscript, and no result should be drawn from it.
* Deconvolution uses locally implemented MuSiC/NNLS procedures. **CIBERSORTx was not run**, so
  deconvolution-method sensitivity is not fully characterized; method agreement between
  BayesPrism, NNLS and MuSiC is quantified in `step17`.
* No germline causal layer: the glioma GWAS summary statistics were not accessible, so
  `step08` rests on transcriptomic and pharmacological evidence only.
* The H3 result is **context-dependent** (IDH-mutant / lower-grade), and the H2 spatial evidence
  derives from regionally dissected bulk profiles (IvyGAP) rather than single-cell-resolution
  spatial transcriptomics.
* Derived result tables and large intermediate objects are **not** distributed in this
  repository; `DATA_MANIFEST.md` records what each step produces and where inputs come from.

## Citation

If you use this code, please cite the manuscript above once published. Author and journal
metadata will be added here upon acceptance.

## License

MIT — see [`LICENSE`](LICENSE).
