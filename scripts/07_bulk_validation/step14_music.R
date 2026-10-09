#!/usr/bin/env Rscript
# H3 rigor (b alt) — MuSiC deconvolution of CGGA_325/693 with the scRNA reference.
# Third deconvolution method (alongside BayesPrism and NNLS) to test §8.5 method agreement.
# Input:  results/bayesprism/{ref.mtx, ref_*.txt, bulk_*.tsv}
# Output: results/bayesprism/fraction_music_<tag>.csv
#
# NB: MuSiC 1.0.0's music_basis() does `sapply(unique(samples), ...)`; with a SINGLE sample
# that collapses to a vector and colMeans() errors ("'x' must be an array of at least two
# dimensions"). Fix: assign each reference cell to one of N_PS random PSEUDO-samples per cell
# type, so unique(samples) has >1 level and the within-cell-type variance can be estimated.
suppressMessages({ library(MuSiC); library(SingleCellExperiment); library(Matrix) })

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA ---
root <- Sys.getenv("GLIOMA_PROJ")
if (!nzchar(root)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .file <- sub("^--file=", "", .args[grep("^--file=", .args)])
  root <- if (length(.file)) normalizePath(file.path(dirname(.file), "..", ".."), mustWork = FALSE) else getwd()
}
DATASETS <- Sys.getenv("GLIOMA_DATA", unset = file.path(dirname(root), "Datasets"))
bp <- file.path(root, "results/bayesprism")
set.seed(42)
N_PS <- 10

ref <- readMM(file.path(bp, "ref.mtx"))            # genes x cells
genes <- readLines(file.path(bp, "ref_genes.txt"))
cells <- readLines(file.path(bp, "ref_cells.txt"))
labels <- readLines(file.path(bp, "ref_labels.txt"))
rownames(ref) <- genes; colnames(ref) <- cells
ref[is.na(ref)] <- 0; ref[!is.finite(ref)] <- 0

ps <- character(length(cells))
for (ct in unique(labels)) {
  idx <- which(labels == ct)
  ps[idx] <- paste0(ct, "_ps", sample(rep(seq_len(N_PS), length.out = length(idx))))
}
sce <- SingleCellExperiment(assays = list(counts = ref),
                            colData = data.frame(cellType = labels, samples = ps, row.names = cells))
cat("sc.sce:", dim(ref), "| pseudo-samples:", length(unique(ps)), "\n")
print(table(labels))

for (tag in c("325", "693")) {
  bulk <- as.matrix(read.table(file.path(bp, paste0("bulk_", tag, ".tsv")),
                               header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
  bulk[is.na(bulk)] <- 0; bulk[!is.finite(bulk)] <- 0
  est <- music_prop(bulk.mtx = bulk, sc.sce = sce, clusters = "cellType",
                    samples = "samples", select.ct = unique(labels), verbose = FALSE)
  write.csv(est$Est.prop.weighted, file.path(bp, paste0("fraction_music_", tag, ".csv")))
  cat("MuSiC", tag, "->", dim(est$Est.prop.weighted), "\n")
}
cat("DONE\n")
