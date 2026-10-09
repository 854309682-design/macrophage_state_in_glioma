#!/usr/bin/env Rscript
# H3 rigor (b) — BayesPrism deconvolution of CGGA_325/693 with the scRNA reference.
# Input: results/bayesprism/{ref.mtx, ref_genes.txt, ref_cells.txt, ref_labels.txt, bulk_*.tsv}
# Output: results/bayesprism/fraction_325.csv, fraction_693.csv
suppressMessages({ library(BayesPrism); library(Matrix) })

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA ---
root <- Sys.getenv("GLIOMA_PROJ")
if (!nzchar(root)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .file <- sub("^--file=", "", .args[grep("^--file=", .args)])
  root <- if (length(.file)) normalizePath(file.path(dirname(.file), "..", ".."), mustWork = FALSE) else getwd()
}
DATASETS <- Sys.getenv("GLIOMA_DATA", unset = file.path(dirname(root), "Datasets"))
bp <- file.path(root, "results/bayesprism")

ref <- as.matrix(t(readMM(file.path(bp, "ref.mtx"))))  # cells x genes, dense (BayesPrism needs matrix)
genes <- readLines(file.path(bp, "ref_genes.txt"))
cells <- readLines(file.path(bp, "ref_cells.txt"))
labels <- readLines(file.path(bp, "ref_labels.txt"))
rownames(ref) <- cells; colnames(ref) <- genes
cat("reference:", dim(ref), "\n")
print(table(labels))

for (tag in c("325", "693")) {
  bulk <- as.matrix(read.table(file.path(bp, paste0("bulk_", tag, ".tsv")),
                               header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
  bulkg <- t(bulk)                                 # -> samples x genes (mixture orientation)
  bulkg[is.na(bulkg)] <- 0; ref[is.na(ref)] <- 0
  cat("mixture", tag, ":", dim(bulkg), "range:", range(bulkg), "\n")
  prism <- new.prism(reference = ref, mixture = bulkg, input.type = "count.matrix",
                     cell.type.labels = labels, cell.state.labels = labels, key = "malignant")
  prism <- run.prism(prism, n.cores = 8)
  saveRDS(prism, file.path(bp, paste0("prism_", tag, ".rds")))   # cache Gibbs result
  # BayesPrism >=2.0 changed the signature: get.fraction(bp, which.theta, state.or.type).
  # The 1st attempt failed here ("argument \"which.theta\" is missing") AFTER the ~1 h Gibbs run,
  # so the prism is now cached above to make a downstream error cheap to recover from.
  theta <- get.fraction(prism, which.theta = "final", state.or.type = "type")  # cell types x samples
  write.csv(theta, file.path(bp, paste0("fraction_", tag, ".csv")))
  cat("wrote fraction_", tag, ".csv (", dim(theta), ")\n", sep = "")
}
cat("DONE\n")
