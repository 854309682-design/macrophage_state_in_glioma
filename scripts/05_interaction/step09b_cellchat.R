#!/usr/bin/env Rscript
# Step 5b — CellChat v2 cell-cell communication on the annotated TME (genes x cells).
# Input: results/cellchat/{expr.mtx, genes.txt, cells.txt, labels.csv}
# Output: results/tables/step09_cellchat_interactions.csv, results/tables/step09_cellchat_summary.csv,
#         objects/step09_cellchat.rds
suppressMessages({
  library(CellChat); library(Matrix)
})

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA ---
root <- Sys.getenv("GLIOMA_PROJ")
if (!nzchar(root)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .file <- sub("^--file=", "", .args[grep("^--file=", .args)])
  root <- if (length(.file)) normalizePath(file.path(dirname(.file), "..", ".."), mustWork = FALSE) else getwd()
}
DATASETS <- Sys.getenv("GLIOMA_DATA", unset = file.path(dirname(root), "Datasets"))
cc <- file.path(root, "results/cellchat")

X <- readMM(file.path(cc, "expr.mtx"))           # genes x cells
genes <- readLines(file.path(cc, "genes.txt"))
cells <- readLines(file.path(cc, "cells.txt"))
lab <- read.csv(file.path(cc, "labels.csv"), row.names = 1, stringsAsFactors = FALSE)
rownames(X) <- genes; colnames(X) <- cells
meta <- data.frame(labels = lab$lineage, row.names = rownames(lab), stringsAsFactors = FALSE)
cat("matrix:", dim(X), "\n")

cellchat <- createCellChat(object = as.matrix(X), meta = meta, group.by = "labels")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)

df.net <- subsetCommunication(cellchat)
write.csv(df.net, file.path(root, "results/tables/step09_cellchat_interactions.csv"), row.names = FALSE)

# pathway-level summary: total communication probability per pathway
p <- cellchat@netP$prob
if (length(dim(p)) == 3) {
  path_sum <- apply(p, 3, function(m) sum(m, na.rm = TRUE))
  write.csv(data.frame(pathway = names(path_sum), prob = as.numeric(path_sum)),
            file.path(root, "results/tables/step09_cellchat_summary.csv"), row.names = FALSE)
}
saveRDS(cellchat, file.path(root, "objects/step09_cellchat.rds"))
cat("n LR interactions:", nrow(df.net), "\n")
cat("top pathways (by total prob):\n")
print(head(sort(path_sum, decreasing = TRUE), 10))
cat("DONE\n")
