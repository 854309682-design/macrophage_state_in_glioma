#!/usr/bin/env Rscript
# Step 5d — NicheNet: which malignant-niche ligands drive the TAM-IS transcriptional program.
# Sender   = malignant hypoxic-niche pseudobulk (step09c_sender_expr.csv)
# Receiver = TAM-IS pseudobulk (step09c_receiver_expr.csv)
# Reference: human ligand_target_matrix (nsga2r final) + lr_network (2021-12-21),
#            Datasets/nichenet/ ; matrix is targets (rows) x ligands (cols).
# Output:
#   results/tables/step09d_nichenet_ligand_activities.tsv   ranked ligand activity (pearson)
#   results/tables/step09d_nichenet_ligand_target.tsv       best ligands x TAM-IS genes
#   results/tables/step09d_nichenet_lr_pairs.tsv            expressed ligand-receptor pairs
#   results/figures/supplementary/step09d_nichenet_ligand_target.{pdf,jpeg}
suppressMessages({
  library(nichenetr); library(dplyr); library(tidyr); library(tibble); library(ggplot2)
})

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA ---
root <- Sys.getenv("GLIOMA_PROJ")
if (!nzchar(root)) {
  .args <- commandArgs(trailingOnly = FALSE)
  .file <- sub("^--file=", "", .args[grep("^--file=", .args)])
  root <- if (length(.file)) normalizePath(file.path(dirname(.file), "..", ".."), mustWork = FALSE) else getwd()
}
DATASETS <- Sys.getenv("GLIOMA_DATA", unset = file.path(dirname(root), "Datasets"))
nich <- file.path(DATASETS, "nichenet")
inp <- file.path(root, "results/nichenet")
tbl <- file.path(root, "results/tables")
fig <- file.path(root, "results/figures/supplementary")
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

DETECT <- 0.10   # a gene is "expressed" if detected in >= 10% of cells within a group
TOPN <- 20       # best upstream ligands to report

sender <- as.matrix(read.csv(file.path(inp, "step09c_sender_expr.csv"), row.names = 1, check.names = FALSE))
receiver <- as.matrix(read.csv(file.path(inp, "step09c_receiver_expr.csv"), row.names = 1, check.names = FALSE))
sen_det <- as.matrix(read.csv(file.path(inp, "step09c_sender_detect.csv"), row.names = 1, check.names = FALSE))
rec_det <- as.matrix(read.csv(file.path(inp, "step09c_receiver_detect.csv"), row.names = 1, check.names = FALSE))

ligand_target_matrix <- readRDS(file.path(nich, "ligand_target_matrix_nsga2r_final.rds"))
lr_network <- readRDS(file.path(nich, "lr_network_human_21122021.rds"))

geneset <- intersect(readLines(file.path(inp, "step09c_geneset_tamis.txt")), rownames(ligand_target_matrix))
stopifnot("all" %in% colnames(sen_det), "all" %in% colnames(rec_det))
background <- intersect(rownames(receiver)[rec_det[, "all"] >= DETECT], rownames(ligand_target_matrix))
expressed_ligands <- intersect(unique(lr_network$from), rownames(sen_det)[sen_det[, "all"] >= DETECT])
expressed_receptors <- intersect(unique(lr_network$to), rownames(rec_det)[rec_det[, "all"] >= DETECT])
potential_ligands <- intersect(expressed_ligands, colnames(ligand_target_matrix))

cat("geneset:", length(geneset), "| background:", length(background),
    "| expressed ligands:", length(expressed_ligands),
    "| potential ligands:", length(potential_ligands),
    "| expressed receptors:", length(expressed_receptors), "\n")

ligand_activities <- predict_ligand_activities(
  geneset = geneset,
  background_expressed_genes = background,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
) %>% arrange(desc(pearson))
write.table(ligand_activities, file.path(tbl, "step09d_nichenet_ligand_activities.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

best_ligands <- head(ligand_activities$test_ligand, TOPN)

# ligand -> target regulatory potential restricted to the TAM-IS gene set
lt <- ligand_target_matrix[geneset, best_ligands, drop = FALSE]
lt_long <- as.data.frame(lt) %>%
  tibble::rownames_to_column("target") %>%
  pivot_longer(-target, names_to = "ligand", values_to = "weight")
write.table(lt_long, file.path(tbl, "step09d_nichenet_ligand_target.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

lr_pairs <- lr_network %>%
  filter(from %in% best_ligands, to %in% expressed_receptors) %>%
  group_by(from, to) %>%
  summarise(database = paste(sort(unique(database)), collapse = ";"), .groups = "drop") %>%
  arrange(from, to)
write.table(lr_pairs, file.path(tbl, "step09d_nichenet_lr_pairs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- figure: ligand x target heatmap (top 40 TAM-IS genes by max regulatory potential) ---
top_targets <- lt_long %>%
  group_by(target) %>% summarise(m = max(weight), .groups = "drop") %>%
  arrange(desc(m)) %>% head(40) %>% pull(target)
sub <- lt_long %>% filter(target %in% top_targets)
sub$ligand <- factor(sub$ligand, levels = best_ligands)

p <- ggplot(sub, aes(x = ligand, y = target, fill = weight)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#b2182b", name = "regulatory\npotential") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank()) +
  labs(title = "NicheNet: malignant hypoxic-niche ligands -> TAM-IS program",
       subtitle = paste0("sender = malignant hypoxic-high, receiver = TAM-IS (",
                         length(best_ligands), " ligands x ", length(top_targets), " targets)"))

ggsave(file.path(fig, "step09d_nichenet_ligand_target.pdf"), p, width = 9, height = 8)
ggsave(file.path(fig, "step09d_nichenet_ligand_target.jpeg"), p, width = 9, height = 8, dpi = 300)

cat("top ligands:\n")
print(head(as.data.frame(ligand_activities), 15))
cat("DONE\n")
