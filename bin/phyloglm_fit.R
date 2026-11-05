# phyloglm_fit.R
# Fits one phylogenetic logistic regression (phyloglm) per feature of ONE gini chunk.
#
# Usage:
#   Rscript phyloglm_fit.R <phyloglm_dataset.tsv> <tree.nwk> <gini_chunk.tsv> <output_partial.csv>
#
#   phyloglm_dataset.tsv : gene | <target 0/1 column> | feature columns ...   (from create_dataset_phyloglm.py)
#   tree.nwk             : FastTree tree of the same cluster (tip labels = gene IDs)
#   gini_chunk.tsv       : feature | gini   (one chunk from SPLIT_GINI; decides WHICH features are fitted)
#   output_partial.csv   : results for ALL fitted features of this chunk, unfiltered.
#                          FDR correction + significance filtering happen once in phyloglm_merge.R.

suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  library(dplyr)
  library(parallel)
})

# ---------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript phyloglm_fit.R <phyloglm_dataset.tsv> <tree.nwk> <gini_chunk.tsv> <output_partial.csv>")
}

data_file   <- args[1]
tree_file   <- args[2]
gini_file   <- args[3]
output_file <- args[4]

# ---------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------
# check.names = FALSE keeps feature names like "243_L" unchanged (R would otherwise make them "X243_L")
df_100 <- read.csv(data_file, header = TRUE, sep = '\t', check.names = FALSE)
rownames(df_100) <- df_100[, 1]   # gene IDs
df_100[, 1] <- NULL               # first remaining column = target (0/1), the rest = features

tree <- read.tree(tree_file)

df_gini <- read.csv(gini_file, header = TRUE, sep = '\t', check.names = FALSE)

# ---------------------------------------------------------------
# Match genes between dataset and tree
# ---------------------------------------------------------------
common <- intersect(rownames(df_100), tree$tip.label)
cat(sprintf("Genes in dataset: %d | tips in tree: %d | shared: %d\n",
            nrow(df_100), length(tree$tip.label), length(common)))

if (length(common) < 4) {
  stop("Fewer than 4 genes shared between dataset and tree. ",
       "Dataset IDs: ", paste(head(rownames(df_100), 3), collapse = ", "),
       " | Tree tips: ", paste(head(tree$tip.label, 3), collapse = ", "))
}

# keep only shared genes, in both objects
df_100 <- df_100[common, , drop = FALSE]
if (length(common) < length(tree$tip.label)) {
  tree <- drop.tip(tree, setdiff(tree$tip.label, common))
}

target_col <- names(df_100)[1]
cat(sprintf("Target column: %s | class counts: %s\n", target_col,
            paste(names(table(df_100[[1]])), table(df_100[[1]]), sep = "=", collapse = ", ")))

# ---------------------------------------------------------------
# Features to test: only those in THIS gini chunk (and present in the dataset)
# ---------------------------------------------------------------
features <- intersect(as.character(df_gini$feature), names(df_100)[-1])
if (length(features) == 0) {
  stop("None of the chunk's features are columns of the dataset. ",
       "Chunk: ", paste(head(df_gini$feature), collapse = ", "),
       " | Dataset: ", paste(head(names(df_100)[-1]), collapse = ", "))
}
n_missing <- nrow(df_gini) - length(features)
if (n_missing > 0) {
  cat(sprintf("WARNING: %d chunk feature(s) not found in the dataset\n", n_missing))
}
cat(sprintf("Fitting %d feature(s)\n", length(features)))

# ---------------------------------------------------------------
# Fit phyloglm for a single feature
# ---------------------------------------------------------------
fit_phyloglm <- function(feature) {

  dat <- df_100[, c(target_col, feature)]
  names(dat) <- c("y", "x")

  m1 <- tryCatch(
    phyloglm(
      formula = y ~ x,
      data = dat,
      phy = tree,
      method = "logistic_IG10",
      boot = 2
    ),
    error = function(e) list(coefficients = NA)
  )

  if (is.na(coef(m1)[1])) {
    return(data.frame(
      feature  = feature,
      Estimate = NA,
      SE       = NA,
      z.value  = NA,
      p.value  = NA,
      alpha    = NA
    ))
  } else {
    m1.sum <- summary(m1)
    return(data.frame(
      feature  = feature,
      Estimate = m1.sum$coefficients["x", 1],
      SE       = m1.sum$coefficients["x", 2],
      z.value  = m1.sum$coefficients["x", 3],
      p.value  = m1.sum$coefficients["x", 6],
      alpha    = m1$alpha
    ))
  }
}

# ---------------------------------------------------------------
# Run in parallel on the cores SLURM gave this task
# ---------------------------------------------------------------
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores)) n_cores <- max(1, detectCores() - 1)

Res_100_pruned <- mclapply(features, fit_phyloglm, mc.cores = n_cores) %>% bind_rows()

n_failed <- sum(is.na(Res_100_pruned$Estimate))
cat(sprintf("Fits done: %d ok, %d failed\n", nrow(Res_100_pruned) - n_failed, n_failed))

# ---------------------------------------------------------------
# Add gini values and write ALL results (no FDR / no filtering here)
# ---------------------------------------------------------------
Res_100_pruned <- left_join(Res_100_pruned, df_gini, by = "feature")

write.csv(Res_100_pruned, output_file, row.names = FALSE)