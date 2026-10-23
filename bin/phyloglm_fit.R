# install.packages("phylolm")
# install.packages("ggtree")
library(ape)
library(phylolm)
library(dplyr)
library(ggplot2)
library(parallel)


# setwd('OneDrive - Friedrich-Schiller-Universität Jena/MCP_struct/Microviridae_analysis/multiple_cluster_analysis/combined_model_phyloglm/')

# df_100 <- read.csv("phyloglm_test_data_all.tsv",header = TRUE, sep = '\t', check.names = FALSE)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript phyloglm.R <data_df.tsv> <tree_file.nwk> [output_file]")
}

data_file <- args[1]
tree_file <- args[2]
gini_file <- args[3]
output_file <- args[4]
plot_file <- args[5]

df_100 <- read.csv(data_file, header = TRUE, sep = '\t', check.names = FALSE)
rownames(df_100) <- df_100[,1]
df_100[,1] <- NULL

tree <- read.tree(tree_file)
#class(tree)           # should be "phylo"
# tree$tip.label        # shows the tip names

df_gini <- read.csv(gini_file, header = TRUE, sep = '\t', check.names = FALSE)

# Features to test: only the ones in THIS chunk's gini file, not every column
# in the full data_df — that's the whole point of splitting into chunks.
features <- intersect(df_gini$feature, names(df_100)[-1])
missing_features <- setdiff(df_gini$feature, names(df_100)[-1])
if (length(missing_features) > 0) {
  warning(sprintf(
    "%d feature(s) from %s not found in %s and will be skipped: %s",
    length(missing_features), gini_file, data_file,
    paste(head(missing_features, 5), collapse = ", ")
  ))
}
if (length(features) == 0) {
  stop(sprintf("No features from %s were found in %s — nothing to fit.", gini_file, data_file))
}

# Function to fit phyloglm for a single feature
fit_phyloglm <- function(feature) {
  #print(feature)  # track progress
  
  # Prepare data
  dat <- df_100[, c(1, which(names(df_100) == feature))]
  names(dat) <- c("y", "x")
  
  # Fit model with error handling
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
  
  # Extract results
  if (is.na(coef(m1)[1])) {
    return(data.frame(
      feature = feature,
      Estimate = NA,
      SE = NA,
      z.value = NA,
      p.value = NA,
      alpha = NA
    ))
  } else {
    m1.sum <- summary(m1)
    return(data.frame(
      feature = feature,
      Estimate = m1.sum$coefficients["x", 1],
      SE = m1.sum$coefficients["x", 2],
      z.value = m1.sum$coefficients["x", 3],
      p.value = m1.sum$coefficients["x", 6],
      alpha = m1$alpha
    ))
  }
}



# Apply to all features
# Res_100_pruned <- lapply(features, fit_phyloglm) %>% bind_rows()
# Detect cores from Slurm allocation (falls back to detectCores() if not set)
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores)) n_cores <- max(1, detectCores() - 1)

# Apply to all features in parallel
Res_100_pruned <- mclapply(features, fit_phyloglm, mc.cores = n_cores) %>% bind_rows()

# View results
#head(Res_100_pruned, n=10)

# join with gini values
Res_100_pruned <- left_join(Res_100_pruned, df_gini, by = "feature")

# NOTE: no padj/filtering here — this is only ONE chunk's p-values, and BH
# correction needs the full pooled set across all chunks to be valid. FDR
# correction happens once, in phyloglm_merge.R, on the combined results.

write.csv(Res_100_pruned,output_file,row.names = FALSE)