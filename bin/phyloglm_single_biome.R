# phyloglm_single_biome.R
# One phylogenetic logistic regression per single-biome cluster, on the GENOME tree:
#
#     biome of the genome (target biome = 1, other biome = 0)  ~  genome carries the cluster (0/1)
#
# Usage:
#   Rscript phyloglm_single_biome.R <cluster_chunk.tsv> <dice_tree.nwk> <metadata_file> <target_biome> <output_partial.tsv>
#
#   cluster_chunk.tsv : chunk of single_biome_clusters.tsv (cluster | representative | biome | ... | genomes)
#   dice_tree.nwk     : genome tree from dice_tree.py (tip labels = genome IDs)
#   metadata_file     : genome <whitespace> biome, with a header line
#   target_biome      : biome coded as 1 (full label or the part after the last ';')
#   output_partial.tsv: ALL fitted clusters of this chunk, unfiltered. FDR + filtering happen
#                       once in phyloglm_single_biome_report.R.
#
# Same direction as phyloglm_fit.R (biome ~ feature). Method is logistic_MPLE (Firth-penalised)
# instead of logistic_IG10: a single-biome cluster never occurs in the other biome, so the data
# are perfectly separated and the unpenalised likelihood has no finite maximum.

suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript phyloglm_single_biome.R <cluster_chunk.tsv> <tree.nwk> <metadata_file> <target_biome> <output_partial.tsv>")
}
chunk_file <- args[1]; tree_file <- args[2]; meta_file <- args[3]
target_arg <- args[4]; output_file <- args[5]

short <- function(x) trimws(sub(".*;", "", x))

# ---------------------------------------------------------------
# Metadata: split at the first whitespace only (biome names contain spaces)
# ---------------------------------------------------------------
lines <- trimws(readLines(meta_file)[-1])
lines <- lines[grepl("\\s", lines)]
meta  <- data.frame(genome = sub("\\s.*$", "", lines),
                    biome  = trimws(sub("^\\S+\\s+", "", lines)), stringsAsFactors = FALSE)
meta  <- meta[!duplicated(meta$genome), ]
biome_of <- setNames(meta$biome, meta$genome)

# ---------------------------------------------------------------
# Tree: tip labels -> genome IDs as used in the metadata (part before the first '|')
# ---------------------------------------------------------------
tree <- read.tree(tree_file)
tree$tip.label <- sub("\\|.*$", "", tree$tip.label)
drop <- tree$tip.label[duplicated(tree$tip.label) | !(tree$tip.label %in% names(biome_of))]
cat(sprintf("Tips in tree: %d | with biome in metadata: %d\n",
            length(tree$tip.label), length(tree$tip.label) - length(drop)))
if (length(drop) > 0) tree <- drop.tip(tree, drop)
if (length(tree$tip.label) < 4) stop("Fewer than 4 genomes shared between tree and metadata")

# neighbour joining can give negative or zero branch lengths; phyloglm needs positive ones
el <- tree$edge.length
el[is.na(el) | el < 0] <- 0
el[el == 0] <- max(1e-8, 1e-6 * mean(el[el > 0]))
tree$edge.length <- el

tips   <- tree$tip.label
biomes <- sort(unique(biome_of[tips]))
if (length(biomes) != 2) {
  stop("phyloglm needs exactly 2 biomes among the genomes of the tree, found ", length(biomes), ": ",
       paste(biomes, collapse = " | "))
}
target <- biomes[biomes == target_arg | short(biomes) == short(target_arg)]
if (length(target) != 1) target <- biomes[2]      # default: alphabetically last, as in create_dataset_phyloglm.py
y <- as.integer(biome_of[tips] == target)
cat(sprintf("Target (1): %s n=%d | Reference (0): %s n=%d\n",
            target, sum(y), setdiff(biomes, target), sum(1 - y)))

# ---------------------------------------------------------------
# Fit one cluster
# ---------------------------------------------------------------
clusters <- read.delim(chunk_file, header = TRUE, quote = "", check.names = FALSE,
                       colClasses = "character")
cat(sprintf("Fitting %d cluster(s)\n", nrow(clusters)))

fit_cluster <- function(i) {
  carriers <- strsplit(clusters$genomes[i], ",", fixed = TRUE)[[1]]
  x   <- as.integer(tips %in% carriers)
  res <- data.frame(cluster = clusters$cluster[i], representative = clusters$representative[i],
                    biome = clusters$biome[i], annotation = clusters$annotation[i],
                    n_genomes = sum(x), n_genomes_biome = sum(biome_of[tips] == clusters$biome[i]),
                    Estimate = NA_real_, SE = NA_real_, z.value = NA_real_, p.value = NA_real_,
                    alpha = NA_real_, stringsAsFactors = FALSE)
  if (sum(x) < 2 || sum(x) == length(x)) return(res)

  dat <- data.frame(y = y, x = x, row.names = tips)
  m1  <- tryCatch(
    suppressWarnings(phyloglm(y ~ x, data = dat, phy = tree, method = "logistic_MPLE")),
    error = function(e) NULL)
  if (is.null(m1)) return(res)

  co <- summary(m1)$coefficients
  res$Estimate <- co["x", "Estimate"]
  res$SE       <- co["x", "StdErr"]
  res$z.value  <- co["x", "z.value"]
  res$p.value  <- co["x", "p.value"]
  res$alpha    <- m1$alpha
  res
}

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)

results <- do.call(rbind, parallel::mclapply(seq_len(nrow(clusters)), fit_cluster, mc.cores = n_cores))

n_failed <- sum(is.na(results$Estimate))
cat(sprintf("Fits done: %d ok, %d failed or skipped\n", nrow(results) - n_failed, n_failed))

write.table(results, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
