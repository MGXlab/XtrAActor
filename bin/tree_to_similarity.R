# tree_to_similarity.R
# Genome tree -> similarity matrix for pyseer's linear mixed model (--lmm --similarity).
#
# Similarity of two genomes = branch length they share from the root down to their most recent
# common ancestor (ape::vcv), the same quantity as pyseer's scripts/phylogeny_distance.py --lmm.
# An unrooted tree is midpoint-rooted first (needs the R package phangorn).
#
# Usage:
#   Rscript tree_to_similarity.R <tree.nwk> <genome_ids.txt> <similarity.tsv>
#
#   genome_ids.txt : genomes to keep (one per line, e.g. the samples of the phenotype file);
#                    other tips are dropped, genomes missing from the tree are reported.

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript tree_to_similarity.R <tree.nwk> <genome_ids.txt> <similarity.tsv>")
tree <- read.tree(args[1])
keep <- trimws(readLines(args[2]))
keep <- keep[keep != ""]

if (is.null(tree$edge.length)) stop("The tree has no branch lengths")

missing <- setdiff(keep, tree$tip.label)
cat(sprintf("Tips in tree: %d | genomes to keep: %d | missing from the tree: %d\n",
            Ntip(tree), length(keep), length(missing)))
if (length(missing) > 0)
  cat("  e.g.", paste(head(missing, 5), collapse = ", "),
      "\n  (these genomes are left out of the association test by pyseer)\n")
drop <- setdiff(tree$tip.label, keep)
if (length(drop) > 0) tree <- drop.tip(tree, drop)
if (Ntip(tree) < 10) stop("Fewer than 10 genomes shared between the tree and the phenotype file")

el <- tree$edge.length
el[is.na(el) | el < 0] <- 0              # neighbour joining can give negative branch lengths
tree$edge.length <- el

if (!is.rooted(tree)) {
  if (!requireNamespace("phangorn", quietly = TRUE))
    stop("The tree is unrooted. Root it first, or install the R package 'phangorn' for midpoint rooting.")
  tree <- phangorn::midpoint(tree)
  cat("Tree was unrooted: midpoint-rooted\n")
}

C <- vcv(tree)
out <- data.frame(C, check.names = FALSE)
con <- file(args[3], "w")
writeLines(paste(c("", colnames(C)), collapse = "\t"), con)
write.table(round(out, 8), con, sep = "\t", quote = FALSE, col.names = FALSE)
close(con)
cat(sprintf("Similarity matrix for %d genomes written to %s\n", ncol(C), args[3]))
