# phyloglm_single_biome_report.R
# Merge the per-chunk phyloglm results of the single-biome clusters, FDR-correct ONCE over all
# clusters, and report the clusters with padj < threshold.
#
# Usage:
#   Rscript phyloglm_single_biome_report.R <partials_dir> <padj_threshold> <all_results.tsv> <report.tsv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript phyloglm_single_biome_report.R <partials_dir> <padj_threshold> <all_results.tsv> <report.tsv>")
}
partials_dir <- args[1]; padj_max <- as.numeric(args[2])
all_file <- args[3]; report_file <- args[4]

files <- sort(list.files(partials_dir, pattern = "_partial\\.tsv$", full.names = TRUE))
res <- do.call(rbind, lapply(files, read.delim, header = TRUE, quote = "",
                             check.names = FALSE, stringsAsFactors = FALSE))

res$padj <- p.adjust(res$p.value, method = "fdr")
res <- res[order(res$padj, res$p.value), ]
num <- c("Estimate", "SE", "z.value", "p.value", "alpha", "padj")
res[num] <- lapply(res[num], signif, 4)

write.table(res, all_file, sep = "\t", quote = FALSE, row.names = FALSE)

report <- res[!is.na(res$padj) & res$padj < padj_max, ]
write.table(report, report_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Single-biome clusters: %d | fitted: %d | padj < %g: %d\n",
            nrow(res), sum(!is.na(res$p.value)), padj_max, nrow(report)))
if (nrow(report) > 0) print(table(report$biome))
