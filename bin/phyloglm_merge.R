library(dplyr)
library(ggplot2)

# 1. Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

# All arguments except the last two are the partial CSV files
n_args <- length(args)
partial_files <- args[1:(n_args - 2)]
output_csv    <- args[n_args - 1]
output_plot   <- args[n_args]

# 2. Read all chunk results and combine them
Res_combined <- lapply(partial_files, read.csv) %>% bind_rows()

# 3. FDR correction ONCE on the complete, combined dataset
Res_combined$padj <- p.adjust(Res_combined$p.value, method = "fdr")

# Significance threshold: padj < 0.05
padj_cutoff <- 0.05
Res_combined$significant <- !is.na(Res_combined$padj) & Res_combined$padj < padj_cutoff

# z-value at which padj is exactly 0.05 (Benjamini-Hochberg):
# with m tested features and k of them significant, a feature needs p <= k * 0.05 / m,
# i.e. a feature with exactly this p-value gets padj = 0.05. This p-value is converted to |z|
# (phyloglm p-values are two-sided normal: p = 2 * pnorm(-|z|)).
# If no feature is significant, k = 1 is used: the z-value a single feature would need.
m <- sum(!is.na(Res_combined$p.value))
k <- max(sum(Res_combined$significant), 1)
p_cutoff <- k * padj_cutoff / m
z_cutoff <- if (m > 0) qnorm(1 - p_cutoff / 2) else NA_real_
cat(sprintf("Features tested: %d | padj < %g: %d | p cutoff: %.3g | |z| cutoff: %.3f\n",
            m, padj_cutoff, sum(Res_combined$significant), p_cutoff, z_cutoff))

# 4. Write the final, complete CSV
write.csv(Res_combined, output_csv, row.names = FALSE)

# plot
p <- ggplot(Res_combined, aes(x = Estimate, y = z.value, label=feature)) +
  geom_point(aes(color = alpha), alpha = 0.2,size=1.5) +   # base layer
  geom_point(
    data = subset(Res_combined, significant),
    aes(color = alpha), size=1.5) +# highlight significant ones
  geom_text(
    data =subset(Res_combined, z.value <= -8),
    size = 3, hjust=1.2)+
  geom_text(
      data =subset(Res_combined, z.value >= 5),
      size = 3, hjust=1.2)+
  { if (!is.na(z_cutoff)) geom_hline(yintercept = c(-z_cutoff, z_cutoff),
                                     linetype = "dashed", color = "black") } +
  #scale_color_gradientn(colours = c("#364B9AFF", "#4A7BB7FF", "#6EA6CDFF", "#98CAE1FF", "#C2E4EFFF", "#EAECCCFF", "#FEDA8BFF", "#FDB366FF", "#F67E4BFF", "#DD3D2DFF", "#A50026FF"),name = "alpha") +
  scale_color_gradient(name = "Phylogenetic\nsignal (alpha)", low = "#82e9e9ff", high = "#0d4d4dff")+
  theme_minimal()+
  theme(
    axis.text = element_text(size=14),
    axis.title = element_text(size=14),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14)
    
  )

ggsave(output_plot, plot = p, width = 10, height = 8, dpi = 300)