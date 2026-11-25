library(dplyr)
library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript plot_zvalue_featimps.R <phyloglm.csv> <Feature_importance.tsv> <output.png> <biome1> <biome2>")
}

phyloglm_results <- args[1]
biome_file <- args[2]
plot_file <- args[3]

# biome2 is the target biome (coded 1 in create_dataset_phyloglm.py --target_biome), so a
# positive z-value means the feature goes with biome2, a negative one with biome1.
# Labels: the part after the last ';' (e.g. "Environmental;Aquatic;Freshwater;Lake" -> "Lake")
short_label <- function(x) trimws(sub(".*;", "", x))
biome1_label <- short_label(args[4])
biome2_label <- short_label(args[5])

df_phylo_res <- read.csv(phyloglm_results)
df_feat_imp <- read.csv(biome_file, sep='\t')
df_merged <- left_join(df_feat_imp, df_phylo_res, by = "feature")

aa_groups <- list(
  "Positive charge" = list(aas = c("R", "H", "K"),          color = "#ef4806"),
  "Negative charge" = list(aas = c("D", "E"),                color = "#234abb"),
  "Polar uncharged"  = list(aas = c("S", "T", "N", "Q"),     color = "#1eae83"),
  "Hydrophobic"      = list(aas = c("A", "V", "I", "L", "M"), color = "#f4b713"),
  "Aromatic"         = list(aas = c("F", "Y", "W"),           color = "#bf21a2"),
  "Other"            = list(aas = c("C", "U", "G", "P"),      color = "#4a4444")
)

# Build lookup vectors from aa_groups
aa_to_group <- unlist(lapply(names(aa_groups), function(g) {
  setNames(rep(g, length(aa_groups[[g]]$aas)), aa_groups[[g]]$aas)
}))
group_colors <- setNames(
  sapply(aa_groups, `[[`, "color"),
  names(aa_groups)
)

# Extract the single-letter AA code from the index (assumes format like "A123" or "pos_A")
# Adjust regex to extract the AA letter - this handles formats like "pos5_A", "A123", "123_A"
df_merged <- df_merged %>%
  mutate(
    aa_letter = toupper(regmatches(feature, regexpr("[A-Za-z]+$", feature))),  # grabs trailing letters
    aa_group  = ifelse(aa_letter %in% names(aa_to_group),
                       aa_to_group[aa_letter], "Other"),
    biome   = ifelse(z.value >= 0, biome2_label, biome1_label)
  )

# Significance threshold: padj < 0.05
padj_cutoff <- 0.05

# z-value at which padj is exactly 0.05 (Benjamini-Hochberg), same as in phyloglm_merge.R:
# with m tested features and k of them significant, a feature needs p <= k * 0.05 / m,
# i.e. a feature with exactly this p-value gets padj = 0.05. This p-value is converted to |z|
# (phyloglm p-values are two-sided normal: p = 2 * pnorm(-|z|)).
# If no feature is significant, k = 1 is used: the z-value a single feature would need.
m <- sum(!is.na(df_phylo_res$p.value))
k <- max(sum(!is.na(df_phylo_res$padj) & df_phylo_res$padj < padj_cutoff), 1)
z_cutoff <- if (m > 0) qnorm(1 - (k * padj_cutoff / m) / 2) else NA_real_

p <- ggplot(df_merged, aes(x = z.value, y =gini_imp , label = feature,
                           fill = aa_group, color = aa_group, shape = biome)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_point(
    data = subset(df_merged, gini_imp >= 0.0005 & padj < padj_cutoff), size = 4,alpha=0.9,key_glyph="rect") +
  geom_text_repel(
    data = subset(df_merged, gini_imp >= 0.006 & padj < padj_cutoff),
    size = 5,
    color = "black",
    box.padding = 0.5,      # padding around each label
    point.padding = 0.3,    # padding between label and point
    max.overlaps = Inf,     # show all labels even if they overlap
    segment.color = "grey50" # color of the line connecting label to point
  )+
  geom_hline(yintercept = 0.001, linetype = "dashed", color = "black") +
  { if (!is.na(z_cutoff)) geom_vline(xintercept = c(-z_cutoff, z_cutoff),
                                     linetype = "dashed", color = "black") } +
  scale_fill_manual(values = group_colors, name = "AA Group", 
                    breaks = names(aa_groups), labels = names(aa_groups)) +
  scale_color_manual(values = group_colors,
                     guide=guide_legend(theme = theme(
                       legend.position = "None"
                       ))) +
  scale_shape_manual(values = setNames(c(24, 25), c(biome2_label, biome1_label)),
                     name = "Biome Indicator") +
  theme_minimal() +
  theme(
    axis.text    = element_text(size = 22,colour = "black"),
    legend.position = "none",
    # legend.text  = element_text(size = 16,colour = "black"),
    # legend.title = element_text(size = 16,colour = "black"),
    axis.minor.ticks.x.bottom = element_blank()
  )


p

ggsave(plot_file)
# #eom_text(
# data = subset(df_merged, gini_imp >= 0.005),
# size = 7, hjust = 1.2, color = "black",) +