library(dplyr)
library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript plot_zvalue_featimps.R <phyloglm_results.csv> <feature_importance.tsv> <plot_file.png>")
}

phyloglm_results <- args[1]
biome_file <- args[2]
plot_file <- args[3]

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
    biome   = ifelse(z.value >= 0, "Ocean", "Lake")
  )

get_z_threshold <- function(z, padj, side = c("pos", "neg")) {
  side <- match.arg(side)
  keep <- !is.na(padj) & !is.na(z) & (if (side == "pos") z > 0 else z < 0)
  d <- data.frame(z = z[keep], padj = padj[keep])
  d <- d[order(d$padj), ]           # padj decreasing as |z| increases -> sort by padj ascending
  if (nrow(d) < 2) return(NA_real_)
  approx(x = d$padj, y = d$z, xout = 0.01, rule = 2)$y
}

z_thresh_pos <- get_z_threshold(df_merged$z.value, df_merged$padj, "pos")
z_thresh_neg <- get_z_threshold(df_merged$z.value, df_merged$padj, "neg")

p <- ggplot(df_merged, aes(x = z.value, y =gini_imp , label = feature,
                           fill = aa_group, color = aa_group, shape = biome)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_point(
    data = subset(df_merged, gini_imp >= 0.0005 & padj <= 0.01), size = 4,alpha=0.9,key_glyph="rect") +
  geom_text_repel(
    data = subset(df_merged, gini_imp >= 0.006 & padj <= 0.01),
    size = 5,
    color = "black",
    box.padding = 0.5,      # padding around each label
    point.padding = 0.3,    # padding between label and point
    max.overlaps = Inf,     # show all labels even if they overlap
    segment.color = "grey50" # color of the line connecting label to point
  )+
  geom_hline(yintercept = 0.001, linetype = "dashed", color = "black") +
  geom_vline(xintercept =  z_thresh_pos, linetype = "dashed", color = "black") +
  geom_vline(xintercept = z_thresh_neg, linetype = "dashed", color = "black") +
  scale_fill_manual(values = group_colors, name = "AA Group", 
                    breaks=c("Positive charged", "Negative charged","Polar uncharged","Hydrophobic","Aromatic","Other"),
                    labels=c("Positive charged", "Negative charged","Polar uncharged","Hydrophobic","Aromatic","Other") ) +
  scale_color_manual(values = group_colors,
                     guide=guide_legend(theme = theme(
                       legend.position = "None"
                       ))) +
  scale_shape_manual(values = c("Ocean" = 24, "Lake" = 25),
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