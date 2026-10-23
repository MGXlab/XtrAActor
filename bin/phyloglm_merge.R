library(dplyr)
library(ggplot2)

# 1. Kommandozeilenargumente einlesen
args <- commandArgs(trailingOnly = TRUE)

# Alle Argumente außer den letzten beiden sind die partial-CSV-Dateien
n_args <- length(args)
partial_files <- args[1:(n_args - 2)]
output_csv    <- args[n_args - 1]
output_plot   <- args[n_args]

# 2. Alle Chunk-Ergebnisse einlesen und zusammenführen
# colClasses forces the column types explicitly — without this, read.csv()
# infers types from the data, and a chunk whose fit_phyloglm() returned only
# NA rows (or, before the phyloglm_fit.R fix, an empty post-filter table)
# comes back with e.g. `feature` as logical instead of character, which then
# makes bind_rows() fail when combining it with the other (correctly typed)
# chunks.
col_spec <- c(feature = "character", Estimate = "numeric", SE = "numeric",
              z.value = "numeric", p.value = "numeric", alpha = "numeric",
              gini = "numeric")

Res_combined <- lapply(partial_files, function(f) {
  read.csv(f, colClasses = col_spec)
}) %>% bind_rows()

# 3. FDR-Korrektur EINMAL auf dem vollständigen, zusammengeführten Datensatz
Res_combined$padj <- p.adjust(Res_combined$p.value, method = "fdr")

# 4. Finale, vollständige CSV schreiben
write.csv(Res_combined, output_csv, row.names = FALSE)

# plot
p <- ggplot(Res_combined, aes(x = Estimate, y = z.value, label=feature)) +
  geom_point(aes(color = alpha), alpha = 0.2,size=1.5) +   # base layer
  geom_point(
    data = subset(Res_combined, -log10(padj) > 3),
    aes(color = alpha), size=1.5) +# highlight significant ones
  geom_text(
    data =subset(Res_combined, z.value <= -8),
    size = 3, hjust=1.2)+
  geom_text(
      data =subset(Res_combined, z.value >= 5),
      size = 3, hjust=1.2)+
  geom_hline(yintercept = 3.9, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -3.9, linetype = "dashed", color = "black") +
  #scale_color_gradientn(colours = c("#364B9AFF", "#4A7BB7FF", "#6EA6CDFF", "#98CAE1FF", "#C2E4EFFF", "#EAECCCFF", "#FEDA8BFF", "#FDB366FF", "#F67E4BFF", "#DD3D2DFF", "#A50026FF"),name = "alpha") +
  scale_color_gradient(name = "Phylogenetic\nsignal (alpha)", low = "#82e9e9ff", high = "#0d4d4dff")+
  theme_minimal()+
  theme(
    axis.text = element_text(size=14),
    axis.title = element_text(size=14),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14)
    
  )



p
ggsave(output_plot)
#print(-log10(0.005))