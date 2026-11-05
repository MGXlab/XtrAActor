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
Res_combined <- lapply(partial_files, read.csv) %>% bind_rows()

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

ggsave(output_plot, plot = p, width = 10, height = 8, dpi = 300)