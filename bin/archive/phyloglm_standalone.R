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
# Features to test (columns 2 onward)
features <- names(df_100)[-c(1)]

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

# adjust the pvalues (Benjamini-Hochberg)
Res_100_pruned$padj <- p.adjust(Res_100_pruned$p.value, method = "fdr")

Res_100_pruned_filt <- filter(Res_100_pruned, padj < 0.01)

write.csv(Res_100_pruned_filt,output_file,row.names = FALSE)

# plot
p <- ggplot(Res_100_pruned, aes(x = Estimate, y = z.value, label=feature)) +
  geom_point(aes(color = alpha), alpha = 0.2,size=1.5) +   # base layer
  geom_point(
    data = subset(Res_100_pruned, -log10(padj) > 3),
    aes(color = alpha), size=1.5) +# highlight significant ones
  geom_text(
    data =subset(Res_100_pruned, z.value <= -8),
    size = 3, hjust=1.2)+
  geom_text(
      data =subset(Res_100_pruned, z.value >= 5),
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
ggsave(plot_file)
#print(-log10(0.005))

