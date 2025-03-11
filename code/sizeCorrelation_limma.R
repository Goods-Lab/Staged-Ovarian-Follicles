## Analyze correlations between continuous variable and gene expression across all genes to 
## identify genes with strongest positive or negative correlations

rm(list=ls())
library(limma)
library(edgeR)

# load normalized count data and metadata
counts <- read.csv("02_staged_follicles_counts_annotated_filtered.csv")
metadata <- read.csv("03.2_staged_follicles_metadata_filtered.csv")


metadata <- metadata[-1]
row.names(metadata) <- metadata$sample_id
counts <- counts[!duplicated(counts$external_gene_name), ]
row.names(counts) <- counts$external_gene_name
# remove Gm genes
rows_to_remove <- grepl("^Gm", rownames(counts))
# Remove those rows
counts <- counts[!rows_to_remove, ]
counts <- counts[, sapply(counts, is.numeric)]
counts <- counts[-1]
same_samples <-intersect(rownames(metadata), colnames(counts))
metadata <- metadata[rownames(metadata) %in% same_samples, ]
counts <- counts[, colnames(counts) %in% same_samples]

## FOLLICLE SIZE

# create DGELst object
dge <- DGEList(counts = counts)

# create design matrix for follicle size
design <- model.matrix(~ follicle_size, data = metadata)

# filter rows that have zero or low counts
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# apply scale normalization to read counts using TMM normaliation method
dge <- calcNormFactors(dge)

# convert counts to logCPM values
logCPM <- cpm(dge, log = TRUE, prior.count=3)

# fit linear model to each gene
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend = TRUE)

# extract result for follicle size
tb_folliclesize <- topTable(fit, coef = "follicle_size", number = nrow(logCPM))
sig_genes_fsize <- subset(tb_folliclesize, adj.P.Val < 0.05)

write.csv(sig_genes_fsize, "folliclesize_genecorrelations.csv")

fs_genes <- rownames(sig_genes_fsize)
write.table(fs_genes, file = "follicle_size_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

## OOCYTE SIZE

# create DGELst object
dge <- DGEList(counts = counts)

# create design matrix for ooyte size
design <- model.matrix(~ oocyte_size, data = metadata)

# filter rows that have zero or low counts
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# apply scale normalization to read counts using TMM normaliation method
dge <- calcNormFactors(dge)

# convert counts to logCPM values
logCPM <- cpm(dge, log = TRUE, prior.count=3)

# fit linear model to each gene
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend = TRUE)

# extract result for follicle size
tb_oocytesize <- topTable(fit, coef = "oocyte_size", number = nrow(logCPM))
sig_genes_osize <- subset(tb_oocytesize, adj.P.Val < 0.05)

write.csv(sig_genes_osize, "oocytesize_genecorrelations.csv")

os_genes <- rownames(sig_genes_osize)
write.table(os_genes, file = "oocyte_size_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

## CREATE PLOTS OF TOP 20 GENES IN FOR EACH VARIABLE

# pull out gene lists
fsize_genes <- row.names(sig_genes_fsize[1:20,])
osize_genes <- row.names(sig_genes_osize[1:20,])

# get tpm data to plot for sig genes

library(tidyverse)
library(broom)

tpm <- read.csv("follicle_DE_tpms.csv")
metadata <- read.csv("03.2_staged_follicles_metadata_filtered.csv")

# convert tpm to dataframe
tpm_df <- as.data.frame(t(tpm))
colnames(tpm_df) <- tpm_df[1,]
tpm_df <- tpm_df[-1,]

# merge metadata with counts dataframe
rownames(metadata) <- metadata$sample_id
metadata <- metadata[,-1]

# take only relevent columns from metadata
metadata <- metadata[, c("stage", "follicle_size","oocyte_size")]

# merge normalized count data with metadata by row name
comb <- merge(tpm_df, metadata, by = "row.names", all = FALSE)
rownames(comb) <- comb$Row.names
comb <- comb[,-1]

# make count matrices corresponding to each stage/cell type gene list
fsize_df <- comb[, fsize_genes]
fsize_df <- merge(fsize_df, metadata, by = "row.names", all = FALSE)
rownames(fsize_df) <- fsize_df$Row.names
fsize_df <- fsize_df[,-1]
fsize_df[, 1:20] <- lapply(fsize_df[, 1:20], function(x) as.numeric(as.character(x)))

osize_df <- comb[, osize_genes]
osize_df <- merge(osize_df, metadata, by = "row.names", all = FALSE)
rownames(osize_df) <- osize_df$Row.names
osize_df <- osize_df[,-1]
osize_df[, 1:20] <- lapply(osize_df[, 1:20], function(x) as.numeric(as.character(x)))

# perform linear regression for each gene against follicle size and oocyte size & plot

library(ggplot2)

# Define a function to create scatter plots for each gene vs follicle size
fs_scatter_plot <- function(data, gene_name) {
  # Fit linear regression
  lm_fit <- lm(data[[gene_name]] ~ follicle_size, data = data)
  
  # Extract coefficients and R-squared
  coef <- coef(lm_fit)
  r_squared <- summary(lm_fit)$r.squared
  
  # Scatter plot with regression line
  p <- ggplot(data, aes(x = follicle_size, y = .data[[gene_name]])) +
    geom_point(color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.5) +
    labs(x = "Follicle Size", y = paste("Expression of", gene_name)) +
    ggtitle(paste(gene_name)) +
    theme_bw() +
    theme(panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # Add equation and R-squared value as annotation
  p + annotate("text", x = Inf, y = Inf, label = paste("R-squared:", round(r_squared, 2)),
               hjust = 1.05, vjust = 1.5, size = 3) +
    theme(plot.margin = margin(10, 10, 10, 10, "pt"))
}


# Create a list of scatter plots for each primordial gene
scatter_plots <- lapply(names(fsize_df)[1:20], function(gene_name) {
  fs_scatter_plot(fsize_df, gene_name)
})

# Arrange scatter plots in a 4 by 5 grid
plot_grid <- cowplot::plot_grid(plotlist = scatter_plots, nrow = 4)

# Open the PDF device with landscape orientation
pdf(file = "FollicleSizeCorrelation.pdf", width = 11, height = 8.5)

# Print the grid of scatter plots
print(plot_grid)

dev.off()

# Create a list of scatter plots for each primary gene
scatter_plots <- lapply(names(osize_df)[1:20], function(gene_name) {
  fs_scatter_plot(osize_df, gene_name)
})

# Arrange scatter plots in a 2 by 4 grid
plot_grid <- cowplot::plot_grid(plotlist = scatter_plots, nrow = 4)

# Open the PDF device with landscape orientation
pdf(file = "OocyteSizeCorrelation.pdf", width = 11, height = 8.5)

# Print the grid of scatter plots
print(plot_grid)

dev.off()

