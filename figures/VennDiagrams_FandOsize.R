## Generate plots to show overview of follicle and oocyte size correlation data (genes from EdgeR analysis)
## Figure 4

if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
if (!requireNamespace("grid", quietly = TRUE)) {
  install.packages("grid")
}

library(VennDiagram)
library(grid)

library(dplyr)

rm(list=ls())

# create violin plots that describe follicle and oocyte size by stage

metadata <- read.csv("03.2_staged_follicles_metadata_filtered.csv")

# follicle size violin 

fs <- metadata %>%
  dplyr::select(stage, follicle_size)

fs$stage <- factor(fs$stage, levels = c('primordial', 'primary', 'secondary'))

stage_colors <- c("primordial" = "#8DB600", "primary" = "#A347FF", "secondary" = "#E30B5D")

fs_plot <- ggplot(fs, aes(factor(stage), follicle_size, fill=stage)) +
  geom_violin(scale = "width", adjust = 1, trim = FALSE, alpha = 0.5) +
  geom_jitter(aes(color = stage), size = 1, width = 0.2) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  scale_fill_manual(values = stage_colors) +  # Apply custom colors to violins
  scale_color_manual(values = stage_colors) +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0)) +
  ggtitle("Follicle size distribution") + xlab("Follicle Stage") + ylab("Size (um)")

fs_plot

# oocyte size violin 

os <- metadata %>%
  dplyr::select(stage, oocyte_size)

os$stage <- factor(os$stage, levels = c('primordial', 'primary', 'secondary'))

stage_colors <- c("primordial" = "#8DB600", "primary" = "#A347FF", "secondary" = "#E30B5D")

os_plot <- ggplot(os, aes(factor(stage), oocyte_size, fill=stage)) +
  geom_violin(scale = "width", adjust = 1, trim = FALSE, alpha = 0.5) +
  geom_jitter(aes(color = stage), size = 1, width = 0.2) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  scale_fill_manual(values = stage_colors) +  # Apply custom colors to violins
  scale_color_manual(values = stage_colors) +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0)) +
  ggtitle("Oocyte size distribution") + xlab("Follicle Stage") + ylab("Size (um)")

os_plot

plot_combined <- plot_grid(fs_plot, os_plot, align = "h", nrow = 1)

# Display the combined plot
plot_combined

# Venn diagram of genes correlated to follicle size and oocyte size 

fs_genes <- read.csv("folliclesize_genecorrelations.csv")
os_genes <- read.csv("oocytesize_genecorrelations.csv")

fs_genelist <- fs_genes$X
os_genelist <- os_genes$X

venn_plot <- venn.diagram(
  x = list("Follicle Size" = fs_genelist, "Oocyte Size" = os_genelist),
  category.names = c("Follicle Size", "Oocyte Size"),
  filename = NULL,
  output = TRUE,
  col = c("black", "black"),
  fill = c("lightpink", "#FFD700"),,
  alpha = 0.5,
  label.col = "black",
  cex = 1.5,
  fontfamily = "Arial",
  fontface = "plain",
  cat.col = c("black", "black"),
  cat.cex = 1.5,
  cat.fontfamily = "Arial",
  cat.pos = c(-15, 15), # Adjust positions if needed
  cat.dist = c(0.03, 0.03),
  cat.just = list(c(0.5, 0.5), c(0.5, 0.5))
)

# Draw the Venn diagram
grid.newpage()
grid.draw(venn_plot)

# create gene lists of overalapping and distinct genes

fs_only <- setdiff(fs_genelist, os_genelist)
os_only <- setdiff(os_genelist, fs_genelist)
both <- intersect(fs_genelist, os_genelist)

write.table(fs_only, "FollicleSize_uniquegenes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(os_only, "OocyteSize_uniquegenes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# determine how many genes are positively and negatively correlated with follicle and oocyte size
fsonly_genes <- fs_genes %>%
  dplyr::filter(X %in% fs_only)
fs_positive <- fsonly_genes %>%
  dplyr::filter(t>0)
fs_negative <- fsonly_genes %>%
  dplyr::filter(t<0)

osonly_genes <- os_genes %>%
  dplyr::filter(X %in% os_only)
os_positive <- osonly_genes %>%
  dplyr::filter(t>0)
os_negative <- osonly_genes %>%
  dplyr::filter(t<0)


# compare oocyte and somatic genes from single cell data

ss <- read.csv("SS_counts.csv")
ss_meta <- read.csv("SS_metadata.csv")

# clean up data
row.names(ss_meta) <- ss_meta$SampleID

ss <- ss[!duplicated(ss$Gene), ]
ss <- ss[complete.cases(ss$Gene), ]
row.names(ss) <- ss$Gene

same_samples <-intersect(rownames(ss_meta), colnames(ss))
ss_meta <- ss_meta[rownames(ss_meta) %in% same_samples, ]
ss <- ss[, colnames(ss) %in% same_samples]

# make counts data in the same order as metadata
ss <- ss %>%
  dplyr::select(rownames(ss_meta))

# check if it's in the same order
all.equal(colnames(ss), rownames(ss_meta))

# Count matrix input
cts <- as.matrix(apply(ss, c(1,2), round))
cts <- cts + 1
coldata <- ss_meta[,c("Sample.type", "Stage_Factor")]

coldata$Sample.type <- factor(coldata$Sample.type, levels = c("Oocyte", "GC"))
coldata$Stage_Factor <- factor(coldata$Stage_Factor, levels = c("0_primordial", "1_trans_to_primary", "2_primary", "3_trans_to_secondary", "4_secondary"))

# Design matrix for DE analysis
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Sample.type)

dds <- dds[ rowSums(counts(dds)) >= 10, ]

# DE analysis with primordial follicles as ref
dds$Sample.type <- relevel(dds$Sample.type, ref = "Oocyte")
dds <- DESeq(dds)

# Output DE results
resultsNames(dds) # lists the coefficients
res1 <- results(dds, name="Sample.type_GC_vs_Oocyte", alpha=0.01) 
df_res1 <- as.data.frame(res1)
df_res1$gene <- rownames(df_res1)

oocyte_genes <- df_res1 %>%
  filter(log2FoldChange <= -1.5 & padj < 0.01)

somatic_genes <- df_res1 %>%
  filter(log2FoldChange >= 1.5 & padj < 0.01)

oocyte_genes <- rownames(oocyte_genes)
somatic_genes <- rownames(somatic_genes)

# somatic genes vs follicle size comparison

venn_plot <- venn.diagram(
  x = list("Follicle Size" = fs_only, "Oocyte Size" = os_only, "Cell Type: Somatic" = somatic_genes, "Cell Type: Oocyte" = oocyte_genes),
  category.names = c("Follicle Size", "Oocyte Size", "Cell Type: Somatic", "Cell Type: Oocyte"),
  filename = NULL,
  output = TRUE,
  col = c("black", "black", "black", "black"),
  fill = c("lightpink", "#FFD700" ,"lightgreen", "lightblue"),
  alpha = 0.5,
  label.col = "black",
  cex = 1.5,
  fontfamily = "Arial",
  fontface = "plain",
  cat.col = c("black", "black", "black", "black"),
  cat.cex = 1.5,
  cat.fontfamily = "Arial",
  #cat.pos = c(-15, 15), # Adjust positions if needed
  #cat.dist = c(0.03, 0.03),
  #cat.just = list(c(0.5, 0.5), c(0.5, 0.5))
)

# Draw the Venn diagram
grid.newpage()
grid.draw(venn_plot)

# intersection of top follicle/oocyte size genes and cell-specific genes

fs_sc <- intersect(fs_only, somatic_genes)
os_oc <- intersect(os_only, oocyte_genes)
fs_oc <- intersect(fs_only, oocyte_genes)
os_sc <- intersect(os_only, somatic_genes)

write.table(fs_sc, file = "folliclesize_somaticcell.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(os_oc, file = "oocytesize_oocyte.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(fs_oc, file = "folliclesize_oocyte.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(os_sc, file = "oocytesize_somaticcell.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# positive vs negative correlation of fs_somatic and os_oocyte specific genes
fs_somatic <- fs_genes %>%
  dplyr::filter(X %in% fs_sc)

fs_somatic_positive <- fs_somatic %>%
  dplyr::filter(t > 0)

fs_somatic_negative <- fs_somatic %>%
  dplyr::filter(t < 0)

fs_oocyte <- fs_genes %>%
  dplyr::filter(X %in% fs_oc)

fs_oocyte_positive <- fs_oocyte %>%
  dplyr::filter(t > 0)

fs_oocyte_negative <- fs_oocyte %>%
  dplyr::filter(t < 0)

os_oocyte <- os_genes %>%
  dplyr::filter(X %in% os_oc)

os_oocyte_positive <- os_oocyte %>%
  dplyr::filter(t > 0)

os_oocyte_negative <- os_oocyte %>%
  dplyr::filter(t < 0)

os_somatic <- os_genes %>%
  dplyr::filter(X %in% os_sc)

os_somatic_positive <- os_somatic %>%
  dplyr::filter(t > 0)

os_somatic_negative <- os_somatic %>%
  dplyr::filter(t < 0)


top_fsgenes <- fsonly_genes[1:47,]
top_fsgenes <- top_fsgenes$X

top_osgenes <- osonly_genes[1:15,]
top_osgenes <- top_osgenes$X

fs_sc <- intersect(top_fsgenes, somatic_genes)
os_oc <- intersect(top_osgenes, oocyte_genes)

# plot top 4 genes from follicle size/somatic cell and oocyte size/oocyte intersections

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

# take only relevant columns from metadata
metadata <- metadata[, c("stage", "follicle_size","oocyte_size")]

# merge normalized count data with metadata by row name
comb <- merge(tpm_df, metadata, by = "row.names", all = FALSE)
rownames(comb) <- comb$Row.names
comb <- comb[,-1]

# make count matrices corresponding to each stage/cell type gene list
fsize_df <- comb[, fs_sc]
fsize_df <- merge(fsize_df, metadata, by = "row.names", all = FALSE)
rownames(fsize_df) <- fsize_df$Row.names
fsize_df <- fsize_df[,-1]
fsize_df[, 1:4] <- lapply(fsize_df[, 1:4], function(x) as.numeric(as.character(x)))

osize_df <- comb[, os_oc]
osize_df <- merge(osize_df, metadata, by = "row.names", all = FALSE)
rownames(osize_df) <- osize_df$Row.names
osize_df <- osize_df[,-1]
osize_df[, 1:4] <- lapply(osize_df[, 1:4], function(x) as.numeric(as.character(x)))

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
scatter_plots <- lapply(names(fsize_df)[1:4], function(gene_name) {
  fs_scatter_plot(fsize_df, gene_name)
})

# Arrange scatter plots in a 2 by 2 grid
plot_grid <- cowplot::plot_grid(plotlist = scatter_plots, nrow = 2)

# Open the PDF device with landscape orientation
#pdf(file = "FollicleSizeCorrelation.pdf", width = 11, height = 8.5)

# Print the grid of scatter plots
print(plot_grid)

#dev.off()

# Define a function to create scatter plots for each gene vs oocyte size
os_scatter_plot <- function(data, gene_name) {
  # Fit linear regression
  lm_fit <- lm(data[[gene_name]] ~ oocyte_size, data = data)
  
  # Extract coefficients and R-squared
  coef <- coef(lm_fit)
  r_squared <- summary(lm_fit)$r.squared
  
  # Scatter plot with regression line
  p <- ggplot(data, aes(x = oocyte_size, y = .data[[gene_name]])) +
    geom_point(color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "red", size = 0.5) +
    labs(x = "Oocyte Size", y = paste("Expression of", gene_name)) +
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


# Create a list of scatter plots for each primary gene
scatter_plots <- lapply(names(osize_df)[1:4], function(gene_name) {
  os_scatter_plot(osize_df, gene_name)
})

# Arrange scatter plots in a 2 by 4 grid
plot_grid <- cowplot::plot_grid(plotlist = scatter_plots, nrow = 2)

# Open the PDF device with landscape orientation
#pdf(file = "OocyteSizeCorrelation.pdf", width = 11, height = 8.5)

# Print the grid of scatter plots
print(plot_grid)

#dev.off()
