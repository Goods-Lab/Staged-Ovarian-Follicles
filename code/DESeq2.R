############################ DESeq2 analysis ################################

## Figure 2A&B

rm(list = ls())

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(limma)

# Load the data
metadata <- read.csv("03.2_staged_follicles_metadata_filtered.csv")
counts <- read.csv("combat-seq_filtereddata_staged_follicles.csv")

metadata_mod <- metadata[-1]
row.names(metadata_mod) <- metadata$sample_id
metadata_mod

counts_mod <- counts[-1]

# remove replicate genes
counts_mod <- counts[!duplicated(counts$X), ]

#remove NA's in gene column
counts_mod <- counts_mod[complete.cases(counts_mod$X), ]

#set row names to unique gene names
row.names(counts_mod) <- counts_mod$X

# extract only numeric data from tpm_mod
counts_mod <- counts_mod[, sapply(counts_mod, is.numeric)]

# remove Gm genes
rows_to_remove <- grepl("^Gm", rownames(counts_mod))
# Remove those rows
counts_mod <- counts_mod[!rows_to_remove, ]

# remove outlier
outlier <- c("D20_250010")
metadata_mod <- metadata_mod[!(rownames(metadata_mod) %in% outlier),]

# take only the intersection of counts and metadata samples
same_samples <-intersect(rownames(metadata_mod), colnames(counts_mod))
metadata_mod <- metadata_mod[rownames(metadata_mod) %in% same_samples, ]
counts_mod <- counts_mod[, colnames(counts_mod) %in% same_samples]

# Count matrix input
cts <- as.matrix(apply(counts_mod, c(1,2), round))
coldata <- metadata_mod[,c("stage", "follicle_size", "oocyte_size")]

coldata$stage <- factor(coldata$stage)

all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

# Design matrix for DE analysis
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ stage)

dds <- dds[ rowSums(counts(dds)) >= 10, ]

# DE analysis
dds$stage <- relevel(dds$stage, ref = "primordial")
dds <- DESeq(dds)

# Make TPM counts
sizeFactors <- sizeFactors(dds)
normalizedCounts <- counts(dds, normalized=TRUE)
librarySize <- colSums(normalizedCounts)
tpmLikeCounts <- sweep(normalizedCounts, 2, librarySize, "/") * 1e6
# Log-transform (log base e of TPM + 1)
logTpmPlus1 <- log(tpmLikeCounts + 1)
# Output the result
write.csv(logTpmPlus1, "follicle_DE_tpms.csv")

# transform TPM counts to z-scores
zscore <- function(x) {
  (x - rowMeans(x)) / apply(x, 1, sd)
}
z_score <- zscore(logTpmPlus1)
z_score <- as.data.frame(z_score)

# To output normalized counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, "normalized_counts.csv")

# Output DE results
resultsNames(dds) # lists the coefficients
res <- results(dds, name="stage_primary_vs_primordial", alpha=0.01) 
write.csv(res, "stage_secondary_vs_primordial.csv")
df_res <- as.data.frame(res)
df_res$gene <- rownames(df_res)

df_res$Significance <- 'Not significant'
df_res$Significance[df_res$padj < 0.05 & df_res$log2FoldChange > 1] <- 'Upregulated'
df_res$Significance[df_res$padj < 0.05 & df_res$log2FoldChange < -1] <- 'Downregulated'

top_upregulated <- df_res %>%
  filter(Significance == 'Upregulated') %>%
  top_n(10, wt=log2FoldChange)

top_downregulated <- df_res %>%
  filter(Significance == 'Downregulated') %>%
  top_n(10, wt=-log2FoldChange)

top_genes <- bind_rows(top_upregulated, top_downregulated)

volcano_plot <- ggplot(df_res, aes(x=log2FoldChange, y=-log10(padj), color=Significance)) +
  geom_point(alpha=0.4) +
  scale_color_manual(values=c('grey'='grey', 'Upregulated'='red', 'Downregulated'='blue')) +
  labs(title="Volcano plot of primary vs primordial", x="Log2 Fold Change", y="-Log10 Adjusted p-value") +
  theme_bw()

# Add labels for the top genes
volcano_plot <- volcano_plot + 
  geom_text_repel(data=top_genes, 
                  aes(label=gene), 
                  size=3.5, 
                  box.padding=0.5, 
                  point.padding=0.3,
                  max.overlaps = 20)

# Print the plot
print(volcano_plot)

# Data transformation and visualization
# Blind dispersion estimation
# Extracting transformed values
vsd <- vst(dds, blind=TRUE)
# rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

# PCA with ggplot function
library(ggplot2)
pcaData <- plotPCA(vsd, intgroup=c("stage"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
stage_colors <- c("primordial" = "#8DB600", "primary" = "#A347FF", "secondary" = "#E30B5D")
ggplot(pcaData, aes(PC1, PC2, color=stage)) +
  geom_point(size=3) +
  theme_minimal() + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values=stage_colors) +
  theme(legend.position="top",
        plot.title=element_text(hjust=0.5),
        panel.border=element_rect(color="black", fill=NA, linewidth=1)) +
  coord_fixed() 

# Extract genes associated with PC1 and PC2
# Extracting transformed values
vst_mat <- assay(vsd)

# Compute PCA
pca_res <- prcomp(t(vst_mat))

# Get the principal components (scores)
pc1 <- pca_res$x[, 1]
pc2 <- pca_res$x[, 2]

# Extract loadings
loadings <- pca_res$rotation

# Sorting and selecting genes based on PC1 and PC2 loadings
# Get indices of top 10 positive and top 10 negative loadings for PC1
top_10_positive_PC1 <- order(pca_res$rotation[, "PC1"], decreasing = TRUE)[1:10]
top_10_negative_PC1 <- order(pca_res$rotation[, "PC1"], decreasing = FALSE)[1:10]

# Get indices of top 10 positive and top 10 negative loadings for PC2
top_10_positive_PC2 <- order(pca_res$rotation[, "PC2"], decreasing = TRUE)[1:10]
top_10_negative_PC2 <- order(pca_res$rotation[, "PC2"], decreasing = FALSE)[1:10]

# Combine indices and extract the actual loadings for plotting
selected_genes_PC1 <- c(top_10_positive_PC1, top_10_negative_PC1)
selected_genes_PC2 <- c(top_10_positive_PC2, top_10_negative_PC2)
loadings_selected_PC1 <- pca_res$rotation[selected_genes_PC1, "PC1"]
loadings_selected_PC2 <- pca_res$rotation[selected_genes_PC2, "PC2"]

# Bar plot for selected genes in PC1
ggplot(data = data.frame(Gene = names(loadings_selected_PC1), Loading = loadings_selected_PC1), aes(x = reorder(Gene, Loading), y = Loading)) +
  geom_bar(stat = "identity", fill = ifelse(loadings_selected_PC1 > 0, "darkgrey", "darkgrey")) +
  theme_minimal() +
  labs(title = "Top 10 Positive and Negative Genes in PC1", x = "Gene", y = "Loading") +
  coord_flip()  # Flips the axes for better visualization of gene names

# Bar plot for selected genes in PC2
ggplot(data = data.frame(Gene = names(loadings_selected_PC2), Loading = loadings_selected_PC2), aes(x = reorder(Gene, Loading), y = Loading)) +
  geom_bar(stat = "identity", fill = ifelse(loadings_selected_PC2 > 0, "darkgreen", "orange")) +
  theme_minimal() +
  labs(title = "Top 10 Positive and Negative Genes in PC2", x = "Gene", y = "Loading") +
  coord_flip()  # Flips the axes for better visualization of gene names

# PCA plot using Generalized PCA
library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$stage <- dds$stage

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = stage)) +
  geom_point(size =3) + coord_fixed() # + ggtitle("glmpca - Generalized PCA")

# Distance between samples
library("pheatmap")
library("RColorBrewer")

sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$stage, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# MDS plot
library(dplyr)
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = stage)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS")

# Gene clustering
library("genefilter")
library(SummarizedExperiment)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 60)

mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("stage", "oocyte_size")])
rownames(anno) <- colnames(mat)
colnames(anno) <- c("stage", "oocyte_size")
anno$stage <- factor(anno$stage, levels = c('primordial', 'primary', 'secondary'))
anno <- anno[order(anno$oocyte_size),]
mat <- mat[,rownames(anno)]
pheatmap(mat, annotation_col = anno[,"stage", drop = FALSE], cluster_cols = FALSE, clustering_distance_rows = "euclidean",
         fontsize=6, breaks = seq(-5, 5, length.out = 101), labels_col = anno$oocyte_size)

# MA plot
plotMA(res, alpha = 0.01, main = "secondary_vs_primary")
