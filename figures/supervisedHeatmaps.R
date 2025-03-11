############################ Supervised Heatmaps ################################

## Figure 2D

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
write.csv(normalized_counts, "follicle_DE_normalized_counts.csv")


# generate heatmap from supervised list of DE genes
# load gene list - generated from R code DEG_List.R

library(ComplexHeatmap)
library(tidyverse)

genes <- read.csv("Staged_Follicles_DEG_list.csv")
genes <- genes[-1]
row.names(genes) <- genes$gene

# remove genes from z_score not present in DEGgene list
DEgenes <- rownames(genes)
z_score_DE <- z_score[DEgenes, , drop = FALSE]

z_mat <- as.matrix(z_score_DE)

dend1 = cluster_between_groups(z_mat, metadata_mod$stage)
metadata_mod$stage <- factor(metadata_mod$stage, levels = c('primordial', 'primary', 'secondary'))

ht = Heatmap(z_mat,
             name = "z",
             cluster_columns=dend1,
             clustering_distance_rows ="pearson",
             clustering_method_rows = "average",
             show_row_dend = FALSE,
             show_row_names = FALSE,
             show_column_names = FALSE,
             row_names_side = "left",
             top_annotation = HeatmapAnnotation(stage = metadata_mod$stage, col = list(stage = c("primordial" = "#8DB600", "primary" = "#A347FF", "secondary" = "#E30B5D"))))

draw(ht, padding = unit(c(0.5, 2, 0.5, 0.5), "cm")) ## see right heatmap in following

# generate heatmap with columns ordered by follicle size
# Sort the metadata by follicle size
metadata_sort <- metadata_mod[order(metadata_mod$oocyte_size), ]

# Reorder the columns of the matrix based on the sorted metadata
z_mat <- z_mat[, rownames(metadata_sort)]

# Define the annotation for the heatmap
ha <- HeatmapAnnotation(stage = metadata_sort$stage,
                        col = list(stage = c("primordial" = "#8DB600",
                                             "primary" = "#A347FF",
                                             "secondary" = "#E30B5D")))

# Create the heatmap with columns ordered by follicle size or oocyte size
ht <- Heatmap(z_mat,
              name = "z",
              cluster_columns = FALSE,
              clustering_distance_rows = "pearson",
              clustering_method_rows = "average",
              show_row_dend = FALSE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              row_names_side = "left",
              top_annotation = ha)

# Draw the heatmap
draw(ht, padding = unit(c(0.5, 2, 0.5, 0.5), "cm"))

# generate heatmap from top & bottom 10 genes in each comparison

# load gene lists
PMDvsPM <- read.csv("stage_primary_vs_primordial_siggenes.csv")
PMDvsSC <- read.csv("stage_secondary_vs_primordial_siggenes.csv")
PMvsSC <- read.csv("stage_secondary_vs_primary_siggenes.csv")

# sort gene lists in descending or ascending order and select 20 most upregulated and 20 most downregulated genes
sorted_PMDvsPM <- PMDvsPM[order(-PMDvsPM$log2FoldChange),]
up_PMDvsPM <- head(sorted_PMDvsPM, 20)
unsorted_PMDvsPM <-PMDvsPM[order(PMDvsPM$log2FoldChange),]
down_PMDvsPM <- head(unsorted_PMDvsPM, 20)
top_PMDvsPM <- rbind(up_PMDvsPM, down_PMDvsPM)

sorted_PMDvsSC <- PMDvsSC[order(-PMDvsSC$log2FoldChange),]
up_PMDvsSC <- head(sorted_PMDvsSC, 20)
unsorted_PMDvsSC <-PMDvsSC[order(PMDvsSC$log2FoldChange),]
down_PMDvsSC <- head(unsorted_PMDvsSC, 20)
top_PMDvsSC <- rbind(up_PMDvsSC, down_PMDvsSC)

# there are only 5 genes in the PM vs SC list, so all are included
top_PMvsSC <- PMvsSC

top_genes <- bind_rows(top_PMDvsPM, top_PMDvsSC, top_PMvsSC) %>%
  distinct(gene, .keep_all = TRUE)

write.csv(top_genes, "Staged_Follicles_topDEG_list.csv")

# remove genes from z_score not present in DEGgene list
topDEgenes <- top_genes$gene
z_score_DE <- z_score[topDEgenes, , drop = FALSE]

z_mat <- as.matrix(z_score_DE)

dend1 = cluster_between_groups(z_mat, metadata_mod$stage)
metadata_mod$stage <- factor(metadata_mod$stage, levels = c('primordial', 'primary', 'secondary'))

ht = Heatmap(z_mat,
             name = "z",
             cluster_columns=dend1,
             clustering_distance_rows ="pearson",
             clustering_method_rows = "average",
             show_row_dend = FALSE,
             show_row_names = TRUE,
             show_column_names = FALSE,
             row_names_side = "left",
             top_annotation = HeatmapAnnotation(stage = metadata_mod$stage,col = list(stage = c("primordial" = "#8DB600", "primary" = "#A347FF", "secondary" = "#E30B5D"))),
             row_names_gp = gpar(fontsize = 6, fontfamily = "sans", fontface = "bold"))

draw(ht, padding = unit(c(0.5, 2, 0.5, 0.5), "cm")) ## see right heatmap in following

## Heatmap of staged follicles vs known primordial/primary single cell signatures from Shikanov et al

# load gene lists

granulosa <- read.csv("granulosa_shikanov.csv")
oocyte <- read.csv("oocyte_shikanov.csv")

# add cell type info and combine lists

granulosa <- granulosa %>%
  mutate(CellType = "granulosa")

oocyte <- oocyte %>%
  mutate(CellType = "oocyte")

allcell <- bind_rows(granulosa, oocyte) %>%
  distinct(Gene.Name, .keep_all = TRUE)

# filter out genes with a score of 1

allcell <- allcell %>%
  filter(Score == 1)

# filter out the z-score matrix to contain only genes contained in allcell
rownames(z_score) <- tolower(rownames(z_score))
allcell$Gene.Name <- tolower(allcell$Gene.Name)

z_score_comp <- z_score %>%
  filter(rownames(z_score) %in% allcell$Gene.Name)

allcell <- allcell %>%
  filter(allcell$Gene.Name %in% rownames(z_score))

z_score_comp <- z_score_comp[match(allcell$Gene.Name, rownames(z_score_comp)), ]
CellType <- allcell$CellType
z_score_ct <- cbind(z_score_comp, CellType)

z_mat <- as.matrix(z_score_comp)

metadata_mod$stage <- factor(metadata_mod$stage, levels = c('primordial', 'primary', 'secondary'))
dend1 = cluster_between_groups(z_mat, metadata_mod$stage)
z_score_ct$CellType <- factor(z_score_ct$CellType, levels = c("oocyte", "granulosa"))

ht = Heatmap(z_mat,
             name = "z",
             cluster_columns=dend1,
             #clustering_distance_rows ="pearson",
             #clustering_method_rows = "average"
             clustering_method_columns = "complete",
             #clustering_method_rows = "complete",
             cluster_rows = FALSE,
             show_row_dend = FALSE,
             show_column_names = FALSE,
             show_row_names = TRUE,
             split = z_score_ct$CellType,
             row_names_side = "left",
             top_annotation = HeatmapAnnotation(stage = metadata_mod$stage, col = list(stage = c("primordial" = "#8DB600", "primary" = "#A347FF", "secondary" = "#E30B5D"))),
             left_annotation = rowAnnotation(celltype = allcell$CellType, col = list(celltype = c("oocyte" = "blue", "granulosa" = "green")), show_legend = FALSE),
             row_names_gp = gpar(fontsize = 6, fontfamily = "sans", fontface = "bold"))

draw(ht, padding = unit(c(0.5, 2, 0.5, 0.5), "cm")) ## see right heatmap in following

## Heatmap of known oocyte/somatic markers from canon

oocyte <- c("Dppa3", "Npm2", "Oog1", "Nobox", "Gdf9", "Bmp15", "Dazl", "Ddx4", "Lhx8", "Pou5f1", "Figla", "Zp2")
somatic <- c("Amh", "Foxl2", "Zeb2", "Nr2f2", "Amhr2", "Fst", "Gatm", "Hmgcs2", "Inha", "Kitl", "Wnt6")
genes <- c(oocyte, somatic)
celltype <- c(rep("oocyte", length(oocyte)), rep("somatic", length(somatic)))
df <- data.frame(gene = genes, cell_type = celltype)

z_score <- z_score[genes, , drop = FALSE]
z_score <- z_score[complete.cases(z_score), ]
z_mat <- as.matrix(z_score)

metadata_mod$stage <- factor(metadata_mod$stage, levels = c('primordial', 'primary', 'secondary'))
dend1 = cluster_between_groups(z_mat, metadata_mod$stage)
df$cell_type <- factor(df$cell_type, levels = c("oocyte", "somatic"))

ht = Heatmap(z_mat,
             name = "z",
             cluster_columns=dend1,
             #clustering_distance_rows ="pearson",
             #clustering_method_rows = "average"
             clustering_method_columns = "complete",
             #clustering_method_rows = "complete",
             cluster_rows = FALSE,
             show_row_dend = FALSE,
             show_column_names = FALSE,
             show_row_names = TRUE,
             split = df$cell_type,
             row_names_side = "left",
             top_annotation = HeatmapAnnotation(stage = metadata_mod$stage, col = list(stage = c("primordial" = "#8DB600", "primary" = "#A347FF", "secondary" = "#E30B5D"))),
             left_annotation = rowAnnotation(celltype = df$cell_type, col = list(celltype = c("oocyte" = "lightpink", "somatic" = "yellow")), show_legend = FALSE),
             row_names_gp = gpar(fontsize = 6, fontfamily = "sans", fontface = "bold"))

draw(ht, padding = unit(c(0.5, 2, 0.5, 0.5), "cm")) ## see right heatmap in following

## Heatmap against known genes by follucilogenesis stage sourced from Pelosi et al 

staged_genes <- read.csv("Pelosi_genelist.csv")

# run this if you lower-cased gene names in z-score df in the previous heatmap
staged_genes$gene <- tolower(staged_genes$gene)

z_score_comp <- z_score %>%
  filter(rownames(z_score) %in% staged_genes$gene)

staged_genes <- staged_genes %>%
  filter(staged_genes$gene %in% rownames(z_score))

z_mat <- as.matrix(z_score_comp)

dend1 = cluster_between_groups(z_mat, metadata_mod$stage)
metadata_mod$stage <- factor(metadata_mod$stage, levels = c('primordial', 'primary', 'secondary'))
staged_genes$phase <- factor(staged_genes$phase, levels = c("primordial", "primary", "secondary", "post-secondary"))
staged_genes$cell.type <- factor(staged_genes$cell.type, levels = c("Somatic", "Germline", "Both/Unknown"))

ht = Heatmap(z_mat,
             name = "z",
             cluster_columns=dend1,
             clustering_distance_rows ="pearson",
             #clustering_method_rows = "average"
             clustering_method_columns = "complete",
             clustering_method_rows = "complete",
             cluster_rows = FALSE,
             show_row_dend = FALSE,
             show_column_names = FALSE,
             show_row_names = TRUE,
             split = staged_genes$phase,
             row_names_side = "left",
             top_annotation = HeatmapAnnotation(stage = metadata_mod$stage, col = list(stage = c("primordial" = "#8DB600", "primary" = "#A347FF", "secondary" = "#E30B5D"))),
             left_annotation = rowAnnotation(stage = staged_genes$phase, col = list(stage = c("primordial" = "blue", "primary" = "green", "secondary" = "red", "post-secondary" = "purple"))),
             right_annotation = rowAnnotation(cell = staged_genes$cell.type, col = list(cell = c("Somatic" = "cyan", "Germline" = "lightblue", "Both/Unknown" = "gray"))),
             row_names_gp = gpar(fontsize = 6, fontfamily = "sans", fontface = "bold"))

draw(ht, padding = unit(c(0.5, 2, 0.5, 0.5), "cm")) ## see right heatmap in following