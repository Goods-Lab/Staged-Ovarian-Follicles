## script to generate lists of DE genes across six comparisons of follicle stage

rm(list = ls())

# load data

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Load the data
metadata <- read.csv("03.2_staged_follicles_metadata_filtered.csv")
counts <- read.csv("comBat-seq_filtereddata_staged_follicles.csv")

metadata_mod <- metadata[-1]
row.names(metadata_mod) <- metadata$sample_id

row.names(counts) <- counts$X

# remove Gm genes
rows_to_remove <- grepl("^Gm", rownames(counts))
# Remove those rows
counts <- counts[!rows_to_remove, ]

# remove outlier
outliers <- c("D20_250010", "S166")
metadata_mod <- metadata_mod[!(rownames(metadata_mod) %in% outliers),]

same_samples <-intersect(rownames(metadata_mod), colnames(counts))
metadata_mod <- metadata_mod[rownames(metadata_mod) %in% same_samples, ]
counts <- counts[, colnames(counts) %in% same_samples]

# Check if this is in the same order
all.equal(colnames(data), metadata$Sample_Name)

# Count matrix input
cts <- as.matrix(apply(counts, c(1,2), round))
coldata <- metadata_mod[,c("stage", "sample_type")]

coldata$stage <- factor(coldata$stage)

# Design matrix for DE analysis
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ stage)

dds <- dds[ rowSums(counts(dds)) >= 10, ]

# DE analysis with primary follicles as ref
dds$stage <- relevel(dds$stage, ref = "primordial")
dds <- DESeq(dds)
# Output DE results
resultsNames(dds) # lists the coefficients

res1 <- results(dds, name="stage_primary_vs_primordial", alpha=0.01) 
df_res1 <- as.data.frame(res1)
df_res1$gene <- rownames(df_res1)

upreg_primary_1 <- df_res1 %>%
  filter(log2FoldChange >=1 & padj < 0.01)

upreg_primordial_1 <- df_res1 %>%
  filter(log2FoldChange <= -1 & padj < 0.01)

df_res1_sig <- df_res1 %>%
  filter((log2FoldChange >= 1 | log2FoldChange <= - 1) & padj < 0.01)
write.csv(df_res1_sig, "stage_primary_vs_primordial_siggenes.csv")

res2 <- results(dds, name="stage_secondary_vs_primordial", alpha=0.01) 
df_res2 <- as.data.frame(res2)
df_res2$gene <- rownames(df_res2)

upreg_secondary_1 <- df_res2 %>%
  filter(log2FoldChange >=1 & padj < 0.01)

upreg_primordial_2 <- df_res2 %>%
  filter(log2FoldChange <= -1 & padj < 0.01)

df_res2_sig <- df_res2 %>%
  filter((log2FoldChange >= 1 | log2FoldChange <= - 1) & padj < 0.01)
write.csv(df_res2_sig, "stage_secondary_vs_primordial_siggenes.csv")

# DE analysis with primary follicles as ref
dds$stage <- relevel(dds$stage, ref = "primary")
dds <- DESeq(dds)
# Output DE results
resultsNames(dds) # lists the coefficients

res3 <- results(dds, name="stage_secondary_vs_primary", alpha=0.01) 
df_res3 <- as.data.frame(res3)
df_res3$gene <- rownames(df_res3)

upreg_primary_2 <- df_res3 %>%
  filter(log2FoldChange <= -1 & padj < 0.01)

upreg_secondary_2 <- df_res3 %>%
  filter(log2FoldChange >= 1 & padj < 0.01)

df_res3_sig <- df_res3 %>%
  filter((log2FoldChange >= 1 | log2FoldChange <= - 1) & padj < 0.01)
write.csv(df_res3_sig, "stage_secondary_vs_primary_siggenes.csv")


# take the union of all sig genes
union_deg <- bind_rows(df_res1_sig, df_res2_sig, df_res3_sig) %>%
  distinct(gene, .keep_all = TRUE)

write.csv(union_deg, "Staged_Follicles_DEG_list.csv")

# take union of each upregulated gene list per stage
union_primordial <- bind_rows(upreg_primordial_1, upreg_primordial_2) %>%
  distinct(gene, .keep_all = TRUE)

union_primary <- bind_rows(upreg_primary_1, upreg_primary_2) %>%
  distinct(gene, .keep_all = TRUE)

union_secondary <- bind_rows(upreg_secondary_1, upreg_secondary_2) %>%
  distinct(gene, .keep_all = TRUE)

write.csv(union_primordial, "upregulated_Primordial_genes.csv")
write.csv(union_primary, "upregulated_Primary_genes.csv")
write.csv(union_secondary, "upregulated_Secondary_genes.csv")

# merge primary and secondary follicle genes

union_primsec <- bind_rows(union_primary, union_secondary) %>%
  distinct(gene, .keep_all = TRUE)

write.csv(union_primsec, "upregulated_PrimarySecondary_genes.csv")

# generate a count matrix with only genes from DEG list

filtered_genes <- rownames(union_deg)
counts_mod <- counts[filtered_genes, , drop = FALSE]

write.csv(counts_mod, "Staged_Follicles_counts_DEGlist.csv")
  


