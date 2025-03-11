############################ TPM global PCA ################################

rm(list = ls())
# Load the data
tpm <- read.csv("~/Dartmouth College Dropbox/Hannah VanBenschoten/Staged follicles/data/follicle_DE_tpms.csv")
filtered_data <- read.csv("03.2_staged_follicles_metadata_filtered.csv")

# remove replicate genes
tpm_mod <- tpm[!duplicated(tpm$X), ]

#remove NA's in gene column
tpm_mod <- tpm_mod[complete.cases(tpm_mod$X), ]

#set row names to unique gene names
row.names(tpm_mod) <- tpm_mod$X

# extract only numeric data from tpm_mod
numeric_tpm <- tpm_mod[, sapply(tpm_mod, is.numeric)]

# set row names on filtered metadata
metadata_mod <- filtered_data
row.names(metadata_mod) <- filtered_data$sample_id

# take only the intersection of tpm and metadata samples
same_samples <-intersect(rownames(metadata_mod), colnames(numeric_tpm))
metadata_mod <- metadata_mod[rownames(metadata_mod) %in% same_samples, ]
numeric_tpm <- numeric_tpm[, colnames(numeric_tpm) %in% same_samples]


# Make eSet
library(Biobase)
tpm.coding.unique.matrix.log <- as.matrix(numeric_tpm)
eSet.log2TPM <- ExpressionSet(assayData = tpm.coding.unique.matrix.log, phenoData = AnnotatedDataFrame(data = metadata_mod))

### PCA analysis
tpm.trim.genes.values <- t(exprs(eSet.log2TPM))

#set infinite values to 0
tpm.trim.genes.values[is.infinite(tpm.trim.genes.values)] <- 0

# Remove columns that all 0
tpm.trim.genes.values <- tpm.trim.genes.values[, which(colSums(tpm.trim.genes.values) != 0)]

# Run PCA
pca.tpm <- prcomp((tpm.trim.genes.values), scale=TRUE)

# PCA plots
# By stage code
library(ggbiplot)
library(ggplot2)
        
ggbiplot(pca.tpm, obs.scale=1, var.scale=1,
         groups=pData(eSet.log2TPM)$stage,
         var.axes=FALSE, labels.size = 4) +
  ggtitle("PCA by follicle stage") +
  ggeasy::easy_center_title()


# By batch
ggbiplot(pca.tpm, obs.scale = 1, var.scale = 1,
         groups = pData(eSet.log2TPM)$batch,
         var.axes = FALSE, labels.size = 4) +
  ggtitle("PCA by batch") +
  ggeasy::easy_center_title()

# By pool
ggbiplot(pca.tpm, obs.scale=1, var.scale=1,
         groups=pData(eSet.log2TPM)$pool,
         var.axes=FALSE, labels.size = 4) +
  ggtitle("PCA by pool") +
  ggeasy::easy_center_title()

# By follicle size
ggbiplot(pca.tpm, obs.scale=1, var.scale=1,
         groups=pData(eSet.log2TPM)$follicle_size,
         var.axes=FALSE, labels.size = 4) +
  ggtitle("PCA by follicle size") +
  ggeasy::easy_center_title()
# 
# # By oocyte size
# ggbiplot(pca.tpm, obs.scale=1, var.scale=1, 
#          groups=pData(eSet.log2TPM)$oocyte_size, 
#          var.axes=FALSE, labels.size = 4) +
#   ggtitle("PCA by oocyte size") +
#   ggeasy::easy_center_title() 
# 
# By ovary side
ggbiplot(pca.tpm, obs.scale=1, var.scale=1,
         groups=as.factor(pData(eSet.log2TPM)$ovary_side),
         var.axes=FALSE, labels.size = 4) +
  ggtitle("PCA by ovary side") +
  ggeasy::easy_center_title()

############ Pull out genes that corr with PC 1 #########
# var.loadings.tpm <- pca.tpm$rotation
# gene.loading.PC1 <- as.data.frame(pca.tpm$rotation[,1])
# colnames(gene.loading.PC1) <- "PC1"
# gene.loading.PC1$genes <- row.names(gene.loading.PC1)
# gene.loading.PC1.ordered <- arrange(gene.loading.PC1,gene.loading.PC1$PC1)
# View(gene.loading.PC1.ordered)
# dim(gene.loading.PC1.ordered)
# 
# PC1.top <- gene.loading.PC1.ordered[1:20,]
# PC1.bottom <- gene.loading.PC1.ordered[39238:39257,]
# PC1.plot <- rbind.data.frame(PC1.top,PC1.bottom)
# row.names(PC1.plot) <- PC1.plot$genes
# PC1.plot$genes <- NULL
# write.table(PC1.plot,file="PC1_plot_20_genes.txt",sep="\t",quote=FALSE,col.names=NA)
# 
# pdf("GeneLoading_PC1.pdf",height = 7,width = 7)
# barplot(PC1.plot[,1],names.arg = row.names(PC1.plot),main ="PC1 Genes",ylab="PC loading",
#         cex.names = 0.65, col="grey",las=2) 
# dev.off()

