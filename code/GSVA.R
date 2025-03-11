# GSVA for staged follicles
rm(list = ls())
library(GSVA)
library(Biobase)
library(GSEABase)
library(GSVAdata)
library(biomaRt)

# load metadata and count matrix
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

all(rownames(metadata_mod)==colnames(counts_mod))
metadata_mod$stage<- factor(metadata_mod$stage)


write.table(rownames(counts_mod), "gene_symbols.txt", sep = "\t",
            row.names = FALSE,  # Do not include row names
            col.names = FALSE,  # Do not include column headers
            quote = FALSE) 

# put table into DAVID and converted to Human orthologue with Entrez ID
# read DAVID output
human_genes <- read.delim("~/Downloads/conv_A8C92A1B60701724085561374.txt")

# 1. Convert mouse gene names to uppercase
mouse_genes_upper <- toupper(rownames(counts_mod))

# 2. Map mouse genes to human genes and then to Entrez IDs
# Assuming that the mouse gene names correspond to human gene names after uppercase
human_genes$MappedMouseGene <- mouse_genes_upper[match(human_genes$From, mouse_genes_upper)]

# 3. Replace Mouse Gene Names with Entrez IDs
# Create a vector to map the mouse genes to their corresponding Entrez IDs
entrez_map <- setNames(human_genes$To, human_genes$MappedMouseGene)

# Map mouse gene names to their corresponding Entrez IDs
mapped_entrez_ids <- entrez_map[toupper(rownames(counts_mod))]

# Filter out NA values in mapped_entrez_ids
valid_indices <- !is.na(mapped_entrez_ids)

# Get the non-NA Entrez IDs and corresponding mouse genes
valid_entrez_ids <- mapped_entrez_ids[valid_indices]
valid_mouse_genes <- rownames(counts_mod)[valid_indices]

# Subset counts_mod to keep only rows with valid mappings
counts_mod_filtered <- counts_mod[valid_mouse_genes, ]

# Replace the row names in counts_mod_filtered with the corresponding non-NA Entrez IDs
rownames(counts_mod_filtered) <- valid_entrez_ids

# Create AnnotatedDataFrame for phenoData in ExpressionSet
phenoData <- AnnotatedDataFrame(data = metadata_mod)

# Create ExpressionSet
eSet <- ExpressionSet(assayData = as.matrix(counts_mod_filtered),
                      phenoData = phenoData)

# import pathway annotation list from KEGG, REACTOME, and BIOCARTA
data(c2BroadSets)
canonicalC2 <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
                                      grep("^REACTOME", names(c2BroadSets)),
                                      grep("^BIOCARTA", names(c2BroadSets)))]

# calculate GSVA enrichment scores

Params <- gsvaParam(eSet,
                         canonicalC2, minSize=5, maxSize=500,
                         kcdf="Poisson")

gsva <- gsva(Params)

# DE expression of pathways
library(limma)

mod <- model.matrix(~ factor(gsva$stage))
colnames(mod) <- c("primary", "primaryVSprimordial", "primaryVSsecondary")
fit <- lmFit(gsva, mod)
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.01)
summary(res)

tt <- topTable(fit, coef=2, n=Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.01]

DEpwys_es <- exprs(gsva[DEpwys, ])
colorLegend <- c("#A347FF","#8DB600","#E30B5D")
names(colorLegend) <- c("primary", "primordial", "secondary")
sample.color.map <- colorLegend[pData(gsva)[, "stage"]]
names(sample.color.map) <- colnames(DEpwys_es)
sampleClustering <- hclust(as.dist(1-cor(DEpwys_es, method="spearman")),
                           method="complete")
geneSetClustering <- hclust(as.dist(1-cor(t(DEpwys_es), method="pearson")),
                            method="complete")

wrap_text_by_words <- function(text, word_limit) {
  # Replace underscores with spaces
  text <- gsub("_", " ", text)
  words <- strsplit(text, " ")[[1]]
  wrapped_lines <- sapply(seq(1, length(words), by = word_limit), function(i) {
    paste(words[i:min(i + word_limit - 1, length(words))], collapse = " ")
  })
  paste(wrapped_lines, collapse = "\n")
}

word_limit <- 4  # Adjust this number to control the number of words per line

# Apply custom wrapping function to row labels
formatted_row_labels <- sapply(rownames(DEpwys_es), function(x) wrap_text_by_words(x, word_limit))

# Create the heatmap
heatmap(DEpwys_es, 
        ColSideColors = sample.color.map,
        xlab = NA,
        ylab = NA,
        margins = c(1, 3),  # Increase margins for better label visibility
        labRow = formatted_row_labels,  # Use formatted row labels
        labCol = NA,  # Remove column labels
        scale = "row",
        Colv = as.dendrogram(sampleClustering),
        Rowv=as.dendrogram(geneSetClustering),
        cexRow = 0.8,  # Size of row labels
        cex.main = 1.5,            # Size of the main title
        cex.lab = 1.2              # Size of x and y axis labels
)

# import pathway annotation list from KEGG, REACTOME, and BIOCARTA
data(c2BroadSets)
canonicalC2 <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
                             grep("^REACTOME", names(c2BroadSets)),
                             grep("^BIOCARTA", names(c2BroadSets)))]

# calculate GSVA enrichment scores

Params <- gsvaParam(eSet,
                    canonicalC2, minSize=5, maxSize=500,
                    kcdf="Poisson")

gsva <- gsva(Params)

# DE expression of pathways
library(limma)

mod <- model.matrix(~ factor(gsva$stage))
colnames(mod) <- c("primary", "primaryVSprimordial", "primaryVSsecondary")
fit <- lmFit(gsva, mod)
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.05)
summary(res)

tt <- topTable(fit, coef=2, n=Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.01]

DEpwys_es <- exprs(gsva[DEpwys, ])
colorLegend <- c("#A347FF","#8DB600","#E30B5D")
names(colorLegend) <- c("primary", "primordial", "secondary")
sample.color.map <- colorLegend[pData(gsva)[, "stage"]]
names(sample.color.map) <- colnames(DEpwys_es)
sampleClustering <- hclust(as.dist(1-cor(DEpwys_es, method="spearman")),
                           method="complete")
geneSetClustering <- hclust(as.dist(1-cor(t(DEpwys_es), method="pearson")),
                            method="complete")

wrap_text_by_words <- function(text, word_limit) {
  # Replace underscores with spaces
  text <- gsub("_", " ", text)
  words <- strsplit(text, " ")[[1]]
  wrapped_lines <- sapply(seq(1, length(words), by = word_limit), function(i) {
    paste(words[i:min(i + word_limit - 1, length(words))], collapse = " ")
  })
  paste(wrapped_lines, collapse = "\n")
}

word_limit <- 4  # Adjust this number to control the number of words per line

# Apply custom wrapping function to row labels
formatted_row_labels <- sapply(rownames(DEpwys_es), function(x) wrap_text_by_words(x, word_limit))

# Create the heatmap
heatmap(DEpwys_es, 
        ColSideColors = sample.color.map,
        xlab = NA,
        ylab = NA,
        margins = c(1, 3),  # Increase margins for better label visibility
        labRow = formatted_row_labels,  # Use formatted row labels
        labCol = NA,  # Remove column labels
        scale = "row",
        Colv = as.dendrogram(sampleClustering),
        Rowv=as.dendrogram(geneSetClustering),
        cexRow = 0.8,  # Size of row labels
        cex.main = 1.5,            # Size of the main title
        cex.lab = 1.2              # Size of x and y axis labels
)


# import GO ontology annotations
library(org.Hs.eg.db)

goannot <- select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns="GO")
genesbygo <- split(goannot$ENTREZID, goannot$GO)

# calculate GSVA enrichment scores

Params <- gsvaParam(eSet,
                    genesbygo, minSize=10, maxSize=500,
                    kcdf="Poisson")

gsva_es <- gsva(Params)

# DE expression of pathways
mod <- model.matrix(~ factor(gsva_es$stage))
colnames(mod) <- c("primary", "primaryVSprimordial", "primaryVSsecondary")
fit <- lmFit(gsva_es, mod)
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.01)
summary(res)

tt <- topTable(fit, coef=2, n=Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.01]

DEpwys_es <- exprs(gsva_es[DEpwys, ])
colorLegend <- c("#A347FF","#8DB600","#E30B5D")
names(colorLegend) <- c("primary", "primordial", "secondary")
sample.color.map <- colorLegend[pData(gsva_es)[, "stage"]]
names(sample.color.map) <- colnames(DEpwys_es)
sampleClustering <- hclust(as.dist(1-cor(DEpwys_es, method="spearman")),
                           method="complete")
geneSetClustering <- hclust(as.dist(1-cor(t(DEpwys_es), method="pearson")),
                            method="complete")

# convert DEpwys_es labels to GO annotations
library(GO.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
go_terms <- Term(GOTERM[DEpwys])

wrap_text_by_words <- function(text, word_limit) {
  # Replace underscores with spaces
  text <- gsub("_", " ", text)
  words <- strsplit(text, " ")[[1]]
  wrapped_lines <- sapply(seq(1, length(words), by = word_limit), function(i) {
    paste(words[i:min(i + word_limit - 1, length(words))], collapse = " ")
  })
  paste(wrapped_lines, collapse = "\n")
}

word_limit <- 6  # Adjust this number to control the number of words per line

# Apply custom wrapping function to row labels
formatted_row_labels <- sapply(go_terms, function(x) wrap_text_by_words(x, word_limit))

# Create the heatmap
heatmap(DEpwys_es, 
        ColSideColors = sample.color.map,
        xlab = NA,
        ylab = NA,
        margins = c(1, 3),  # Increase margins for better label visibility
        labRow = formatted_row_labels,  # Use formatted row labels
        labCol = NA,  # Remove column labels
        scale = "row",
        Colv = as.dendrogram(sampleClustering),
        Rowv=as.dendrogram(geneSetClustering),
        cexRow = 0.8,  # Size of row labels
        cex.main = 1.5,            # Size of the main title
        cex.lab = 1.2              # Size of x and y axis labels
)



