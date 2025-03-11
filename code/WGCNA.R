# generate heatmaps from individual modules

rm(list=ls())

library(tidyverse)
library(magrittr) 
library(DESeq2)
library(WGCNA)

# Load the data
data <- read.csv("Staged_Follicles_counts_DEGlist.csv", row.names=1)
metadata <- read.csv("03.2_staged_follicles_metadata_filtered.csv")

# take only the intersection of counts and metadata samples
same_samples <-intersect(metadata$sample_id, colnames(data))
metadata <- metadata[metadata$sample_id %in% same_samples, ]

# Make the data in the order of the metadata
data <- data %>%
  dplyr::select(metadata$sample_id)

# Check if this is in the same order
all.equal(colnames(data), metadata$sample_id)

# Remove rows with less than 10 counts
data <- data[rowSums(data[, -1]) > 10, ]

# Transform condition as a factor
metadata$stage <- factor(metadata$stage)
levels(metadata$stage)

# Create a `DESeqDataSet` object
dds <- DESeqDataSetFromMatrix(
  countData = round(data), # Our prepped data frame with counts
  colData = metadata, # Data frame with annotation for our samples
  design = ~stage
)

# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
dds_norm <- vst(dds)

# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(dds_norm) %>%
  t() # Transpose this data

# Determine parameters for WGCNA
sft <- pickSoftThreshold(normalized_counts,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)

sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

# Plot the model fitting by the power soft threshold so we can decide on a soft-threshold for power.
the_plot <- ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()
the_plot 

# Choose power parameter close to 0.8.
# Run WGCNA
cor <- WGCNA::cor
bwnet <- blockwiseModules(normalized_counts,
                          maxBlockSize = 20000, # What size chunks (how many genes) the calculations should be run in
                          TOMType = "signed", # topological overlap matrix
                          power = 9, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          randomSeed = 1234, # there's some randomness associated with this calculation
                          # so we should set a seed
)

# Write main WGCNA results object to file
#readr::write_rds(bwnet,
#                 file = file.path("results", "follicle_wgcna_results.RDS")
#)

# Explore results
module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)

all.equal(metadata$sample_id, rownames(module_eigengenes))

# Create the design matrix from the `condition` variable
des_mat <- model.matrix(~ metadata$stage)

# lmFit() needs a transposed version of the matrix
library(limma)
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

head(stats_df)

# Module 3 seems to be the most differentially expressed across condition groups. 
# Now we can do some investigation into this module.
module_3_df <- module_eigengenes %>%
  tibble::rownames_to_column("Sample") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metadata %>%
                      dplyr::select(sample_id, stage),
                    by = c("Sample" = "sample_id")
  )

ggplot(
  module_3_df,
  aes(
    x = stage,
    y = ME8,
    color = stage
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()

# What genes are a part of module 1?
gene_module_key <- tibble::enframe(bwnet$colors, name = "gene_name", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))

# gene_module_key %>%
#   dplyr::filter(module == "ME3")

# Save results
# readr::write_tsv(gene_module_key,
#                  file = file.path("ME_wgcna_gene_to_module_ECTO.tsv"))

# Make a custom heatmap function
make_module_heatmap <- function(module_name,
                                expression_mat = normalized_counts,
                                metadata_df = metadata,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  
  # Set up the module eigengene with its Sample_ID
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("sample_id")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the condition and sample ID columns
    dplyr::select(sample_id, stage) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "sample_id") %>%
    # Arrange by patient and time point
    dplyr::arrange(stage) %>%
    # Store sample
    tibble::column_to_rownames("sample_id")
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply condition labels
    stage = col_annot_df$stage,
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in time_point
    col = list(stage = c("primordial" = "#8DB600", "primary" = "#A347FF",
                         "secondary" = "#E30B5D"))
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene_name)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )
  
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = FALSE,
                                     show_column_names = FALSE)
  
  # Return heatmap
  return(heatmap)
}

mod_1_heatmap <- make_module_heatmap(module_name = "ME8")

# Print out the plot
mod_1_heatmap

# Write genes from current module into the file
write.csv(mod_5_heatmap@matrix, "ME11_genes.csv")
