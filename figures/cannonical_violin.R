
## Figure 2C & 3B - violin plot of canonical follicle development genes and DE genes

rm(list=ls())

# load TMP and normalized count data
tpm <- read.csv("follicle_DE_tpms.csv")
metadata <- read.csv("03.2_staged_follicles_metadata_filtered.csv")

metadata_mod <- metadata[-1]
row.names(metadata_mod) <- metadata$sample_id
row.names(tpm) <- tpm$X
tpm <- tpm[-1]
same_samples <-intersect(rownames(metadata_mod), colnames(tpm))
metadata_mod <- metadata_mod[rownames(metadata_mod) %in% same_samples, ]

# generate heatmap from supervised list of genes

library(ComplexHeatmap)
library(tidyverse)
library(cowplot)

# make violin plots of canonical follicle development genes

# ordered by expression
# # canon <- c("Zp3",
#            "Gdf9",
#            "Sohlh1",
#                 "Wnt6",
#            "Foxl2",
#                 "Amh",
#                 "Dazl",
#                 "Fgf8",
#                 "Umodl1",
#                 "Bmp4"
# )

# grouped by expression pattern
canon <- c("Zp3",
           "Gdf9",
           "Bmp15",
           "Dazl",
           "Fgf8",
           "Wnt6",
           "Foxl2",
           "Amh"
)


# transpose count matrix to have samples in rows and genes in columns
tpm_t <- as.data.frame(t(as.matrix(tpm)))

# make row names into a sample ID column
tpm_t <- rownames_to_column(tpm_t, var = "SampleID")

# add stage column to transposed matrix
tpm_t$stage <- metadata_mod$stage
tpm_t$stage <- factor(tpm_t$stage, levels = c('primordial', 'primary', 'secondary'))

tpm_canon <- reshape2::melt(tpm_t, id.vars = c("SampleID", "stage"), measure.vars = canon,
                                variable.name = "Gene", value.name = "Expr")

stage_colors <- c("primordial" = "#8DB600", "primary" = "#A347FF", "secondary" = "#E30B5D")

ggplot(tpm_canon, aes(factor(stage), Expr, fill=stage)) +
  geom_violin(scale = "width", adjust = 1, trim = FALSE, alpha = 0.5) +
  geom_jitter(aes(color = stage), size = 1, width = 0.2) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(Gene), scales = "free", switch = "y") +
  scale_fill_manual(values = stage_colors) +  # Apply custom colors to violins
  scale_color_manual(values = stage_colors) +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0)) +
  ggtitle("Canonical Follicle Development Genes") + xlab("Follicle Stage") + ylab("Expression Level")


## Top variable genes per stage

primordial <- c("Stpg1", "Camk1g", "Tex13b", "Card10", "Ttc22", "Colec11", "Sptb")
primary <- c("Myo7a", "Ank1", "Eya2", "Pabpc1l")
secondary <- c("Shisa8", "Fbxw20", "Mrap")

top_genes <- unique(c(primordial, primary, secondary))

tpm_top <- tpm %>%
  filter(rownames(tpm) %in% top_genes)

# transpose count matrix to have samples in rows and genes in columns
tpm_t <- as.data.frame(t(as.matrix(tpm_top)))

# make row names into a sample ID column
tpm_t <- rownames_to_column(tpm_t, var = "SampleID")

# add stage column to transposed matrix
tpm_t$stage <- metadata_mod$stage
tpm_t$stage <- factor(tpm_t$stage, levels = c('primordial', 'primary', 'secondary'))

tpm_de <- reshape2::melt(tpm_t, id.vars = c("SampleID", "stage"), measure.vars = top_genes,
                            variable.name = "Gene", value.name = "Expr")

stage_colors <- c("primordial" = "#8DB600", "primary" = "#A347FF", "secondary" = "#E30B5D")

ggplot(tpm_de, aes(factor(stage), Expr, fill=stage)) +
  geom_violin(scale = "width", adjust = 1, trim = FALSE, alpha = 0.5) +
  geom_jitter(aes(color = stage), size = 1, width = 0.2) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(Gene), scales = "free", switch = "y") +
  scale_fill_manual(values = stage_colors) +  # Apply custom colors to violins
  scale_color_manual(values = stage_colors) +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0)) +
  ggtitle("Differentially Expressed Genes") + xlab("Follicle Stage") + ylab("Expression Level")

