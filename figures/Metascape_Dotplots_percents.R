## Generating dot plots of Metascape pathways
# load Metascape file for given follicle stage

rm(list = ls())

library(readxl)
library(ggplot2)
library(stringr)
library(cowplot)

# Primordial

primordial_df <- read_excel("~/Downloads/metascape_result.t_ug56k9h.xlsx", sheet = "Enrichment")

# filter out everything but member terms in relevant categories
primordial_df <- primordial_df %>%
  filter(GroupID == "1_Member" | GroupID == "2_Member" | GroupID == "3_Member" | GroupID == "4_Member"
         | GroupID == "5_Member"| GroupID == "6_Member" | GroupID == "8_Member"
         | GroupID == "11_Member" | GroupID == "13_Member" | GroupID == "15_Member"| GroupID == "16_Member"
         | GroupID == "17_Member" | GroupID == "18_Member" | GroupID == "19_Member"| GroupID == "20_Member")

# convert InTerm_InList column to a percentage
fraction_to_percentage <- function(fraction) {
  parts <- strsplit(fraction, "/")[[1]]
  numerator <- as.numeric(parts[1])
  denominator <- as.numeric(parts[2])
  percentage <- (numerator / denominator) * 100
  return(percentage)
}

percentages <- sapply(primordial_df$InTerm_InList, fraction_to_percentage)

primordial_df$Genes_In_Set <- percentages

# filter all terms with less than 30% genes in set
primordial_df <- primordial_df %>%
  filter(Genes_In_Set > 30)

primordial_df$Description <- str_wrap(primordial_df$Description, width = 35)

plot_primordial <- ggplot(primordial_df, aes(x = Genes_In_Set, y = Description)) +
  geom_segment(aes(x = 0, xend = Genes_In_Set, y = Description, yend = Description), color = "grey") +
  geom_point(aes(size = Genes_In_Set, color = LogP)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = -2) +
  scale_size_continuous(range = c(2, 9)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Pathways enriched in primordial follicles",
       x = "Genes in set (%)",
       y = "",
       color = "-log10(P)") +
  geom_point(aes(size = Genes_In_Set), shape = 21, color = "black", fill = NA, stroke = 0.5)

plot_primordial

# Primary and Secondary 

df <- read_excel("~/Downloads/metascape_result.tqq5l3gkn.xlsx", sheet = "Enrichment")

# filter out everything but member terms in relevant categories
df <- df %>%
  filter(GroupID == "1_Member" | GroupID == "2_Member" | GroupID == "3_Member"
         | GroupID == "5_Member"| GroupID == "6_Member" | GroupID == "7_Member" | GroupID == "8_Member" | GroupID == "10_Member"
         | GroupID == "11_Member" | GroupID == "12_Member" |GroupID == "14_Member" |GroupID == "15_Member"| GroupID == "16_Member"
         | GroupID == "17_Member" | GroupID == "18_Member" | GroupID == "19_Member"| GroupID == "20_Member")


percentages <- sapply(df$InTerm_InList, fraction_to_percentage)

df$Genes_In_Set <- percentages

# filter all terms with less than 20% genes in set
df <- df %>%
  filter(Genes_In_Set > 20)

df$Description <- str_wrap(df$Description, width = 35)

plot_secondary <- ggplot(df, aes(x = Genes_In_Set, y = Description)) +
  geom_segment(aes(x = 0, xend = Genes_In_Set, y = Description, yend = Description), color = "grey") +
  geom_point(aes(size = Genes_In_Set, color = LogP)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = -2) +
  scale_size_continuous(range = c(2, 9)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Pathways enriched in primary and secondary follicles",
       x = "Genes in set (%)",
       y = "",
       color = "-log10(P)") +
  geom_point(aes(size = Genes_In_Set), shape = 21, color = "black", fill = NA, stroke = 0.5)

plot_combined <- plot_grid(plot_primordial, plot_secondary, align = "h", nrow = 1)

# Display the combined plot
plot_combined

write.csv(df$Description, "primarysecondary_terms.csv")
write.csv(primordial_df$Description, "primordial_terms.csv")

