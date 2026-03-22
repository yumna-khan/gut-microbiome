# Genomics Assignment 3

# This R script analyses shotgun metagenomic data comparing vegan and omnivore gut microbiomes. It imports Bracken species-level abundances into a phyloseq object, performs taxonomic visualization, calculates alpha and beta diversity, and conducts differential abundance testing using ANCOM-BC2.


# Import libraries
library(BiocManager)
library(phyloseq)
library(tidyverse)
library(vegan)
library(microbiome)
library(ANCOMBC)

#------------------------ Relative Abundance -------------------------

# Create metadata for the samples
metadata = data.frame(
    sample = c("SRR8146963", "SRR8146968", "SRR8146973",
               "SRR8146935", "SRR8146936", "SRR8146938"),
    diet = c("vegan", "vegan", "vegan",
             "omnivore", "omnivore", "omnivore")
  ) %>% column_to_rownames("sample")

# Read in bracken reports
read_bracken <- function(file) {
  sample_name <- gsub("_bracken.report", "", basename(file))
  read.delim(file, header = FALSE,
             col.names = c("percent", "reads_covered", "reads_assigned", "rank", "taxid", "name")) %>%
    filter(rank == "S") %>%
    select(name, reads_assigned) %>%
    mutate(name = trimws(name), 
           sample = sample_name)
}

# Load all files and make OTU table
otu <- list.files(pattern = "_bracken.report") %>%
  map_dfr(read_bracken) %>%
  pivot_wider(names_from = sample, values_from = reads_assigned, values_fill = 0) %>%
  column_to_rownames("name") %>%
  as.matrix()

# Create phyloseq object
physeq <- phyloseq(
  otu_table(otu, taxa_are_rows = TRUE),
  sample_data(metadata)
)

physeq

# Convert to relative abundance
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# Get top 10 taxa
top10 <- names(sort(taxa_sums(physeq), decreasing = TRUE)[1:10])
physeq_top10 <- prune_taxa(top10, physeq_rel)

# Add diet label to sample data
sample_data(physeq_top10)$label <- paste0(
  sample_data(physeq_top10)$diet, 
  "\n(", 
  rownames(sample_data(physeq_top10)), 
  ")"
)

# Plot
plot_bar(physeq_top10, fill = "OTU") +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = sample_data(physeq_top10)$label) +
  labs(title = "Top 10 Most Abundant Species by Diet Group",
       x = "Sample",
       y = "Relative Abundance",
       fill = "Species") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 7))


#------------------------- Alpha Diversity --------------------------
# Calculate alpha diversity
alpha_div <- estimate_richness(physeq, measures = c("Observed", "Shannon"))

# Add diet group to the results
alpha_div$diet <- sample_data(physeq)$diet
alpha_div$sample <- rownames(alpha_div)

# Plot
alpha_div %>%
  pivot_longer(cols = c("Observed", "Shannon"), 
               names_to = "metric", 
               values_to = "value") %>%
  ggplot(aes(x = diet, y = value, fill = diet)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Alpha Diversity by Diet Group",
       x = "Diet",
       y = "Value") +
  theme_classic() +
  theme(legend.position = "none")

# Summary of values from boxplot
alpha_div %>%
  group_by(diet) %>%
  summarise(
    median_observed = median(Observed),
    range_observed = paste(min(Observed), max(Observed), sep = "-"),
    median_shannon = median(Shannon),
    range_shannon = paste(round(min(Shannon),2), round(max(Shannon),2), sep = "-")
  )


# Wilcoxon test for each metric
wilcox.test(Observed ~ diet, data = alpha_div)
wilcox.test(Shannon ~ diet, data = alpha_div)


#------------------------- Beta Diversity --------------------------
# Calculate Bray-Curtis dissimilarity
bray <- phyloseq::distance(physeq_rel, method = "bray")

# PCoA ordination
pcoa <- ordinate(physeq_rel, method = "PCoA", distance = bray)

# Extract PCoA coordinates
pcoa_df <- as.data.frame(pcoa$vectors[, 1:2])

# Add diet based on rownames
pcoa_df$diet <- sample_data(physeq_rel)$diet

# Check if it looks right
print(pcoa_df)

# Plot
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = diet)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("omnivore" = "salmon", "vegan" = "blue")) +
  labs(title = "Beta Diversity - Bray-Curtis PCoA",
       x = "PC1 [48.3%]",
       y = "PC2 [23.6%]",
       color = "Diet") +
  theme_classic()

# PERMANOVA
metadata_df <- data.frame(diet = sample_data(physeq_rel)$diet)
otu_df <- as.data.frame(t(otu_table(physeq_rel)))

permanova <- adonis2(otu_df ~ diet, data = metadata_df,
                     permutations = 999,
                     method = "bray")
print(permanova)

#----------------------- Differential Abundance -----------------------
# Differential Abundance
# Run ANCOM-BC2
ancom_result <- ancombc2(
  data = physeq,
  fix_formula = "diet",
  group = "diet",
  p_adj_method = "BH",
  struc_zero = TRUE,
  neg_lb = TRUE
)

# Extract results
res <- ancom_result$res

# Get top 20 taxa by log fold change
top20 <- res %>%
  arrange(p_dietvegan) %>%
  head(20) %>%
  mutate(direction = ifelse(lfc_dietvegan > 0, "Higher in Vegan", "Higher in Omnivore"))

# Plot
ggplot(top20, aes(x = reorder(taxon, lfc_dietvegan),
                  y = lfc_dietvegan,
                  fill = direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Higher in Vegan" = "#2ca25f",
                               "Higher in Omnivore" = "#de2d26")) +
  labs(title = "Differential Abundance (ANCOM-BC2): Omnivore vs Vegan",
       x = "Species",
       y = "Log Fold Change (Vegan vs Omnivore)",
       fill = "Direction") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom")

