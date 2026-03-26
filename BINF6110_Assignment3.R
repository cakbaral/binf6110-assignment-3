##**************************
## BINF*6110 - Assignment 3
##
## Student Name: Cyrus Akbarally
##
## Student Number: 1099054
##
## 2026-03-27
##
##**************************

rm(list = ls())

## _ Install packages----
## install.packages("tidyverse")
## install.packages("viridis")
## install.packages("vegan")
## install.packages("dplyr")
## install.packages("ggplot2")
## install.packages("BiocManager")
## BiocManager::install("phyloseq")
## BiocManager::install("biomformat")
## BiocManager::install("microbiome")
## BiocManager::install("rhdf5")
## BiocManager::install("ANCOMBC")


#--


## _ Packages used----
library(tidyverse)
library(readr)
library(viridis)
library(vegan)
library(dplyr)
library(ggplot2)
library(ANCOMBC)
library(phyloseq)
library(rhdf5)
library(biomformat)
library(microbiome)
## setwd("~/BINF6110 Genomic Methods for Bioinformatics/Assignment_3")


#--


### 1. IMPORTING AND PREPARING KRAKEN2/BRACKEN TAXON ABUNDANCE TABLES----

# Import BIOM file with Kraken2/Bracken results & metadata
human_ps <- import_biom("human_gut_table.biom")

human_diet <- as.data.frame(read_csv("human_gut_diet.csv", col_types = cols(Sample = col_character())))


# Check sample names on phyloseq object, replace with names from metadata table, and attach metadata
sample_names(human_ps)
sample_names(human_ps) <- human_diet$Sample

rownames(human_diet) <- human_diet$Sample
human_diet$Sample <- NULL

metadata <- sample_data(human_diet)
metadata$Sample <- rownames(metadata)
sample_data(human_ps) <- metadata

## Taxonomy reference table
taxonomy_reference <- as.data.frame(tax_table(human_ps))


#--


### 2. VISUALIZING RELATIVE TAXONOMIC ABUNDANCE OF SAMPLES----


# Select for top 20 taxa across all samples
top_20_taxa <- top_taxa(human_ps, n = 20)
human_ps_filtered <- prune_taxa(top_20_taxa, human_ps)

# Convert to relative abundance
human_ps_rel <- transform_sample_counts(human_ps_filtered, function(x) x / sum(x))

# Plot at species level
rank_names(human_ps)
human_ps_sp <- tax_glom(human_ps_rel, taxrank = "Rank7")

# Convert to df
human_ps_df <- psmelt(human_ps_sp)

# Plot relative abundance of samples
ggplot(human_ps_df, aes(x = Sample, y = Abundance, fill = Rank7)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Diet, scales = "free_x") +
  labs(title = "Relative Abundance of Human Gut Microbiome for Italian Vegan and Omnivore Diets (Top 20 Species)", y = "Relative Abundance", x = "Sample", fill = "Species")


#--


### 3. ALPHA DIVERSITY (SHANNON AND SIMPSON INDEX)----


# Alpha diversity and visualization with phyloseq

# Calculate Shannon, Simpson, and Chao1 indices for all samples and attach metadata
human_alpha <- as.data.frame(estimate_richness(human_ps, measures = c("Shannon", "Simpson", "Chao1")))
human_alpha$Diet <- human_diet$Diet

# Calculate average indices and standard deviation for omnivore and vegan diet groups
human_alpha_mean_sd <- as.data.frame(human_alpha %>%
  group_by(Diet) %>%
  summarise(
    Mean_Shannon = mean(Shannon),
    SD_Shannon = sd(Shannon),
    Mean_Simpson = mean(Simpson),
    SD_Simpson = sd(Simpson),
    Mean_Chao1 = mean(Chao1),
    SD_Chao1 = sd(Chao1)
  ))

# Plot alpha diversity

# Shannon index - Taxa evenness
plot_richness(human_ps, measures = "Shannon", color = "variable") +
  facet_wrap(~Diet, scales = "free_x") +
  labs(title = "Alpha Diversity of Human Gut Microbiome for Italian Vegan and Omnivore Diets - Shannon Index", x = "Sample", y = "Shannon Index", color = "Index")

# Simpson index - Taxa dominance
plot_richness(human_ps, measures = "Simpson", color = "variable") +
  facet_wrap(~Diet, scales = "free_x") +
  labs(title = "Alpha Diversity of Human Gut Microbiome for Italian Vegan and Omnivore Diets - Simpson Index", x = "Sample", y = "Simpson Index", color = "Index")

# Chao1 index - Taxa richness
plot_richness(human_ps, measures = "Chao1", color = "variable") +
  facet_wrap(~Diet, scales = "free_x") +
  labs(title = "Alpha Diversity of Human Gut Microbiome for Italian Vegan and Omnivore Diets - Chao1 Index", x = "Sample", y = "Chao1 Index", color = "Index")


#--


### 4. BETA DIVERSITY (JACCARD AND BRAY-CURTIS SIMILARITY)----


# Determine Jaccard and Bray-Curtis distance matrices
human_jaccard <- phyloseq::distance(human_ps, method = "jaccard")
human_bray_curtis <- phyloseq::distance(human_ps, method = "bray")

# PCoA ordination
human_pcoa_jaccard <- ordinate(human_ps, method = "PCoA", distance = "jaccard")
human_pcoa_bray_curtis <- ordinate(human_ps, method = "PCoA", distance = "bray")


# PCoA Plots for Jaccard and Bray-Curtis Beta Diversity Measures
plot_ordination(human_ps, human_pcoa_jaccard, color = "Diet", title = "Jaccard Beta Diversity PCoA of Italian Vegan and Omnivore Human Gut Microbiomes") +
  geom_point(size = 4) +
  stat_ellipse(aes(color = Diet, group = Diet), linetype = 2)

plot_ordination(human_ps, human_pcoa_bray_curtis, color = "Diet", title = "Bray-Curtis Beta Diversity PCoA of Italian Vegan and Omnivore Human Gut Microbiomes") +
  geom_point(size = 4) +
  stat_ellipse(aes(color = Diet, group = Diet), linetype = 2)


# Run PERMANOVA

# Metadata for statistical test
meta_test <- data.frame(sample_data(human_ps))

# Conduct statistical test and view results
jaccard_permanova <- adonis2(human_jaccard ~ Diet, data = meta_test, permutations = 999)
bray_curtis_permanova <- adonis2(human_bray_curtis ~ Diet, data = meta_test, permutation = 999)

jaccard_permanova
bray_curtis_permanova


#--


### 5. DIFFERENTIAL ABUNDANCE (ANCOM-BC2)----


# Run ANCOM-BC2 on the human gut samples
human_ancombc.out <- ancombc2(
  data = human_ps, tax_level = "Rank7",
  fix_formula = "Diet", rand_formula = NULL,
  p_adj_method = "holm", pseudo_sens = TRUE,
  prv_cut = 0, lib_cut = 1000, s0_perc = 0.05,
  group = "Diet", struc_zero = TRUE, neg_lb = TRUE
)

human_ancombc_results <- as.data.frame(human_ancombc.out$res)


# Subset results by p-value < 0.05 for the significant taxa
human_ancombc_top <- subset(human_ancombc.out$res, p_DietVegan < 0.05)
human_ancombc_top


# Barplot of log fold changes of taxa from differential abundance
ggplot(human_ancombc_top, aes(x = lfc_DietVegan, y = reorder(taxon, lfc_DietVegan))) +
  geom_point(aes(color = p_DietVegan < 0.05), size = 3) +
  geom_errorbar(aes(
    xmin = lfc_DietVegan - se_DietVegan,
    xmax = lfc_DietVegan + se_DietVegan
  )) +
  geom_vline(xintercept = 0, color = "red") +
  labs(
    title = "Differential Abundance of Taxa in Italian Vegan and Omnivore Human Gut Microbiomes - ANCOM-BC2", x = "Log Fold Change (DietVegan)",
    y = "Species"
  )
