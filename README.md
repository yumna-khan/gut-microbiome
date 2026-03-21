# Gut Microbiome Metagenomics: Omnivore vs Vegan Dietary Comparison

## General Overview
This project examines differences in gut microbiome composition between vegan and omnivore diets using publicly available shotgun metagenomic sequencing data. Reads were quality controlled with fastp, taxonomically classified using Kraken2, and refined with Bracken to generate species-level abundance profiles. Downstream analyses in R (phyloseq, vegan, ANCOM-BC2) included taxonomic visualization, alpha and beta diversity, and differential abundance testing. While no statistically significant differences were detected, this project demonstrates a reproducible workflow for microbiome analysis.

## Table of Contents
- [Introduction](#introduction)
- [Methods](#methods)
  - [1. Data Description](#1-data-description)
  - [2. Quality Control](#2-quality-control)
  - [3. Genome Assembly](#3-genome-assembly)
  - [4. Alignment](#4-alignment)
  - [5. Variant Calling](#5-variant-calling)
  - [6. Visualization](#6-visualization)
- [Results](#results)
- [Discussion](#discussion)
- [References](#references)

## Introduction
The gut microbiome performs essential metabolic functions including vitamin biosynthesis, breakdown of indigestible compounds, and production of beneficial or harmful metabolites that interact with the host (Filippis et al., 2019). Diet is among the most influential external factors shaping the gut microbiome. The shift towards Western dietary patterns, which is characterized by high fat and protein consumption with low fibre intake, has been associated with a loss of microbial diversity with downstream consequences for human health (Filippis et al., 2019). One genus of particular interest is _Prevotella_, which is typically associated with agrarian diets rich in fruits and vegetables. _Prevotella_ has controversial impacts on the human gut: it is associated with the production of short chain fatty acids (health promoting), improved glucose metabolism, and contains anti-inflammatory effects, yet has also been implicated in inflammatory conditions, insulin resistance, and impaired glucose tolerance (Filippis et al., 2019). This inconsistency suggests that different _Prevotella_ strains may drive opposing metabolic responses, making taxonomic classification at the species level critical to understanding its true role in gut health.

To investigate this, shotgun metagenomics is employed, as it enables genus and species level taxonomic profiling of the entire gut microbial community, and provides a more complete picture of microbiome composition (Filippis et al., 2019). The metagenomic data used in this study are derived from the De Filippis et al. (2019) cohort of Italian omnivores and vegans. While omnivore, vegetarian, and vegan diets each produce distinct gut microbiome signatures, the contrast between omnivores and vegans represents the most pronounced divergence. Fackelmann et al. (2025) demonstrated across 21,561 individuals that microbial profiles most accurately distinguished vegans from omnivores (mean Area Under the Curve (AUC) = 0.90) compared to vegetarians versus omnivores (AUC = 0.82) or vegetarians versus vegans (AUC = 0.84), making these two groups the most informative comparison for investigating dietary differences in _Prevotella_ abundance.

To enable accurate downstream taxonomic classification, raw metagenomic reads must first undergo quality control and preprocessing. Tools such as fastp facilitate this step by integrating quality assessment and adapter trimming within a single pipeline, effectively replacing the traditional combination of FastQC and Trimmomatic while operating faster and producing high-quality outputs (Chen et al., 2018). Following preprocessing, taxonomic classification can be conducted using Kraken2, which has demonstrated higher precision, recall, and F1 scores than MetaPhlAn, with diversity estimates that more closely reflect known community composition (Wood et al., 2019). However, Kraken2 relies on reference database completeness, meaning novel or poorly represented organisms may remain undetected. Species-level abundance can then be refined using Bracken, which redistributes reads assigned to higher taxonomic levels to the species level, enabling accurate abundance estimation even when closely related species are present (Lu et al., 2017).

Once taxonomic profiles are established, they provide a foundation for characterizing microbial diversity within and between samples. Within-sample diversity is commonly described using alpha diversity metrics such as Observed features and the Shannon index. Observed features capture species richness, whereas the Shannon index incorporates both richness and evenness, providing complementary perspectives on community structure. However, a limitation of the Shannon index is that variation in its value cannot be attributed to richness or evenness independently (Cassol et al., 2025). In contrast, differences in community composition across samples can be evaluated through beta diversity. Bray–Curtis dissimilarity is frequently applied because it incorporates taxon abundance, unlike presence–absence metrics such as the Jaccard index. These dissimilarities can then be visualized using principal coordinates analysis (PCoA), with statistical differences between groups often evaluated using PERMANOVA. Nevertheless, PERMANOVA assumes equal dispersion among groups, an assumption that may not hold when sample sizes are small.

To identify taxa that differ between dietary groups, differential abundance analysis can be performed using ANCOM-BC2, which corrects for compositional bias and sampling fraction differences through log-ratio transformations (Lin & Peddada, 2020). Unlike DESeq2, which was originally developed for RNA-seq data and does not fully account for the compositional nature of microbiome datasets, and unlike methods such as LEfSe and edgeR that have been reported to produce elevated false discovery rates; ANCOM-BC2 has demonstrated more consistent and reliable performance across microbiome analyses (Nearing et al., 2022). Additionally, ANCOM-BC2 explicitly accounts for structural zeros (taxa that are entirely absent from one group) which may be particularly relevant in dietary comparisons where certain taxa, such as _Prevotella_ strains, may exhibit diet-specific distributions (Nearing et al., 2022).

Overall, this study aims to characterize and compare the gut microbiome composition of Italian omnivores and vegans in order to investigate how dietary patterns shape microbial diversity and species-level _Prevotella_ abundance.

## Methods
### 1. Data Description


### 2. 

### 3. 

### 4. 

### 5. 

### 6. Visualization


## Results




## Discussion




## References


Wick, R. R., P Howden, B., & P Stinear, T. (2025b, August 28). Autocycler: long-read consensus assembly for bacterial genomes. Oxford Academic. https://academic.oup.com/bioinformatics/article/41/9/btaf474/8242761
