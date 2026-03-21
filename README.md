# Gut Microbiome Metagenomics: Omnivore vs Vegan Dietary Comparison

## General Overview
This project examines differences in gut microbiome composition between vegan and omnivore diets using publicly available shotgun metagenomic sequencing data. Reads were quality controlled with fastp, taxonomically classified using Kraken2, and refined with Bracken to generate species-level abundance profiles. Downstream analyses in R (phyloseq, vegan, ANCOM-BC2) included taxonomic visualization, alpha and beta diversity, and differential abundance testing. While no statistically significant differences were detected, this project demonstrates a reproducible workflow for microbiome analysis.

## Table of Contents
- [Introduction](#introduction)
- [Methods](#methods)
  - [1. Data Description](#1-data-description)
  - [2. Quality Control](#2-quality-control)
  - [3. Taxanomic Classification](#3-taxaonmic-classification)
  - [4. Taxanomic Abundance](#4-taxanomic-abundance)
  - [5. Alpha Diversity](#5-alpha-diversity)
  - [6. Beta Diversity](#6-beta-diversity)
  - [7. Differential Abundance](#7-differential-abundance)
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
Raw Illumina shotgun metagenomic sequencing data were obtained from the NCBI Sequence Read Archive (SRA). The dataset consists of six human gut microbiome samples from individuals in Turin, including three vegan and three omnivore samples.

**Vegan Samples:**
- SRR8146963
- SRR8146968
- SRR8146973

**Omnivore samples:**
- SRR8146935
- SRR8146936
- SRR8146938

Data retrieval was performed on a high-performance computing (HPC) system using the SRA Toolkit (v3.0.9). Raw sequencing files were downloaded using the `prefetch` command and subsequently converted to paired-end FASTQ format using `fasterq-dump`.

### 2. Quality Control
Raw paired-end Illumina reads were quality controlled and adapter-trimmed using fastp (v1.0.1) (Chen et al., 2018), which performs quality assessment and filtering in a single step. Adapter sequences for paired-end reads were automatically detected using the `--detect_adapter_for_pe option`. Bases with a Phred quality score below Q20 were removed, corresponding to a base call accuracy of 99%, and reads shorter than 50 bp after trimming were discarded (Wang et al., 2025). Filtered reads were stored in a dedicated directory, with log files and HTML reports saved separately for record-keeping.

### 3. Taxanomic Classification
Taxonomic composition of each sample was determined using Kraken2 (v2.1.6) (Wood et al., 2019), a k-mer–based classifier that assigns reads to taxa with high precision. To provide a database for classification, the Kraken2 standard database (8 Gb, 2023-10-09) was downloaded and extracted:

```wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20231009.tar.gz```

```tar -xzvf k2_standard_08gb_20231009.tar.gz```

Trimmed paired-end reads from fastp were classified against this database. To reduce false positives, the confidence threshold was increased from the default 0 to the recommended 0.15. 

Moreover, to improve species-level abundance estimates, particularly for _Prevotella_ species, Bracken (v3.0) was applied to re-estimate read counts (Lu et al., 2017). Re-estimation was performed at the species level (`-l S`) using a read length of 150 bp (`-r 150`) with a minimum threshold of 10 reads (`-t 10`) to ensure reliable abundance estimation (Lu et al., 2017). 

The resulting Bracken reports were then imported directly into R for downstream analysis; conversion to BIOM format was not required.

### 4. Taxanomic Abundance
Taxonomic abundance was visualized in R (v4.5.2) using the `phyloseq` and `ggplot2` packages. Relative abundance was calculated by normalizing read counts per sample, and the top 10 most abundant species were visualized as stacked bar charts to summarize community composition across vegan and omnivore samples.

### 5. Alpha Diversity
Within-sample diversity was quantified using Observed features and the Shannon index, which capture species richness and the integration of richness and evenness respectively, together providing complementary perspectives on within-sample diversity (Cassol et al., 2025). Alpha diversity metrics were visualized as boxplots in R and statistical differences between dietary groups were evaluated using the Wilcoxon rank-sum test, which was selected as a non-parametric alternative appropriate for small sample sizes and data that cannot be assumed to follow a normal distribution.

### 6. Beta Diversity
Differences in community composition between dietary groups were assessed using Bray-Curtis dissimilarity, which accounts for taxon abundance. Ordination was performed using Principal Coordinates Analysis (PCoA) and statistical significance between groups was assessed with PERMANOVA, implemented via the `vegan` package in R. However, PERMANOVA assumes equal dispersion between groups, which may not hold with small sample sizes and results should therefore be interpreted cautiously.

### 7. Differential Abundance
Differential abundance between vegans and omnivores was evaluated using ANCOM-BC2 (v1.0), which accounts for compositional bias, sampling fraction differences, and structural zeros (Lin & Peddada, 2020; Nearing et al., 2022). Benjamini-Hochberg correction was applied for multiple testing (`p_adj_method = "BH"`), as it controls the false discovery rate and is more appropriate than the more conservative Holm correction when testing many species simultaneously with limited sample sizes. In addition, structural zeros detection was enabled (`struc_zero = TRUE, neg_lb = TRUE`) to identify taxa entirely absent from one dietary group, which is particularly relevant for diet-specific _Prevotella_ species.

## Results




## Discussion




## References


Wick, R. R., P Howden, B., & P Stinear, T. (2025b, August 28). Autocycler: long-read consensus assembly for bacterial genomes. Oxford Academic. https://academic.oup.com/bioinformatics/article/41/9/btaf474/8242761
