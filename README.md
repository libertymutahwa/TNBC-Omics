# TNBC-Omics: Differential Gene Expression Analysis of Triple-Negative Breast Cancer Using DESeq2

## 📖 Overview
**TNBC-Omics** is a bioinformatics pipeline for analyzing differential gene expression (DGE) in Triple-Negative Breast Cancer (TNBC) compared to normal tissue samples. This project integrates five publicly available GEO datasets (GSE52194, GSE63582, GSE142731, GSE171957, GSE206998) to identify differentially expressed genes (DEGs) using the DESeq2 package in R. It employs robust data preprocessing, SQL-like metadata management, and advanced visualizations (e.g., volcano plots, heatmaps, PCA) to uncover molecular signatures of TNBC. Additionally, it performs protein-protein interaction (PPI) network analysis using STRINGdb and Gene Ontology (GO) and KEGG pathway enrichment to provide biological insights into TNBC mechanisms.

### Key Features
- **Multi-Dataset Integration**: Combines five GEO datasets for robust TNBC vs. normal DGE analysis.
- **Differential Expression**: Identifies DEGs using DESeq2 with log2 fold change (LFC) shrinkage for reliable results.
- **Visualizations**: Generates boxplots, PCA plots, heatmaps (top 10, 20, 50 DEGs), volcano plots, and dispersion plots.
- **Enrichment Analysis**: Performs GO (Biological Process) and KEGG pathway analysis for top DEGs.
- **PPI Networks**: Constructs protein interaction networks for top DEGs using STRINGdb.
- **Data Science Ready**: Outputs structured CSV files for DEGs and enrichment results, suitable for downstream machine learning or statistical analysis.

## 🛠️ Installation
Follow these steps to set up and run the project locally.

### Prerequisites
- **R**: Version 4.0+ (tested with 4.2.3 or later recommended).
- **R Libraries**: Install required packages via CRAN and Bioconductor:
  ```R
  # Install CRAN packages
  install.packages(c("pheatmap", "ggplot2", "RColorBrewer", "dplyr", "tibble", "readr", "grid", "stringr"))

  # Install Bioconductor packages
  if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install(c("DESeq2", "EnhancedVolcano", "STRINGdb", "AnnotationDbi", "org.Hs.eg.db", "vsn", "clusterProfiler", "enrichplot", "apeglm"))
  ```
