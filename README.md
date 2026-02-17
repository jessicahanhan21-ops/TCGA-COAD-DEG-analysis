# TCGA-COAD-DEG-analysis
RNA-seq analysis of TCGA-COAD: identifying DEGs, pathway enrichment, and visualization with R.

This repository contains R code for analyzing RNA-seq data from The Cancer Genome Atlas Colon Adenocarcinoma (TCGA-COAD) project. The goal is to identify differentially expressed genes (DEGs) between primary tumor and normal tissue samples.

## Data Source

- **Datasets**: GDC TCGA-COAD (Colon Adenocarcinoma)
- **Samples**: 471 primary tumor, 41 solid tissue normal
- **Data type**: STAR‑Counts (gene-level raw counts) and phenotype files downloaded from ([https://xenabrowser.net/datapages/])

## Methods

1. **Data preprocessing**  1. **Data preprocessing**
   - Counts were filtered (≥10 counts in ≥10% of samples)  
   - DESeq2 was used for normalization and differential expression testing  
2. **Differential expression**  2. **Differential expression**
   - Threshold: adjusted p‑value < 0.05 and |log2 fold change| > 1  
3. **Gene ranking**  
   - Combined score = -log10(padj) × |log2FC|  
   - Genes with baseMean > 100 were prioritized  
4. **Functional enrichment**  
   - Top 500 ranked genes submitted to Enrichr (GO Biological Process 2023)  
5. **Validation**  
   - Checked expression of 17 canonical cell cycle genes (e.g., CDK1, MKI67, MCM2‑7)

## Key Results

- **5,981 significant DEGs** (3,252 upregulated, 2,729 downregulated)
- Top upregulated genes: *MMP7*, *KRT80*, *CDH3*, *FOXQ1* – many previously linked to colorectal cancer
- Enriched pathways: extracellular matrix organization, cell cycle regulation, antimicrobial humoral response
- All 17 cell cycle genes showed expected upregulation, confirming analysis accuracy

## Repository Contents

- `TCGA_COAD_analysis.R` – Main R script with all steps (fully commented)
- `figures/` – Contains volcano plot, heatmap, and enrichment bar plot
- `results_summary.csv` – Top genes ranked by combined score

## Requirements

- R (version ≥ 4.0)
- Packages: `DESeq2`, `ggplot2`, `pheatmap`, `enrichR`, `org.Hs.eg.db`, `dplyr`, `readr`

## How to Run

1. Clone this repository  
2. Download the TCGA-COAD STAR‑Counts and clinical files from the(https://xenabrowser.net/datapages/)
3. Update the file paths in the script to point to your downloaded files  
4. Run the script in R or RStudio  

## Feedback

This project was a learning exercise in applying a standard differential expression workflow to public cancer data. Suggestions and improvements are always welcome!
