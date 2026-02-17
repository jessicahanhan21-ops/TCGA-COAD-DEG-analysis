
# TCGA-COAD Analysis 

library(readr)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(enrichR)

# ================================================================================
# LOAD DATA
# ================================================================================

# 1. Load the expression counts
counts <- read_tsv("C:/Users/jessi/Downloads/TCGA-COAD.star_counts.tsv.gz")

# 2. Load the phenotype 
pheno <- read_tsv ("C:/Users/jessi/Downloads/TCGA-COAD.clinical.tsv.gz")

# ================================================================================
# PREPARE METADATA
# ================================================================================

metadata <- pheno[, c("sample", "sample_type.samples")]
colnames(metadata) <- c("sample_id", "condition")

# Keep only Primary Tumor and Normal
metadata <- metadata[metadata$condition %in% c("Primary Tumor", "Solid Tissue Normal"), ]

# Set factor levels (Normal as reference)
metadata$condition <- factor(metadata$condition, 
                             levels = c("Solid Tissue Normal", "Primary Tumor"))

print(table(metadata$condition))
# ================================================================================
# PREPARE COUNT MATRIX AND REVERSE LOG TRANSFORMATION
# ================================================================================

counts_df <- as.data.frame(counts)
gene_ids <- counts_df[[1]]

# Extract numeric data
count_matrix <- counts_df[, -1]
count_matrix <- apply(count_matrix, 2, as.numeric)
rownames(count_matrix) <- gene_ids

# Check data format
sample_vals <- count_matrix[1, 1:5]
cat("Sample values (BEFORE transformation):\n")
cat(paste(round(sample_vals, 2), collapse=", "), "\n")
cat("Max value:", round(max(count_matrix), 2), "\n")
cat("Median value:", round(median(count_matrix), 2), "\n\n")

# data is log2(count+1) - REVERSE IT!

raw_counts <- round(2^count_matrix - 1)

# Verify transformation
cat("Sample values (AFTER transformation):\n")
cat(paste(raw_counts[1, 1:5], collapse=", "), "\n")
cat("Max value:", max(raw_counts), "\n")
cat("Median value:", median(raw_counts), "\n\n")

if (max(raw_counts) < 100) {
  cat("WARNING: Max count is still very low!\n")
  cat(" Something may be wrong with the transformation.\n\n")
} else {
  cat("✓ Transformation successful! Values now look like raw counts.\n\n")
}
# ================================================================================
# MATCH SAMPLES
# ================================================================================

metadata$sample_id <- gsub("-", ".", metadata$sample_id)
colnames(raw_counts) <- gsub("-", ".", colnames(raw_counts))

common_samples <- intersect(metadata$sample_id, colnames(raw_counts))

cat("Matching samples:", length(common_samples), "\n")

metadata <- metadata[metadata$sample_id %in% common_samples, ]
metadata <- metadata[order(metadata$sample_id), ]

raw_counts <- raw_counts[, common_samples]
raw_counts <- raw_counts[, order(colnames(raw_counts))]

# Align
raw_counts <- raw_counts[, metadata$sample_id]

print(table(metadata$condition))
# ================================================================================
# QUALITY CONTROL
# ================================================================================

genes_before <- nrow(raw_counts)

# Keep genes with at least 10 counts in at least 10% of samples
min_samples <- ceiling(0.1 * ncol(raw_counts))
keep_genes <- rowSums(raw_counts >= 10) >= min_samples
raw_counts_filtered <- raw_counts[keep_genes, ]

genes_after <- nrow(raw_counts_filtered)
cat("Genes before filtering:", genes_before, "\n")
cat("Genes after filtering:", genes_after, "\n")
cat("Removed:", genes_before - genes_after, "\n\n")

# ================================================================================
# RUN DESEQ2
# ================================================================================
# Ensure integer counts
raw_counts_filtered <- round(raw_counts_filtered)

dds <- DESeqDataSetFromMatrix(
  countData = raw_counts_filtered,
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("condition", "Primary Tumor", "Solid Tissue Normal"))

summary(res)
# ================================================================================
# FILTER SIGNIFICANT GENES
# ================================================================================

res_sorted <- res[order(res$padj), ]
sig_genes <- subset(res_sorted, padj < 0.05 & abs(log2FoldChange) > 1)

cat("Significant genes (padj < 0.05, |log2FC| > 1):", nrow(sig_genes), "\n")
cat("  - Upregulated:", sum(sig_genes$log2FoldChange > 1), "\n")
cat("  - Downregulated:", sum(sig_genes$log2FoldChange < -1), "\n\n")

if (nrow(sig_genes) < 1000) {
  cat(" WARNING: Expected 10,000+ genes, got", nrow(sig_genes), "\n")
  cat("   Check if transformation was done correctly.\n\n")
} else {
  cat("✓ Good! Gene count looks normal.\n\n")
}
# ================================================================================
# ANNOTATE GENES
# ================================================================================

sig_genes_df <- as.data.frame(sig_genes)
sig_genes_df$Ensembl_Full <- rownames(sig_genes_df)

clean_ids <- gsub("\\..*", "", sig_genes_df$Ensembl_Full)

sig_genes_df$Symbol <- mapIds(
  org.Hs.eg.db,
  keys = clean_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

sig_genes_df <- sig_genes_df[!is.na(sig_genes_df$Symbol), ]

cat("Genes with symbols:", nrow(sig_genes_df), "\n\n")

# ================================================================================
# VALIDATION - CELL CYCLE GENES
# ================================================================================

cell_cycle_genes <- c("CDK1", "CCNB1", "CCNB2", "CCNA2", "AURKA", "AURKB",
                      "BUB1", "BUB1B", "TOP2A", "PCNA", "MKI67", 
                      "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7")

all_res_df <- as.data.frame(res)
all_res_df$Ensembl_Full <- rownames(all_res_df)
clean_ids_all <- gsub("\\..*", "", all_res_df$Ensembl_Full)
all_res_df$Symbol <- mapIds(org.Hs.eg.db, keys = clean_ids_all, 
                            column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

cell_cycle_check <- all_res_df[all_res_df$Symbol %in% cell_cycle_genes, 
                               c("Symbol", "log2FoldChange", "baseMean", "padj")]
cell_cycle_check <- cell_cycle_check[!is.na(cell_cycle_check$Symbol), ]
cell_cycle_check <- cell_cycle_check[order(-cell_cycle_check$log2FoldChange), ]

print(cell_cycle_check)

positive_count <- sum(cell_cycle_check$log2FoldChange > 0, na.rm=TRUE)
significant_count <- sum(cell_cycle_check$padj < 0.05, na.rm=TRUE)
total_count <- nrow(cell_cycle_check)

cat("  - Positive log2FC:", positive_count, "/", total_count, "\n")
cat("  - Significant (padj < 0.05):", significant_count, "/", total_count, "\n\n")

if (positive_count >= 0.8 * total_count && significant_count >= 10) {
  cat("✓✓✓ VALIDATION PASSED ✓✓✓\n")
  cat("Cell cycle genes are upregulated and significant!\n")
  cat("Analysis is CORRECT!\n")
} else {
  cat("VALIDATION FAILED ️\n")
  cat("Something is still wrong with the analysis!\n")
}
# ================================================================================
# CREATE RANKINGS
# ================================================================================

# By significance
sig_by_padj <- sig_genes_df[order(sig_genes_df$padj), ]

# By fold change
sig_by_fc_up <- sig_genes_df[sig_genes_df$log2FoldChange > 0, ]
sig_by_fc_up <- sig_by_fc_up[order(-sig_by_fc_up$log2FoldChange), ]

sig_by_fc_down <- sig_genes_df[sig_genes_df$log2FoldChange < 0, ]
sig_by_fc_down <- sig_by_fc_down[order(sig_by_fc_down$log2FoldChange), ]

# Combined score
reliable_genes <- sig_genes_df[sig_genes_df$baseMean > 100, ]
reliable_genes$combined_score <- -log10(reliable_genes$padj) * abs(reliable_genes$log2FoldChange)
reliable_genes_sorted <- reliable_genes[order(-reliable_genes$combined_score), ]

print(reliable_genes_sorted[1:min(20, nrow(reliable_genes_sorted)), 
                            c("Symbol", "log2FoldChange", "baseMean", "padj", "combined_score")])

# ================================================================================
# VISUALIZATIONS
# ================================================================================

# Volcano plot

res_plot <- as.data.frame(res)
res_plot$group <- "Not Significant"
res_plot$group[res_plot$padj < 0.05 & res_plot$log2FoldChange > 1] <- "Upregulated"
res_plot$group[res_plot$padj < 0.05 & res_plot$log2FoldChange < -1] <- "Downregulated"

ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = c("Downregulated" = "blue", 
                                "Not Significant" = "grey", 
                                "Upregulated" = "red")) +
  theme_minimal() +
  labs(title = "Differential Gene Expression in Colorectal Cancer",
       subtitle = paste("TCGA-COAD:", nrow(metadata), "samples"),
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-value") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")
---------------------------------------------------------------------------------

# Heat map
  
if (nrow(reliable_genes_sorted) >= 20) {
    
    # 2. Extract Variance Stabilized Data
    vsd <- vst(dds, blind = FALSE)
    
    # 3. Pick the top 20 genes
    top_20 <- reliable_genes_sorted[1:20, ]
    
    # 4. Extract data and ensure Sample Alignment
    heatmap_data <- assay(vsd)[rownames(top_20), ]
    heatmap_data <- heatmap_data[, metadata$sample_id]
    
    # 5. Set row labels to Gene Symbols
    gene_labels <- ifelse(is.na(top_20$Symbol), top_20$Ensembl_Full, top_20$Symbol)
    rownames(heatmap_data) <- gene_labels
    
    # 6. Create the Sample Annotation (Top Color Bar)
    annotation_col <- data.frame(SampleType = metadata$condition)
    rownames(annotation_col) <- metadata$sample_id
    
    # 7. Generate the Heat map
    pheatmap(heatmap_data,
             annotation_col = annotation_col,
             show_colnames = FALSE,       
             cluster_cols = TRUE,         
             scale = "row",               
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
             main = "Top 20 Genes by Combined Score (TCGA-COAD)")
  }
  
# ================================================================================
# ENRICHMENT ANALYSIS
# ================================================================================
# 1. Point the enrichment tool to your BALANCED list

top_balanced_genes <- reliable_genes_sorted$Symbol[1:500]

# 2. Run Enrichment
dbs <- c("GO_Biological_Process_2023")
enriched_balanced <- enrichr(top_balanced_genes, dbs)

# 3. Generate the Bar Plot
# This will now show the pathways associated with your best genes

plotEnrich(enriched_balanced[["GO_Biological_Process_2023"]], 
           showTerms = 15, 
           numChar = 50, 
           y = "Count", 
           orderBy = "P.value") +
  labs(title = "Core Biological Drivers of Colorectal Cancer") +
  theme_minimal()

# Save the top 20 genes list to your working directory
write.csv(top_20, "results_summary.csv", row.names = FALSE)
