library(SoupX)
library(Matrix)
library(Seurat) 
# library(hdf5r)
library(DropletUtils)

# Define paths to H5 files
setwd("/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/fastq_file_Dyslexia_r2_cellRanger/AD_AG_2658")
list.files()

# ============================================================================
# Configuration
# ============================================================================
raw_h5      <- "AD_AG_2658_raw_feature_bc_matrix.h5"
filtered_h5 <- "AD_AG_2658_sample_filtered_feature_bc_matrix.h5"
sample_name <- "AD_AG_2658"

# Create output directories
plot_dir <- "plots_soupx"
out_dir  <- paste0(sample_name, "_soupx_corrected")
dir.create(plot_dir, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)

# ============================================================================
# Load 10X Data from H5 files
# ============================================================================
cat("Loading 10X data from H5 files...\n")
tod_sce <- read10xCounts(raw_h5)           # Table of droplets (all barcodes)
toc_sce <- read10xCounts(filtered_h5)      # Table of counts (cell-containing barcodes)

cat("Raw matrix dimensions:", nrow(counts(tod_sce)), "genes x", ncol(counts(tod_sce)), "droplets\n")
cat("Filtered matrix dimensions:", nrow(counts(toc_sce)), "genes x", ncol(counts(toc_sce)), "cells\n")

# ============================================================================
# Check if rownames are ENSEMBL IDs
# ============================================================================
cat("\n=== Checking TOD (raw) rowname format ===\n")
cat("First 10 rownames:\n")
print(head(rownames(tod_sce), 10))

# Check for ENSEMBL pattern in TOD
is_ensembl_tod <- grepl("^ENS", rownames(tod_sce))
cat("\nENSEMBL IDs detected:", sum(is_ensembl_tod), "out of", length(is_ensembl_tod), "genes\n")
cat("Percentage ENSEMBL:", round(100 * sum(is_ensembl_tod) / length(is_ensembl_tod), 1), "%\n")

# Show examples of non-ENSEMBL if any exist
if(sum(!is_ensembl_tod) > 0) {
  cat("\nExamples of non-ENSEMBL rownames:\n")
  print(head(rownames(tod_sce)[!is_ensembl_tod], 10))
}

cat("\n=== Checking TOC (filtered) rowname format ===\n")
cat("First 10 rownames:\n")
print(head(rownames(toc_sce), 10))

# Check for ENSEMBL pattern in TOC
is_ensembl_toc <- grepl("^ENS", rownames(toc_sce))
cat("\nENSEMBL IDs detected:", sum(is_ensembl_toc), "out of", length(is_ensembl_toc), "genes\n")
cat("Percentage ENSEMBL:", round(100 * sum(is_ensembl_toc) / length(is_ensembl_toc), 1), "%\n")

# Show examples of non-ENSEMBL if any exist
if(sum(!is_ensembl_toc) > 0) {
  cat("\nExamples of non-ENSEMBL rownames:\n")
  print(head(rownames(toc_sce)[!is_ensembl_toc], 10))
}

# Check if both have same rownames
cat("\n=== Comparing TOD and TOC rownames ===\n")
cat("Same rownames in both matrices:", identical(rownames(tod_sce), rownames(toc_sce)), "\n")

# ============================================================================
# Convert DelayedMatrix to sparse dgCMatrix
# ============================================================================
cat("Converting DelayedMatrix to sparse matrices (this may take a moment)...\n")

# Simple conversion using as() - this works with DelayedMatrix
tod <- as(counts(tod_sce), "dgCMatrix")
toc <- as(counts(toc_sce), "dgCMatrix")

# Add row names (genes)
rownames(tod) <- rownames(tod_sce)
rownames(toc) <- rownames(toc_sce)

# Add column names (barcodes)
if ("Barcode" %in% colnames(colData(tod_sce))) {
  colnames(tod) <- colData(tod_sce)$Barcode
  colnames(toc) <- colData(toc_sce)$Barcode
} else {
  colnames(tod) <- paste0("DROPLET-", seq_len(ncol(tod)))
  colnames(toc) <- paste0("CELL-", seq_len(ncol(toc)))
}

cat("Conversion complete!\n")

# ============================================================================
# Create SoupChannel
# ============================================================================
cat("Creating SoupChannel object...\n")
sc <- SoupChannel(tod = tod, toc = toc, calcSoupProfile = TRUE)

# ============================================================================
# Quick clustering for autoEstCont
# ============================================================================
cat("\nPerforming quick clustering for contamination estimation...\n")

# Create a quick Seurat object for clustering
cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(counts = toc, project = sample_name)

# Standard preprocessing
cat("Normalizing data...\n")
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

cat("Finding variable features...\n")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", 
                                    nfeatures = 2000, verbose = FALSE)

cat("Scaling data...\n")
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

cat("Running PCA...\n")
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)

cat("Finding neighbors...\n")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20, verbose = FALSE)

cat("Finding clusters...\n")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)

# Extract clusters
clusters <- seurat_obj$seurat_clusters
names(clusters) <- colnames(seurat_obj)

cat("Clustering complete! Found", length(unique(clusters)), "clusters\n")
print(table(clusters))

# Add clusters to SoupChannel
sc <- setClusters(sc, clusters)

sc

str(sc)

sc$toc[1: 5, 1:5]
rownames(sc$toc[1: 5, 1:5])
dim(sc$toc)

# Check for mixed content (ENSEMBL vs gene symbols)
has_ensembl <- grepl("^ENS", rownames(sc$toc))
table(has_ensembl)  # Shows count of ENSEMBL vs non-ENSEMBL

head(sc$metaData, 2)

head(sc$soupProfile, 2)

# ============================================================================
# Estimate Contamination with autoEstCont
# ============================================================================
cat("\nAuto-estimating contamination fraction...\n")
png(file.path(plot_dir, paste0(sample_name, "_SoupX_autoEstCont.png")),
    width = 7, height = 6, units = "in", res = 300)
sc <- autoEstCont(sc, doPlot = TRUE)
dev.off()

# Extract rho estimate
rho_estimate <- sc$fit$rhoEst

# Check if estimation succeeded
if(is.null(rho_estimate)) {
  stop("autoEstCont failed - you may need to manually estimate contamination")
}

cat("Estimated contamination fraction:", round(rho_estimate, 4), "\n")

# Also print the confidence interval
cat("95% FWHM interval: [", round(sc$fit$rhoFWHM[1], 4), ", ", 
    round(sc$fit$rhoFWHM[2], 4), "]\n", sep="")

writeLines(paste("Contamination fraction (rho):", rho_estimate,
                 "\n95% FWHM interval: [", sc$fit$rhoFWHM[1], ", ", sc$fit$rhoFWHM[2], "]"),
           file.path(plot_dir, paste0(sample_name, "_contamination_estimate.txt")))

str(sc)

str(sc$fit)

head(sc$fit$dd, 2)

sc$fit$priorRho

head(sc$fit$posterior) 

sc$fit$rhoEst

sc$fit$rhoFWHM 

head(sc$fit$markersUsed, 2)

head(sc$soupProfile, 2)
dim(sc$soupProfile)

# ============================================================================
# Map ALL genes from ENSEMBL to Gene Symbols
# ============================================================================
cat("Mapping ALL genes from ENSEMBL to gene symbols...\n")

library(org.Hs.eg.db)

# Get all genes from the soup profile (it's already a dataframe)
all_genes <- sc$soupProfile
all_genes$ensembl <- rownames(all_genes)

cat("Total genes to map:", nrow(all_genes), "\n")

# Map ENSEMBL to SYMBOL for ALL genes
symbols <- mapIds(org.Hs.eg.db, 
                  keys = all_genes$ensembl,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")

# Add to dataframe
all_genes$gene_symbol <- symbols

# Use symbol if available, otherwise keep ENSEMBL ID
all_genes$gene_name <- ifelse(is.na(all_genes$gene_symbol), 
                               all_genes$ensembl, 
                               all_genes$gene_symbol)

# Check mapping success rate
cat("Successfully mapped:", sum(!is.na(all_genes$gene_symbol)), "genes\n")
cat("Unmapped (kept as ENSEMBL):", sum(is.na(all_genes$gene_symbol)), "genes\n")

head(all_genes)
# Replace rownames with gene_name
# rownames(all_genes) <- all_genes$gene_name
# head(all_genes)

# how many unique symbols and duplicates :

# Check for duplicates in gene_name
cat("\nChecking for duplicate gene names...\n")

# Count duplicates
duplicate_counts <- table(all_genes$gene_name)
duplicates <- duplicate_counts[duplicate_counts > 1]

cat("Total unique gene names:", length(unique(all_genes$gene_name)), "\n")
cat("Total genes:", nrow(all_genes), "\n")
cat("Number of duplicated gene names:", length(duplicates), "\n")
cat("Total duplicate entries:", sum(duplicate_counts[duplicate_counts > 1]), "\n")

# Show examples of duplicated genes
cat("\nTop duplicated gene names:\n")
print(head(sort(duplicates, decreasing = TRUE), 20))

# Show which ENSEMBL IDs map to the same gene (example)
if(length(duplicates) > 0) {
  cat("\nExample: Multiple ENSEMBL IDs for", names(duplicates)[1], ":\n")
  example_gene <- names(duplicates)[1]
  print(all_genes[all_genes$gene_name == example_gene, c("ensembl", "gene_symbol", "est", "counts")])
}

cat("")
cat("printing duplicated genes")
print(duplicates)

# ============================================================================
# Solution: Make rownames unique by appending ENSEMBL ID to duplicates
#           Or keep ENSEMBL for duplicates (alternative approach)
# ============================================================================
cat("\nCreating unique rownames...\n")

# Create unique gene names by adding suffix for duplicates
# all_genes$unique_gene_name <- make.unique(all_genes$gene_name, sep = "_")

# Or keep ENSEMBL for duplicates (alternative approach)
all_genes$unique_gene_name <- ifelse(duplicated(all_genes$gene_name) | duplicated(all_genes$gene_name, fromLast = TRUE),
                                      all_genes$ensembl,
                                      all_genes$gene_name)
# Set unique rownames
rownames(all_genes) <- all_genes$unique_gene_name

cat("Rownames are now unique!\n")

# Show preview
cat("\nPreview with unique names:\n")
print(head(all_genes[, c("ensembl", "gene_symbol", "gene_name", "unique_gene_name", "est", "counts")], 10))

# ============================================================================
# Verify duplicate names after making them unique
# ============================================================================
cat("\n=== Verifying rownames after making them unique ===\n")

# Check if rownames are now unique
cat("Total rownames:", length(rownames(all_genes)), "\n")
cat("Unique rownames:", length(unique(rownames(all_genes))), "\n")
cat("Are all rownames unique?", length(rownames(all_genes)) == length(unique(rownames(all_genes))), "\n")

# Check for any remaining duplicates
remaining_duplicates <- sum(duplicated(rownames(all_genes)))
cat("Number of remaining duplicate rownames:", remaining_duplicates, "\n")

if(remaining_duplicates > 0) {
  cat("\nWARNING: Still have duplicates!\n")
  dup_names <- rownames(all_genes)[duplicated(rownames(all_genes))]
  cat("Duplicated names:", paste(head(dup_names, 10), collapse = ", "), "\n")
} else {
  cat("\n✓ SUCCESS: All rownames are now unique!\n")
}

# Show examples of how duplicates were handled
cat("\n=== Examples of duplicate handling ===\n")
# Find genes that got suffixes added
has_suffix <- grepl("_[0-9]+$", rownames(all_genes))
if(sum(has_suffix) > 0) {
  cat("Examples of genes with suffixes added:\n")
  print(head(all_genes[has_suffix, c("ensembl", "gene_symbol", "gene_name", "unique_gene_name")], 10))
}

# Summary
cat("\n=== Summary ===\n")
cat("Original gene names (before unique):", length(unique(all_genes$gene_name)), "\n")
cat("Unique rownames (after make.unique):", length(unique(rownames(all_genes))), "\n")
cat("Genes with suffixes added:", sum(has_suffix), "\n")

# Show a preview
cat("\nPreview of mapped genes:\n")
print(head(all_genes[, c("ensembl", "gene_symbol", "gene_name", "est", "counts")], 10))

# Save ALL genes with mappings (sorted by contamination level)
all_genes_sorted <- all_genes[order(-all_genes$est), ]
write.csv(all_genes_sorted,
          file.path(plot_dir, paste0(sample_name, "_all_genes_mapped.csv")),
          row.names = FALSE)

# Also save TOP 100 soup genes separately
top_soup_genes <- head(all_genes_sorted, 100)
write.csv(top_soup_genes,
          file.path(plot_dir, paste0(sample_name, "_top_soup_genes.csv")),
          row.names = FALSE)

cat("\nSaved all genes to:", file.path(plot_dir, paste0(sample_name, "_all_genes_mapped.csv")), "\n")
cat("Saved top 100 soup genes to:", file.path(plot_dir, paste0(sample_name, "_top_soup_genes.csv")), "\n")



# ============================================================================
# Additional Diagnostic Plots
# ============================================================================
cat("Generating diagnostic plots...\n")

# 1. Contamination fraction distribution across cells
if(!is.null(sc$metaData$rho)) {
  png(file.path(plot_dir, paste0(sample_name, "_contamination_dist.png")),
      width = 6, height = 5, units = "in", res = 300)
  hist(sc$metaData$rho, breaks = 50,
       main = paste(sample_name, "- Contamination per Cell"),
       xlab = "Contamination Fraction (rho)",
       col = "skyblue", border = "white")
  abline(v = rho_estimate, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = paste("Mean =", round(rho_estimate, 4)), 
         col = "red", lty = 2, lwd = 2)
  dev.off()
}

# 2. Top contaminating genes (ambient RNA profile)
cat("Identifying top soup genes...\n")
if(is.data.frame(sc$soupProfile)) {
  top_soup_genes <- head(sc$soupProfile[order(-sc$soupProfile$est), ], 100)
} else {
  # If soupProfile is a vector
  top_soup_genes <- data.frame(
    gene = names(head(sort(sc$soupProfile, decreasing = TRUE), 100)),
    fraction = head(sort(sc$soupProfile, decreasing = TRUE), 100)
  )
}

head(top_soup_genes)

top_soup_genes$ensembl = rownames(top_soup_genes)

library(org.Hs.eg.db)

# Map ENSEMBL to SYMBOL
symbols <- mapIds(org.Hs.eg.db, 
                  keys = top_soup_genes$ensembl,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")

# Add to dataframe
top_soup_genes$gene_symbol <- symbols

# Use symbol if available, otherwise keep ENSEMBL ID
top_soup_genes$gene_name <- ifelse(is.na(top_soup_genes$gene_symbol), 
                                    top_soup_genes$ensembl, 
                                    top_soup_genes$gene_symbol)

head(top_soup_genes)

write.csv(top_soup_genes,
          file.path(plot_dir, paste0(sample_name, "_top_soup_genes.csv")),
          row.names = TRUE)

# ============================================================================
# Correct Counts
# ============================================================================
cat("Adjusting counts for ambient RNA contamination...\n")
adj <- adjustCounts(sc, roundToInt = FALSE)

# ============================================================================
# Save Corrected Data
# ============================================================================
cat("Saving corrected count matrix...\n")

# Option 1: Save as RDS (recommended - preserves sparse matrix format)
saveRDS(adj, file.path(out_dir, paste0(sample_name, "_soupx_corrected.rds")))

# Option 2: Save in 10X format (for compatibility with other tools)
writeMM(adj, file.path(out_dir, "matrix.mtx"))
write.table(rownames(adj), file.path(out_dir, "features.tsv"),
            quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(colnames(adj), file.path(out_dir, "barcodes.tsv"),
            quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

# Compress the matrix file
system(paste0("gzip -f ", file.path(out_dir, "matrix.mtx")))

# Option 3: Save back as H5 using DropletUtils
cat("Saving as H5 format...\n")
write10xCounts(path = file.path(out_dir, paste0(sample_name, "_soupx_corrected.h5")),
               x = adj, 
               type = "HDF5", 
               overwrite = TRUE,
               gene.id = rownames(adj),
               gene.symbol = rownames(adj))

# Save the SoupChannel object for later inspection
saveRDS(sc, file.path(out_dir, paste0(sample_name, "_soupChannel.rds")))

# ============================================================================
# Save Session Info
# ============================================================================
writeLines(capture.output(sessionInfo()),
           file.path(plot_dir, paste0(sample_name, "_sessionInfo.txt")))

cat("\n=== SoupX correction complete ===\n")
cat("Corrected matrix saved to:", out_dir, "\n")
cat("Diagnostic plots saved to:", plot_dir, "\n")
cat("\nOutput files:\n")
cat("  - RDS format:", file.path(out_dir, paste0(sample_name, "_soupx_corrected.rds")), "\n")
cat("  - H5 format:", file.path(out_dir, paste0(sample_name, "_soupx_corrected.h5")), "\n")
cat("  - 10X MTX format:", file.path(out_dir, "matrix.mtx.gz"), "\n")
cat("\nTo load corrected matrix:\n")
cat("  adj <- readRDS('", file.path(out_dir, paste0(sample_name, "_soupx_corrected.rds")), "')\n", sep = "")

# The runs were on 10X cloud :

# What we do not do in this script :
# we do not : Update rownames of corrected matrix with gene symbols
# we do not : Fix features.tsv to include both ENSEMBL IDs and gene symbols

# We will change the ENSEMBL genes into gene symbols on the adjusted counts matrix




