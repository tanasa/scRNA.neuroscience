library(SoupX)
library(Matrix)
library(DropletUtils)
library(Seurat)

# ============================================================================
# Configuration
# ============================================================================
# Set working directory
setwd("/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/fastq_file_Dyslexia_r2_cellRanger/AD_AG_2658")

# Define paths to H5 files
# raw_h5      <- "AD_AG_2658_raw_feature_bc_matrix.h5"
# filtered_h5 <- "AD_AG_2658_sample_filtered_feature_bc_matrix.h5"
# sample_name <- "AD_AG_2658"

# Define paths to H5 files
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
