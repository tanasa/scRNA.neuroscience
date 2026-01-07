##########################################################################################################
# RNA VELOCITY ANALYSIS - FLEXIBLE VERSION
# Adapted for Kallisto/La Manno or any velocity matrices
# 
# USAGE: Update the velocyto_dir path below to point to YOUR velocity_matrices_output folder
##########################################################################################################

library(Seurat)
library(SeuratWrappers)
library(velocyto.R)
library(Matrix)
library(ggplot2)
library(cowplot)
library(dplyr)

# ============================================================================
# EDIT THIS PATH - Point to YOUR velocity_matrices_output folder
# ============================================================================

# OPTION 1: If files are in velocity_matrices_output subfolder
velocyto_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/JK4488_01_combined/pooled/JK4488_01_kallisto_lamanno/velocity_matrices_output"

# OPTION 2: If files are directly in JK4488_01_kallisto_lamanno folder
# velocyto_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/JK4488_01_combined/pooled/JK4488_01_kallisto_lamanno"

# Output directory (will be created in same parent folder)
output_dir <- file.path(dirname(velocyto_dir), "velocity_analysis_output")

# ============================================================================
# ============================================================================

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║         RNA VELOCITY ANALYSIS - STARTING                  ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n\n")

cat("=== Directory Setup ===\n")
cat("Input directory:", velocyto_dir, "\n")
cat("Output directory:", output_dir, "\n\n")

# Check if directory exists
if (!dir.exists(velocyto_dir)) {
  cat("❌ ERROR: Input directory does not exist!\n")
  cat("   Path:", velocyto_dir, "\n\n")
  cat("Please update the 'velocyto_dir' variable at the top of this script\n")
  cat("to point to the folder containing your velocity matrices.\n\n")
  stop("Directory not found")
}

setwd(velocyto_dir)

##########################################################################################################
# STEP 1: VERIFY AND LOAD MATRICES
##########################################################################################################

cat("=== Checking for Required Files ===\n")

# Check if files exist
required_files <- c("spliced.mtx", "unspliced.mtx", "features.tsv", "barcodes.tsv")

cat("Looking for files in:", getwd(), "\n\n")

file_status <- sapply(required_files, file.exists)
for (i in seq_along(required_files)) {
  status <- ifelse(file_status[i], "✅", "❌")
  cat(status, required_files[i], "\n")
}

missing_files <- required_files[!file_status]

if (length(missing_files) > 0) {
  cat("\n❌ ERROR: Missing required files:", paste(missing_files, collapse=", "), "\n")
  cat("\nFiles in current directory:\n")
  print(list.files())
  stop("Required files not found!")
}

cat("\n✅ All required files present!\n\n")

# Load matrices
cat("=== Loading Matrices ===\n")

spliced <- readMM("spliced.mtx")
unspliced <- readMM("unspliced.mtx")

cat("✓ Spliced matrix loaded:", dim(spliced)[1], "genes x", dim(spliced)[2], "cells\n")
cat("✓ Unspliced matrix loaded:", dim(unspliced)[1], "genes x", dim(unspliced)[2], "cells\n")

# Check if dimensions match
if (!all(dim(spliced) == dim(unspliced))) {
  stop("ERROR: Spliced and unspliced matrices have different dimensions!")
}

# Load barcodes and features
barcodes <- read.table("barcodes.tsv", header = FALSE, stringsAsFactors = FALSE)
features <- read.table("features.tsv", header = FALSE, stringsAsFactors = FALSE)

cat("✓ Barcodes loaded:", nrow(barcodes), "\n")
cat("✓ Features loaded:", nrow(features), "rows x", ncol(features), "columns\n\n")

# Verify dimensions match
if (nrow(features) != nrow(spliced)) {
  cat("⚠️  WARNING: Number of genes in features.tsv (", nrow(features), 
      ") doesn't match matrix (", nrow(spliced), ")\n")
}
if (nrow(barcodes) != ncol(spliced)) {
  cat("⚠️  WARNING: Number of barcodes (", nrow(barcodes), 
      ") doesn't match matrix (", ncol(spliced), ")\n")
}

# ⭐ Use gene symbols (V2) instead of Ensembl IDs (V1)
cat("=== Assigning Gene Names ===\n")

if (ncol(features) >= 2) {
  gene_names <- make.unique(features$V2)
  cat("✓ Using gene symbols from column 2 (V2)\n")
  cat("  Example genes:", paste(head(gene_names, 3), collapse=", "), "\n")
} else if (ncol(features) == 1) {
  gene_names <- make.unique(features$V1)
  cat("⚠️  Only 1 column in features.tsv, using V1\n")
  cat("  Example genes:", paste(head(gene_names, 3), collapse=", "), "\n")
} else {
  stop("ERROR: features.tsv has unexpected format!")
}

# Assign rownames and colnames
rownames(spliced) <- gene_names
colnames(spliced) <- barcodes$V1

rownames(unspliced) <- gene_names
colnames(unspliced) <- barcodes$V1

cat("✓ Dimension names assigned\n\n")

##########################################################################################################
# STEP 2: QC & FILTERING
##########################################################################################################

cat("=== Quality Control ===\n")

nGene <- Matrix::colSums(spliced > 0)
nUMI <- Matrix::colSums(spliced)

cat("Basic statistics:\n")
cat("  Median genes/cell:", median(nGene), "\n")
cat("  Median UMI/cell:", median(nUMI), "\n")
cat("  Mean genes/cell:", round(mean(nGene), 1), "\n")
cat("  Mean UMI/cell:", round(mean(nUMI), 1), "\n\n")

# Calculate MT%
mt_genes <- grep("^MT-", rownames(spliced), value = TRUE)
if (length(mt_genes) > 0) {
  percent_mito <- Matrix::colSums(spliced[mt_genes, ]) / nUMI * 100
  cat("✓ Found", length(mt_genes), "mitochondrial genes\n")
  cat("  MT genes:", paste(mt_genes, collapse=", "), "\n")
  cat("  Median MT%:", round(median(percent_mito), 2), "%\n")
} else {
  percent_mito <- rep(0, ncol(spliced))
  cat("⚠️  No MT genes found (genes don't start with 'MT-')\n")
  cat("  This is OK if your gene names don't follow this convention\n")
}

# QC filtering
cat("\nQC Filters:\n")
cat("  • Genes/cell: 200 - 6000\n")
cat("  • UMI/cell: ≥ 500\n")
cat("  • MT%: < 10%\n\n")

cells_pass <- which(nGene >= 200 & nGene <= 6000 & nUMI >= 500 & percent_mito < 10)

cat("Cells before QC:", ncol(spliced), "\n")
cat("Cells after QC:", length(cells_pass), 
    sprintf(" (%.1f%% retained)\n\n", 100 * length(cells_pass) / ncol(spliced)))

spliced <- spliced[, cells_pass]
unspliced <- unspliced[, cells_pass]

##########################################################################################################
# STEP 3: CREATE SEURAT OBJECT
##########################################################################################################

cat("=== Creating Seurat Object ===\n")

seurat_obj <- CreateSeuratObject(
  counts = spliced,
  project = "Velocity_Analysis",
  min.cells = 3,
  min.features = 200
)

seurat_obj$percent_mito <- percent_mito[cells_pass]
cat("✓ Seurat object created:", ncol(seurat_obj), "cells x", nrow(seurat_obj), "genes\n\n")

##########################################################################################################
# STEP 4: SEURAT PROCESSING
##########################################################################################################

cat("=== Seurat Workflow ===\n")

seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
cat("✓ Normalized\n")

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
cat("✓ Variable features (", length(VariableFeatures(seurat_obj)), " genes)\n", sep="")

seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
cat("✓ Scaled\n")

seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)
cat("✓ PCA\n")

seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
cat("✓ UMAP\n")

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)
cat("✓ Clustering (", length(unique(seurat_obj$seurat_clusters)), " clusters)\n\n", sep="")

# Save UMAP
p_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP Clustering") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "UMAP_clusters.png"), p_umap, width = 8, height = 6, dpi = 300)
cat("✓ UMAP plot saved\n\n")

##########################################################################################################
# STEP 5: PREPARE VELOCITY MATRICES
##########################################################################################################

cat("=== Preparing Velocity Matrices ===\n")

filtered_cells <- colnames(seurat_obj)
spliced_filtered <- spliced[, filtered_cells]
unspliced_filtered <- unspliced[, filtered_cells]

cat("✓ Subset to", ncol(spliced_filtered), "Seurat-filtered cells\n")

# Filter genes
spliced_gene_counts <- Matrix::rowSums(spliced_filtered > 0)
unspliced_gene_counts <- Matrix::rowSums(unspliced_filtered > 0)
genes_keep <- names(which(spliced_gene_counts >= 10 & unspliced_gene_counts >= 10))

cat("✓ Genes before filtering:", nrow(spliced_filtered), "\n")
cat("✓ Genes after filtering:", length(genes_keep), "(≥10 cells in both matrices)\n")

spliced_filtered <- spliced_filtered[genes_keep, ]
unspliced_filtered <- unspliced_filtered[genes_keep, ]

cat("✓ Final dimensions:", dim(spliced_filtered)[1], "genes x", 
    dim(spliced_filtered)[2], "cells\n\n")

##########################################################################################################
# STEP 6: ESTIMATE VELOCITY
##########################################################################################################

cat("=== Estimating RNA Velocity ===\n")
cat("Parameters: kCells=50, fit.quantile=0.05\n")
cat("This may take several minutes...\n\n")

set.seed(42)

# ⭐ APPLIED
rvel <- gene.relative.velocity.estimates(
  emat = spliced_filtered,      # ✅ Sparse dgCMatrix
  nmat = unspliced_filtered,    # ✅ Sparse dgCMatrix
  deltaT = 1,
  kCells = 50,
  fit.quantile = 0.05,
  verbose = TRUE                # ✅ Comma present
)

cat("\n✓ Velocity estimation complete!\n")
cat("  Genes in velocity:", nrow(rvel$current), "\n")
cat("  Cells in velocity:", ncol(rvel$current), "\n\n")

# glimpse(rvel, max.level=1)
# List of 11
# $ cellKNN       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
# $ conv.nmat.norm:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
# $ conv.emat.norm:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
# $ gamma         : Named num [1:8315] 0.0774 0.334 0.5239 0.1539 0.0793 ...
#  ..- attr(*, "names")= chr [1:8315] "ENSG00000241860" "ENSG00000292994" "ENSG00000237491" "ENSG00000228794" ...
# $ projected     :Formal class 'dgeMatrix' [package "Matrix"] with 4 slots
# $ current       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
# $ deltaE        :Formal class 'dgeMatrix' [package "Matrix"] with 4 slots
# $ deltaT        : num 1
# $ ko            :'data.frame': 14791 obs. of  4 variables:
# $ mult          : num 1000
# $ kCells        : num 50

##########################################################################################################
# STEP 7: PROJECT ON UMAP
##########################################################################################################

cat("=== Projecting Velocity on UMAP ===\n")

umap_coords <- Embeddings(seurat_obj, "umap")
common_cells <- intersect(rownames(umap_coords), colnames(rvel$current))
umap_coords <- umap_coords[common_cells, ]

cat("✓ Matched", length(common_cells), "cells between UMAP and velocity\n")

# head(umap_coords)
#                    umap_1    umap_2
# GGTGATTAGGTACATA -2.338967 -1.803730
# CATGCCTTCGTTCTGC  5.222305 -4.557078
# GAACACTAGAGATGCC  4.288650  4.042195
# TTCATGTCAATAACCC -9.680163  1.380516
# CCCGAAGGTGAGAACC -3.373738  1.354640

# Create velocity plot
cat("Creating velocity visualization...\n")

png(file.path(output_dir, "velocity_UMAP_arrows.png"),
    width = 10, height = 10, units = "in", res = 300)

show.velocity.on.embedding.cor(
  emb = umap_coords,
  vel = rvel,
  n = 200,
  scale = 'sqrt',
  cell.colors = ac(x = as.numeric(seurat_obj$seurat_clusters[common_cells]), alpha = 0.5),
  cex = 0.8,
  arrow.scale = 3,
  show.grid.flow = TRUE,
  min.grid.cell.mass = 0.5,
  grid.n = 40,
  arrow.lwd = 1,
  do.par = TRUE,
  cell.border.alpha = 0.1
)

dev.off()
cat("✓ Velocity plot saved\n\n")

##########################################################################################################
# STEP 8: SAVE RESULTS
##########################################################################################################

cat("=== Saving Results ===\n")

# saveRDS(seurat_obj, file.path(output_dir, "seurat_object.rds"))
# cat("✓ seurat_object.rds\n")

saveRDS(rvel, file.path(output_dir, "velocity_results.rds"))
cat("✓ velocity_results.rds\n")

writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
cat("✓ session_info.txt\n\n")

##########################################################################################################
# SUMMARY
##########################################################################################################

cat("\n")
cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║                                                           ║\n")
cat("║        ✓✓✓  ANALYSIS COMPLETE  ✓✓✓                       ║\n")
cat("║                                                           ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("=== Summary ===\n")
cat("Cells analyzed:", ncol(seurat_obj), "\n")
cat("Genes:", nrow(seurat_obj), "\n")
cat("Clusters:", length(unique(seurat_obj$seurat_clusters)), "\n")
cat("Velocity genes:", nrow(spliced_filtered), "\n\n")

cat("=== Output Files ===\n")
cat("Location:", output_dir, "\n")
cat("  • seurat_object.rds\n")
cat("  • velocity_results.rds\n")
cat("  • UMAP_clusters.png\n")
cat("  • velocity_UMAP_arrows.png\n")
cat("  • session_info.txt\n\n")

cat("✅ All done! Check your output directory for results.\n\n")
