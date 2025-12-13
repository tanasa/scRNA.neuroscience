library(Seurat)
library(Matrix)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(data.table)
library(vioplot)
library(harmony)
library(cowplot)
library(patchwork)
library(future)
library(foreach)
library(doParallel)
library(fs)
library(MAST)
library(DoubletFinder)
library("destiny")

# options(future.globals.maxSize = 20 * 1024^3) 
# future::plan(multicore, workers = parallel::detectCores() - 1)

# Replace the future setup with:
# options(future.globals.maxSize = 20 * 1024^3)
# future::plan(sequential)  # Use sequential instead of parallel for stability
# (options(future.seed = TRUE))

packageVersion("Seurat")
set.seed(1337)

# --- Jupyter inline rendering presets ---
# Sharper inline figures
options(repr.plot.width = 8,   # inches
        repr.plot.height = 6,
        repr.plot.res = 160)   # dpi

# Crisper text/lines on some servers
options(bitmapType = "cairo")

# Clean ggplot look & readable text
library(ggplot2)
theme_set(theme_classic(base_size = 14))

# Optional: default point/text sizes in ggplot
update_geom_defaults("point", list(size = 1.2, alpha = 0.8))
update_geom_defaults("text",  list(size = 4))

getwd()

# Path to RDS file

mtgsea <- "/mnt/nfs/CX000008_DS1/projects/btanasa/brain_refs_MTG_SEA/brain_refs_MTG_SEA.from.matrix.metadata.sce.rds"
mtg_sea <- readRDS(mtgsea)
# Inspect its structure
class(mtg_sea)
# str(mtg_sea)

# .$ cluster_label              : chr [1:76533] "Inh L1-2 SST CCNJL" "Exc L5-6 FEZF2 IFNG-AS1" "Exc L3-5 RORB LINC01202" "Exc L2 LINC00507 GLRA3" ...
# .$ class_label                : chr [1:76533] "GABAergic" "Glutamatergic" "Glutamatergic" "Glutamatergic" ...
# .$ subclass_label             : chr [1:76533] "Sst" "L5/6 NP" "L5 IT" "L2/3 IT" ...

# Check available columns in colData
colnames(colData(mtg_sea))

# Get unique values for each label
unique_clusters    <- unique(mtg_sea$cluster_label)
unique_classes     <- unique(mtg_sea$class_label)
unique_subclasses  <- unique(mtg_sea$subclass_label)

# Print counts + first few
cat("cluster_label:", length(unique_clusters), "unique\n")
print(unique_clusters)

cat("subclass_label:", length(unique_subclasses), "unique\n")
print(unique_subclasses)

cat("class_label:", length(unique_classes), "unique\n")
print(unique_classes)

# there are no ENSEMBL IDs

# Check what dimensionality reductions are available
reducedDimNames(mtg_sea)

# If UMAP exists, extract it
if ("UMAP" %in% reducedDimNames(mtg_sea)) {
    umap_coords <- reducedDim(mtg_sea, "UMAP")
    print("UMAP found!")
    print(dim(umap_coords))
} else {
    print("No UMAP found. Available reductions:")
    print(reducedDimNames(mtg_sea))
}

library(scater)
library(scran)

# Compute PCA first (if not present)
if (!"PCA" %in% reducedDimNames(mtg_sea)) {
    mtg_sea <- runPCA(mtg_sea, ncomponents = 50)
}

# Compute UMAP
mtg_sea <- runUMAP(mtg_sea, dimred = "PCA")

# Create the plot and assign to variable
p <- plotReducedDim(mtg_sea, dimred = "UMAP", colour_by = "class_label")

# Save with ggsave
ggsave("mtg_sea_UMAP_by_class.png", plot = p, width = 10, height = 8, dpi = 300)

options(repr.plot.width = 14, repr.plot.height = 10)
p



library(ggplot2)
library(dplyr)

# Get UMAP coordinates and labels
umap_coords <- reducedDim(mtg_sea, "UMAP")
umap_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  subclass = mtg_sea$subclass_label
)

# Calculate centroids for each subclass
centroids <- umap_df %>%
  group_by(subclass) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2)
  )

# Create plot from scratch
p_subclass <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = subclass)) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_text(data = centroids, 
            aes(x = UMAP1, y = UMAP2, label = subclass),
            color = "black",
            size = 6,
            fontface = "bold",
            check_overlap = TRUE) +
  theme_classic(base_size = 12) +
  labs(title = "MTG-SEA UMAP by Subclass",
       color = "Subclass")

ggsave("mtg_sea_UMAP_by_subclass_labeled.png", plot = p_subclass, 
       width = 14, height = 10, dpi = 300)

p_subclass 

library(ggrepel)

# p_subclass <- plotReducedDim(mtg_sea, dimred = "UMAP", colour_by = "subclass_label") +
#  geom_text_repel(data = centroids, 
#                  aes(x = UMAP1, y = UMAP2, label = subclass),
#                  color = "black",
#                  size = 3,
#                  fontface = "bold",
#                  bg.color = "white",
#                  bg.r = 0.1,
#                  max.overlaps = 20)

# ggsave("mtg_sea_UMAP_by_subclass_labeled.png", plot = p_subclass, 
#         width = 14, height = 10, dpi = 300)

# p_subclass library(scater)
library(scran)
library(ggplot2)
library(dplyr)

# Get UMAP coordinates and CLASS labels
umap_coords <- reducedDim(mtg_sea, "UMAP")
umap_df_class <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  class = mtg_sea$class_label
)

# Calculate centroids for CLASS labels (not subclass)
centroids_class <- umap_df_class %>%
  group_by(class) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2)
  )

# Create the plot with labels overlaid
p <- plotReducedDim(mtg_sea, dimred = "UMAP", colour_by = "class_label") +
  geom_text(data = centroids_class, 
            aes(x = UMAP1, y = UMAP2, label = class),
            color = "black",
            size = 8,
            fontface = "bold",
            check_overlap = TRUE, 
            inherit.aes = FALSE)  # ADD THIS!

# Save with ggsave
ggsave("mtg_sea_UMAP_by_class_labeled.png", plot = p, width = 10, height = 8, dpi = 300)

options(repr.plot.width = 14, repr.plot.height = 10)
p



# read the object to annotate :

# sample_id <- "JK4488_01"

# base_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan"
# path <- file.path(
#  base_dir,
#  paste0(sample_id, "_combined"),
#  "pooled",
#  paste0(sample_id, "_cellranger_unfilteredMatrix")
# )

# Show the target path
# cat("Target path:\n", path, "\n")

# Ensure the parent path exists before setwd (fail fast if not)
# if (!dir.exists(path)) {
#  stop("Path does not exist: ", path)
# }


# Move into the directory and list files
# setwd(path)
# cat("Working directory set to:\n", getwd(), "\n\n")
# print(list.files())

# Create <sample_id>.soupX inside this path
# singleR_dir <- file.path(path, paste0(sample_id, ".singleR"))

# if (!dir.exists(singleR_dir)) {
#  dir.create(singleR_dir, showWarnings = FALSE)
#  cat("\nCreated folder:\n", singleR_dir, "\n")
# } else {
#  cat("\nFolder already exists:\n", singleR_dir, "\n")
# }

# output_dir = singleR_dir
# outfn <- function(name) file.path(output_dir, paste0(sample_id, "_", name))

# === Path and sample setup ===

sample_id <- "AD_AG_2658"

# Base directory (up to zanalysis_bogdan)
base_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan"

soupX_dir <- file.path(base_dir, "fastq_file_Dyslexia_r2_cellRanger", sample_id, 
                       paste0(sample_id, ".soupX.gs.analysis"))

cat("Target path:\n", soupX_dir, "\n")
if (!dir.exists(soupX_dir)) stop("Path does not exist: ", soupX_dir)

setwd(soupX_dir)
cat("Working directory set to:\n", getwd(), "\n\n")

# === Create SingleR output folder ===
singleR_dir <- file.path(soupX_dir, paste0(sample_id, ".singleR"))
if (!dir.exists(singleR_dir)) {
  dir.create(singleR_dir, showWarnings = FALSE)
  cat("Created folder:\n", singleR_dir, "\n")
} else {
  cat("Folder already exists:\n", singleR_dir, "\n")
}

# === Output helper ===
output_dir <- singleR_dir
outfn <- function(name) file.path(output_dir, paste0(sample_id, "_", name))

# === Load Seurat object ===
seurat_rds <- file.path(soupX_dir, paste0(sample_id, "_soupX.gs.rds"))
seu <- readRDS(seurat_rds)
cat("\nLoaded Seurat:", ncol(seu), "cells,", nrow(seu), "genes\n")

########################################################################## new code in Seurat, use resolution = 0.6 
##########################################################################

library(Seurat)
library(ggplot2)

process_seurat <- function(obj,
                           assay = "RNA",
                           nfeatures = 3000,
                           npcs = 50,
                           dims_use = 1:30,
                           resolution = 0.6,
                           plot_dir = "plots",
                           seed = 42) {
  set.seed(seed)
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  DefaultAssay(obj) <- assay
  
  # Calculate mitochondrial percentage
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # Standard Seurat workflow
  obj <- NormalizeData(obj, verbose = FALSE) %>%
    FindVariableFeatures(nfeatures = nfeatures, verbose = FALSE) %>%
    ScaleData(features = VariableFeatures(.), verbose = FALSE) %>%
    RunPCA(features = VariableFeatures(.), npcs = npcs, verbose = FALSE) %>%
    FindNeighbors(dims = dims_use, verbose = FALSE) %>%
    FindClusters(resolution = resolution, verbose = FALSE) %>%
    RunUMAP(dims = dims_use, verbose = FALSE)
  
  # Generate plots
  generate_standard_plots(obj, plot_dir, npcs)
  
  obj
}

# Helper function for plotting
generate_standard_plots <- function(obj, plot_dir, npcs) {
  save_plot <- function(filename, plot_expr) {
    png(file.path(plot_dir, filename), width = 7, height = 6, units = "in", res = 300)
    print(plot_expr)
    dev.off()
  }
  
  save_plot("01_VariableFeatures.png", VariableFeaturePlot(obj))
  save_plot("02_PCA_Elbow.png", ElbowPlot(obj, ndims = npcs))
  save_plot("03_PCA_Loadings.png", VizDimLoadings(obj, dims = 1:2, reduction = "pca"))
  save_plot("04_UMAP_by_cluster.png", 
            DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", 
                    label = TRUE, repel = TRUE))
  
  # Feature plot for QC metrics
  qc_features <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                           colnames(obj@meta.data))
  if (length(qc_features) > 0) {
    save_plot("05_UMAP_QC_metrics.png", 
              FeaturePlot(obj, reduction = "umap", features = qc_features))
  }
}

# Usage
seu <- process_seurat(seu,
                      assay = "RNA",
                      nfeatures = 3000,
                      npcs = 50,
                      dims_use = 1:30,
                      resolution = 0.5,
                      plot_dir = "plots",
                      seed = 123)

##########################################################################
##########################################################################
##########################################################################
##########################################################################

# === Convert to SingleCellExperiment ===
# sce <- as.SingleCellExperiment(seu)
# cat("Converted to SCE:", ncol(sce), "cells,", nrow(sce), "genes\n")

########################################################################## a new function : seurat -- > seu
##########################################################################
##########################################################################
##########################################################################

library(Seurat)
library(SingleCellExperiment)

seurat_to_sce_complete <- function(seurat_obj, 
                                    assay = "RNA") {
  
  cat("Converting Seurat to SingleCellExperiment...\n")
  
  # ========================================
  # 0. Check available layers
  # ========================================
  available_layers <- Layers(seurat_obj[[assay]])
  cat("Available layers:", paste(available_layers, collapse = ", "), "\n\n")
  
  # ========================================
  # 1. Extract count and normalized data using layers
  # ========================================
  counts_matrix <- GetAssayData(seurat_obj, assay = assay, layer = "counts")
  data_matrix <- GetAssayData(seurat_obj, assay = assay, layer = "data")
  
  # Create SCE with counts and logcounts
  sce <- SingleCellExperiment(
    assays = list(
      counts = counts_matrix,
      logcounts = data_matrix
    )
  )
  
  cat("✓ Transferred counts and logcounts\n")
  
  # ========================================
  # 2. Transfer all metadata
  # ========================================
  colData(sce) <- DataFrame(seurat_obj@meta.data)
  cat("✓ Transferred cell metadata (", ncol(colData(sce)), " columns)\n", sep = "")
  
  # ========================================
  # 3. Transfer gene metadata (if any)
  # ========================================
  
  # Handle both Assay5 (v5) and Assay (v4) objects
  assay_obj <- seurat_obj[[assay]]
  
  # Try to get meta.features (v4) or meta.data (v5)
  gene_meta <- tryCatch({
    if (inherits(assay_obj, "Assay5")) {
      # Seurat v5: meta.data slot
      if (length(slot(assay_obj, "meta.data")) > 0) {
        slot(assay_obj, "meta.data")
      } else {
        data.frame(row.names = rownames(sce))
      }
    } else {
      # Seurat v4: meta.features slot
      assay_obj@meta.features
    }
  }, error = function(e) {
    data.frame(row.names = rownames(sce))
  })
  
  # Add variable feature information
  var_features <- VariableFeatures(seurat_obj, assay = assay)
  
  if (nrow(gene_meta) > 0) {
    gene_meta$highly_variable <- rownames(gene_meta) %in% var_features
    
    # Ensure same order as SCE genes
    gene_meta_ordered <- gene_meta[match(rownames(sce), rownames(gene_meta)), , drop = FALSE]
    rowData(sce) <- DataFrame(gene_meta_ordered)
    
    cat("✓ Transferred gene metadata\n")
  } else {
    # Create minimal gene metadata with variable feature flag
    gene_meta <- data.frame(
      highly_variable = rownames(sce) %in% var_features,
      row.names = rownames(sce)
    )
    rowData(sce) <- DataFrame(gene_meta)
    cat("✓ Added variable feature information\n")
  }
  
  # ========================================
  # 4. Transfer dimensionality reductions
  # ========================================
  reductions <- names(seurat_obj@reductions)
  
  for (reduction in reductions) {
    reduction_data <- Embeddings(seurat_obj, reduction = reduction)
    
    # Standardize reduction names for SCE
    sce_reduction_name <- switch(
      tolower(reduction),
      "pca" = "PCA",
      "umap" = "UMAP",
      "tsne" = "TSNE",
      "harmony" = "HARMONY",
      toupper(reduction)  # Default: uppercase
    )
    
    reducedDim(sce, sce_reduction_name) <- reduction_data
    cat("✓ Transferred", reduction, "→", sce_reduction_name, 
        "(", ncol(reduction_data), "dimensions)\n")
  }
  
  # ========================================
  # 5. Transfer additional information
  # ========================================
  
  # Store Seurat clustering info
  if ("seurat_clusters" %in% colnames(colData(sce))) {
    sce$clusters <- sce$seurat_clusters
  }
  
  # Store variable features as metadata
  metadata(sce)$variable_features <- var_features
  metadata(sce)$n_variable_features <- length(var_features)
  
  # Try to get HVG method
  hvg_method <- tryCatch({
    if (inherits(assay_obj, "Assay5")) {
      # Try to extract from misc slot or commands
      "vst"  # Default assumption
    } else {
      assay_obj@var.features.method
    }
  }, error = function(e) "unknown")
  
  metadata(sce)$hvg_method <- hvg_method
  
  # Store Seurat parameters if available
  if (length(seurat_obj@commands) > 0) {
    metadata(sce)$seurat_commands <- names(seurat_obj@commands)
  }
  
  cat("✓ Transferred variable features and parameters\n")
  
  # ========================================
  # Summary
  # ========================================
  cat("\n=== CONVERSION SUMMARY ===\n")
  cat("Dimensions:", nrow(sce), "genes x", ncol(sce), "cells\n")
  cat("Assays:", paste(assayNames(sce), collapse = ", "), "\n")
  cat("Reductions:", paste(reducedDimNames(sce), collapse = ", "), "\n")
  cat("Metadata columns:", ncol(colData(sce)), "\n")
  cat("Variable features:", length(metadata(sce)$variable_features), "\n\n")
  
  return(sce)
}

# ========================================
# Usage Example
# ========================================

# Convert your Seurat object
sce_obj <- seurat_to_sce_complete(seu, assay = "RNA")

# ========================================
# Verify the conversion
# ========================================

cat("=== VERIFICATION ===\n")
cat("Assays available:", paste(assayNames(sce_obj), collapse = ", "), "\n")
cat("Reductions available:", paste(reducedDimNames(sce_obj), collapse = ", "), "\n")
cat("Cell metadata columns:", ncol(colData(sce_obj)), "\n")
cat("Gene metadata columns:", ncol(rowData(sce_obj)), "\n")
cat("Variable features stored:", length(metadata(sce_obj)$variable_features), "\n")

# Preview first few cells
cat("\nFirst 5 cells metadata preview:\n")
print(head(as.data.frame(colData(sce_obj)[, 1:min(5, ncol(colData(sce_obj)))])))

# ========================================
# Optional: Quick visualization to confirm
# ========================================

library(scater)

if ("UMAP" %in% reducedDimNames(sce_obj)) {
  p <- plotReducedDim(sce_obj, dimred = "UMAP", 
                      colour_by = "seurat_clusters",
                      text_by = "seurat_clusters")
  print(p)
}

########################################################################## continuing the previous script
##########################################################################
##########################################################################
##########################################################################

# === Inspect ===
# sce_obj <- sce

cat("\nClass:", class(sce_obj), "\n")
cat("Dimensions:", dim(sce_obj), "\n")
cat("Assays:", paste(assayNames(sce_obj), collapse = ", "), "\n\n")

if ("counts" %in% assayNames(sce_obj)) {
  cat("Counts preview:\n")
  print(assay(sce_obj, "counts")[1:5, 1:5])
}

# --- New section: Read the .soupX RDS object ---

# Construct the path to the SoupX subfolder
# soupX_dir <- file.path(path, paste0(sample_id, ".soupX"))

# cat("\nSoupX folder path:\n", soupX_dir, "\n")

# Check that the subfolder exists
# if (!dir.exists(soupX_dir)) {
#  stop("SoupX folder does not exist: ", soupX_dir)
# }

# Build full filename for the RDS object
# sce_filename <- paste0(sample_id, "_sce_obj.filtered.soupX.rds")
# sce_file <- file.path(soupX_dir, sce_filename)

# cat("\nLooking for RDS file:\n", sce_file, "\n")

# Check file presence and read it
# if (!file.exists(sce_file)) {
#  stop("SoupX-corrected SCE object not found: ", sce_file)
# } else {
#  cat("\nReading SoupX-corrected SingleCellExperiment object...\n")
#  sce_obj <- readRDS(sce_file)
#  cat("Loaded successfully ✅\n\n")
# }

# Inspect basic info
cat("Object class:", class(sce_obj), "\n")
cat("Dimensions (genes x cells):", dim(sce_obj), "\n\n")

assays(sce_obj)

library(SingleCellExperiment)

# --- Check assays available in the object ---
cat("=== Available assays ===\n")
print(assayNames(sce_obj))   # List the assay layers (e.g., counts, logcounts, corrected)

# If you want to see the structure of one assay:
cat("\nAssay 'counts' preview (first 5 genes, first 5 cells):\n")
print(assay(sce_obj, "counts")[1:5, 1:5])

# --- Check if UMAP embedding exists ---
cat("\n=== Reduced dimensions (dimensionality reductions) ===\n")
print(reducedDimNames(sce_obj))   # e.g., "PCA", "UMAP", "TSNE"

# If UMAP exists, preview the first few coordinates
if ("UMAP" %in% reducedDimNames(sce_obj)) {
  cat("\nUMAP coordinates (first 5 cells):\n")
  print(head(reducedDim(sce_obj, "UMAP")))
} else {
  cat("\nNo UMAP embedding found in sce_obj.\n")
}

# --- Show cell-level metadata ---
cat("\n=== Cell metadata (colData) ===\n")
print(colnames(colData(sce_obj)))

cat("\nFirst 5 rows of cell metadata:\n")
print(head(as.data.frame(colData(sce_obj))))

reducedDimNames(sce_obj)

library(SingleCellExperiment)
library(ggplot2)

# 1) Grab UMAP coordinates
stopifnot("UMAP" %in% reducedDimNames(sce_obj))
umap_df <- as.data.frame(reducedDim(sce_obj, "UMAP"))

# 2) Inspect and standardize column names to UMAP1 / UMAP2
print(colnames(umap_df))
if (ncol(umap_df) >= 2) {
  # If names are missing or different, set them
  if (!all(c("UMAP1","UMAP2") %in% colnames(umap_df))) {
    colnames(umap_df)[1:2] <- c("UMAP1","UMAP2")
  }
} else {
  stop("UMAP reduction does not have at least two columns.")
}

# 3) Attach a grouping variable from colData (choose what you want to color by)
umap_df$cluster <- sce_obj$seurat_clusters        # or sce_obj$ident, Phase, etc.

# 4) Plot
p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = as.factor(cluster))) +
  geom_point(size = 0.4, alpha = 0.7) +
  labs(color = "Cluster", title = paste(sample_id, "- UMAP by cluster")) +
  theme_classic()

print(p)

# Optional: save to your singleR output dir if you defined `outfn()`
if (exists("outfn")) ggsave(filename = outfn("UMAP_by_cluster.png"), plot = p, width = 6, height = 5, dpi = 300)


##########################################################################
##########################################################################
########################################################################## annotations with MTG_SEA

# SingleR annotations with MTG_SEA

suppressPackageStartupMessages({
  library(SingleR)
  library(SingleCellExperiment)
  library(scuttle)
})


## Ensure logcounts exist (SingleR prefers log-normalized assays)
if (!"logcounts" %in% assayNames(sce_obj)) {
  sce_obj <- scuttle::logNormCounts(sce_obj)
}
if (!"logcounts" %in% assayNames(mtg_sea)) {
  mtg_sea <- scuttle::logNormCounts(mtg_sea)
}

# --- normalize IDs and find common genes ---
sce_ids <- sub("\\.\\d+$", "", rownames(sce_obj))
mtg_ids <- sub("\\.\\d+$", "", rownames(mtg_sea))
common_genes <- intersect(sce_ids, mtg_ids)

cat("sce_obj genes:", length(sce_ids), "\n")
cat("mtg_sea genes:", length(mtg_ids), "\n")
cat("Common genes:", length(common_genes), "\n")

# --- subset AND align order identically ---
sce_sub <- sce_obj[match(common_genes, sce_ids), ]
mtg_sub <- mtg_sea[match(common_genes, mtg_ids), ]

# make rownames identical to avoid surprises
rownames(sce_sub) <- common_genes
rownames(mtg_sub) <- common_genes
stopifnot(identical(rownames(sce_sub), rownames(mtg_sub)))


dim(sce_obj)
dim(mtg_sea)

dim(sce_sub)
dim(mtg_sub)

# plotScoreHeatmap() displays the scores for all cells across all reference labels, which allows users to inspect the confidence of the predicted labels across the dataset. 
# Ideally, each cell (i.e., column of the heatmap) should have one score that is obviously larger than the rest, 
# indicating that it is unambiguously assigned to a single label. 
# A spread of similar scores for a given cell indicates that the assignment is uncertain, 
# though this may be acceptable if the uncertainty is distributed across similar cell types that cannot be easily resolved.

cat("delta distribution")

# Another diagnostic is based on the per-cell “deltas”, i.e., the difference between the score for the assigned label and the median across all labels for each cell. 
# Low deltas indicate that the assignment is uncertain, which is especially relevant if the cell’s true label does not exist in the reference. 
# We can inspect these deltas across cells for each label using the plotDeltaDistribution() function

# https://biostatsquid.com/singler-tutorial/

# Well, by itself, the SingleR algorithm will always assign a label to every cell. But what if the cell’s true label isn’t in the reference dataset? 
# It will still assign it a label. For this and other reasons, SingleR can assign incorrect labels.

# The developers tried to mitigate this  by removing poor-quality predictions with “low” scores. 
# They basically compute a “delta” value for each cell, which is the gap, or the difference between the score for the assigned label and the median score across all labels. 
# If the delta is small, this indicates that the cell matches all labels with the same confidence, so the assigned label is not very meaningful.

# This way, SingleR can discard cells with low delta values caused by (i) ambiguous assignments 
# with closely related reference labels and (ii) incorrect assignments that match poorly to all reference labels – 
# so in the pruned_labels column you will find ‘cleaner’ or more reliable labels.



# --- run SingleR (choose the label field you want) ---
# pred1 <- SingleR(
#  test   = sce_sub,
#  ref    = mtg_sub,
#  labels = colData(mtg_sub)$cluster_label   # or subclass_label / class_label,
#  de.method = 'wilcox'
# )

# Save the complete predictions as RDS
# saveRDS(pred1, outfn("singleR_pred1.mtg_sea.cluster_labels.rds"))

# Save predictions as CSV (main results)
# pred1_df <- data.frame(
#  cell_id = rownames(pred1),
#  labels = pred1$labels,
#  pruned_labels = pred1$pruned.labels,
#  delta_next = pred1$delta.next,
#  stringsAsFactors = FALSE
# )
# head(pred1_df)
# write.csv(pred1_df, outfn("singleR_pred1.mtg_sea.cluster_labels.csv"), row.names = FALSE)

# inspect
# cat("n cells:", ncol(sce_sub), "\n")
# print(unique(pred1$labels))
# print(table(pred1$labels))

# str(pred1)

#  ..@ listData       :List of 4
#  .. ..$ scores       : num [1:15103, 1:127] 0.295 0.115 0.166 0.129 0.169 ...
#  .. .. ..- attr(*, "dimnames")=List of 2
#  .. .. .. ..$ : NULL
#  .. .. .. ..$ : chr [1:127] "Astro L1 FGFR3 SERPINI2" "Astro L1-6 FGFR3 AQP1" "Astro L1-6 FGFR3 PLCG1" "Endo L2-5 NOSTRIN SRGN" ...
#  .. ..$ labels       : chr [1:15103] "OPC L1-6 PDGFRA COL20A1" "Micro L1-6 TYROBP CD74" "Inh L3-5 SST GGTLC3" "Oligo L5-6 OPALIN LDLRAP1" ...
#  .. ..$ delta.next   : num [1:15103] 0.2445 0.2661 0.1594 0.0477 0.2012 ...
#  .. ..$ pruned.labels: chr [1:15103] "OPC L1-6 PDGFRA COL20A1" "Micro L1-6 TYROBP CD74" "Inh L3-5 SST GGTLC3" "Oligo L5-6 OPALIN LDLRAP1" ...

# length(pred1$labels)
# length(pred1$pruned.labels)

# print(unique(pred1$labels))

# print(table(pred1$labels))

# summary(is.na(pred1$labels))

# print(unique(pred1$pruned.labels))

# print(table(pred1$pruned.labels))

# summary(is.na(pred1$pruned.labels))

# options(repr.plot.width = 12, repr.plot.height = 10)  # adjust to your liking
# plotScoreHeatmap(pred1)

# options(repr.plot.width = 12, repr.plot.height = 10)

# Create the plot
# p_heatmap <- plotScoreHeatmap(pred1)

# Save it
# ggsave(outfn("singleR_pred1.mtg_sea.score_heatmap.png"), 
#       plot = p_heatmap, width = 12, height = 10, dpi = 300)

# options(repr.plot.width = 12, repr.plot.height = 60)  # adjust to your liking
# plotDeltaDistribution(pred1, ncol = 4, dots.on.top = FALSE)

# options(repr.plot.width = 12, repr.plot.height = 60)

# Save as PNG
# png(outfn("singleR_pred1.mtg_sea.delta_distribution.png"), 
#    width = 12, height = 60, units = "in", res = 300)
# plotDeltaDistribution(pred1, ncol = 4, dots.on.top = FALSE)
# dev.off()

# Save as PDF (better for very tall plots)
# pdf(outfn("singleR_pred1.mtg_sea.delta_distribution.pdf"), 
#    width = 12, height = 60)
# plotDeltaDistribution(pred1, ncol = 4, dots.on.top = FALSE)
# dev.off()

# cat("Delta distribution plot saved!\n")

# Keep any types that have more than 10 cells to the label, and put those labels back on our Seurat object and plot our on our umap.

# lbls.keep <- table(pred1$labels)>10
# sce_sub$SingleR.labels <- ifelse(lbls.keep[pred1$labels], pred1$labels, 'Other')
# DimPlot(sce_sub, reduction='umap', group.by='SingleR.labels')

# dim(sce_sub)
# dim(pred1)

# --- Transfer SingleR predictions to sce_sub ---
# sce_sub$singleR_labels.mtg_sea <- pred1$labels
# sce_sub$singleR_scores.mtg_sea <- pred1$scores  # optional: confidence scores

# Verify the transfer
# cat("Labels transferred successfully!\n")
# table(sce_sub$singleR_labels.mtg_sea)

# Optional: add pruned labels (removes low-confidence predictions)
# sce_sub$singleR_labels_pruned.mtg_sea <- pred1$pruned.labels

# library(scater)

# If you have UMAP computed
# plotReducedDim(sce_sub, dimred = "UMAP", colour_by = "singleR_labels.mtg_sea")

# If you have UMAP computed
# plotReducedDim(sce_sub, dimred = "UMAP", colour_by = "singleR_labels_pruned.mtg_sea")

# If you want to transfer to the FULL sce_obj (not just subset):

# Create vectors for both label types
# all_labels <- rep(NA, ncol(sce_obj))
# all_labels_pruned <- rep(NA, ncol(sce_obj))
# names(all_labels) <- colnames(sce_obj)
# names(all_labels_pruned) <- colnames(sce_obj)

# Fill in the predicted labels
# all_labels[colnames(sce_sub)] <- pred1$labels
# all_labels_pruned[colnames(sce_sub)] <- pred1$pruned.labels

# Add to full object
# sce_obj$singleR_labels.mtg_sea <- all_labels
# sce_obj$singleR_labels_pruned.mtg_sea <- all_labels_pruned

# Summary
# cat("=== SingleR Label Transfer Summary ===\n")
# cat("Total cells in sce_obj:", ncol(sce_obj), "\n")
# cat("Cells with regular labels:", sum(!is.na(all_labels)), "\n")
# cat("Cells with pruned labels:", sum(!is.na(all_labels_pruned)), "\n")
# cat("Cells removed by pruning:", sum(!is.na(all_labels)) - sum(!is.na(all_labels_pruned)), "\n")

# For many cell types, use a better color palette

# library(viridis)

# n_types <- length(unique(na.omit(sce_obj$singleR_labels.mtg_sea)))
# mycols <- viridis(n_types, option = "turbo")

# options(repr.plot.width = 12, repr.plot.height = 12)

# p_umap_fancy <- plotReducedDim(
#  sce_obj, 
#  dimred = "UMAP", 
#  colour_by = "singleR_labels.mtg_sea"
# ) +
#  scale_color_manual(values = mycols, na.value = "grey80") +
#  theme_minimal(base_size = 12) +
#  theme(
#    legend.position = "bottom",
#    legend.box = "horizontal",
#    legend.text = element_text(size = 8),
#    legend.key.size = unit(0.5, "lines"),
#    plot.margin = margin(10, 10, 10, 10)
#  ) +
#  guides(
#    colour = guide_legend(
#      nrow = 4,              # Number of legend rows
#      byrow = TRUE,
#      override.aes = list(size = 3),
#      title.position = "top"
#    )
#  ) +
#  labs(title = "SingleR Cell Type Annotation (MTG SEA Reference)",
#       colour = "Cell Type")

# print(p_umap_fancy)

# ggsave(outfn("UMAP_singleR_labels.mtg_sea.fancy.png"), 
#        plot = p_umap_fancy, width = 12, height = 12, dpi = 300)

# For many cell types, use a better color palette - PRUNED LABELS
# library(viridis)

# n_types_pruned <- length(unique(na.omit(sce_obj$singleR_labels_pruned.mtg_sea)))
# mycols_pruned <- viridis(n_types_pruned, option = "turbo")

# options(repr.plot.width = 12, repr.plot.height = 12)

# p_umap_fancy_pruned <- plotReducedDim(
#  sce_obj, 
#  dimred = "UMAP", 
#  colour_by = "singleR_labels_pruned.mtg_sea"
# ) +
#  scale_color_manual(values = mycols_pruned, na.value = "grey80") +
#  theme_minimal(base_size = 12) +
#  theme(
#    legend.position = "bottom",
#    legend.box = "horizontal",
#    legend.text = element_text(size = 8),
#    legend.key.size = unit(0.5, "lines"),
#    plot.margin = margin(10, 10, 10, 10)
#  ) +
#  guides(
#    colour = guide_legend(
#      nrow = 4,              # Number of legend rows
#      byrow = TRUE,
#      override.aes = list(size = 3),
#      title.position = "top"
#    )
#  ) +
#  labs(title = "SingleR Cell Type Annotation - Pruned (High Confidence)",
#       colour = "Cell Type")

# print(p_umap_fancy_pruned)

# ggsave(outfn("UMAP_singleR_labels_pruned.mtg_sea.fancy.png"), 
#        plot = p_umap_fancy_pruned, width = 12, height = 12, dpi = 300)

# Create a summary table
# label_counts <- table(sce_obj$singleR_labels.mtg_sea)
# label_counts_df <- data.frame(
#  cell_type = names(label_counts),
#  count = as.numeric(label_counts)
# )
# label_counts_df <- label_counts_df[order(label_counts_df$count, decreasing = TRUE), ]

# cat("\nTop 30 cell types by count:\n")
# print(head(label_counts_df, 30))

# Save the summary
# write.csv(label_counts_df, outfn("singleR_label_counts.mtg_sea.csv"), row.names = FALSE)



# --- run SingleR (choose the label field you want) ---
pred2 <- SingleR(
  test   = sce_sub,
  ref    = mtg_sub,
  labels = colData(mtg_sub)$subclass_label,
  de.method = 'wilcox'
)

# inspect
cat("n cells:", ncol(sce_sub), "\n")
print(head(table(pred2$labels), 10))

# Save the complete predictions as RDS
saveRDS(pred2, outfn("singleR_pred2.mtg_sea.subclass_labels.rds"))

# Save predictions as CSV (main results)
pred2_df <- data.frame(
  cell_id = rownames(pred2),
  labels = pred2$labels,
  pruned_labels = pred2$pruned.labels,
  delta_next = pred2$delta.next,
  stringsAsFactors = FALSE
)

head(pred2_df)
write.csv(pred2_df, outfn("singleR_pred2.mtg_sea.subclass_labels.csv"), row.names = FALSE)

# inspect
cat("n cells:", ncol(sce_sub), "\n")
print(unique(pred2$labels))
print(table(pred2$labels))
summary(is.na(pred2$labels))
print(unique(pred2$pruned.labels))
print(table(pred2$pruned.labels))
summary(is.na(pred2$pruned.labels))

# Plot and save score heatmap
options(repr.plot.width = 12, repr.plot.height = 10)
p_heatmap2 <- plotScoreHeatmap(pred2)
ggsave(outfn("singleR_pred2.mtg_sea.score_heatmap.png"), 
       plot = p_heatmap2, width = 12, height = 10, dpi = 300)

# Plot and save delta distribution
options(repr.plot.width = 12, repr.plot.height = 40)
png(outfn("singleR_pred2.mtg_sea.delta_distribution.png"), 
    width = 12, height = 40, units = "in", res = 300)
plotDeltaDistribution(pred2, ncol = 4, dots.on.top = FALSE)
dev.off()

pdf(outfn("singleR_pred2.mtg_sea.delta_distribution.pdf"), 
    width = 12, height = 40)
plotDeltaDistribution(pred2, ncol = 4, dots.on.top = FALSE)
dev.off()

# Transfer to sce_sub
sce_sub$singleR_labels.mtg_sea.subclass <- pred2$labels
sce_sub$singleR_labels_pruned.mtg_sea.subclass <- pred2$pruned.labels

# Transfer to full sce_obj
all_labels2 <- rep(NA, ncol(sce_obj))
all_labels2_pruned <- rep(NA, ncol(sce_obj))
names(all_labels2) <- colnames(sce_obj)
names(all_labels2_pruned) <- colnames(sce_obj)

all_labels2[colnames(sce_sub)] <- pred2$labels
all_labels2_pruned[colnames(sce_sub)] <- pred2$pruned.labels

sce_obj$singleR_labels.mtg_sea.subclass <- all_labels2
sce_obj$singleR_labels_pruned.mtg_sea.subclass <- all_labels2_pruned

cat("=== SingleR pred2 (Subclass) Transfer Summary ===\n")
cat("Total cells in sce_obj:", ncol(sce_obj), "\n")
cat("Cells with regular labels:", sum(!is.na(all_labels2)), "\n")
cat("Cells with pruned labels:", sum(!is.na(all_labels2_pruned)), "\n")
cat("Cells removed by pruning:", sum(!is.na(all_labels2)) - sum(!is.na(all_labels2_pruned)), "\n")

# ========================================
# PHASE 1: UMAP with regular labels (pred2 - subclass)
# ========================================


library(viridis)
library(ggplot2)
library(dplyr)

# Ensure UMAP exists
if (!"UMAP" %in% reducedDimNames(sce_obj)) {
  cat("Computing UMAP...\n")
  if (!"PCA" %in% reducedDimNames(sce_obj)) {
    sce_obj <- runPCA(sce_obj, ncomponents = 50)
  }
  sce_obj <- runUMAP(sce_obj, dimred = "PCA")
}

# Extract UMAP coordinates and labels
umap_coords <- reducedDim(sce_obj, "UMAP")
umap_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  subclass = colData(sce_obj)$singleR_labels.mtg_sea.subclass
)

# Calculate centroids for each subclass
centroids_subclass <- umap_df %>%
  filter(!is.na(subclass)) %>%
  group_by(subclass) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2)
  )

# Create color palette
n_types2 <- length(unique(na.omit(colData(sce_obj)$singleR_labels.mtg_sea.subclass)))
mycols2 <- viridis(n_types2, option = "turbo")

options(repr.plot.width = 12, repr.plot.height = 12)

p_umap_fancy2 <- plotReducedDim(
  sce_obj, 
  dimred = "UMAP", 
  colour_by = "singleR_labels.mtg_sea.subclass"
) +
  scale_color_manual(values = mycols2, na.value = "grey80") +
  geom_text(data = centroids_subclass, 
            aes(x = UMAP1, y = UMAP2, label = subclass),
            color = "black",
            size = 6,
            fontface = "bold",
            check_overlap = TRUE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 12),  # INCREASE THIS (was 8)
    legend.key.size = unit(0.8, "lines"),   # Also increase key size to match
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(
    colour = guide_legend(
      nrow = 4,
      byrow = TRUE,
      override.aes = list(size = 3),
      title.position = "top"
    )
  ) +
  labs(title = "SingleR Subclass Annotation (MTG SEA Reference)",
       colour = "Cell Subclass")


print(p_umap_fancy2)

ggsave(outfn("UMAP_singleR_labels.mtg_sea.subclass.fancy_with_labels.png"), 
       plot = p_umap_fancy2, width = 12, height = 12, dpi = 300)

options(repr.plot.width = 12, repr.plot.height = 12)

p_umap_fancy2 <- plotReducedDim(
  sce_obj, 
  dimred = "UMAP", 
  colour_by = "singleR_labels.mtg_sea.subclass"
) +
  scale_color_manual(values = mycols2, na.value = "grey80") +
  geom_text(data = centroids_subclass, 
            aes(x = UMAP1, y = UMAP2, label = subclass),
            color = "black",
            size = 10,
            fontface = "bold",
            check_overlap = TRUE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 12),  # INCREASE THIS (was 8)
    legend.key.size = unit(0.8, "lines"),   # Also increase key size to match
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(
    colour = guide_legend(
      nrow = 8,
      byrow = TRUE,
      override.aes = list(size = 3),
      title.position = "top"
    )
  ) +
  labs(title = "SingleR Subclass Annotation (MTG SEA Reference)",
       colour = "Cell Subclass")


print(p_umap_fancy2)

ggsave(outfn("UMAP_singleR_labels.mtg_sea.subclass.fancy_with_labels2.png"), 
       plot = p_umap_fancy2, width = 12, height = 12, dpi = 300)



# ========================================
# PHASE 2: UMAP with pruned labels (pred2 - subclass)
# ========================================

# library(viridis)

# n_types2_pruned <- length(unique(na.omit(sce_obj$singleR_labels_pruned.mtg_sea.subclass)))
# mycols2_pruned <- viridis(n_types2_pruned, option = "turbo")

# options(repr.plot.width = 12, repr.plot.height = 12)

# p_umap_fancy2_pruned <- plotReducedDim(
#  sce_obj, 
#  dimred = "UMAP", 
#  colour_by = "singleR_labels_pruned.mtg_sea.subclass"
# ) +
#  scale_color_manual(values = mycols2_pruned, na.value = "grey80") +
#  theme_minimal(base_size = 12) +
#  theme(
#    legend.position = "bottom",
#    legend.box = "horizontal",
#    legend.text = element_text(size = 8),
#    legend.key.size = unit(0.5, "lines"),
#    plot.margin = margin(10, 10, 10, 10)
#  ) +
#  guides(
#    colour = guide_legend(
#      nrow = 4,
#      byrow = TRUE,
#      override.aes = list(size = 3),
#      title.position = "top"
#    )
#  ) +
#  labs(title = "SingleR Subclass Annotation - Pruned (High Confidence)",
#       colour = "Cell Subclass")

# print(p_umap_fancy2_pruned)

# ggsave(outfn("UMAP_singleR_labels_pruned.mtg_sea.subclass.fancy.png"), 
#        plot = p_umap_fancy2_pruned, width = 12, height = 12, dpi = 300)

# ========================================
# PHASE 2: UMAP with pruned labels (pred2 - subclass) + LABEL OVERLAY
# ========================================
library(viridis)
library(ggplot2)
library(dplyr)

# Extract UMAP coordinates and PRUNED labels
umap_coords <- reducedDim(sce_obj, "UMAP")
umap_df_pruned <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  subclass_pruned = colData(sce_obj)$singleR_labels_pruned.mtg_sea.subclass
)

# Calculate centroids for each PRUNED subclass
centroids_subclass_pruned <- umap_df_pruned %>%
  filter(!is.na(subclass_pruned)) %>%
  group_by(subclass_pruned) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2)
  )

# Create color palette
n_types2_pruned <- length(unique(na.omit(colData(sce_obj)$singleR_labels_pruned.mtg_sea.subclass)))
mycols2_pruned <- viridis(n_types2_pruned, option = "turbo")

options(repr.plot.width = 12, repr.plot.height = 12)

p_umap_fancy2_pruned <- plotReducedDim(
  sce_obj, 
  dimred = "UMAP", 
  colour_by = "singleR_labels_pruned.mtg_sea.subclass"
) +
  scale_color_manual(values = mycols2_pruned, na.value = "grey80") +
  # ADD CLUSTER LABELS ON TOP
  geom_text(data = centroids_subclass_pruned, 
            aes(x = UMAP1, y = UMAP2, label = subclass_pruned),
            color = "black",
            size = 6,
            fontface = "bold",
            check_overlap = TRUE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "lines"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(
    colour = guide_legend(
      nrow = 4,
      byrow = TRUE,
      override.aes = list(size = 3),
      title.position = "top"
    )
  ) +
  labs(title = "SingleR Subclass Annotation - Pruned (High Confidence)",
       colour = "Cell Subclass")

print(p_umap_fancy2_pruned)

ggsave(outfn("UMAP_singleR_labels_pruned.mtg_sea.subclass.fancy_with_labels.png"), 
       plot = p_umap_fancy2_pruned, width = 12, height = 12, dpi = 300)

# ========================================
# PHASE 2: UMAP with pruned labels (pred2 - subclass) + LABEL OVERLAY
# ========================================
library(viridis)
library(ggplot2)
library(dplyr)

# Extract UMAP coordinates and PRUNED labels
umap_coords <- reducedDim(sce_obj, "UMAP")
umap_df_pruned <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  subclass_pruned = colData(sce_obj)$singleR_labels_pruned.mtg_sea.subclass
)

# Calculate centroids for each PRUNED subclass
centroids_subclass_pruned <- umap_df_pruned %>%
  filter(!is.na(subclass_pruned)) %>%
  group_by(subclass_pruned) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2)
  )

# Create color palette
n_types2_pruned <- length(unique(na.omit(colData(sce_obj)$singleR_labels_pruned.mtg_sea.subclass)))
mycols2_pruned <- viridis(n_types2_pruned, option = "turbo")

options(repr.plot.width = 12, repr.plot.height = 12)

p_umap_fancy2_pruned <- plotReducedDim(
  sce_obj, 
  dimred = "UMAP", 
  colour_by = "singleR_labels_pruned.mtg_sea.subclass"
) +
  scale_color_manual(values = mycols2_pruned, na.value = "grey80") +
  # ADD CLUSTER LABELS ON TOP
  geom_text(data = centroids_subclass_pruned, 
            aes(x = UMAP1, y = UMAP2, label = subclass_pruned),
            color = "black",
            size = 8,
            fontface = "bold",
            check_overlap = TRUE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "lines"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(
    colour = guide_legend(
      nrow = 4,
      byrow = TRUE,
      override.aes = list(size = 3),
      title.position = "top"
    )
  ) +
  labs(title = "SingleR Subclass Annotation - Pruned (High Confidence)",
       colour = "Cell Subclass")

print(p_umap_fancy2_pruned)

ggsave(outfn("UMAP_singleR_labels_pruned.mtg_sea.subclass.fancy_with_labels2.png"), 
       plot = p_umap_fancy2_pruned, width = 12, height = 12, dpi = 300)



# ========================================
# PHASE 3: Summary table (pred2 - subclass)
# ========================================

# Create a summary table
label_counts2 <- table(sce_obj$singleR_labels.mtg_sea.subclass)
label_counts2_df <- data.frame(
  cell_subclass = names(label_counts2),
  count = as.numeric(label_counts2)
)
label_counts2_df <- label_counts2_df[order(label_counts2_df$count, decreasing = TRUE), ]

cat("\nTop 30 cell subclasses by count:\n")
print(head(label_counts2_df, 30))

# Save the summary
write.csv(label_counts2_df, outfn("singleR_label_counts.mtg_sea.subclass.csv"), row.names = FALSE)



# ========================================
# COMPREHENSIVE SCE OBJECT STRUCTURE
# ========================================

cat("=== BASIC INFORMATION ===\n")
cat("Object class:", class(sce_obj), "\n")
cat("Dimensions (genes x cells):", nrow(sce_obj), "genes x", ncol(sce_obj), "cells\n")
cat("Memory size:", format(object.size(sce_obj), units = "auto"), "\n\n")

# ========================================
# ASSAYS (count matrices)
# ========================================
cat("=== ASSAYS (Expression Matrices) ===\n")
print(assayNames(sce_obj))
cat("\n")

# Preview first assay
if (length(assayNames(sce_obj)) > 0) {
  cat("Preview of first assay '", assayNames(sce_obj)[1], "':\n", sep = "")
  print(assay(sce_obj, assayNames(sce_obj)[1])[1:5, 1:5])
  cat("\n")
}

# ========================================
# DIMENSIONALITY REDUCTIONS
# ========================================
cat("=== DIMENSIONALITY REDUCTIONS ===\n")
print(reducedDimNames(sce_obj))

if (length(reducedDimNames(sce_obj)) > 0) {
  for (dimred in reducedDimNames(sce_obj)) {
    dims <- dim(reducedDim(sce_obj, dimred))
    cat(sprintf("  %s: %d cells x %d dimensions\n", dimred, dims[1], dims[2]))
  }
}
cat("\n")

# ========================================
# COLUMN DATA (Cell metadata)
# ========================================
cat("=== COLUMN DATA (Cell Metadata) ===\n")
cat("Number of metadata columns:", ncol(colData(sce_obj)), "\n")
cat("\nColumn names:\n")
print(colnames(colData(sce_obj)))
cat("\n")

# Show summary of each metadata column
cat("Summary of metadata columns:\n")
for (col in colnames(colData(sce_obj))) {
  col_data <- colData(sce_obj)[[col]]
  
  if (is.numeric(col_data)) {
    cat(sprintf("  %-50s: numeric, range [%.2f, %.2f]\n", 
                col, min(col_data, na.rm = TRUE), max(col_data, na.rm = TRUE)))
  } else if (is.factor(col_data) || is.character(col_data)) {
    n_unique <- length(unique(na.omit(col_data)))
    n_na <- sum(is.na(col_data))
    cat(sprintf("  %-50s: %d unique values, %d NAs\n", col, n_unique, n_na))
  } else {
    cat(sprintf("  %-50s: %s\n", col, class(col_data)))
  }
}
cat("\n")

# ========================================
# ROW DATA (Gene metadata)
# ========================================
cat("=== ROW DATA (Gene Metadata) ===\n")
cat("Number of gene metadata columns:", ncol(rowData(sce_obj)), "\n")

if (ncol(rowData(sce_obj)) > 0) {
  cat("\nGene metadata column names:\n")
  print(colnames(rowData(sce_obj)))
  cat("\n")
  
  # Preview first few rows
  cat("Preview of gene metadata (first 5 genes):\n")
  print(head(rowData(sce_obj), 5))
} else {
  cat("No gene metadata available\n")
}
cat("\n")

# ========================================
# SINGLER ANNOTATIONS (if present)
# ========================================
cat("=== SINGLER ANNOTATIONS ===\n")
singleR_cols <- grep("singleR", colnames(colData(sce_obj)), value = TRUE)

if (length(singleR_cols) > 0) {
  for (col in singleR_cols) {
    n_annotated <- sum(!is.na(colData(sce_obj)[[col]]))
    n_unique <- length(unique(na.omit(colData(sce_obj)[[col]])))
    cat(sprintf("  %-50s: %5d cells annotated, %3d unique types\n", 
                col, n_annotated, n_unique))
  }
  
  cat("\n  Unique cell types in singleR_labels.mtg_sea.subclass:\n")
  if ("singleR_labels.mtg_sea.subclass" %in% colnames(colData(sce_obj))) {
    labels_table <- table(colData(sce_obj)$singleR_labels.mtg_sea.subclass)
    labels_sorted <- sort(labels_table, decreasing = TRUE)
    print(head(labels_sorted, 20))
  }
} else {
  cat("No SingleR annotations found\n")
}
cat("\n")

# ========================================
# SAMPLE PREVIEW
# ========================================
cat("=== SAMPLE PREVIEW (First 5 cells) ===\n")
sample_cols <- c("singleR_labels.mtg_sea.subclass", 
                 "singleR_labels_pruned.mtg_sea.subclass")
sample_cols <- sample_cols[sample_cols %in% colnames(colData(sce_obj))]

if (length(sample_cols) > 0) {
  preview_df <- as.data.frame(colData(sce_obj)[1:min(5, ncol(sce_obj)), sample_cols])
  rownames(preview_df) <- colnames(sce_obj)[1:min(5, ncol(sce_obj))]
  print(preview_df)
}



# ========================================
# Save the SCE object with all annotations
# ========================================

cat("\n=== Saving SingleCellExperiment Object ===\n")

# Save the complete sce_obj with all SingleR annotations
output_file <- outfn("sce_obj.filtered.soupX.annotated.singleR.rds")

saveRDS(sce_obj, file = output_file)

cat("SCE object saved successfully to:\n")
cat(output_file, "\n")
cat("File size:", format(file.size(output_file), units = "auto"), "\n\n")

# Verify the save
cat("Verifying save...\n")
test_load <- readRDS(output_file)
cat("Verification successful! Object can be loaded.\n")
cat("Loaded object class:", class(test_load), "\n")
cat("Loaded object dimensions:", dim(test_load), "\n\n")



# ========================================
# Convert SCE to Seurat and Save
# ========================================

library(Seurat)
library(SingleCellExperiment)

cat("\n=== Converting SCE to Seurat Object ===\n")

# Convert SCE to Seurat
seurat_obj <- as.Seurat(sce_obj, counts = "counts", data = "logcounts")

cat("Conversion successful!\n")
cat("Seurat object dimensions:", nrow(seurat_obj), "genes x", ncol(seurat_obj), "cells\n\n")

# Transfer all metadata from colData to Seurat metadata
cat("Transferring metadata...\n")
for (col in colnames(colData(sce_obj))) {
  if (!col %in% colnames(seurat_obj@meta.data)) {
    seurat_obj@meta.data[[col]] <- colData(sce_obj)[[col]]
  }
}

# Transfer dimensionality reductions
cat("Transferring dimensionality reductions...\n")
for (dimred in reducedDimNames(sce_obj)) {
  dimred_data <- reducedDim(sce_obj, dimred)
  
  if (dimred == "PCA") {
    seurat_obj[["pca"]] <- CreateDimReducObject(
      embeddings = dimred_data,
      key = "PC_",
      assay = DefaultAssay(seurat_obj)
    )
  } else if (dimred == "UMAP") {
    # Ensure UMAP has exactly 2 dimensions
    if (ncol(dimred_data) == 2) {
      colnames(dimred_data) <- c("UMAP_1", "UMAP_2")
      seurat_obj[["umap"]] <- CreateDimReducObject(
        embeddings = dimred_data,
        key = "UMAP_",
        assay = DefaultAssay(seurat_obj)
      )
    }
  } else if (dimred == "TSNE") {
    seurat_obj[["tsne"]] <- CreateDimReducObject(
      embeddings = dimred_data,
      key = "tSNE_",
      assay = DefaultAssay(seurat_obj)
    )
  }
}

cat("\n=== Seurat Object Summary (Before Saving) ===\n")
cat("Assays:", paste(names(seurat_obj@assays), collapse = ", "), "\n")
cat("Reductions:", paste(names(seurat_obj@reductions), collapse = ", "), "\n")
cat("Metadata columns:", ncol(seurat_obj@meta.data), "\n")
cat("SingleR columns:", length(grep("singleR", colnames(seurat_obj@meta.data), value = TRUE)), "\n\n")

# Save the Seurat object
cat("=== Saving Seurat Object ===\n")
output_file <- outfn("seurat_obj.filtered.soupX.annotated.singleR.rds")

saveRDS(seurat_obj, file = output_file)

cat("Seurat object saved successfully to:\n")
cat(output_file, "\n")
cat("File size:", format(file.size(output_file), units = "auto"), "\n\n")

# Verify the save
cat("Verifying save...\n")
test_load <- readRDS(output_file)
cat("Verification successful!\n")
cat("Loaded object class:", class(test_load), "\n")
cat("Loaded object dimensions:", nrow(test_load), "genes x", ncol(test_load), "cells\n")
cat("SingleR columns present:", 
    length(grep("singleR", colnames(test_load@meta.data), value = TRUE)), "\n\n")



# ========================================
# Read SCE object
# ========================================

cat("\n=== Loading SingleCellExperiment Object ===\n")
sce_obj <- readRDS(outfn("sce_obj.filtered.soupX.annotated.singleR.rds"))

cat("SCE Object loaded successfully!\n\n")

# SCE Object Summary
cat("=== SCE OBJECT SUMMARY ===\n")
cat("Class:", class(sce_obj), "\n")
cat("Dimensions:", nrow(sce_obj), "genes x", ncol(sce_obj), "cells\n")
cat("Memory size:", format(object.size(sce_obj), units = "auto"), "\n\n")

# Assays
cat("Assays:\n")
print(assayNames(sce_obj))
cat("\n")

# Dimensionality reductions
cat("Dimensionality Reductions:\n")
print(reducedDimNames(sce_obj))
cat("\n")

# Metadata columns
cat("Metadata Columns (", ncol(colData(sce_obj)), " total):\n", sep = "")
print(colnames(colData(sce_obj)))
cat("\n")

# SingleR annotations summary
cat("=== SingleR Annotations (SCE) ===\n")
singleR_cols_sce <- grep("singleR_labels", colnames(colData(sce_obj)), value = TRUE)
for (col in singleR_cols_sce) {
  n_annotated <- sum(!is.na(colData(sce_obj)[[col]]))
  n_unique <- length(unique(na.omit(colData(sce_obj)[[col]])))
  cat(sprintf("%-45s: %5d cells, %3d unique types\n", col, n_annotated, n_unique))
}
cat("\n")





# ========================================
# Read Seurat object
cat("\n=== Loading Seurat Object ===\n")
seurat_obj <- readRDS(outfn("seurat_obj.filtered.soupX.annotated.singleR.rds"))

cat("Seurat Object loaded successfully!\n\n")

# Seurat Object Summary
cat("=== SEURAT OBJECT SUMMARY ===\n")
cat("Class:", class(seurat_obj), "\n")
cat("Dimensions:", nrow(seurat_obj), "genes x", ncol(seurat_obj), "cells\n")
cat("Memory size:", format(object.size(seurat_obj), units = "auto"), "\n\n")

# Assays
cat("Assays:\n")
print(names(seurat_obj@assays))
cat("\n")

# Dimensionality reductions
cat("Dimensionality Reductions:\n")
print(names(seurat_obj@reductions))
cat("\n")

# Metadata columns
cat("Metadata Columns (", ncol(seurat_obj@meta.data), " total):\n", sep = "")
print(colnames(seurat_obj@meta.data))
cat("\n")

# SingleR annotations summary
cat("=== SingleR Annotations (Seurat) ===\n")
singleR_cols_seurat <- grep("singleR_labels", colnames(seurat_obj@meta.data), value = TRUE)
for (col in singleR_cols_seurat) {
  n_annotated <- sum(!is.na(seurat_obj@meta.data[[col]]))
  n_unique <- length(unique(na.omit(seurat_obj@meta.data[[col]])))
  cat(sprintf("%-45s: %5d cells, %3d unique types\n", col, n_annotated, n_unique))
}
cat("\n")



# Show first few rows of SingleR labels in Seurat
cat("=== Sample of Seurat SingleR Labels (first 5 cells) ===\n")

# Find SingleR columns
singleR_cols <- grep("singleR_labels", colnames(seurat_obj@meta.data), value = TRUE)

if (length(singleR_cols) > 0) {
  n_rows_to_show <- min(5, nrow(seurat_obj@meta.data))
  
  # Create a data frame with cell IDs and SingleR labels
  labels_df <- data.frame(
    cell_id = rownames(seurat_obj@meta.data)[1:n_rows_to_show]
  )
  
  for (col in singleR_cols) {
    labels_df[[col]] <- seurat_obj@meta.data[[col]][1:n_rows_to_show]
  }
  
  print(labels_df)
} else {
  cat("No SingleR label columns found\n")
}

cat("\n")

