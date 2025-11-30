# ============================================================================
# Script: integration.18samples.packageSeuratIntegrate.v15.R
# Description: Seurat integration analysis for 18 samples using SeuratIntegrate (INTEGRATED WITH CLUSTERING METRICS)
# Methods: Harmony, Scanorama, scVI (Combat, BBKNN, trVAE excluded)
# Output: PNG files only (PDF disabled)
# 
# NEW IN v15:
# - Integrated clustering comparison metrics (ARI, NMI, Silhouette, LISI, Purity, Batch Correction Scores)
# 
# FIXES APPLIED (v14):
# 1. CRITICAL: Fixed scVI reduction name inconsistency - checks actual reduction name
# 2. MODERATE: Removed redundant initial clustering (now only done in post-processing)
# 3. MODERATE: Fixed graph name consistency (all use post-processing graphs)
# 4. MINOR: Added explicit reduction = "pca" parameters for clarity
# 5. MINOR: Kept initial UMAP for visualization but clarified its purpose
#
# Previous fixes:
# Fix v13: Improved layer detection to check ALL layers (not just counts) using grepl("\\.", layers)
#          to properly detect split layers like 'counts.1', 'data.2', 'scale.data.3'
# Fix v12: Corrected cluster visualization to use method-specific cluster references
# Features: UMAP plots with correct cluster labels for all integration methods
# Features: Error handling, memory monitoring (85% threshold), optimization, comprehensive logging
# ============================================================================

# Record start time
script_start_time <- Sys.time()
cat("\n")
cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("Script started at:", format(script_start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

Sys.setenv(TMPDIR = "/mnt/nfs/CX000008_DS1/projects/btanasa/tmp")
tempdir()

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
library(SeuratIntegrate)
library(lisi)
library(distances)
library(aricode)  # For ARI, NMI metrics
library(cluster)  # For silhouette score

# To increase the limit (e.g., to 70 GiB)
options(future.globals.maxSize = 700 * 1024^3)

# ============================================================================
# Function Definitions
# ============================================================================

# Summary function
show_summary <- function(obj, name) {
  cat(sprintf("\n%s:\n", name))
  cat(sprintf("  Genes: %d\n", nrow(obj)))
  cat(sprintf("  Cells: %d\n", ncol(obj)))
  cat(sprintf("  UMAP: %s\n", ifelse("umap" %in% names(obj@reductions), "Yes", "No")))
}

# Memory monitoring function
check_memory_usage <- function(threshold = 0.85, stop_on_exceed = TRUE) {
  # Get system memory information (Linux)
  tryCatch({
    # Get total and available memory from /proc/meminfo
    meminfo <- readLines("/proc/meminfo", n = 3)
    mem_total_kb <- as.numeric(gsub("[^0-9]", "", meminfo[1]))
    mem_available_kb <- as.numeric(gsub("[^0-9]", "", meminfo[3]))
    
    if (is.na(mem_available_kb)) {
      # Fallback: use MemFree if MemAvailable not available
      mem_free_kb <- as.numeric(gsub("[^0-9]", "", meminfo[2]))
      mem_used_kb <- mem_total_kb - mem_free_kb
    } else {
      mem_used_kb <- mem_total_kb - mem_available_kb
    }
    
    mem_used_gb <- mem_used_kb / 1024^2
    mem_total_gb <- mem_total_kb / 1024^2
    mem_percent <- mem_used_kb / mem_total_kb
    
    cat(sprintf("Memory Status: %.2f GB used / %.2f GB total (%.1f%%)\n", 
                mem_used_gb, mem_total_gb, mem_percent * 100))
    
    if (mem_percent >= threshold) {
      warning_msg <- sprintf(
        "WARNING: Memory usage (%.1f%%) exceeds threshold (%.1f%%)!\n",
        mem_percent * 100, threshold * 100
      )
      cat(warning_msg)
      
      if (stop_on_exceed) {
        stop("Script stopped due to high memory usage. Please free memory or reduce data size.\n")
      }
      return(FALSE)
    }
    return(TRUE)
  }, error = function(e) {
    # Fallback: Use R's gc() to estimate memory
    cat("Note: Could not read system memory info, using R memory stats\n")
    mem_stats <- gc(verbose = FALSE)
    # Get max memory used by R
    max_mem_mb <- sum(mem_stats[, "used"]) / 1024
    cat(sprintf("R Memory: %.2f MB used\n", max_mem_mb))
    # Can't check system total, so just warn
    if (max_mem_mb > 100000) {  # > 100 GB
      warning("R is using a large amount of memory (>100 GB). Consider monitoring system memory.\n")
    }
    return(TRUE)
  })
}

# ============================================================================
# Setup and Configuration
# ============================================================================

# Check and set working directory
work_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_seurat.BT.analysis.sep30.nov21.18samples"
if (!dir.exists(work_dir)) {
  stop("Working directory does not exist: ", work_dir)
}
setwd(work_dir)
cat("Working directory set to:", getwd(), "\n")

# https://academic.oup.com/bioinformatics/article/41/6/btaf358/8171806
# a new package SEURAT INTEGRATE
# https://github.com/cbib/Seurat-Integrate/blob/main/vignettes/SeuratIntegrate.Rmd

# Point methods to your *path-based* env
UpdateEnvCache("scvi",
  conda.env = "/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/SeuratIntegrate",
  conda.env.is.path = TRUE
)
UpdateEnvCache("scanorama",
  conda.env = "/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/SeuratIntegrate",
  conda.env.is.path = TRUE
)

# Verify
getCache()

library(reticulate)

# Use your specific conda environment
use_condaenv("/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/SeuratIntegrate", required = TRUE)

# Verify it's active
py_config()

# ---- Load Individual Seurat Objects ----

# To increase the limit (e.g., to 70 GiB)
options(future.globals.maxSize = 700 * 1024^3)

# ---- Load Individual Seurat Objects ----
tryCatch({
  # Check if files exist before loading
  sample_files <- c("AD_AG_2658_seurat_obj.filtered.soupX.rds",
                    "AD_AG_2822_seurat_obj.filtered.soupX.rds",
                    "AD_AG_3184_seurat_obj.filtered.soupX.rds",
                    "AD_AG_3188_seurat_obj.filtered.soupX.rds",
                    "AD_CTL_2273_S2_seurat_obj.filtered.soupX.rds",
                    "AD_CTL_2354_S3_seurat_obj.filtered.soupX.rds",
                    "AD_CTL_2812_S1_seurat_obj.filtered.soupX.rds",
                    "Dy_AD_2328_S5_seurat_obj.filtered.soupX.rds",
                    "Dy_AD_2368_S6_seurat_obj.filtered.soupX.rds",
                    "Dy_AD_2660_S4_seurat_obj.filtered.soupX.rds",
                    "Dy_AG_2659_seurat_obj.filtered.soupX.rds",
                    "Dy_AG_2922_seurat_obj.filtered.soupX.rds",
                    "Dy_AG_3090_seurat_obj.filtered.soupX.rds",
                    "Dy_AG_3296_seurat_obj.filtered.soupX.rds",
                    "Healthy_AG_2394_seurat_obj.filtered.soupX.rds",
                    "Healthy_AG_2594_seurat_obj.filtered.soupX.rds",
                    "Healthy_AG_3041_seurat_obj.filtered.soupX.rds",
                    "Healthy_AG_3267_seurat_obj.filtered.soupX.rds")
  
  missing_files <- sample_files[!file.exists(sample_files)]
  if (length(missing_files) > 0) {
    stop("Missing files: ", paste(missing_files, collapse = ", "))
  }
  
  AD_AG_2658_seurat <- readRDS("AD_AG_2658_seurat_obj.filtered.soupX.rds")
  AD_AG_2822_seurat <- readRDS("AD_AG_2822_seurat_obj.filtered.soupX.rds")
  AD_AG_3184_seurat <- readRDS("AD_AG_3184_seurat_obj.filtered.soupX.rds")
  AD_AG_3188_seurat <- readRDS("AD_AG_3188_seurat_obj.filtered.soupX.rds")
  AD_CTL_2273_S2_seurat <- readRDS("AD_CTL_2273_S2_seurat_obj.filtered.soupX.rds")
  AD_CTL_2354_S3_seurat <- readRDS("AD_CTL_2354_S3_seurat_obj.filtered.soupX.rds")
  AD_CTL_2812_S1_seurat <- readRDS("AD_CTL_2812_S1_seurat_obj.filtered.soupX.rds")
  Dy_AD_2328_S5_seurat <- readRDS("Dy_AD_2328_S5_seurat_obj.filtered.soupX.rds")
  Dy_AD_2368_S6_seurat <- readRDS("Dy_AD_2368_S6_seurat_obj.filtered.soupX.rds")
  Dy_AD_2660_S4_seurat <- readRDS("Dy_AD_2660_S4_seurat_obj.filtered.soupX.rds")
  Dy_AG_2659_seurat <- readRDS("Dy_AG_2659_seurat_obj.filtered.soupX.rds")
  Dy_AG_2922_seurat <- readRDS("Dy_AG_2922_seurat_obj.filtered.soupX.rds")
  Dy_AG_3090_seurat <- readRDS("Dy_AG_3090_seurat_obj.filtered.soupX.rds")
  Dy_AG_3296_seurat <- readRDS("Dy_AG_3296_seurat_obj.filtered.soupX.rds")
  Healthy_AG_2394_seurat <- readRDS("Healthy_AG_2394_seurat_obj.filtered.soupX.rds")
  Healthy_AG_2594_seurat <- readRDS("Healthy_AG_2594_seurat_obj.filtered.soupX.rds")
  Healthy_AG_3041_seurat <- readRDS("Healthy_AG_3041_seurat_obj.filtered.soupX.rds")
  Healthy_AG_3267_seurat <- readRDS("Healthy_AG_3267_seurat_obj.filtered.soupX.rds")
  cat("All 18 Seurat objects loaded successfully\n")
  # Check memory after loading
  check_memory_usage(threshold = 0.85, stop_on_exceed = TRUE)
}, error = function(e) {
  cat("Error loading Seurat objects:", conditionMessage(e), "\n")
  stop("Failed to load input files. Exiting.\n")
})

# Create sample summary table
tryCatch({
  # Create output directory for plots early
  dir.create("SeuratIntegrate_integration_plots", showWarnings = FALSE)
  
  summary_df <- data.frame(
    Sample = c("AD_AG_2658", "AD_AG_2822", "AD_AG_3184", "AD_AG_3188", 
               "AD_CTL_2273_S2", "AD_CTL_2354_S3", "AD_CTL_2812_S1",
               "Dy_AD_2328_S5", "Dy_AD_2368_S6", "Dy_AD_2660_S4",
               "Dy_AG_2659", "Dy_AG_2922", "Dy_AG_3090", "Dy_AG_3296",
               "Healthy_AG_2394", "Healthy_AG_2594", "Healthy_AG_3041", "Healthy_AG_3267"),
    Genes = c(nrow(AD_AG_2658_seurat), nrow(AD_AG_2822_seurat), nrow(AD_AG_3184_seurat), nrow(AD_AG_3188_seurat),
              nrow(AD_CTL_2273_S2_seurat), nrow(AD_CTL_2354_S3_seurat), nrow(AD_CTL_2812_S1_seurat),
              nrow(Dy_AD_2328_S5_seurat), nrow(Dy_AD_2368_S6_seurat), nrow(Dy_AD_2660_S4_seurat),
              nrow(Dy_AG_2659_seurat), nrow(Dy_AG_2922_seurat), nrow(Dy_AG_3090_seurat), nrow(Dy_AG_3296_seurat),
              nrow(Healthy_AG_2394_seurat), nrow(Healthy_AG_2594_seurat), nrow(Healthy_AG_3041_seurat), nrow(Healthy_AG_3267_seurat)),
    Cells = c(ncol(AD_AG_2658_seurat), ncol(AD_AG_2822_seurat), ncol(AD_AG_3184_seurat), ncol(AD_AG_3188_seurat),
              ncol(AD_CTL_2273_S2_seurat), ncol(AD_CTL_2354_S3_seurat), ncol(AD_CTL_2812_S1_seurat),
              ncol(Dy_AD_2328_S5_seurat), ncol(Dy_AD_2368_S6_seurat), ncol(Dy_AD_2660_S4_seurat),
              ncol(Dy_AG_2659_seurat), ncol(Dy_AG_2922_seurat), ncol(Dy_AG_3090_seurat), ncol(Dy_AG_3296_seurat),
              ncol(Healthy_AG_2394_seurat), ncol(Healthy_AG_2594_seurat), ncol(Healthy_AG_3041_seurat), ncol(Healthy_AG_3267_seurat))
  )
  
  write.csv(summary_df, file = "SeuratIntegrate_integration_plots/sample_summary_table.csv", row.names = FALSE)
  # Note: We keep summary_df for now as it may be referenced later, but clean up if needed
  gc(verbose = FALSE)  # Garbage collection after saving summary
  cat("Sample summary table saved to SeuratIntegrate_integration_plots/sample_summary_table.csv\n")
}, error = function(e) {
  cat("Error saving summary table:", conditionMessage(e), "\n")
})

# ---- Setup conda environments ----
# (Already loaded at the top, no need to reload)

# ---- Add sample IDs and merge ----
tryCatch({
  AD_AG_2658_seurat$sample_id <- "AD_AG_2658"
  AD_AG_2822_seurat$sample_id <- "AD_AG_2822"
  AD_AG_3184_seurat$sample_id <- "AD_AG_3184"
  AD_AG_3188_seurat$sample_id <- "AD_AG_3188"
  AD_CTL_2273_S2_seurat$sample_id <- "AD_CTL_2273_S2"
  AD_CTL_2354_S3_seurat$sample_id <- "AD_CTL_2354_S3"
  AD_CTL_2812_S1_seurat$sample_id <- "AD_CTL_2812_S1"
  Dy_AD_2328_S5_seurat$sample_id <- "Dy_AD_2328_S5"
  Dy_AD_2368_S6_seurat$sample_id <- "Dy_AD_2368_S6"
  Dy_AD_2660_S4_seurat$sample_id <- "Dy_AD_2660_S4"
  Dy_AG_2659_seurat$sample_id <- "Dy_AG_2659"
  Dy_AG_2922_seurat$sample_id <- "Dy_AG_2922"
  Dy_AG_3090_seurat$sample_id <- "Dy_AG_3090"
  Dy_AG_3296_seurat$sample_id <- "Dy_AG_3296"
  Healthy_AG_2394_seurat$sample_id <- "Healthy_AG_2394"
  Healthy_AG_2594_seurat$sample_id <- "Healthy_AG_2594"
  Healthy_AG_3041_seurat$sample_id <- "Healthy_AG_3041"
  Healthy_AG_3267_seurat$sample_id <- "Healthy_AG_3267"
  cat("Sample IDs assigned successfully\n")
}, error = function(e) {
  cat("Error assigning sample IDs:", conditionMessage(e), "\n")
  stop("Failed to assign sample IDs. Exiting.\n")
})

# Merge all samples
tryCatch({
  merged_seurat <- merge(AD_AG_2658_seurat,
                         y = c(AD_AG_2822_seurat, AD_AG_3184_seurat, AD_AG_3188_seurat,
                               AD_CTL_2273_S2_seurat, AD_CTL_2354_S3_seurat, AD_CTL_2812_S1_seurat,
                               Dy_AD_2328_S5_seurat, Dy_AD_2368_S6_seurat, Dy_AD_2660_S4_seurat,
                               Dy_AG_2659_seurat, Dy_AG_2922_seurat, Dy_AG_3090_seurat, Dy_AG_3296_seurat,
                               Healthy_AG_2394_seurat, Healthy_AG_2594_seurat, 
                               Healthy_AG_3041_seurat, Healthy_AG_3267_seurat),
                         add.cell.ids = c("AD_AG_2658", "AD_AG_2822", "AD_AG_3184", "AD_AG_3188",
                                          "AD_CTL_2273_S2", "AD_CTL_2354_S3", "AD_CTL_2812_S1",
                                          "Dy_AD_2328_S5", "Dy_AD_2368_S6", "Dy_AD_2660_S4",
                                          "Dy_AG_2659", "Dy_AG_2922", "Dy_AG_3090", "Dy_AG_3296",
                                          "Healthy_AG_2394", "Healthy_AG_2594", "Healthy_AG_3041", "Healthy_AG_3267"))
  cat("Merged Seurat object created successfully\n")
  cat(sprintf("  Total genes: %d\n", nrow(merged_seurat)))
  cat(sprintf("  Total cells: %d\n", ncol(merged_seurat)))
  
  # Memory optimization: Remove individual seurat objects after merging
  rm(AD_AG_2658_seurat, AD_AG_2822_seurat, AD_AG_3184_seurat, AD_AG_3188_seurat,
     AD_CTL_2273_S2_seurat, AD_CTL_2354_S3_seurat, AD_CTL_2812_S1_seurat,
     Dy_AD_2328_S5_seurat, Dy_AD_2368_S6_seurat, Dy_AD_2660_S4_seurat,
     Dy_AG_2659_seurat, Dy_AG_2922_seurat, Dy_AG_3090_seurat, Dy_AG_3296_seurat,
     Healthy_AG_2394_seurat, Healthy_AG_2594_seurat, Healthy_AG_3041_seurat, Healthy_AG_3267_seurat)
  gc(verbose = FALSE)  # Garbage collection to free memory
  cat("  Individual Seurat objects removed from memory\n")
  
  # Check memory after merging
  check_memory_usage(threshold = 0.85, stop_on_exceed = TRUE)
}, error = function(e) {
  cat("Error merging Seurat objects:", conditionMessage(e), "\n")
  stop("Failed to merge objects. Exiting.\n")
})

# Join layers if already split (to avoid error on re-runs)
tryCatch({
  # Check ALL layers in RNA assay (not just counts)
  all_rna_layers <- tryCatch({
    Layers(merged_seurat[["RNA"]])
  }, error = function(e) {
    NULL
  })
  
  # Check if any layer names contain a period (indicating split layers like counts.1, data.2, etc.)
  if (!is.null(all_rna_layers) && any(grepl("\\.", all_rna_layers))) {
    cat(sprintf("RNA assay has split layers detected (%d total layers) - joining layers before re-splitting\n", length(all_rna_layers)))
    cat("Split layer names:", paste(all_rna_layers[grepl("\\.", all_rna_layers)], collapse = ", "), "\n")
    merged_seurat <- JoinLayers(merged_seurat, assay = "RNA")
    cat("RNA layers joined successfully\n")
  } else {
    cat("No split layers detected in RNA assay\n")
  }
  
  # Also check SCT assay if it exists
  if ("SCT" %in% Assays(merged_seurat)) {
    all_sct_layers <- tryCatch({
      Layers(merged_seurat[["SCT"]])
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(all_sct_layers) && any(grepl("\\.", all_sct_layers))) {
      cat(sprintf("SCT assay has split layers detected - joining layers\n"))
      merged_seurat <- JoinLayers(merged_seurat, assay = "SCT")
      cat("SCT layers joined successfully\n")
    }
  }
  
  gc(verbose = FALSE)
}, error = function(e) {
  cat("Warning: Could not check/join layers:", conditionMessage(e), "\n")
  cat("Attempting to proceed with split...\n")
})

# Define batch variable
batch.var <- "sample_id"

# ---- Preprocessing with Seurat ----
# Split by sample for independent normalization
tryCatch({
  # Check if layers are already split before attempting to split
  rna_layers <- tryCatch({
    Layers(merged_seurat[["RNA"]], search = "counts")
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(rna_layers) && length(rna_layers) > 1) {
    cat(sprintf("Layers are already split (%d layers detected) - skipping split step\n", length(rna_layers)))
    cat("Proceeding with SCTransform on already-split layers\n")
  } else {
    # Layers are not split, proceed with splitting
    merged_seurat <- split(merged_seurat, f = merged_seurat$sample_id)
    cat("Split merged object by sample_id\n")
  }
}, error = function(e) {
  error_msg <- conditionMessage(e)
  if (grepl("already split", error_msg, ignore.case = TRUE)) {
    cat("Warning: Layers are already split (detected from error message)\n")
    cat("Proceeding with SCTransform on existing split layers\n")
  } else {
    cat("Error splitting by sample_id:", error_msg, "\n")
    stop("Failed to split object. Exiting.\n")
  }
})

# SCTransform normalization
tryCatch({
  merged_seurat <- SCTransform(merged_seurat)
  gc(verbose = FALSE, full = TRUE)  # Full garbage collection after SCTransform
  cat("SCTransform normalization completed successfully\n")
}, error = function(e) {
  cat("Error in SCTransform:", conditionMessage(e), "\n")
  stop("SCTransform failed. Exiting.\n")
})

# PCA
tryCatch({
  merged_seurat <- RunPCA(merged_seurat, npcs = 50)
  gc(verbose = FALSE)  # Garbage collection after PCA
  cat("PCA completed successfully\n")
  # Check memory after PCA
  check_memory_usage(threshold = 0.85, stop_on_exceed = TRUE)
}, error = function(e) {
  cat("Error in RunPCA:", conditionMessage(e), "\n")
  stop("PCA failed. Exiting.\n")
})

# Find neighbors on unintegrated data (FIX: Added explicit reduction parameter)
tryCatch({
  merged_seurat <- FindNeighbors(merged_seurat, reduction = "pca", dims = 1:30, k.param = 20)
  gc(verbose = FALSE)  # Garbage collection after creating graph
  cat("FindNeighbors completed successfully\n")
}, error = function(e) {
  cat("Error in FindNeighbors:", conditionMessage(e), "\n")
  stop("FindNeighbors failed. Exiting.\n")
})

# Run UMAP on unintegrated data for initial visualization
# NOTE: This is kept for visualization purposes only. Method-specific UMAPs will be created
# after integration and post-processing in the sections below.
tryCatch({
  merged_seurat <- RunUMAP(merged_seurat, reduction = "pca", dims = 1:30, n.neighbors = 20)
  gc(verbose = FALSE)  # Garbage collection after UMAP
  cat("RunUMAP (unintegrated, for initial visualization) completed successfully\n")
}, error = function(e) {
  cat("Error in RunUMAP (unintegrated):", conditionMessage(e), "\n")
  stop("RunUMAP failed. Exiting.\n")
})

# FIX: Removed redundant initial clustering section (lines 416-451 in original)
# Clustering will be done properly in post-processing sections below with method-specific graphs
# This eliminates redundancy and ensures consistency

# Visualize unintegrated data (sample_id only, no clusters yet)
tryCatch({
  DimPlot(merged_seurat, group.by = "sample_id", reduction = "umap", label = FALSE)
  cat("Unintegrated data visualization plot created (sample_id only)\n")
}, error = function(e) {
  cat("Warning: Error creating visualization plot:", conditionMessage(e), "\n")
  # Don't stop, just warn - visualization is not critical
})

# Save unintegrated UMAP with sample_id (clusters will be added after post-processing)
tryCatch({
  p_sample_unintegrated <- DimPlot(merged_seurat, group.by = "sample_id", 
                                   reduction = "umap", label = FALSE, pt.size = 0.5) + 
    ggtitle("Unintegrated - UMAP colored by sample_id")
  ggsave("SeuratIntegrate_integration_plots/SeuratIntegrate_umap_unintegrated_samples.png", 
         plot = p_sample_unintegrated, width = 10, height = 8, dpi = 300)
  rm(p_sample_unintegrated)
  gc(verbose = FALSE)
  cat("Unintegrated UMAP with sample_id saved\n")
}, error = function(e) {
  cat("Error saving unintegrated UMAP with sample_id:", conditionMessage(e), "\n")
})

# ---- Integration with SeuratIntegrate ----
tryCatch({
  merged_seurat <- DoIntegrate(
    object = merged_seurat,
    
    # integrations
    SeuratIntegrate::HarmonyIntegration(orig = "pca", dims = 1:30),
    SeuratIntegrate::ScanoramaIntegration(),
    SeuratIntegrate::scVIIntegration(),
    
    # additional parameters
    use.hvg = TRUE,
    use.future = TRUE
  )
  cat("Integration with SeuratIntegrate completed successfully\n")
  cat(sprintf("  Total genes: %d\n", nrow(merged_seurat)))
  cat(sprintf("  Total cells: %d\n", ncol(merged_seurat)))
  cat("  Available reductions:", paste(names(merged_seurat@reductions), collapse = ", "), "\n")
  gc(verbose = FALSE, full = TRUE)  # Full garbage collection after integration
  # Check memory after integration
  check_memory_usage(threshold = 0.85, stop_on_exceed = TRUE)
}, error = function(e) {
  cat("Error in DoIntegrate:", conditionMessage(e), "\n")
  stop("Integration failed. Exiting.\n")
})


# Save the object immediately after integration
tryCatch({
  saveRDS(merged_seurat, file = "18samples_merged_after_DoIntegrate.v15.rds")
  gc(verbose = FALSE)  # Garbage collection after saving
  cat("Object saved after DoIntegrate: 18samples_merged_after_DoIntegrate.v15.rds\n")
}, error = function(e) {
  cat("Error saving intermediate object:", conditionMessage(e), "\n")
})


# ---- Post-processing: Unintegrated (dimension reduction) ----
# FIX: This section now properly creates clusters for unintegrated data
# using a dedicated graph, eliminating the redundancy from the removed initial clustering
try({
  DefaultAssay(merged_seurat) <- "SCT"
  merged_seurat <- FindNeighbors(merged_seurat, reduction = "pca", dims = 1:30, k.param = 20,
                                 return.neighbor = TRUE, graph.name = "knn.unintegrated")
  gc(verbose = FALSE)  # Cleanup after FindNeighbors
  merged_seurat <- SymmetrizeKnn(merged_seurat, graph.name = "knn.unintegrated")
  # Create clusters on unintegrated data for reference
  merged_seurat <- FindClusters(merged_seurat, graph.name = "knn.unintegrated_symmetric", 
                                resolution = 0.5, verbose = FALSE)
  merged_seurat$unintegrated_clusters_ref <- merged_seurat$seurat_clusters
  merged_seurat <- FindOptimalClusters(merged_seurat, graph.name = "knn.unintegrated_symmetric",
                                       cluster.name = "unintegrated_clusters_{metric}",
                                       cell.var = "unintegrated_clusters_ref",
                                       optimisation.metric = c("nmi", "ari"))
  gc(verbose = FALSE)  # Cleanup after clustering
  cat("Post-processing: Unintegrated completed successfully\n")
}, silent = FALSE)

# ---- Post-processing: Harmony (dimension reduction) ----
try({
  merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30,
                                 return.neighbor = TRUE, graph.name = "knn.harmony")
  gc(verbose = FALSE)  # Cleanup after FindNeighbors
  merged_seurat <- SymmetrizeKnn(merged_seurat, graph.name = "knn.harmony")
  # Create clusters on Harmony data for reference
  merged_seurat <- FindClusters(merged_seurat, graph.name = "knn.harmony_symmetric", 
                                resolution = 0.5, verbose = FALSE)
  merged_seurat$harmony_clusters_ref <- merged_seurat$seurat_clusters
  merged_seurat <- FindOptimalClusters(merged_seurat, graph.name = "knn.harmony_symmetric",
                                       cluster.name = "harmony_clusters_{metric}",
                                       cell.var = "harmony_clusters_ref",
                                       optimisation.metric = c("nmi", "ari"))
  gc(verbose = FALSE)  # Cleanup after clustering
  cat("Post-processing: Harmony completed successfully\n")
}, silent = FALSE)

# ---- Post-processing: Scanorama (dimension reduction) ----
try({
  merged_seurat <- FindNeighbors(merged_seurat, reduction = "integrated.scanorama", dims = 1:30,
                                 return.neighbor = TRUE, graph.name = "knn.scanorama")
  gc(verbose = FALSE)  # Cleanup after FindNeighbors
  merged_seurat <- SymmetrizeKnn(merged_seurat, graph.name = "knn.scanorama")
  # Create clusters on Scanorama data for reference
  merged_seurat <- FindClusters(merged_seurat, graph.name = "knn.scanorama_symmetric", 
                                resolution = 0.5, verbose = FALSE)
  merged_seurat$scanorama_clusters_ref <- merged_seurat$seurat_clusters
  merged_seurat <- FindOptimalClusters(merged_seurat, graph.name = "knn.scanorama_symmetric",
                                       cluster.name = "scanorama_clusters_{metric}",
                                       cell.var = "scanorama_clusters_ref",
                                       optimisation.metric = c("nmi", "ari"))
  gc(verbose = FALSE)  # Cleanup after clustering
  cat("Post-processing: Scanorama completed successfully\n")
}, silent = FALSE)

# ---- Post-processing: scVI (dimension reduction) ----
# FIX: Check actual scVI reduction name to handle case sensitivity
try({
  # Check available reductions to find scVI reduction (case-insensitive)
  available_reductions <- names(merged_seurat@reductions)
  scvi_reduction <- NULL
  
  # Try different possible names
  if ("integrated.scVI" %in% available_reductions) {
    scvi_reduction <- "integrated.scVI"
  } else if ("integrated.scvi" %in% available_reductions) {
    scvi_reduction <- "integrated.scvi"
  } else {
    # Try case-insensitive search
    scvi_matches <- grep("scvi", available_reductions, ignore.case = TRUE, value = TRUE)
    if (length(scvi_matches) > 0) {
      scvi_reduction <- scvi_matches[1]
      cat(sprintf("  Note: Using scVI reduction '%s' (found via case-insensitive search)\n", scvi_reduction))
    }
  }
  
  if (is.null(scvi_reduction)) {
    stop("scVI reduction not found. Available reductions: ", paste(available_reductions, collapse = ", "))
  }
  
  cat(sprintf("  Using scVI reduction: %s\n", scvi_reduction))
  
  merged_seurat <- FindNeighbors(merged_seurat, reduction = scvi_reduction, dims = 1:30,
                                 return.neighbor = TRUE, graph.name = "knn.scvi")
  gc(verbose = FALSE)  # Cleanup after FindNeighbors
  merged_seurat <- SymmetrizeKnn(merged_seurat, graph.name = "knn.scvi")
  # Create clusters on scVI data for reference
  merged_seurat <- FindClusters(merged_seurat, graph.name = "knn.scvi_symmetric", 
                                resolution = 0.5, verbose = FALSE)
  merged_seurat$scvi_clusters_ref <- merged_seurat$seurat_clusters
  merged_seurat <- FindOptimalClusters(merged_seurat, graph.name = "knn.scvi_symmetric",
                                       cluster.name = "scvi_clusters_{metric}",
                                       cell.var = "scvi_clusters_ref",
                                       optimisation.metric = c("nmi", "ari"))
  gc(verbose = FALSE)  # Cleanup after clustering
  cat("Post-processing: scVI completed successfully\n")
}, silent = FALSE)

# Create UMAPs for each integration method
tryCatch({
  merged_seurat <- RunUMAP(merged_seurat, dims = 1:30, reduction = "harmony",
                           reduction.name = "umap.harmony")
  gc(verbose = FALSE)  # Garbage collection after UMAP
  cat("UMAP for Harmony created successfully\n")
}, error = function(e) {
  cat("Error creating UMAP for Harmony:", conditionMessage(e), "\n")
})

tryCatch({
  merged_seurat <- RunUMAP(merged_seurat, dims = 1:30, reduction = "integrated.scanorama",
                           reduction.name = "umap.scanorama")
  gc(verbose = FALSE)  # Garbage collection after UMAP
  cat("UMAP for Scanorama created successfully\n")
}, error = function(e) {
  cat("Error creating UMAP for Scanorama:", conditionMessage(e), "\n")
})

# FIX: Use the detected scVI reduction name for UMAP creation
tryCatch({
  # Use the same scvi_reduction variable if available, otherwise detect again
  if (!exists("scvi_reduction")) {
    available_reductions <- names(merged_seurat@reductions)
    if ("integrated.scVI" %in% available_reductions) {
      scvi_reduction <- "integrated.scVI"
    } else if ("integrated.scvi" %in% available_reductions) {
      scvi_reduction <- "integrated.scvi"
    } else {
      scvi_matches <- grep("scvi", available_reductions, ignore.case = TRUE, value = TRUE)
      if (length(scvi_matches) > 0) {
        scvi_reduction <- scvi_matches[1]
      } else {
        stop("scVI reduction not found for UMAP creation")
      }
    }
  }
  
  merged_seurat <- RunUMAP(merged_seurat, dims = 1:30, reduction = scvi_reduction,
                           reduction.name = "umap.scvi")
  gc(verbose = FALSE)  # Garbage collection after UMAP
  cat(sprintf("UMAP for scVI created successfully (using reduction: %s)\n", scvi_reduction))
}, error = function(e) {
  cat("Error creating UMAP for scVI:", conditionMessage(e), "\n")
})



# ---- Save individual UMAP plots ----
# Note: plots directory already created earlier when saving summary table

# Unintegrated
tryCatch({
  p1 <- DimPlot(merged_seurat, group.by = batch.var, reduction = "umap")
  # ggsave("SeuratIntegrate_integration_plots/umap_unintegrated.pdf", plot = p1, width = 8, height = 6)  # PDF disabled - saving PNG only
  ggsave("SeuratIntegrate_integration_plots/umap_unintegrated.png", plot = p1, width = 8, height = 6, dpi = 300)
  rm(p1)
  gc(verbose = FALSE)
  cat("Unintegrated UMAP plot saved\n")
}, error = function(e) {
  cat("Error saving unintegrated UMAP plot:", conditionMessage(e), "\n")
})

# Unintegrated - clusters (with method-specific reference clusters)
tryCatch({
  p1_clusters <- DimPlot(merged_seurat, group.by = "unintegrated_clusters_ref", 
                        reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.5) + 
    ggtitle("Unintegrated - UMAP colored by unintegrated clusters")
  ggsave("SeuratIntegrate_integration_plots/SeuratIntegrate_umap_unintegrated_clusters.png", 
         plot = p1_clusters, width = 10, height = 8, dpi = 300)
  rm(p1_clusters)
  gc(verbose = FALSE)
  cat("Unintegrated UMAP with method-specific clusters saved\n")
}, error = function(e) {
  cat("Error saving unintegrated UMAP with method-specific clusters:", conditionMessage(e), "\n")
})


# Harmony
tryCatch({
  p3 <- DimPlot(merged_seurat, group.by = batch.var, reduction = "umap.harmony")
  # ggsave("SeuratIntegrate_integration_plots/umap_harmony.pdf", plot = p3, width = 8, height = 6)  # PDF disabled - saving PNG only
  ggsave("SeuratIntegrate_integration_plots/umap_harmony.png", plot = p3, width = 8, height = 6, dpi = 300)
  rm(p3)
  gc(verbose = FALSE)
  cat("Harmony UMAP plot saved\n")
}, error = function(e) {
  cat("Error saving Harmony UMAP plot:", conditionMessage(e), "\n")
})

# Harmony - clusters
tryCatch({
  p3_clusters <- DimPlot(merged_seurat, group.by = "harmony_clusters_ref", 
                        reduction = "umap.harmony", label = TRUE, label.size = 4, pt.size = 0.5) + 
    ggtitle("Harmony - UMAP colored by Harmony clusters")
  ggsave("SeuratIntegrate_integration_plots/SeuratIntegrate_umap_harmony_clusters.png", 
         plot = p3_clusters, width = 10, height = 8, dpi = 300)
  rm(p3_clusters)
  gc(verbose = FALSE)
  cat("Harmony UMAP with clusters saved\n")
}, error = function(e) {
  cat("Error saving Harmony UMAP with clusters:", conditionMessage(e), "\n")
})

# Scanorama
tryCatch({
  p4 <- DimPlot(merged_seurat, group.by = batch.var, reduction = "umap.scanorama")
  # ggsave("SeuratIntegrate_integration_plots/umap_scanorama.pdf", plot = p4, width = 8, height = 6)  # PDF disabled - saving PNG only
  ggsave("SeuratIntegrate_integration_plots/umap_scanorama.png", plot = p4, width = 8, height = 6, dpi = 300)
  rm(p4)
  gc(verbose = FALSE)
  cat("Scanorama UMAP plot saved\n")
}, error = function(e) {
  cat("Error saving Scanorama UMAP plot:", conditionMessage(e), "\n")
})

# Scanorama - clusters
tryCatch({
  p4_clusters <- DimPlot(merged_seurat, group.by = "scanorama_clusters_ref", 
                         reduction = "umap.scanorama", label = TRUE, label.size = 4, pt.size = 0.5) + 
    ggtitle("Scanorama - UMAP colored by Scanorama clusters")
  ggsave("SeuratIntegrate_integration_plots/SeuratIntegrate_umap_scanorama_clusters.png", 
         plot = p4_clusters, width = 10, height = 8, dpi = 300)
  rm(p4_clusters)
  gc(verbose = FALSE)
  cat("Scanorama UMAP with clusters saved\n")
}, error = function(e) {
  cat("Error saving Scanorama UMAP with clusters:", conditionMessage(e), "\n")
})

# scVI
tryCatch({
  p5 <- DimPlot(merged_seurat, group.by = batch.var, reduction = "umap.scvi")
  # ggsave("SeuratIntegrate_integration_plots/umap_scvi.pdf", plot = p5, width = 8, height = 6)  # PDF disabled - saving PNG only
  ggsave("SeuratIntegrate_integration_plots/umap_scvi.png", plot = p5, width = 8, height = 6, dpi = 300)
  rm(p5)
  gc(verbose = FALSE)
  cat("scVI UMAP plot saved\n")
}, error = function(e) {
  cat("Error saving scVI UMAP plot:", conditionMessage(e), "\n")
})

# scVI - clusters
tryCatch({
  p5_clusters <- DimPlot(merged_seurat, group.by = "scvi_clusters_ref", 
                        reduction = "umap.scvi", label = TRUE, label.size = 4, pt.size = 0.5) + 
    ggtitle("scVI - UMAP colored by scVI clusters")
  ggsave("SeuratIntegrate_integration_plots/SeuratIntegrate_umap_scvi_clusters.png", 
         plot = p5_clusters, width = 10, height = 8, dpi = 300)
  rm(p5_clusters)
  gc(verbose = FALSE)
  cat("scVI UMAP with clusters saved\n")
}, error = function(e) {
  cat("Error saving scVI UMAP with clusters:", conditionMessage(e), "\n")
})



# Create a combined plot with all methods (excluding Combat, BBKNN, trVAE)
tryCatch({
  library(patchwork)
  # Recreate plots for combined visualization
  p1_combined <- DimPlot(merged_seurat, group.by = batch.var, reduction = "umap")
  p3_combined <- DimPlot(merged_seurat, group.by = batch.var, reduction = "umap.harmony")
  p4_combined <- DimPlot(merged_seurat, group.by = batch.var, reduction = "umap.scanorama")
  p5_combined <- DimPlot(merged_seurat, group.by = batch.var, reduction = "umap.scvi")
  
  combined_plot <- p1_combined + p3_combined + p4_combined + p5_combined + 
    plot_layout(ncol = 2, guides = 'collect')
  
  # ggsave("SeuratIntegrate_integration_plots/umap_all_methods_combined.pdf", 
  #        plot = combined_plot, width = 12, height = 12)  # PDF disabled - saving PNG only
  ggsave("SeuratIntegrate_integration_plots/umap_all_methods_combined.png", 
         plot = combined_plot, width = 12, height = 12, dpi = 300)
  rm(p1_combined, p3_combined, p4_combined, p5_combined, combined_plot)
  gc(verbose = FALSE)
  cat("Combined UMAP plot saved\n")
}, error = function(e) {
  cat("Error creating combined plot:", conditionMessage(e), "\n")
})

# Create a combined plot with clusters for all methods
tryCatch({
  library(patchwork)
  # Recreate plots for combined cluster visualization using method-specific clusters
  p1_clusters_combined <- DimPlot(merged_seurat, group.by = "unintegrated_clusters_ref", 
                                  reduction = "umap", label = TRUE, label.size = 3) + 
    ggtitle("Unintegrated")
  p3_clusters_combined <- DimPlot(merged_seurat, group.by = "harmony_clusters_ref", 
                                  reduction = "umap.harmony", label = TRUE, label.size = 3) + 
    ggtitle("Harmony")
  p4_clusters_combined <- DimPlot(merged_seurat, group.by = "scanorama_clusters_ref", 
                                  reduction = "umap.scanorama", label = TRUE, label.size = 3) + 
    ggtitle("Scanorama")
  p5_clusters_combined <- DimPlot(merged_seurat, group.by = "scvi_clusters_ref", 
                                  reduction = "umap.scvi", label = TRUE, label.size = 3) + 
    ggtitle("scVI")
  
  combined_clusters_plot <- p1_clusters_combined + p3_clusters_combined + 
    p4_clusters_combined + p5_clusters_combined + 
    plot_layout(ncol = 2, guides = 'collect')
  
  ggsave("SeuratIntegrate_integration_plots/SeuratIntegrate_umap_all_methods_clusters_combined.png", 
         plot = combined_clusters_plot, width = 12, height = 12, dpi = 300)
  rm(p1_clusters_combined, p3_clusters_combined, p4_clusters_combined, 
     p5_clusters_combined, combined_clusters_plot)
  gc(verbose = FALSE)
  cat("Combined UMAP plot with clusters saved\n")
}, error = function(e) {
  cat("Error creating combined cluster plot:", conditionMessage(e), "\n")
})

cat("All plots saved to 'SeuratIntegrate_integration_plots/' directory\n")

# Memory optimization: Final garbage collection before saving
gc(verbose = FALSE)
# Final memory check before saving
check_memory_usage(threshold = 0.85, stop_on_exceed = TRUE)

# ---- Save the object ----
tryCatch({
  saveRDS(merged_seurat, file = "18_samples_SeuratIntegrate.filtered.soupX.DoIntegrate.v15.end.rds")
  gc(verbose = FALSE)  # Garbage collection after saving
  cat("\n=== Analysis Complete ===\n")
  cat("Seurat object saved: 18_samples_SeuratIntegrate.filtered.soupX.DoIntegrate.v15.end.rds\n")
  cat("Individual UMAP plots: SeuratIntegrate_integration_plots/ directory (PNG only)\n")
}, error = function(e) {
  cat("\n=== ERROR SAVING FINAL OBJECT ===\n")
  cat("Error:", conditionMessage(e), "\n")
  cat("Attempting to save with alternative name...\n")
  tryCatch({
    saveRDS(merged_seurat, file = "18_samples_SeuratIntegrate.filtered.soupX.DoIntegrate.v15.end.backup.rds")
    cat("Saved with backup name: 18_samples_SeuratIntegrate.filtered.soupX.DoIntegrate.v15.end.backup.rds\n")
  }, error = function(e2) {
    cat("Failed to save even with backup name:", conditionMessage(e2), "\n")
  })
})

# ============================================================================
# Clustering Comparison Metrics
# ============================================================================
cat("\n")
cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("CLUSTERING COMPARISON METRICS\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

# Create output directory for metrics
OUTPUT_DIR <- "clustering_comparison_metrics"
dir.create(OUTPUT_DIR, showWarnings = FALSE)
cat(sprintf("Metrics output directory: %s/\n\n", OUTPUT_DIR))

# Define cluster column names for each method (using _ref columns for comparison)
cluster_columns <- list(
  unintegrated = "unintegrated_clusters_ref",
  harmony = "harmony_clusters_ref",
  scanorama = "scanorama_clusters_ref",
  scvi = "scvi_clusters_ref"
)

# Check which columns exist
available_methods <- list()
for (method in names(cluster_columns)) {
  col_name <- cluster_columns[[method]]
  if (col_name %in% colnames(merged_seurat@meta.data)) {
    available_methods[[method]] <- col_name
    cat(sprintf("✓ Found: %s (%s)\n", method, col_name))
  } else {
    cat(sprintf("✗ Missing: %s (%s)\n", method, col_name))
  }
}

if (length(available_methods) < 2) {
  cat("Warning: Need at least 2 clustering methods to compare!\n")
  cat("Skipping clustering comparison metrics.\n\n")
} else {
  cat(sprintf("\nComparing %d methods\n\n", length(available_methods)))
  
  # ============================================================================
  # 1. Basic Cluster Statistics
  # ============================================================================
  cat("=", rep("=", 70), "\n", sep = "")
  cat("1. BASIC CLUSTER STATISTICS\n")
  cat("=", rep("=", 70), "\n", sep = "")
  
  cluster_stats <- data.frame(
    Method = character(),
    N_Clusters = numeric(),
    N_Cells = numeric(),
    Min_Cluster_Size = numeric(),
    Max_Cluster_Size = numeric(),
    Median_Cluster_Size = numeric(),
    Mean_Cluster_Size = numeric(),
    SD_Cluster_Size = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (method in names(available_methods)) {
    col_name <- available_methods[[method]]
    clusters <- merged_seurat@meta.data[[col_name]]
    
    cluster_sizes <- table(clusters)
    n_clust <- length(unique(clusters))
    n_cells <- length(clusters)
    
    stats <- data.frame(
      Method = method,
      N_Clusters = n_clust,
      N_Cells = n_cells,
      Min_Cluster_Size = min(cluster_sizes),
      Max_Cluster_Size = max(cluster_sizes),
      Median_Cluster_Size = median(cluster_sizes),
      Mean_Cluster_Size = round(mean(cluster_sizes), 1),
      SD_Cluster_Size = round(sd(cluster_sizes), 1)
    )
    
    cluster_stats <- rbind(cluster_stats, stats)
    
    cat(sprintf("\n%s:\n", method))
    cat(sprintf("  Clusters: %d\n", n_clust))
    cat(sprintf("  Cells: %d\n", n_cells))
    cat(sprintf("  Cluster size range: %d - %d\n", min(cluster_sizes), max(cluster_sizes)))
    cat(sprintf("  Mean cluster size: %.1f ± %.1f\n", mean(cluster_sizes), sd(cluster_sizes)))
  }
  
  print(cluster_stats)
  write.csv(cluster_stats, 
            file = paste0(OUTPUT_DIR, "/cluster_statistics.csv"),
            row.names = FALSE)
  
  # ============================================================================
  # 2. Adjusted Rand Index (ARI) - Pairwise Comparison
  # ============================================================================
  cat("\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("2. ADJUSTED RAND INDEX (ARI) - Pairwise Comparison\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("ARI ranges from -1 to 1, where 1 = perfect agreement, 0 = random\n\n")
  
  n_methods <- length(available_methods)
  ari_matrix <- matrix(NA, nrow = n_methods, ncol = n_methods)
  rownames(ari_matrix) <- names(available_methods)
  colnames(ari_matrix) <- names(available_methods)
  
  for (i in 1:n_methods) {
    method1 <- names(available_methods)[i]
    col1 <- available_methods[[method1]]
    clusters1 <- merged_seurat@meta.data[[col1]]
    
    for (j in 1:n_methods) {
      method2 <- names(available_methods)[j]
      col2 <- available_methods[[method2]]
      clusters2 <- merged_seurat@meta.data[[col2]]
      
      ari_value <- ARI(clusters1, clusters2)
      ari_matrix[i, j] <- ari_value
      
      if (i != j) {
        cat(sprintf("%s vs %s: ARI = %.4f\n", method1, method2, ari_value))
      }
    }
  }
  
  # Set diagonal to 1
  diag(ari_matrix) <- 1.0
  
  print(ari_matrix)
  write.csv(ari_matrix, 
            file = paste0(OUTPUT_DIR, "/ari_pairwise_comparison.csv"),
            row.names = TRUE)
  
  # ============================================================================
  # 3. Normalized Mutual Information (NMI) - Pairwise Comparison
  # ============================================================================
  cat("\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("3. NORMALIZED MUTUAL INFORMATION (NMI) - Pairwise Comparison\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("NMI ranges from 0 to 1, where 1 = perfect agreement, 0 = independent\n\n")
  
  nmi_matrix <- matrix(NA, nrow = n_methods, ncol = n_methods)
  rownames(nmi_matrix) <- names(available_methods)
  colnames(nmi_matrix) <- names(available_methods)
  
  for (i in 1:n_methods) {
    method1 <- names(available_methods)[i]
    col1 <- available_methods[[method1]]
    clusters1 <- merged_seurat@meta.data[[col1]]
    
    for (j in 1:n_methods) {
      method2 <- names(available_methods)[j]
      col2 <- available_methods[[method2]]
      clusters2 <- merged_seurat@meta.data[[col2]]
      
      nmi_value <- NMI(clusters1, clusters2)
      nmi_matrix[i, j] <- nmi_value
      
      if (i != j) {
        cat(sprintf("%s vs %s: NMI = %.4f\n", method1, method2, nmi_value))
      }
    }
  }
  
  # Set diagonal to 1
  diag(nmi_matrix) <- 1.0
  
  print(nmi_matrix)
  write.csv(nmi_matrix, 
            file = paste0(OUTPUT_DIR, "/nmi_pairwise_comparison.csv"),
            row.names = TRUE)
  
  # ============================================================================
  # 4. Silhouette Score (requires reduction coordinates)
  # ============================================================================
  cat("\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("4. SILHOUETTE SCORE (Cluster Quality)\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("Silhouette ranges from -1 to 1, higher = better separation\n\n")
  
  silhouette_scores <- data.frame(
    Method = character(),
    Reduction = character(),
    Mean_Silhouette = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Define reduction names for each method
  reduction_map <- list(
    unintegrated = "umap",
    harmony = "umap.harmony",
    scanorama = "umap.scanorama",
    scvi = "umap.scvi"
  )
  
  for (method in names(available_methods)) {
    col_name <- available_methods[[method]]
    reduction_name <- reduction_map[[method]]
    
    if (reduction_name %in% names(merged_seurat@reductions)) {
      # Get coordinates
      coords <- merged_seurat@reductions[[reduction_name]]@cell.embeddings
      clusters <- as.numeric(as.factor(merged_seurat@meta.data[[col_name]]))
      
      # Calculate distance matrix (sample if too large)
      n_cells <- nrow(coords)
      if (n_cells > 10000) {
        # Sample for large datasets
        sample_idx <- sample(1:n_cells, 10000)
        coords_sample <- coords[sample_idx, ]
        clusters_sample <- clusters[sample_idx]
        cat(sprintf("%s: Sampling %d cells for silhouette calculation\n", method, length(sample_idx)))
      } else {
        coords_sample <- coords
        clusters_sample <- clusters
      }
      
      # Calculate distance matrix
      dist_matrix <- dist(coords_sample)
      
      # Calculate silhouette
      sil <- silhouette(clusters_sample, dist_matrix)
      mean_sil <- mean(sil[, "sil_width"])
      
      silhouette_scores <- rbind(silhouette_scores, 
                                 data.frame(
                                   Method = method,
                                   Reduction = reduction_name,
                                   Mean_Silhouette = round(mean_sil, 4)
                                 ))
      
      cat(sprintf("%s (%s): Mean Silhouette = %.4f\n", method, reduction_name, mean_sil))
    } else {
      cat(sprintf("%s: Reduction '%s' not found, skipping\n", method, reduction_name))
    }
  }
  
  print(silhouette_scores)
  write.csv(silhouette_scores, 
            file = paste0(OUTPUT_DIR, "/silhouette_scores.csv"),
            row.names = FALSE)
  
  # ============================================================================
  # 5. Integration Quality Metrics (LISI - Local Inverse Simpson Index)
  # ============================================================================
  cat("\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("5. INTEGRATION QUALITY (LISI - Local Inverse Simpson Index)\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("Lower LISI = better integration (cells from different samples mix well)\n\n")
  
  lisi_scores <- data.frame(
    Method = character(),
    Reduction = character(),
    Mean_LISI = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (method in names(available_methods)) {
    reduction_name <- reduction_map[[method]]
    
    if (reduction_name %in% names(merged_seurat@reductions)) {
      coords <- merged_seurat@reductions[[reduction_name]]@cell.embeddings
      
      # Sample if too large
      n_cells <- nrow(coords)
      if (n_cells > 10000) {
        sample_idx <- sample(1:n_cells, 10000)
        coords_sample <- coords[sample_idx, ]
        meta_sample <- merged_seurat@meta.data[sample_idx, batch.var, drop = FALSE]
        cat(sprintf("%s: Sampling %d cells for LISI calculation\n", method, length(sample_idx)))
      } else {
        coords_sample <- coords
        meta_sample <- merged_seurat@meta.data[, batch.var, drop = FALSE]
      }
      
      # Calculate LISI
      lisi_result <- compute_lisi(coords_sample, meta_sample, batch.var)
      mean_lisi <- mean(lisi_result[[batch.var]])
      
      lisi_scores <- rbind(lisi_scores,
                           data.frame(
                             Method = method,
                             Reduction = reduction_name,
                             Mean_LISI = round(mean_lisi, 4)
                           ))
      
      cat(sprintf("%s (%s): Mean LISI = %.4f\n", method, reduction_name, mean_lisi))
    }
  }
  
  print(lisi_scores)
  write.csv(lisi_scores, 
            file = paste0(OUTPUT_DIR, "/lisi_scores.csv"),
            row.names = FALSE)
  
  # ============================================================================
  # 6. Cluster Purity (how well clusters separate by sample)
  # ============================================================================
  cat("\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("6. CLUSTER PURITY (Sample Separation)\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("Higher purity = clusters are more sample-specific (less integration)\n")
  cat("Lower purity = clusters mix samples well (better integration)\n\n")
  
  purity_scores <- data.frame(
    Method = character(),
    Mean_Purity = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (method in names(available_methods)) {
    col_name <- available_methods[[method]]
    clusters <- merged_seurat@meta.data[[col_name]]
    samples <- merged_seurat@meta.data[[batch.var]]
    
    # Calculate purity for each cluster
    cluster_purities <- numeric()
    unique_clusters <- unique(clusters)
    
    for (clust in unique_clusters) {
      clust_cells <- which(clusters == clust)
      clust_samples <- samples[clust_cells]
      
      # Purity = proportion of most common sample in cluster
      sample_counts <- table(clust_samples)
      max_count <- max(sample_counts)
      purity <- max_count / length(clust_cells)
      cluster_purities <- c(cluster_purities, purity)
    }
    
    mean_purity <- mean(cluster_purities)
    
    purity_scores <- rbind(purity_scores,
                          data.frame(
                            Method = method,
                            Mean_Purity = round(mean_purity, 4)
                          ))
    
    cat(sprintf("%s: Mean Cluster Purity = %.4f\n", method, mean_purity))
  }
  
  print(purity_scores)
  write.csv(purity_scores, 
            file = paste0(OUTPUT_DIR, "/cluster_purity.csv"),
            row.names = FALSE)
  
  # ============================================================================
  # 7. Summary Table
  # ============================================================================
  cat("\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("7. SUMMARY TABLE\n")
  cat("=", rep("=", 70), "\n", sep = "")
  
  # Combine all metrics
  summary_table <- cluster_stats %>%
    left_join(silhouette_scores %>% select(Method, Mean_Silhouette), by = "Method") %>%
    left_join(purity_scores, by = "Method") %>%
    left_join(lisi_scores %>% select(Method, Mean_LISI), by = "Method")
  
  print(summary_table)
  write.csv(summary_table, 
            file = paste0(OUTPUT_DIR, "/comparison_summary.csv"),
            row.names = FALSE)
  
  # ============================================================================
  # 8. Batch Correction Quality Scores (AddScoreRegressPC & AddScoreDensityPC)
  # ============================================================================
  cat("\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("8. BATCH CORRECTION QUALITY SCORES\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("Adding regression and density scores for integration quality assessment\n\n")
  
  # Define integration methods and their corresponding reductions
  integration_methods <- list(
    unintegrated = list(name = "unintegrated", reduction = "pca"),
    harmony = list(name = "harmony", reduction = "harmony"),
    scanorama = list(name = "scanorama", reduction = "integrated.scanorama"),
    scvi = list(name = "scvi", reduction = if(exists("scvi_reduction")) scvi_reduction else "integrated.scVI")
  )
  
  # Filter to only methods that have reductions available
  available_integrations <- list()
  for (method in names(integration_methods)) {
    reduction_name <- integration_methods[[method]]$reduction
    if (reduction_name %in% names(merged_seurat@reductions)) {
      available_integrations[[method]] <- integration_methods[[method]]
      cat(sprintf("✓ Will add scores for: %s (reduction: %s)\n", 
                  method, reduction_name))
    } else {
      cat(sprintf("✗ Skipping: %s (reduction '%s' not found)\n", 
                  method, reduction_name))
    }
  }
  
  if (length(available_integrations) > 0) {
    # Add regression scores
    cat("\nAdding regression scores (AddScoreRegressPC)...\n")
    tryCatch({
      for (method in names(available_integrations)) {
        int_name <- available_integrations[[method]]$name
        reduction_name <- available_integrations[[method]]$reduction
        
        merged_seurat <- AddScoreRegressPC(merged_seurat, 
                                          integration = int_name,
                                          batch.var = batch.var, 
                                          reduction = reduction_name)
        cat(sprintf("  ✓ Added regression score for: %s\n", method))
      }
      cat("Regression scores added successfully!\n")
    }, error = function(e) {
      cat(sprintf("Error adding regression scores: %s\n", conditionMessage(e)))
    })
    
    # Add density scores
    cat("\nAdding density scores (AddScoreDensityPC)...\n")
    tryCatch({
      for (method in names(available_integrations)) {
        int_name <- available_integrations[[method]]$name
        reduction_name <- available_integrations[[method]]$reduction
        
        merged_seurat <- AddScoreDensityPC(merged_seurat, 
                                          integration = int_name,
                                          batch.var = batch.var, 
                                          reduction = reduction_name)
        cat(sprintf("  ✓ Added density score for: %s\n", method))
      }
      cat("Density scores added successfully!\n")
    }, error = function(e) {
      cat(sprintf("Error adding density scores: %s\n", conditionMessage(e)))
    })
    
    # Extract and save the scores
    cat("\nExtracting batch correction scores...\n")
    score_columns <- grep("ScoreRegressPC|ScoreDensityPC", colnames(merged_seurat@meta.data), value = TRUE)
    
    if (length(score_columns) > 0) {
      cat(sprintf("Found %d score columns:\n", length(score_columns)))
      for (col in score_columns) {
        cat(sprintf("  - %s\n", col))
      }
      
      # Create summary table of scores
      score_summary <- data.frame(
        Score_Column = score_columns,
        Mean_Score = numeric(length(score_columns)),
        Median_Score = numeric(length(score_columns)),
        SD_Score = numeric(length(score_columns)),
        stringsAsFactors = FALSE
      )
      
      for (i in seq_along(score_columns)) {
        col_name <- score_columns[i]
        scores <- merged_seurat@meta.data[[col_name]]
        score_summary$Mean_Score[i] <- round(mean(scores, na.rm = TRUE), 4)
        score_summary$Median_Score[i] <- round(median(scores, na.rm = TRUE), 4)
        score_summary$SD_Score[i] <- round(sd(scores, na.rm = TRUE), 4)
      }
      
      print(score_summary)
      write.csv(score_summary, 
                file = paste0(OUTPUT_DIR, "/batch_correction_scores_summary.csv"),
                row.names = FALSE)
      
      cat("\nBatch correction scores summary saved!\n")
    } else {
      cat("Warning: No score columns found in metadata\n")
    }
  } else {
    cat("No integration methods with available reductions found.\n")
  }
  
  cat("\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("All clustering comparison metrics saved to:", OUTPUT_DIR, "\n")
  cat("=", rep("=", 70), "\n", sep = "")
  cat("\n")
  
  # Garbage collection after metrics calculation
  gc(verbose = FALSE)
}

# ============================================================================
# Final Summary
# ============================================================================
script_end_time <- Sys.time()
script_duration <- difftime(script_end_time, script_start_time, units = "mins")

cat("\n")
cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("SCRIPT EXECUTION SUMMARY\n")
cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("Start time:", format(script_start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("End time:  ", format(script_end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Duration: ", round(as.numeric(script_duration), 2), "minutes\n")
cat("\n")


# Save full session info to file
tryCatch({
  sink("SeuratIntegrate_integration_plots/session_info.txt")
  print(sessionInfo())
  sink()
  gc(verbose = FALSE)  # Cleanup after sessionInfo
  cat("Full session info saved to SeuratIntegrate_integration_plots/session_info.txt\n")
}, error = function(e) {
  cat("Error saving session info:", conditionMessage(e), "\n")
})

# Final object summary
if (exists("merged_seurat")) {
  cat("Final Seurat Object Summary:\n")
  cat(sprintf("  Genes: %d\n", nrow(merged_seurat)))
  cat(sprintf("  Cells: %d\n", ncol(merged_seurat)))
  cat("  Available reductions:", paste(names(merged_seurat@reductions), collapse = ", "), "\n")
  cat("\n")
}

cat("Output Files:\n")
cat("  Main object: 18_samples_SeuratIntegrate.filtered.soupX.DoIntegrate.v15.end.rds\n")
cat("  Intermediate: 18samples_merged_after_DoIntegrate.v15.rds\n")
cat("  UMAP plots: SeuratIntegrate_integration_plots/ directory (PNG only)\n")
cat("  Clustering metrics: clustering_comparison_metrics/ directory\n")
cat("\n")


cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("Analysis complete!\n")
cat(paste0(rep("=", 70), collapse = ""), "\n")

