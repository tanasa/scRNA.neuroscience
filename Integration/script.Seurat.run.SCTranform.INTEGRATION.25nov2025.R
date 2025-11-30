# ============================================================================
# Script: integration.18samples.SeuratStandard.vX.scTranform.v5.version.claude.v15.R
# Description: Seurat v5 SCTransform integration (OFFICIAL WORKFLOW - v15)
# v15 Changes:
#   - Added JoinLayers section at the end (after line 1567) to join RNA layers
#   - Added save of joined layers object: 18.samples.integrate.SeuratStandard.SCTransform.v5.OFFICIAL.joinedLayers.rds
#   - Fixed log filename from 16samples to 18samples
#   - Fixed output filename from 16.samples to 18.samples
# v14 Features (maintained):
#   - Added comprehensive configuration section (from v6)
#   - Increased from 16 to 18 samples (added AD_AG_2658, Dy_AG_3090)
#   - All plot parameters now configurable at top
#   - JoinLayers section commented out (uncomment if needed for specific analyses)
# v13 Features (maintained):
#   - Improved Traditional Harmony workflow with explicit split layer handling
#   - Uses separate PCA reduction (pca.traditional) to preserve SCTransform PCA
#   - Includes split plots by sample for traditional Harmony
#   - Better error handling and logging
#   - Fixed JoinLayers: Only joins RNA layers (not SCT) to prevent errors (now commented out)
#   - Removed SCT layer joining before scVI (not needed)
# Based on: https://satijalab.org/seurat/articles/seurat5_integration
# Features:
#   - Unintegrated UMAP/clustering analysis
#   - Unintegrated plots (by sample, by clusters at multiple resolutions)
#   - Layer re-split checks (skips if already split)
#   - Cluster resolutions: 0.1, 0.2, 0.3, 0.4, 0.5, 1.0 (removed 2.0)
#   - QC plots (nFeature_RNA, nCount_RNA, percent.mt) with violin plots
#   - Combined plot by sample (all methods comparison)
#   - Resolution comparison grids for all methods
#   - Split plots by sample for better visualization
# ============================================================================

# ============================================================================
# Configuration Section (Comprehensive - all parameters in one place)
# ============================================================================
# Memory and compute settings
MEMORY_THRESHOLD <- 0.85              # Stop script if memory usage exceeds 85%
FUTURE_GLOBALS_MAX_SIZE <- 600 * 1024^3  # 600 GB

# Dimensionality settings
PCA_DIMS <- 1:30     # Dimensions to use for PCA-based operations (RPCA, Harmony)
SCVI_DIMS <- 1:10    # Dimensions to use for scVI (typically fewer dims - scVI learns compressed latent)
PCA_NDIMS <- 50      # Number of PCs to compute initially

# Clustering settings
CLUSTER_RESOLUTIONS <- c(0.1, 0.2, 0.3, 0.4, 0.5, 1.0)  # Clustering resolutions to test
DEFAULT_RESOLUTION <- 0.5  # Default resolution for combined plots and comparisons

# Plot settings
PLOT_WIDTH <- 10          # Default plot width in inches
PLOT_HEIGHT <- 8          # Default plot height in inches
PLOT_DPI <- 300           # Plot resolution (dots per inch)
PLOT_PT_SIZE <- 0.2       # Point size for UMAP plots
PLOT_LABEL_SIZE <- 3      # Label size for cluster labels

# Combined plot settings
COMBINED_PLOT_WIDTH <- 16              # Width for combined multi-panel plots
COMBINED_PLOT_HEIGHT <- 16             # Height for combined multi-panel plots
COMBINED_PLOT_WIDTH_SAMPLE <- 20       # Width for sample comparison plots
COMBINED_PLOT_HEIGHT_SAMPLE <- 16      # Height for sample comparison plots

# Resolution comparison grid settings
RES_GRID_WIDTH <- 18      # Width for resolution comparison grids
RES_GRID_HEIGHT <- 12     # Height for resolution comparison grids
RES_GRID_NCOL <- 3        # Number of columns in resolution grids

# Split plot settings
SPLIT_PLOT_WIDTH <- 20    # Width for split plots by sample
SPLIT_PLOT_HEIGHT <- 16   # Height for split plots by sample
SPLIT_PLOT_NCOL <- 4      # Number of columns for split plots

# Directory settings
WORK_DIR <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_seurat.BT.analysis.sep30.nov22.16samples"
TMP_DIR <- "/mnt/nfs/CX000008_DS1/projects/btanasa/tmp"
OUTPUT_DIR <- "plots_StandardSeurat_SCTransform_v5_OFFICIAL"

# scVI settings
SCVI_CONDA_ENV <- "/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/SeuratIntegrate"

# ============================================================================
# Record start time and setup logging
# ============================================================================
script_start_time <- Sys.time()
log_file <- paste0("integration_18samples_v5_OFFICIAL_", format(script_start_time, "%Y%m%d_%H%M%S"), ".log")
log_con <- file(log_file, open = "wt")

log_message <- function(..., level = "INFO", timestamp = TRUE) {
  msg <- paste(...)
  if (timestamp) {
    time_str <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    full_msg <- paste0("[", time_str, "] [", level, "] ", msg)
  } else {
    full_msg <- paste0("[", level, "] ", msg)
  }
  cat(full_msg, "\n")
  cat(full_msg, "\n", file = log_con, append = TRUE)
  flush(log_con)
}

log_cat <- function(...) {
  msg <- paste(...)
  cat(msg)
  cat(msg, file = log_con, append = TRUE)
  flush(log_con)
}

time_operation <- function(operation_name, expr) {
  log_message(sprintf("Starting: %s", operation_name))
  step_start <- Sys.time()
  result <- tryCatch({
    eval(expr)
  }, error = function(e) {
    step_duration <- difftime(Sys.time(), step_start, units = "mins")
    log_message(sprintf("ERROR in %s after %.2f minutes: %s", 
                       operation_name, as.numeric(step_duration), conditionMessage(e)), 
                level = "ERROR")
    stop(e)
  })
  step_duration <- difftime(Sys.time(), step_start, units = "mins")
  log_message(sprintf("Completed: %s (Duration: %.2f minutes)", 
                     operation_name, as.numeric(step_duration)))
  return(result)
}

sink(file = log_con, append = TRUE, type = "output", split = TRUE)

log_message(paste0(rep("=", 70), collapse = ""), timestamp = FALSE)
log_cat("OFFICIAL Seurat v5 SCTransform Integration Workflow\n")
log_cat("Script started at:", format(script_start_time, "%Y-%m-%d %H:%M:%S"), "\n")
log_cat(paste0(rep("=", 70), collapse = ""), "\n\n")

Sys.setenv(TMPDIR = TMP_DIR)

# ============================================================================
# Load Required Libraries
# ============================================================================
log_message("Loading required libraries...")

library(Seurat)
library(Matrix)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(harmony)
library(cowplot)
library(patchwork)
library(future)
library(SeuratIntegrate)

options(future.globals.maxSize = FUTURE_GLOBALS_MAX_SIZE)
log_message(sprintf("Future globals max size set to: %.0f GB", FUTURE_GLOBALS_MAX_SIZE / 1024^3))

# ============================================================================
# Function Definitions
# ============================================================================

show_summary <- function(obj, name) {
  log_message(sprintf("%s: %d genes, %d cells", name, nrow(obj), ncol(obj)))
}

check_memory_usage <- function(threshold = 0.85, stop_on_exceed = TRUE) {
  if (.Platform$OS.type == "unix") {
    meminfo <- tryCatch({
      readLines("/proc/meminfo", n = 3)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(meminfo)) {
      total_kb <- as.numeric(gsub("[^0-9]", "", meminfo[1]))
      available_kb <- as.numeric(gsub("[^0-9]", "", meminfo[3]))
      
      used_kb <- total_kb - available_kb
      usage_percent <- used_kb / total_kb
      
      log_message(sprintf("Memory usage: %.1f%% (%.1f GB used / %.1f GB total)", 
                         usage_percent * 100, 
                         used_kb / 1024^2, 
                         total_kb / 1024^2))
      
      if (usage_percent > threshold && stop_on_exceed) {
        log_message(sprintf("ERROR: Memory usage (%.1f%%) exceeds threshold (%.0f%%)", 
                           usage_percent * 100, threshold * 100), level = "ERROR")
        stop(sprintf("Memory usage (%.1f%%) exceeds threshold (%.0f%%). Stopping script.", 
                    usage_percent * 100, threshold * 100))
      }
      
      return(usage_percent)
    } else {
      log_message("Could not read memory info", level = "WARNING")
      return(NULL)
    }
  } else {
    log_message("Memory monitoring not available on this platform", level = "WARNING")
    return(NULL)
  }
}

check_and_calculate_percent_mt <- function(obj, sample_name) {
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    log_message(sprintf("  %s: percent.mt not found, calculating...", sample_name))
    mt_genes <- grep("^MT-", rownames(obj), value = TRUE)
    if (length(mt_genes) > 0) {
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
      log_message(sprintf("  %s: Calculated percent.mt using ^MT- pattern (human)", sample_name))
    } else {
      mt_genes <- grep("^mt-", rownames(obj), value = TRUE)
      if (length(mt_genes) > 0) {
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
        log_message(sprintf("  %s: Calculated percent.mt using ^mt- pattern (mouse)", sample_name))
      } else {
        log_message(sprintf("  %s: Warning - No mitochondrial genes found", sample_name), level = "WARNING")
      }
    }
  } else {
    log_message(sprintf("  %s: percent.mt already exists", sample_name))
  }
  return(obj)
}

# ============================================================================
# Check working directory
# ============================================================================
if (!dir.exists(WORK_DIR)) {
  stop("Working directory does not exist: ", WORK_DIR)
}
setwd(WORK_DIR)
log_message("Working directory set to:", getwd())

dir.create(OUTPUT_DIR, showWarnings = FALSE)
log_message(sprintf("Output directory: %s", OUTPUT_DIR))

# ============================================================================
# Load ORIGINAL filtered.soupX.rds files
# ============================================================================
log_message("Loading ORIGINAL filtered.soupX.rds files...")

tryCatch({
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
  
  log_message("All 18 Seurat objects loaded successfully")
  check_memory_usage(threshold = MEMORY_THRESHOLD, stop_on_exceed = TRUE)
}, error = function(e) {
  log_message("Error loading Seurat objects:", conditionMessage(e), level = "ERROR")
  stop("Failed to load input files. Exiting.\n")
})

# ============================================================================
# Check percent.mt
# ============================================================================
log_message("Checking for percent.mt in all samples...")

AD_AG_2658_seurat <- check_and_calculate_percent_mt(AD_AG_2658_seurat, "AD_AG_2658")
AD_AG_2822_seurat <- check_and_calculate_percent_mt(AD_AG_2822_seurat, "AD_AG_2822")
AD_AG_3184_seurat <- check_and_calculate_percent_mt(AD_AG_3184_seurat, "AD_AG_3184")
AD_AG_3188_seurat <- check_and_calculate_percent_mt(AD_AG_3188_seurat, "AD_AG_3188")
AD_CTL_2273_S2_seurat <- check_and_calculate_percent_mt(AD_CTL_2273_S2_seurat, "AD_CTL_2273_S2")
AD_CTL_2354_S3_seurat <- check_and_calculate_percent_mt(AD_CTL_2354_S3_seurat, "AD_CTL_2354_S3")
AD_CTL_2812_S1_seurat <- check_and_calculate_percent_mt(AD_CTL_2812_S1_seurat, "AD_CTL_2812_S1")
Dy_AD_2328_S5_seurat <- check_and_calculate_percent_mt(Dy_AD_2328_S5_seurat, "Dy_AD_2328_S5")
Dy_AD_2368_S6_seurat <- check_and_calculate_percent_mt(Dy_AD_2368_S6_seurat, "Dy_AD_2368_S6")
Dy_AD_2660_S4_seurat <- check_and_calculate_percent_mt(Dy_AD_2660_S4_seurat, "Dy_AD_2660_S4")
Dy_AG_2659_seurat <- check_and_calculate_percent_mt(Dy_AG_2659_seurat, "Dy_AG_2659")
Dy_AG_2922_seurat <- check_and_calculate_percent_mt(Dy_AG_2922_seurat, "Dy_AG_2922")
Dy_AG_3090_seurat <- check_and_calculate_percent_mt(Dy_AG_3090_seurat, "Dy_AG_3090")
Dy_AG_3296_seurat <- check_and_calculate_percent_mt(Dy_AG_3296_seurat, "Dy_AG_3296")
Healthy_AG_2394_seurat <- check_and_calculate_percent_mt(Healthy_AG_2394_seurat, "Healthy_AG_2394")
Healthy_AG_2594_seurat <- check_and_calculate_percent_mt(Healthy_AG_2594_seurat, "Healthy_AG_2594")
Healthy_AG_3041_seurat <- check_and_calculate_percent_mt(Healthy_AG_3041_seurat, "Healthy_AG_3041")
Healthy_AG_3267_seurat <- check_and_calculate_percent_mt(Healthy_AG_3267_seurat, "Healthy_AG_3267")

log_message("percent.mt check completed")

# ============================================================================
# Set orig.ident
# ============================================================================
AD_AG_2658_seurat$orig.ident <- "AD_AG_2658"
AD_AG_2822_seurat$orig.ident <- "AD_AG_2822"
AD_AG_3184_seurat$orig.ident <- "AD_AG_3184"
AD_AG_3188_seurat$orig.ident <- "AD_AG_3188"
AD_CTL_2273_S2_seurat$orig.ident <- "AD_CTL_2273_S2"
AD_CTL_2354_S3_seurat$orig.ident <- "AD_CTL_2354_S3"
AD_CTL_2812_S1_seurat$orig.ident <- "AD_CTL_2812_S1"
Dy_AD_2328_S5_seurat$orig.ident <- "Dy_AD_2328_S5"
Dy_AD_2368_S6_seurat$orig.ident <- "Dy_AD_2368_S6"
Dy_AD_2660_S4_seurat$orig.ident <- "Dy_AD_2660_S4"
Dy_AG_2659_seurat$orig.ident <- "Dy_AG_2659"
Dy_AG_2922_seurat$orig.ident <- "Dy_AG_2922"
Dy_AG_3090_seurat$orig.ident <- "Dy_AG_3090"
Dy_AG_3296_seurat$orig.ident <- "Dy_AG_3296"
Healthy_AG_2394_seurat$orig.ident <- "Healthy_AG_2394"
Healthy_AG_2594_seurat$orig.ident <- "Healthy_AG_2594"
Healthy_AG_3041_seurat$orig.ident <- "Healthy_AG_3041"
Healthy_AG_3267_seurat$orig.ident <- "Healthy_AG_3267"

# ============================================================================
# Create list and merge
# ============================================================================
seurat_list <- list(
  AD_AG_2658 = AD_AG_2658_seurat,
  AD_AG_2822 = AD_AG_2822_seurat,
  AD_AG_3184 = AD_AG_3184_seurat,
  AD_AG_3188 = AD_AG_3188_seurat,
  AD_CTL_2273_S2 = AD_CTL_2273_S2_seurat,
  AD_CTL_2354_S3 = AD_CTL_2354_S3_seurat,
  AD_CTL_2812_S1 = AD_CTL_2812_S1_seurat,
  Dy_AD_2328_S5 = Dy_AD_2328_S5_seurat,
  Dy_AD_2368_S6 = Dy_AD_2368_S6_seurat,
  Dy_AD_2660_S4 = Dy_AD_2660_S4_seurat,
  Dy_AG_2659 = Dy_AG_2659_seurat,
  Dy_AG_2922 = Dy_AG_2922_seurat,
  Dy_AG_3090 = Dy_AG_3090_seurat,
  Dy_AG_3296 = Dy_AG_3296_seurat,
  Healthy_AG_2394 = Healthy_AG_2394_seurat,
  Healthy_AG_2594 = Healthy_AG_2594_seurat,
  Healthy_AG_3041 = Healthy_AG_3041_seurat,
  Healthy_AG_3267 = Healthy_AG_3267_seurat
)

for (i in seq_along(seurat_list)) {
  seurat_list[[i]]$sample <- names(seurat_list)[i]
}

# ============================================================================
# Merge samples
# ============================================================================
obj <- time_operation("Merge samples", {
  result <- merge(seurat_list[[1]], y = seurat_list[-1], 
                  add.cell.ids = names(seurat_list),
                  merge.data = FALSE)
  
  log_cat(sprintf("  Total genes: %d\n", nrow(result)))
  log_cat(sprintf("  Total cells: %d\n", ncol(result)))
  
  check_memory_usage(threshold = MEMORY_THRESHOLD, stop_on_exceed = TRUE)
  
  # Clean up
  rm(AD_AG_2822_seurat, AD_AG_3184_seurat, AD_AG_3188_seurat,
     AD_CTL_2273_S2_seurat, AD_CTL_2354_S3_seurat, AD_CTL_2812_S1_seurat,
     Dy_AD_2328_S5_seurat, Dy_AD_2368_S6_seurat, Dy_AD_2660_S4_seurat,
     Dy_AG_2659_seurat, Dy_AG_2922_seurat, Dy_AG_3296_seurat,
     Healthy_AG_2394_seurat, Healthy_AG_2594_seurat, Healthy_AG_3041_seurat, 
     Healthy_AG_3267_seurat, envir = .GlobalEnv)
  rm(seurat_list, envir = .GlobalEnv)
  gc(verbose = FALSE)
  
  result
})

# ============================================================================
# OFFICIAL WORKFLOW STEP 1: Split by sample (check if already split)
# ============================================================================
log_message("Step 1: Checking if layers are already split...")

tryCatch({
  # Check if layers are already split
  rna_layers <- Layers(obj[["RNA"]], search = "counts")
  if (length(rna_layers) > 1) {
    log_message(sprintf("RNA assay already has %d layers split - skipping split step", length(rna_layers)))
    log_message(sprintf("  Layers found: %s", paste(head(rna_layers, 5), collapse = ", "), ifelse(length(rna_layers) > 5, "...", "")))
  } else {
    # Layers are not split, need to split them
    log_message("RNA assay layers not split - splitting by sample for Seurat v5 workflow...")
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
    log_message(sprintf("Split into %d layers", length(Layers(obj[["RNA"]]))))
  }
  gc(verbose = FALSE)
}, error = function(e) {
  # If error says layers are already split, that's okay - just log and continue
  if (grepl("already split", conditionMessage(e), ignore.case = TRUE)) {
    log_message("Layers are already split (detected from error message) - proceeding with SCTransform", level = "INFO")
    # Verify layers exist
    rna_layers <- tryCatch({
      Layers(obj[["RNA"]], search = "counts")
    }, error = function(e2) {
      NULL
    })
    if (!is.null(rna_layers) && length(rna_layers) > 1) {
      log_message(sprintf("Confirmed: %d layers are already split", length(rna_layers)))
    }
  } else {
    log_message("Error checking/splitting layers:", conditionMessage(e), level = "ERROR")
    stop("Failed to check/split layers. Exiting.\n")
  }
})

# ============================================================================
# OFFICIAL WORKFLOW STEP 2: SCTransform
# ============================================================================
obj <- time_operation("Step 2: SCTransform (OFFICIAL method)", {
  log_message("Running SCTransform as per official Seurat v5 documentation...")
  result <- SCTransform(obj, 
                       vars.to.regress = "percent.mt",
                       verbose = FALSE)
  
  log_cat(sprintf("  SCT assay created\n"))
  gc(verbose = FALSE)
  check_memory_usage(threshold = MEMORY_THRESHOLD, stop_on_exceed = TRUE)
  
  result
})

# ============================================================================
# OFFICIAL WORKFLOW STEP 3: RunPCA
# ============================================================================
obj <- time_operation("Step 3: RunPCA", {
  result <- RunPCA(obj, npcs = PCA_NDIMS, verbose = FALSE)
  gc(verbose = FALSE)
  result
})

# Elbow plot
tryCatch({
  elbow_plot <- ElbowPlot(obj, ndims = PCA_NDIMS, reduction = "pca")
  ggsave(paste0(OUTPUT_DIR, "/pca_elbow_plot.png"), 
         plot = elbow_plot, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
  rm(elbow_plot)
  gc(verbose = FALSE)
  log_message("PCA elbow plot saved")
}, error = function(e) {
  log_message("Error saving elbow plot:", conditionMessage(e), level = "ERROR")
})

# ============================================================================
# UNINTEGRATED ANALYSIS: FindNeighbors, FindClusters, RunUMAP
# ============================================================================
log_message("Performing unintegrated analysis (before integration)...")

# FindNeighbors on unintegrated PCA
tryCatch({
  if (!"pca" %in% names(obj@reductions)) {
    stop("PCA reduction not found. Please run RunPCA first.")
  }
  obj <- FindNeighbors(obj, 
                      dims = PCA_DIMS, 
                      reduction = "pca", 
                      graph.name = "graph_unintegrated",
                      verbose = FALSE)
  log_message("FindNeighbors (unintegrated) completed successfully")
}, error = function(e) {
  log_message("Error in FindNeighbors (unintegrated):", conditionMessage(e), level = "ERROR")
  stop("FindNeighbors (unintegrated) failed. Exiting.\n")
})

# FindClusters at multiple resolutions for unintegrated data
tryCatch({
  for (res in CLUSTER_RESOLUTIONS) {
    cluster_name <- sprintf("clusters_unintegrated_res%.1f", res)
    obj <- FindClusters(obj, 
                       resolution = res, 
                       graph.name = "graph_unintegrated",
                       cluster.name = cluster_name,
                       verbose = FALSE)
  }
  log_message(sprintf("FindClusters (unintegrated) completed for resolutions: %s", 
                     paste(CLUSTER_RESOLUTIONS, collapse = ", ")))
  gc(verbose = FALSE)
}, error = function(e) {
  log_message("Error in FindClusters (unintegrated):", conditionMessage(e), level = "ERROR")
  stop("FindClusters (unintegrated) failed. Exiting.\n")
})

# RunUMAP for unintegrated data
tryCatch({
  if (!"pca" %in% names(obj@reductions)) {
    stop("PCA reduction not found. Please run RunPCA first.")
  }
  obj <- RunUMAP(obj, 
                dims = PCA_DIMS, 
                reduction = "pca", 
                reduction.name = "umap.unintegrated",
                verbose = FALSE)
  log_message("RunUMAP (unintegrated) completed successfully")
}, error = function(e) {
  log_message("Error in RunUMAP (unintegrated):", conditionMessage(e), level = "ERROR")
  stop("RunUMAP (unintegrated) failed. Exiting.\n")
})

# Save unintegrated plots
tryCatch({
  # Ensure output directory exists
  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
    log_message(sprintf("Created output directory: %s", OUTPUT_DIR))
  }
  
  # Check if UMAP reduction exists
  if (!"umap.unintegrated" %in% names(obj@reductions)) {
    log_message("Error: 'umap.unintegrated' reduction not found. Cannot create unintegrated plots.", level = "ERROR")
    stop("umap.unintegrated reduction missing")
  }
  
  # Plot by sample
  log_message("Creating unintegrated UMAP plot by sample...")
  p_unintegrated_sample <- DimPlot(obj, reduction = "umap.unintegrated", group.by = "orig.ident")
  output_file <- paste0(OUTPUT_DIR, "/umap_unintegrated_by_sample.png")
  ggsave(output_file, plot = p_unintegrated_sample, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
  log_message(sprintf("Saved: %s", output_file))
  rm(p_unintegrated_sample)
  
  # Plot by clusters at multiple resolutions
  log_message("Creating unintegrated cluster plots at multiple resolutions...")
  log_message(sprintf("Available cluster columns: %s", 
                     paste(grep("^clusters_unintegrated_res", colnames(obj@meta.data), value = TRUE), collapse = ", ")))
  
  for (res in CLUSTER_RESOLUTIONS) {
    cluster_col <- sprintf("clusters_unintegrated_res%.1f", res)
    if (cluster_col %in% colnames(obj@meta.data)) {
      log_message(sprintf("Creating plot for %s...", cluster_col))
      p <- DimPlot(obj, reduction = "umap.unintegrated", group.by = cluster_col, label = TRUE)
      # Use same format as column name for filename
      filename_res <- sprintf("%.1f", res)
      output_file <- paste0(OUTPUT_DIR, "/umap_unintegrated_clusters_res", filename_res, ".png")
      ggsave(output_file, plot = p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
      log_message(sprintf("Saved: %s", output_file))
      rm(p)
    } else {
      log_message(sprintf("Warning: Cluster column '%s' not found in metadata. Available columns: %s", 
                         cluster_col, paste(grep("clusters_unintegrated", colnames(obj@meta.data), value = TRUE), collapse = ", ")), 
                 level = "WARNING")
    }
  }
  gc(verbose = FALSE)
  log_message("Unintegrated plots saved")
}, error = function(e) {
  log_message("Error saving unintegrated plots:", conditionMessage(e), level = "ERROR")
})

# ============================================================================
# Quality Control (QC) Plots (after UMAP is created)
# ============================================================================
log_message("Creating QC plots...")
tryCatch({
  # Check if QC metrics exist
  if (!"nFeature_RNA" %in% colnames(obj@meta.data)) {
    obj[["nFeature_RNA"]] <- colSums(obj[["RNA"]] > 0)
  }
  if (!"nCount_RNA" %in% colnames(obj@meta.data)) {
    obj[["nCount_RNA"]] <- colSums(obj[["RNA"]])
  }
  
  # QC plots on UMAP (FeaturePlot)
  if ("umap.unintegrated" %in% names(obj@reductions)) {
    p_qc1 <- FeaturePlot(obj, features = "nFeature_RNA", reduction = "umap.unintegrated", 
                         pt.size = PLOT_PT_SIZE, raster = TRUE) + 
      ggtitle("nFeature_RNA") +
      theme(plot.title = element_text(hjust = 0.5))
    
    p_qc2 <- FeaturePlot(obj, features = "nCount_RNA", reduction = "umap.unintegrated", 
                         pt.size = PLOT_PT_SIZE, raster = TRUE) + 
      ggtitle("nCount_RNA") +
      theme(plot.title = element_text(hjust = 0.5))
    
    p_qc3 <- FeaturePlot(obj, features = "percent.mt", reduction = "umap.unintegrated", 
                         pt.size = PLOT_PT_SIZE, raster = TRUE) + 
      ggtitle("percent.mt") +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Combine QC plots
    qc_combined <- p_qc1 + p_qc2 + p_qc3 + plot_layout(ncol = RES_GRID_NCOL, guides = 'collect')
    ggsave(paste0(OUTPUT_DIR, "/qc_metrics_unintegrated.png"), 
           plot = qc_combined, width = RES_GRID_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
    
    rm(p_qc1, p_qc2, p_qc3, qc_combined)
  }
  
  # Violin plots for QC metrics (don't need UMAP)
  p_violin1 <- VlnPlot(obj, features = "nFeature_RNA", group.by = "orig.ident", 
                      pt.size = 0, combine = FALSE)[[1]] + 
    NoLegend() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p_violin2 <- VlnPlot(obj, features = "nCount_RNA", group.by = "orig.ident", 
                      pt.size = 0, combine = FALSE)[[1]] + 
    NoLegend() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p_violin3 <- VlnPlot(obj, features = "percent.mt", group.by = "orig.ident", 
                      pt.size = 0, combine = FALSE)[[1]] + 
    NoLegend() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  violin_combined <- p_violin1 + p_violin2 + p_violin3 + plot_layout(ncol = RES_GRID_NCOL)
  ggsave(paste0(OUTPUT_DIR, "/qc_violin_plots.png"), 
         plot = violin_combined, width = RES_GRID_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
  
  rm(p_violin1, p_violin2, p_violin3, violin_combined)
  gc(verbose = FALSE)
  log_message("QC plots saved")
}, error = function(e) {
  log_message("Error creating QC plots (will continue without them):", conditionMessage(e), level = "WARNING")
})

# Resolution comparison grid for unintegrated data
tryCatch({
  log_message("Creating resolution comparison grid for unintegrated data...")
  plots_res <- list()
  for (res in CLUSTER_RESOLUTIONS) {
    cluster_col <- sprintf("clusters_unintegrated_res%.1f", res)
    if (cluster_col %in% colnames(obj@meta.data)) {
      p <- DimPlot(obj, reduction = "umap.unintegrated", group.by = cluster_col, 
                  label = TRUE, label.size = 2, pt.size = 0.1) +
        ggtitle(sprintf("Res %.1f", res)) +
        NoLegend()
      plots_res[[length(plots_res) + 1]] <- p
    }
  }
  if (length(plots_res) > 0) {
    res_grid <- wrap_plots(plots_res, ncol = RES_GRID_NCOL, guides = 'collect') +
      plot_annotation(title = "Unintegrated: Cluster Resolution Comparison")
    ggsave(paste0(OUTPUT_DIR, "/umap_unintegrated_resolution_comparison.png"), 
           plot = res_grid, width = RES_GRID_WIDTH, height = RES_GRID_HEIGHT, dpi = PLOT_DPI)
    rm(plots_res, res_grid)
    gc(verbose = FALSE)
    log_message("Unintegrated resolution comparison grid saved")
  }
}, error = function(e) {
  log_message("Error creating resolution comparison grid:", conditionMessage(e), level = "WARNING")
})

# Split plots by sample for unintegrated data
tryCatch({
  log_message("Creating split plots by sample for unintegrated data...")
  default_res <- DEFAULT_RESOLUTION
  cluster_col <- sprintf("clusters_unintegrated_res%.1f", default_res)
  if (cluster_col %in% colnames(obj@meta.data)) {
    p_split <- DimPlot(obj, reduction = "umap.unintegrated", group.by = cluster_col, 
                      split.by = "orig.ident", ncol = SPLIT_PLOT_NCOL, label = TRUE, pt.size = 0.1)
    ggsave(paste0(OUTPUT_DIR, "/umap_unintegrated_split_by_sample.png"), 
           plot = p_split, width = SPLIT_PLOT_WIDTH, height = SPLIT_PLOT_HEIGHT, dpi = PLOT_DPI)
    rm(p_split)
    gc(verbose = FALSE)
    log_message("Unintegrated split plots by sample saved")
  }
}, error = function(e) {
  log_message("Error creating split plots:", conditionMessage(e), level = "WARNING")
})

# ============================================================================
# OFFICIAL WORKFLOW STEP 4: IntegrateLayers with normalization.method = "SCT"
# This is the KEY difference from my previous script!
# ============================================================================

# --- RPCA Integration with SCT ---
obj <- time_operation("Step 4a: IntegrateLayers (RPCA + SCT)", {
  log_message("Using normalization.method = 'SCT' as per official documentation")
  result <- IntegrateLayers(
    object = obj,
    method = RPCAIntegration,
    normalization.method = "SCT",  # KEY PARAMETER!
    verbose = FALSE
  )
  log_cat("  Integration complete\n")
  log_cat("  Available reductions:", paste(names(result@reductions), collapse = ", "), "\n")
  gc(verbose = FALSE)
  check_memory_usage(threshold = MEMORY_THRESHOLD, stop_on_exceed = TRUE)
  result
})

# ============================================================================
# OFFICIAL WORKFLOW STEP 5: FindNeighbors using integrated.dr reduction
# ============================================================================
obj <- time_operation("Step 5: FindNeighbors (integrated.dr)", {
  result <- FindNeighbors(obj, 
                         dims = PCA_DIMS, 
                         reduction = "integrated.dr",  # Note: uses integrated.dr
                         verbose = FALSE)
  result
})

# ============================================================================
# OFFICIAL WORKFLOW STEP 6: FindClusters
# ============================================================================
log_message("Step 6: FindClusters at multiple resolutions...")
for (res in CLUSTER_RESOLUTIONS) {
  obj <- FindClusters(obj, 
                     resolution = res,
                     verbose = FALSE)
  # Rename cluster column for clarity
  # Try both formats: SCT_snn_res.0.1 and SCT_snn_res0.1
  cluster_col_old1 <- paste0("SCT_snn_res.", res)
  cluster_col_old2 <- paste0("SCT_snn_res", sprintf("%.1f", res))
  cluster_col_new <- sprintf("rpca_sct_clusters_res%.1f", res)
  
  if (cluster_col_old1 %in% colnames(obj@meta.data)) {
    obj@meta.data[[cluster_col_new]] <- obj@meta.data[[cluster_col_old1]]
    log_message(sprintf("Renamed %s to %s", cluster_col_old1, cluster_col_new))
  } else if (cluster_col_old2 %in% colnames(obj@meta.data)) {
    obj@meta.data[[cluster_col_new]] <- obj@meta.data[[cluster_col_old2]]
    log_message(sprintf("Renamed %s to %s", cluster_col_old2, cluster_col_new))
  } else {
    log_message(sprintf("Warning: Could not find cluster column for resolution %.1f (tried %s and %s)", 
                       res, cluster_col_old1, cluster_col_old2), level = "WARNING")
  }
}
log_message(sprintf("Clustering completed for resolutions: %s", 
                   paste(CLUSTER_RESOLUTIONS, collapse = ", ")))

# ============================================================================
# OFFICIAL WORKFLOW STEP 7: RunUMAP using integrated.dr
# ============================================================================
obj <- time_operation("Step 7: RunUMAP (integrated.dr)", {
  result <- RunUMAP(obj, 
                   dims = PCA_DIMS, 
                   reduction = "integrated.dr",  # Note: uses integrated.dr
                   verbose = FALSE)
  result
})

# Save RPCA+SCT integrated plots
tryCatch({
  # Ensure output directory exists
  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
    log_message(sprintf("Created output directory: %s", OUTPUT_DIR))
  }
  
  # Check if UMAP reduction exists
  if (!"umap" %in% names(obj@reductions)) {
    log_message("Error: 'umap' reduction not found. Cannot create RPCA plots.", level = "ERROR")
    stop("umap reduction missing")
  }
  
  log_message("Creating RPCA+SCT UMAP plot by sample...")
  p1 <- DimPlot(obj, reduction = "umap", group.by = "orig.ident")
  output_file1 <- paste0(OUTPUT_DIR, "/umap_rpca_sct_by_sample.png")
  ggsave(output_file1, plot = p1, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
  log_message(sprintf("Saved: %s", output_file1))
  
  log_message("Creating RPCA cluster plots at multiple resolutions...")
  log_message(sprintf("Available cluster columns: %s", 
                     paste(grep("^rpca_sct_clusters_res", colnames(obj@meta.data), value = TRUE), collapse = ", ")))
  
  for (res in CLUSTER_RESOLUTIONS) {
    cluster_col <- sprintf("rpca_sct_clusters_res%.1f", res)
    if (cluster_col %in% colnames(obj@meta.data)) {
      log_message(sprintf("Creating plot for %s...", cluster_col))
      p <- DimPlot(obj, reduction = "umap", group.by = cluster_col, label = TRUE)
      # Use same format as column name for filename
      filename_res <- sprintf("%.1f", res)
      output_file <- paste0(OUTPUT_DIR, "/umap_rpca_sct_clusters_res", filename_res, ".png")
      ggsave(output_file, plot = p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
      log_message(sprintf("Saved: %s", output_file))
      rm(p)
    } else {
      log_message(sprintf("Warning: Cluster column '%s' not found in metadata. Available columns: %s", 
                         cluster_col, paste(grep("rpca_sct_clusters", colnames(obj@meta.data), value = TRUE), collapse = ", ")), 
                 level = "WARNING")
    }
  }
  rm(p1)
  gc(verbose = FALSE)
  log_message("RPCA+SCT plots saved")
}, error = function(e) {
  log_message("Error saving plots:", conditionMessage(e), level = "ERROR")
})

# Resolution comparison grid for RPCA
tryCatch({
  log_message("Creating resolution comparison grid for RPCA...")
  plots_res <- list()
  for (res in CLUSTER_RESOLUTIONS) {
    cluster_col <- sprintf("rpca_sct_clusters_res%.1f", res)
    if (cluster_col %in% colnames(obj@meta.data)) {
      p <- DimPlot(obj, reduction = "umap", group.by = cluster_col, 
                  label = TRUE, label.size = 2, pt.size = 0.1) +
        ggtitle(sprintf("Res %.1f", res)) +
        NoLegend()
      plots_res[[length(plots_res) + 1]] <- p
    }
  }
  if (length(plots_res) > 0) {
    res_grid <- wrap_plots(plots_res, ncol = RES_GRID_NCOL, guides = 'collect') +
      plot_annotation(title = "RPCA+SCT: Cluster Resolution Comparison")
    ggsave(paste0(OUTPUT_DIR, "/umap_rpca_resolution_comparison.png"), 
           plot = res_grid, width = RES_GRID_WIDTH, height = RES_GRID_HEIGHT, dpi = PLOT_DPI)
    rm(plots_res, res_grid)
    gc(verbose = FALSE)
    log_message("RPCA resolution comparison grid saved")
  }
}, error = function(e) {
  log_message("Error creating RPCA resolution comparison grid:", conditionMessage(e), level = "WARNING")
})

# Split plots by sample for RPCA
tryCatch({
  log_message("Creating split plots by sample for RPCA...")
  default_res <- DEFAULT_RESOLUTION
  cluster_col <- sprintf("rpca_sct_clusters_res%.1f", default_res)
  if (cluster_col %in% colnames(obj@meta.data)) {
    p_split <- DimPlot(obj, reduction = "umap", group.by = cluster_col, 
                      split.by = "orig.ident", ncol = SPLIT_PLOT_NCOL, label = TRUE, pt.size = 0.1)
    ggsave(paste0(OUTPUT_DIR, "/umap_rpca_split_by_sample.png"), 
           plot = p_split, width = SPLIT_PLOT_WIDTH, height = SPLIT_PLOT_HEIGHT, dpi = PLOT_DPI)
    rm(p_split)
    gc(verbose = FALSE)
    log_message("RPCA split plots by sample saved")
  }
}, error = function(e) {
  log_message("Error creating RPCA split plots:", conditionMessage(e), level = "WARNING")
})

# ============================================================================
# Additional Integration: Harmony with SCT
# ============================================================================
# Layers should still be split from RPCA integration - no need to join/re-split
log_message("Preparing for Harmony integration - checking layers...")

tryCatch({
  # Check if layers are already split (they should be from RPCA)
  rna_layers <- Layers(obj[["RNA"]], search = "counts")
  if (length(rna_layers) > 1) {
    log_message(sprintf("RNA assay already has %d layers split - using existing layers for Harmony", length(rna_layers)))
  } else {
    log_message("Warning: Layers not split - splitting for Harmony integration...", level = "WARNING")
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
    log_message(sprintf("Split into %d layers for Harmony", length(Layers(obj[["RNA"]]))))
  }
}, error = function(e) {
  if (grepl("already split", conditionMessage(e), ignore.case = TRUE)) {
    log_message("Layers are already split - proceeding with Harmony integration", level = "INFO")
  } else {
    log_message("Error checking/splitting layers for Harmony:", conditionMessage(e), level = "ERROR")
    stop("Failed to prepare layers for Harmony. Exiting.\n")
  }
})

# https://satijalab.org/seurat/reference/harmonyintegration

obj <- time_operation("Harmony Integration with SCT", {
  log_message("Running Harmony integration with normalization.method = 'SCT'")
  result <- IntegrateLayers(
    object = obj,
    method = HarmonyIntegration,
    normalization.method = "SCT",
    orig.reduction = "pca",
    new.reduction = "harmony",
    verbose = FALSE
  )
  log_cat("  Harmony integration complete\n")
  gc(verbose = FALSE)
  check_memory_usage(threshold = MEMORY_THRESHOLD, stop_on_exceed = TRUE)
  result
})

# Harmony downstream analysis
obj <- FindNeighbors(obj, 
                    dims = PCA_DIMS, 
                    reduction = "harmony",
                    graph.name = "harmony_snn",
                    verbose = FALSE)

for (res in CLUSTER_RESOLUTIONS) {
  obj <- FindClusters(obj, 
                     resolution = res,
                     graph.name = "harmony_snn",
                     cluster.name = sprintf("harmony_sct_clusters_res%.1f", res),
                     verbose = FALSE)
}

obj <- RunUMAP(obj, 
              dims = PCA_DIMS, 
              reduction = "harmony",
              reduction.name = "umap.harmony",
              verbose = FALSE)

# Save Harmony plots
tryCatch({
  p1 <- DimPlot(obj, reduction = "umap.harmony", group.by = "orig.ident")
  ggsave(paste0(OUTPUT_DIR, "/umap_harmony_sct_by_sample.png"), 
         plot = p1, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
  
  for (res in CLUSTER_RESOLUTIONS) {
    cluster_col <- sprintf("harmony_sct_clusters_res%.1f", res)
    if (cluster_col %in% colnames(obj@meta.data)) {
      p <- DimPlot(obj, reduction = "umap.harmony", group.by = cluster_col, label = TRUE)
      ggsave(paste0(OUTPUT_DIR, "/umap_harmony_sct_clusters_res", res, ".png"), 
             plot = p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
      rm(p)
    }
  }
  rm(p1)
  gc(verbose = FALSE)
  log_message("Harmony+SCT plots saved")
}, error = function(e) {
  log_message("Error saving Harmony plots:", conditionMessage(e), level = "ERROR")
})

# Resolution comparison grid for Harmony
tryCatch({
  log_message("Creating resolution comparison grid for Harmony...")
  plots_res <- list()
  for (res in CLUSTER_RESOLUTIONS) {
    cluster_col <- sprintf("harmony_sct_clusters_res%.1f", res)
    if (cluster_col %in% colnames(obj@meta.data)) {
      p <- DimPlot(obj, reduction = "umap.harmony", group.by = cluster_col, 
                  label = TRUE, label.size = 2, pt.size = 0.1) +
        ggtitle(sprintf("Res %.1f", res)) +
        NoLegend()
      plots_res[[length(plots_res) + 1]] <- p
    }
  }
  if (length(plots_res) > 0) {
    res_grid <- wrap_plots(plots_res, ncol = RES_GRID_NCOL, guides = 'collect') +
      plot_annotation(title = "Harmony+SCT: Cluster Resolution Comparison")
    ggsave(paste0(OUTPUT_DIR, "/umap_harmony_resolution_comparison.png"), 
           plot = res_grid, width = RES_GRID_WIDTH, height = RES_GRID_HEIGHT, dpi = PLOT_DPI)
    rm(plots_res, res_grid)
    gc(verbose = FALSE)
    log_message("Harmony resolution comparison grid saved")
  }
}, error = function(e) {
  log_message("Error creating Harmony resolution comparison grid:", conditionMessage(e), level = "WARNING")
})

# Split plots by sample for Harmony
tryCatch({
  log_message("Creating split plots by sample for Harmony...")
  default_res <- DEFAULT_RESOLUTION
  cluster_col <- sprintf("harmony_sct_clusters_res%.1f", default_res)
  if (cluster_col %in% colnames(obj@meta.data)) {
    p_split <- DimPlot(obj, reduction = "umap.harmony", group.by = cluster_col, 
                      split.by = "orig.ident", ncol = SPLIT_PLOT_NCOL, label = TRUE, pt.size = 0.1)
    ggsave(paste0(OUTPUT_DIR, "/umap_harmony_split_by_sample.png"), 
           plot = p_split, width = SPLIT_PLOT_WIDTH, height = SPLIT_PLOT_HEIGHT, dpi = PLOT_DPI)
    rm(p_split)
    gc(verbose = FALSE)
    log_message("Harmony split plots by sample saved")
  }
}, error = function(e) {
  log_message("Error creating Harmony split plots:", conditionMessage(e), level = "WARNING")
})

# ============================================================================
# Additional Integration: scVI with raw counts
# ============================================================================
# NOTE: scVI expects raw counts from RNA assay, not normalized SCT data
# RNA layers must be split for integration
# SCT layers are NOT joined - they remain split (only RNA layers are joined at the end)
log_message("Preparing for scVI integration - checking layers (scVI uses raw counts from RNA)...")

tryCatch({
  # Check if RNA layers are split (needed for integration)
  # scVI uses raw counts from RNA assay, so we need RNA layers split
  rna_layers <- Layers(obj[["RNA"]], search = "counts")
  if (length(rna_layers) > 1) {
    log_message(sprintf("RNA assay has %d layers split - ready for scVI integration", length(rna_layers)))
  } else {
    log_message("Warning: RNA layers not split - splitting for scVI integration...", level = "WARNING")
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
    log_message(sprintf("Split RNA into %d layers for scVI", length(Layers(obj[["RNA"]]))))
  }
  
  # Note: We do NOT join SCT layers before scVI
  # scVI uses raw counts from RNA assay, so only RNA layers need to be split
  # SCT layers can remain split - they will be joined at the end if needed
}, error = function(e) {
  if (grepl("already split", conditionMessage(e), ignore.case = TRUE)) {
    log_message("Layers are already split - proceeding with scVI integration", level = "INFO")
  } else {
    log_message("Error checking/splitting layers for scVI:", conditionMessage(e), level = "ERROR")
    stop("Failed to prepare layers for scVI. Exiting.\n")
  }
})

obj <- time_operation("scVI Integration with raw counts", {
  log_message("Running scVI integration with raw counts from RNA assay (as per scVI documentation)")
  log_message("NOTE: scVI expects raw counts from RNA assay, not normalized SCT data")
  
  # Explicitly set default assay to RNA to ensure scVI uses RNA counts
  original_assay <- DefaultAssay(obj)
  DefaultAssay(obj) <- "RNA"
  log_message(sprintf("Set default assay to RNA (was: %s)", original_assay))
  
  # Check RNA assay layers (scVI needs raw counts from RNA assay)
  log_message("Checking RNA assay layer status before scVI integration...")
  rna_layers <- tryCatch({
    Layers(obj[["RNA"]], search = "counts")
  }, error = function(e) {
    NULL
  })
  if (!is.null(rna_layers)) {
    log_message(sprintf("RNA assay has %d layers with counts", length(rna_layers)))
    rna_layers_str <- paste(head(rna_layers, 5), collapse = ", ")
    if (length(rna_layers) > 5) rna_layers_str <- paste0(rna_layers_str, ", ...")
    log_message(sprintf("RNA layers: %s", rna_layers_str))
  } else {
    log_message("Warning: Could not detect RNA layers", level = "WARNING")
  }
  
  # Verify we're using counts layer
  log_message("Verifying RNA assay has counts layer...")
  if ("counts" %in% Layers(obj[["RNA"]])) {
    log_message("✓ RNA assay has 'counts' layer - ready for scVI")
  } else {
    log_message("Warning: 'counts' layer not found in RNA assay", level = "WARNING")
  }
  
  # scVI integration with raw counts from RNA assay
  # No normalization.method parameter - scVI uses raw counts from default assay (RNA)
  result <- tryCatch({
    IntegrateLayers(
      object = obj,
      method = scVIIntegration,
      # Explicitly using RNA assay (set as default above) with counts layer
      # No normalization.method - scVI uses raw counts from RNA assay
      new.reduction = "integrated.scvi",
      conda_env = SCVI_CONDA_ENV,
      verbose = TRUE
    )
  }, error = function(e) {
    error_msg <- conditionMessage(e)
    log_message(sprintf("Error during scVI integration: %s", error_msg), level = "ERROR")
    # Restore original assay before stopping
    DefaultAssay(obj) <- original_assay
    stop("scVI integration failed")
  })
  
  # Restore original default assay
  DefaultAssay(result) <- original_assay
  log_message(sprintf("Restored default assay to: %s", original_assay))
  
  log_cat("  scVI integration complete (used raw counts from RNA assay)\n")
  gc(verbose = FALSE)
  check_memory_usage(threshold = MEMORY_THRESHOLD, stop_on_exceed = TRUE)
  result
})

# scVI downstream analysis (use fewer dims for scVI)
SCVI_DIMS <- 1:10

obj <- FindNeighbors(obj, 
                    dims = SCVI_DIMS, 
                    reduction = "integrated.scvi",
                    graph.name = "scvi_snn",
                    verbose = FALSE)

for (res in CLUSTER_RESOLUTIONS) {
  obj <- FindClusters(obj, 
                     resolution = res,
                     graph.name = "scvi_snn",
                     cluster.name = sprintf("scvi_sct_clusters_res%.1f", res),
                     verbose = FALSE)
}

obj <- RunUMAP(obj, 
              dims = SCVI_DIMS, 
              reduction = "integrated.scvi",
              reduction.name = "umap.scvi",
              verbose = FALSE)

# Save scVI plots
tryCatch({
  p1 <- DimPlot(obj, reduction = "umap.scvi", group.by = "orig.ident")
  ggsave(paste0(OUTPUT_DIR, "/umap_scvi_sct_by_sample.png"), 
         plot = p1, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
  
  for (res in CLUSTER_RESOLUTIONS) {
    cluster_col <- sprintf("scvi_sct_clusters_res%.1f", res)
    if (cluster_col %in% colnames(obj@meta.data)) {
      p <- DimPlot(obj, reduction = "umap.scvi", group.by = cluster_col, label = TRUE)
      ggsave(paste0(OUTPUT_DIR, "/umap_scvi_sct_clusters_res", res, ".png"), 
             plot = p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
      rm(p)
    }
  }
  rm(p1)
  gc(verbose = FALSE)
  log_message("scVI+SCT plots saved")
}, error = function(e) {
  log_message("Error saving scVI plots:", conditionMessage(e), level = "ERROR")
})

# Resolution comparison grid for scVI
tryCatch({
  log_message("Creating resolution comparison grid for scVI...")
  plots_res <- list()
  for (res in CLUSTER_RESOLUTIONS) {
    cluster_col <- sprintf("scvi_sct_clusters_res%.1f", res)
    if (cluster_col %in% colnames(obj@meta.data)) {
      p <- DimPlot(obj, reduction = "umap.scvi", group.by = cluster_col, 
                  label = TRUE, label.size = 2, pt.size = 0.1) +
        ggtitle(sprintf("Res %.1f", res)) +
        NoLegend()
      plots_res[[length(plots_res) + 1]] <- p
    }
  }
  if (length(plots_res) > 0) {
    res_grid <- wrap_plots(plots_res, ncol = RES_GRID_NCOL, guides = 'collect') +
      plot_annotation(title = "scVI+SCT: Cluster Resolution Comparison")
    ggsave(paste0(OUTPUT_DIR, "/umap_scvi_resolution_comparison.png"), 
           plot = res_grid, width = RES_GRID_WIDTH, height = RES_GRID_HEIGHT, dpi = PLOT_DPI)
    rm(plots_res, res_grid)
    gc(verbose = FALSE)
    log_message("scVI resolution comparison grid saved")
  }
}, error = function(e) {
  log_message("Error creating scVI resolution comparison grid:", conditionMessage(e), level = "WARNING")
})

# Split plots by sample for scVI
tryCatch({
  log_message("Creating split plots by sample for scVI...")
  default_res <- DEFAULT_RESOLUTION
  cluster_col <- sprintf("scvi_sct_clusters_res%.1f", default_res)
  if (cluster_col %in% colnames(obj@meta.data)) {
    p_split <- DimPlot(obj, reduction = "umap.scvi", group.by = cluster_col, 
                      split.by = "orig.ident", ncol = SPLIT_PLOT_NCOL, label = TRUE, pt.size = 0.1)
    ggsave(paste0(OUTPUT_DIR, "/umap_scvi_split_by_sample.png"), 
           plot = p_split, width = SPLIT_PLOT_WIDTH, height = SPLIT_PLOT_HEIGHT, dpi = PLOT_DPI)
    rm(p_split)
    gc(verbose = FALSE)
    log_message("scVI split plots by sample saved")
  }
}, error = function(e) {
  log_message("Error creating scVI split plots:", conditionMessage(e), level = "WARNING")
})

# ============================================================================
# ALTERNATIVE HARMONY WORKFLOW: Traditional log-normalization pipeline
# ============================================================================
# This section demonstrates Harmony integration using the traditional workflow:
# Raw counts → NormalizeData (log-normalize) → FindVariableFeatures → ScaleData → RunPCA → Harmony
# This is an alternative to the SCTransform workflow shown above
# This workflow works with split layers (Seurat v5 structure)
# ============================================================================

log_message("\n========================================")
log_message("HARMONY TRADITIONAL WORKFLOW (with split layers)")
log_message("========================================")
log_message("Pipeline: Raw counts → NormalizeData → FindVariableFeatures → ScaleData → RunPCA → Harmony")

tryCatch({
  log_message("Starting traditional Harmony workflow with split layers...")
  
  # ============================================================================
  # Step 1: Verify RNA assay and layers are split
  # ============================================================================
  DefaultAssay(obj) <- "RNA"
  log_message("Step 1: Using RNA assay with raw counts")
  
  # Check if layers are split
  rna_layers <- Layers(obj[["RNA"]], search = "counts")
  if (length(rna_layers) > 1) {
    log_message(sprintf("✓ RNA assay has %d layers split (ready for per-layer processing)", length(rna_layers)))
  } else {
    log_message("Warning: RNA layers not split - this workflow expects split layers", level = "WARNING")
    log_message("Splitting RNA assay by orig.ident...")
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
    log_message(sprintf("✓ Split into %d layers", length(Layers(obj[["RNA"]]))))
  }
  
  # ============================================================================
  # Step 2: NormalizeData (log-normalize) - works on split layers
  # ============================================================================
  log_message("Step 2: Normalizing data (log-normalization) on split layers...")
  obj <- NormalizeData(obj,
                      normalization.method = "LogNormalize",
                      scale.factor = 10000,
                      verbose = FALSE)
  log_message("✓ Log-normalization complete (applied to each layer)")
  gc(verbose = FALSE)
  
  # ============================================================================
  # Step 3: FindVariableFeatures - works on split layers
  # ============================================================================
  log_message("Step 3: Finding variable features across all layers...")
  obj <- FindVariableFeatures(obj,
                             selection.method = "vst",
                             nfeatures = 2000,
                             verbose = FALSE)
  log_message(sprintf("✓ Found %d variable features", length(VariableFeatures(obj))))
  gc(verbose = FALSE)
  
  # ============================================================================
  # Step 4: ScaleData - works on split layers
  # ============================================================================
  log_message("Step 4: Scaling data (regressing out percent.mt)...")
  obj <- ScaleData(obj,
                  features = VariableFeatures(obj),
                  vars.to.regress = "percent.mt",
                  verbose = FALSE)
  log_message("✓ Data scaling complete (applied to each layer)")
  gc(verbose = FALSE)
  check_memory_usage(threshold = MEMORY_THRESHOLD, stop_on_exceed = TRUE)
  
  # ============================================================================
  # Step 5: RunPCA - works on split layers (using separate reduction name)
  # ============================================================================
  log_message("Step 5: Running PCA on variable features...")
  obj <- RunPCA(obj,
               features = VariableFeatures(obj),
               npcs = PCA_NDIMS,
               reduction.name = "pca.traditional",
               verbose = FALSE)
  log_message("✓ PCA complete (computed per layer, then combined)")
  log_message(sprintf("  PCA reduction created: pca.traditional with %d dimensions", ncol(obj[["pca.traditional"]])))
  gc(verbose = FALSE)
  
  # ============================================================================
  # Step 6: Harmony integration - works with PCA from split layers
  # ============================================================================
  log_message("Step 6: Running Harmony integration...")
  log_message("  Using PCA reduction (from split layers): pca.traditional")
  log_message("  Grouping by: orig.ident")
  
  obj <- RunHarmony(obj,
                   group.by.vars = "orig.ident",
                   reduction = "pca.traditional",
                   reduction.save = "harmony.traditional",
                   verbose = FALSE)
  log_message("✓ Harmony integration complete (traditional workflow)")
  log_message(sprintf("  Created reduction: %s", "harmony.traditional"))
  log_message(sprintf("  Harmony reduction dimensions: %d", ncol(obj[["harmony.traditional"]])))
  gc(verbose = FALSE)
  check_memory_usage(threshold = MEMORY_THRESHOLD, stop_on_exceed = TRUE)
  
  # ============================================================================
  # Downstream analysis with traditional Harmony
  # ============================================================================
  log_message("Performing downstream analysis with traditional Harmony...")
  
  # FindNeighbors using harmony.traditional reduction
  log_message("  Finding neighbors...")
  obj <- FindNeighbors(obj,
                      reduction = "harmony.traditional",
                      dims = PCA_DIMS,
                      graph.name = "harmony_traditional_snn",
                      verbose = FALSE)
  log_message("✓ Neighbors found")
  
  # FindClusters at multiple resolutions
  log_message("  Finding clusters at multiple resolutions...")
  for (res in CLUSTER_RESOLUTIONS) {
    obj <- FindClusters(obj,
                       resolution = res,
                       graph.name = "harmony_traditional_snn",
                       cluster.name = sprintf("harmony_traditional_clusters_res%.1f", res),
                       verbose = FALSE)
  }
  log_message(sprintf("✓ Clustering complete for resolutions: %s", paste(CLUSTER_RESOLUTIONS, collapse = ", ")))
  
  # RunUMAP
  log_message("  Running UMAP...")
  obj <- RunUMAP(obj,
                reduction = "harmony.traditional",
                dims = PCA_DIMS,
                reduction.name = "umap.harmony.traditional",
                verbose = FALSE)
  log_message("✓ UMAP complete (traditional Harmony)")
  
  # ============================================================================
  # Save plots for traditional Harmony
  # ============================================================================
  log_message("Creating plots for traditional Harmony workflow...")
  
  # Plot by sample
  p1 <- DimPlot(obj, reduction = "umap.harmony.traditional", group.by = "orig.ident")
  ggsave(paste0(OUTPUT_DIR, "/umap_harmony_traditional_by_sample.png"),
         plot = p1, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
  log_message("✓ Saved: umap_harmony_traditional_by_sample.png")
  rm(p1)
  
  # Plot by clusters at multiple resolutions
  for (res in CLUSTER_RESOLUTIONS) {
    cluster_col <- sprintf("harmony_traditional_clusters_res%.1f", res)
    if (cluster_col %in% colnames(obj@meta.data)) {
      p <- DimPlot(obj, reduction = "umap.harmony.traditional", group.by = cluster_col, label = TRUE)
      filename_res <- sprintf("%.1f", res)
      ggsave(paste0(OUTPUT_DIR, "/umap_harmony_traditional_clusters_res", filename_res, ".png"),
             plot = p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
      rm(p)
    }
  }
  log_message("✓ Saved cluster plots for all resolutions")
  
  # Resolution comparison grid
  plots_res <- list()
  for (res in CLUSTER_RESOLUTIONS) {
    cluster_col <- sprintf("harmony_traditional_clusters_res%.1f", res)
    if (cluster_col %in% colnames(obj@meta.data)) {
      p <- DimPlot(obj, reduction = "umap.harmony.traditional", group.by = cluster_col,
                  label = TRUE, label.size = 2, pt.size = 0.1) +
        ggtitle(sprintf("Res %.1f", res)) +
        NoLegend()
      plots_res[[length(plots_res) + 1]] <- p
    }
  }
  if (length(plots_res) > 0) {
    res_grid <- wrap_plots(plots_res, ncol = RES_GRID_NCOL, guides = 'collect') +
      plot_annotation(title = "Harmony (Traditional): Cluster Resolution Comparison")
    ggsave(paste0(OUTPUT_DIR, "/umap_harmony_traditional_resolution_comparison.png"),
           plot = res_grid, width = RES_GRID_WIDTH, height = RES_GRID_HEIGHT, dpi = PLOT_DPI)
    log_message("✓ Saved: umap_harmony_traditional_resolution_comparison.png")
    rm(plots_res, res_grid)
  }
  
  # Split plots by sample
  default_res <- DEFAULT_RESOLUTION
  cluster_col <- sprintf("harmony_traditional_clusters_res%.1f", default_res)
  if (cluster_col %in% colnames(obj@meta.data)) {
    p_split <- DimPlot(obj, reduction = "umap.harmony.traditional", group.by = cluster_col,
                      split.by = "orig.ident", ncol = SPLIT_PLOT_NCOL, label = TRUE, pt.size = 0.1)
    ggsave(paste0(OUTPUT_DIR, "/umap_harmony_traditional_split_by_sample.png"),
           plot = p_split, width = SPLIT_PLOT_WIDTH, height = SPLIT_PLOT_HEIGHT, dpi = PLOT_DPI)
    log_message("✓ Saved: umap_harmony_traditional_split_by_sample.png")
    rm(p_split)
  }
  
  gc(verbose = FALSE)
  
  # ============================================================================
  # Summary
  # ============================================================================
  log_message("\n✓ Traditional Harmony workflow complete!")
  log_message("  Reduction: harmony.traditional")
  log_message("  UMAP: umap.harmony.traditional")
  log_message("  Cluster columns: harmony_traditional_clusters_res*")
  log_message("  All plots saved to:", OUTPUT_DIR)
  
}, error = function(e) {
  log_message(sprintf("Error in traditional Harmony workflow: %s", conditionMessage(e)), level = "ERROR")
  log_message("Continuing with main workflow...", level = "WARNING")
})

# ============================================================================
# Join layers for final analysis - ONLY RNA assay
# ============================================================================
# Join layers for final analysis (RNA assay only) - COMMENTED OUT
# Note: Commenting out as it may not be required for all downstream analyses
# Uncomment if you need joined layers for specific operations
# ============================================================================
# log_message("Joining layers for final analysis (RNA assay only)...")
# tryCatch({
#   # Check if RNA assay has split layers
#   rna_layers <- tryCatch({
#     Layers(obj[["RNA"]], search = "counts")
#   }, error = function(e) {
#     NULL
#   })
#   
#   if (!is.null(rna_layers) && length(rna_layers) > 1) {
#     log_message(sprintf("RNA assay has %d layers - joining RNA layers", length(rna_layers)))
#     obj <- JoinLayers(obj, assay = "RNA")
#     log_message("RNA layers joined successfully")
#   } else {
#     log_message("RNA assay has no split layers or already joined")
#   }
# }, error = function(e) {
#   log_message(sprintf("Error joining RNA layers: %s", conditionMessage(e)), level = "ERROR")
#   log_message("Continuing without joining layers...", level = "WARNING")
# })

# ============================================================================
# Create comparison plots
# ============================================================================
tryCatch({
  library(patchwork)
  log_message("Creating method comparison plots...")
  
  default_res <- DEFAULT_RESOLUTION
  
  plots_list <- list()
  
  # Unintegrated
  cluster_col <- sprintf("clusters_unintegrated_res%.1f", default_res)
  if (cluster_col %in% colnames(obj@meta.data) && "umap.unintegrated" %in% names(obj@reductions)) {
    p_unintegrated <- DimPlot(obj, reduction = "umap.unintegrated", 
                             group.by = cluster_col, label = TRUE, label.size = 3, pt.size = PLOT_PT_SIZE) +
      ggtitle("Unintegrated")
    plots_list[[length(plots_list) + 1]] <- p_unintegrated
  }
  
  # RPCA+SCT
  cluster_col <- sprintf("rpca_sct_clusters_res%.1f", default_res)
  if (cluster_col %in% colnames(obj@meta.data) && "umap" %in% names(obj@reductions)) {
    p_rpca <- DimPlot(obj, reduction = "umap", 
                     group.by = cluster_col, label = TRUE, label.size = 3, pt.size = PLOT_PT_SIZE) +
      ggtitle("RPCA + SCT")
    plots_list[[length(plots_list) + 1]] <- p_rpca
  }
  
  # Harmony+SCT
  cluster_col <- sprintf("harmony_sct_clusters_res%.1f", default_res)
  if (cluster_col %in% colnames(obj@meta.data) && "umap.harmony" %in% names(obj@reductions)) {
    p_harmony <- DimPlot(obj, reduction = "umap.harmony", 
                        group.by = cluster_col, label = TRUE, label.size = 3, pt.size = PLOT_PT_SIZE) +
      ggtitle("Harmony + SCT")
    plots_list[[length(plots_list) + 1]] <- p_harmony
  }
  
  # scVI+SCT
  cluster_col <- sprintf("scvi_sct_clusters_res%.1f", default_res)
  if (cluster_col %in% colnames(obj@meta.data) && "umap.scvi" %in% names(obj@reductions)) {
    p_scvi <- DimPlot(obj, reduction = "umap.scvi", 
                     group.by = cluster_col, label = TRUE, label.size = 3, pt.size = PLOT_PT_SIZE) +
      ggtitle("scVI + SCT")
    plots_list[[length(plots_list) + 1]] <- p_scvi
  }
  
  if (length(plots_list) >= 2) {
    combined_plot <- wrap_plots(plots_list, ncol = 2, guides = 'collect') +
      plot_annotation(title = sprintf("SCTransform Integration Comparison (res%.1f)", default_res))
    
    ggsave(paste0(OUTPUT_DIR, "/umap_all_methods_sct_comparison.png"), 
           plot = combined_plot, width = 16, height = 16, dpi = 300)
    
    log_message("Comparison plot saved")
  }
  
  rm(plots_list)
  if (exists("p_unintegrated")) rm(p_unintegrated)
  if (exists("p_rpca")) rm(p_rpca)
  if (exists("p_harmony")) rm(p_harmony)
  if (exists("p_scvi")) rm(p_scvi)
  if (exists("combined_plot")) rm(combined_plot)
  gc(verbose = FALSE)
}, error = function(e) {
  log_message("Error creating comparison plots:", conditionMessage(e), level = "ERROR")
})

# ============================================================================
# Combined plot by sample (all methods) - quick visual comparison
# ============================================================================
tryCatch({
  library(patchwork)
  log_message("Creating combined plot by sample for all methods...")
  
  plots_sample_list <- list()
  
  # Unintegrated by sample
  if ("umap.unintegrated" %in% names(obj@reductions)) {
    p_unintegrated_sample <- DimPlot(obj, reduction = "umap.unintegrated", 
                                    group.by = "orig.ident", pt.size = PLOT_PT_SIZE) +
      ggtitle("Unintegrated")
    plots_sample_list[[length(plots_sample_list) + 1]] <- p_unintegrated_sample
  }
  
  # RPCA by sample
  if ("umap" %in% names(obj@reductions)) {
    p_rpca_sample <- DimPlot(obj, reduction = "umap", 
                            group.by = "orig.ident", pt.size = PLOT_PT_SIZE) +
      ggtitle("RPCA + SCT")
    plots_sample_list[[length(plots_sample_list) + 1]] <- p_rpca_sample
  }
  
  # Harmony by sample
  if ("umap.harmony" %in% names(obj@reductions)) {
    p_harmony_sample <- DimPlot(obj, reduction = "umap.harmony", 
                               group.by = "orig.ident", pt.size = PLOT_PT_SIZE) +
      ggtitle("Harmony + SCT")
    plots_sample_list[[length(plots_sample_list) + 1]] <- p_harmony_sample
  }
  
  # scVI by sample
  if ("umap.scvi" %in% names(obj@reductions)) {
    p_scvi_sample <- DimPlot(obj, reduction = "umap.scvi", 
                            group.by = "orig.ident", pt.size = PLOT_PT_SIZE) +
      ggtitle("scVI + SCT")
    plots_sample_list[[length(plots_sample_list) + 1]] <- p_scvi_sample
  }
  
  if (length(plots_sample_list) >= 2) {
    combined_sample_plot <- wrap_plots(plots_sample_list, ncol = 2, guides = 'collect') +
      plot_annotation(title = "Integration Methods Comparison by Sample")
    
    ggsave(paste0(OUTPUT_DIR, "/umap_all_methods_by_sample_comparison.png"), 
           plot = combined_sample_plot, width = 16, height = 16, dpi = 300)
    
    log_message("Combined plot by sample saved")
  }
  
  rm(plots_sample_list)
  if (exists("p_unintegrated_sample")) rm(p_unintegrated_sample)
  if (exists("p_rpca_sample")) rm(p_rpca_sample)
  if (exists("p_harmony_sample")) rm(p_harmony_sample)
  if (exists("p_scvi_sample")) rm(p_scvi_sample)
  if (exists("combined_sample_plot")) rm(combined_sample_plot)
  gc(verbose = FALSE)
}, error = function(e) {
  log_message("Error creating combined plot by sample:", conditionMessage(e), level = "ERROR")
})

# ============================================================================
# Save final object
# ============================================================================
gc(verbose = FALSE)
check_memory_usage(threshold = MEMORY_THRESHOLD, stop_on_exceed = TRUE)

tryCatch({
  saveRDS(obj, file = "18.samples.integrate.SeuratStandard.SCTransform.v5.OFFICIAL.rds")
  gc(verbose = FALSE)
  log_message("Final object saved: 18.samples.integrate.SeuratStandard.SCTransform.v5.OFFICIAL.rds")
}, error = function(e) {
  log_message("Error saving final object:", conditionMessage(e), level = "ERROR")
  stop("Failed to save final object")
})

# ============================================================================
# Final Summary
# ============================================================================
script_end_time <- Sys.time()
script_duration <- difftime(script_end_time, script_start_time, units = "mins")

log_cat("\n")
log_cat(paste0(rep("=", 70), collapse = ""), "\n")
log_cat("SCRIPT EXECUTION SUMMARY\n")
log_cat(paste0(rep("=", 70), collapse = ""), "\n")
log_cat("Start time:", format(script_start_time, "%Y-%m-%d %H:%M:%S"), "\n")
log_cat("End time:  ", format(script_end_time, "%Y-%m-%d %H:%M:%S"), "\n")
log_cat("Duration: ", round(as.numeric(script_duration), 2), "minutes\n")
log_cat("\n")

if (exists("obj")) {
  log_cat("Final Object Summary:\n")
  log_cat(sprintf("  Genes: %d\n", nrow(obj)))
  log_cat(sprintf("  Cells: %d\n", ncol(obj)))
  log_cat("  Reductions:", paste(names(obj@reductions), collapse = ", "), "\n")
  log_cat("  Assays:", paste(Assays(obj), collapse = ", "), "\n")
}

log_cat("\nOutput:\n")
log_cat("  Object: 18.samples.integrate.SeuratStandard.SCTransform.v5.OFFICIAL.rds\n")
log_cat("  Plots: ", OUTPUT_DIR, "/\n")
log_cat("  Log: ", log_file, "\n")
log_cat(paste0(rep("=", 70), collapse = ""), "\n")
log_message("Analysis complete!", timestamp = FALSE)
log_cat(paste0(rep("=", 70), collapse = ""), "\n")

# Join RNA layers at the end
log_message("Joining RNA layers for final object...")
tryCatch({
  rna_layers <- tryCatch({
    Layers(obj[["RNA"]], search = "counts")
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(rna_layers) && length(rna_layers) > 1) {
    log_message(sprintf("RNA assay has %d layers - joining RNA layers", length(rna_layers)))
    obj <- JoinLayers(obj, assay = "RNA")
    log_message("RNA layers joined successfully")
  } else {
    log_message("RNA assay has no split layers or already joined")
  }
}, error = function(e) {
  log_message(sprintf("Error joining RNA layers: %s", conditionMessage(e)), level = "ERROR")
  log_message("Continuing without joining layers...", level = "WARNING")
})

# Save object with joined layers
log_message("Saving object with joined layers...")
tryCatch({
  saveRDS(obj, file = "18.samples.integrate.SeuratStandard.SCTransform.v5.OFFICIAL.joinedLayers.rds")
  gc(verbose = FALSE)
  log_message("Object with joined layers saved: 18.samples.integrate.SeuratStandard.SCTransform.v5.OFFICIAL.joinedLayers.rds")
}, error = function(e) {
  log_message("Error saving object with joined layers:", conditionMessage(e), level = "ERROR")
  log_message("Continuing without saving joined layers object...", level = "WARNING")
})

sink()
close(log_con)
cat("Log file closed: ", log_file, "\n")

