################################################################################
#                                                                              #
#           COMPREHENSIVE SINGLE-CELL ANALYSIS SCRIPT                          #
#           Combined Script: Visualization, QC, and Table Generation           #
#           HARMONY VERSION                                                    #
#                                                                              #
################################################################################

# ==============================================================================
# SECTION 0: SETUP AND INITIALIZATION
# ==============================================================================

cat("\n==================================================\n")
cat("LOADING REQUIRED LIBRARIES\n")
cat("==================================================\n\n")

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(writexl)

cat("All libraries loaded successfully!\n")

# ==============================================================================
# CREATE OUTPUT DIRECTORY STRUCTURE
# ==============================================================================

cat("\n==================================================\n")
cat("CREATING OUTPUT DIRECTORY STRUCTURE\n")
cat("==================================================\n\n")

# Main output directory
output_dir <- "jk.18.samples.integrate.Seurat5.SCTransform.multiple.cluster.res.normalized_RNA.singleR.harmony.out"
dir.create(output_dir, showWarnings = FALSE)

# Subdirectories
plots_dir <- file.path(output_dir, "plots")
qc_dir <- file.path(output_dir, "qc")
tables_dir <- file.path(output_dir, "tables")

dir.create(plots_dir, showWarnings = FALSE)
dir.create(qc_dir, showWarnings = FALSE)
dir.create(tables_dir, showWarnings = FALSE)

# Additional subdirectories
dir.create(file.path(plots_dir, "by_resolution"), showWarnings = FALSE)
dir.create(file.path(qc_dir, "by_resolution"), showWarnings = FALSE)
dir.create(file.path(qc_dir, "seurat_style"), showWarnings = FALSE)
dir.create(file.path(tables_dir, "comprehensive"), showWarnings = FALSE)
dir.create(file.path(tables_dir, "QC"), showWarnings = FALSE)

cat("Output directory structure created:\n")
cat("  - ", output_dir, "/\n")
cat("    - plots/\n")
cat("      - by_resolution/\n")
cat("    - qc/\n")
cat("      - by_resolution/\n")
cat("      - seurat_style/\n")
cat("    - tables/\n")
cat("      - comprehensive/\n")
cat("      - QC/\n\n")

# ==============================================================================
# LOAD SEURAT OBJECT
# ==============================================================================

cat("\n==================================================\n")
cat("LOADING SEURAT OBJECT\n")
cat("==================================================\n\n")

# a <- readRDS("jk.18.samples.integrate.SeuratIntegrate.28nov.multiple.cluster.res.normalized_RNA.singleR.rds")

cat("Seurat object loaded successfully!\n")
cat("Object name: a\n")
cat("Number of cells:", ncol(a), "\n")
cat("Number of features:", nrow(a), "\n\n")

# ==============================================================================
# VALIDATION CHECKS
# ==============================================================================

cat("Performing validation checks...\n\n")

# Show available metadata columns
meta_cols <- colnames(a@meta.data)
cat("Available metadata columns (first 30):\n")
cat(paste(head(meta_cols, 30), collapse=", "), "\n")
if(length(meta_cols) > 30) {
  cat(sprintf("... and %d more columns\n", length(meta_cols) - 30))
}
cat("\n")

# Check if required columns exist
required_cols <- c("singleR_labels", "singleR_labels_pruned", "orig.ident")
missing_cols <- setdiff(required_cols, meta_cols)
if(length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse=", "))
}
cat("✓ All required metadata columns found\n")

# Check for condition column (with alternatives)
condition_col <- NULL
condition_alternatives <- c("condition", "sample", "group", "sample_group", "treatment", "Group", "Condition", "sample_condition", "Sample")
for (alt in condition_alternatives) {
  if (alt %in% meta_cols) {
    condition_col <- alt
    cat(sprintf("✓ Found condition column: '%s'\n", condition_col))
    break
  }
}

if (is.null(condition_col)) {
  cat("⚠ WARNING: No 'condition' column found (checked: ", paste(condition_alternatives, collapse=", "), ")\n")
  cat("  Creating a dummy 'condition' column with value 'All' to allow script execution.\n")
  cat("  If you have a condition column with a different name, please rename it to 'condition' in your Seurat object.\n")
  a@meta.data$condition <- "All"
  cat("  ✓ Dummy condition column created\n")
} else if (condition_col != "condition") {
  cat(sprintf("  Note: Using '%s' as condition column. Renaming to 'condition' for script compatibility...\n", condition_col))
  a@meta.data$condition <- a@meta.data[[condition_col]]
  cat("  ✓ Renamed successfully\n")
}
cat("\n")

# Check if UMAP exists
if(!"umap.harmony" %in% names(a@reductions)) {
  stop("UMAP reduction 'umap.harmony' not found!")
}
cat("✓ UMAP reduction 'umap.harmony' found\n\n")

# Check which harmony cluster columns exist in meta.data
cat("Checking for harmony cluster columns in meta.data...\n")
resolutions_check <- c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0")
meta_cols <- colnames(a@meta.data)

harmony_cols_found <- list()
for (res in resolutions_check) {
  cluster_col_sct <- paste0("harmony_sct_clusters_res", res)
  cluster_col_std <- paste0("harmony_clusters_res", res)
  
  if (cluster_col_sct %in% meta_cols) {
    harmony_cols_found[[res]] <- cluster_col_sct
    cat(sprintf("  ✓ Found: %s\n", cluster_col_sct))
  } else if (cluster_col_std %in% meta_cols) {
    harmony_cols_found[[res]] <- cluster_col_std
    cat(sprintf("  ✓ Found: %s\n", cluster_col_std))
  } else {
    cat(sprintf("  ✗ WARNING: Neither %s nor %s found for resolution %s\n", cluster_col_sct, cluster_col_std, res))
  }
}

cat(sprintf("\nTotal harmony cluster columns found: %d out of %d resolutions\n\n", 
            length(harmony_cols_found), length(resolutions_check)))

# Set resolutions variable for use in rest of script
resolutions <- resolutions_check

# ==============================================================================
# SECTION 1: SINGLER LABELS PROCESSING
# ==============================================================================

cat("\n==================================================\n")
cat("SECTION 1: PROCESSING SINGLER LABELS\n")
cat("==================================================\n\n")

# First, let's diagnose the issue
cat("Checking the data types and NAs:\n")
cat("Class of singleR_labels_pruned:", class(a@meta.data$singleR_labels_pruned), "\n")
cat("NAs in singleR_labels_pruned:", sum(is.na(a@meta.data$singleR_labels_pruned)), "\n")

# Convert to factors
a@meta.data$singleR_labels <- as.factor(a@meta.data$singleR_labels)
a@meta.data$singleR_labels_pruned <- as.factor(a@meta.data$singleR_labels_pruned)

# CORRECT way to relabel NAs:
a@meta.data$singleR_labels_pruned_labeled <- as.character(a@meta.data$singleR_labels_pruned)
a@meta.data$singleR_labels_pruned_labeled[is.na(a@meta.data$singleR_labels_pruned_labeled)] <- "Low_Confidence"
a@meta.data$singleR_labels_pruned_labeled <- as.factor(a@meta.data$singleR_labels_pruned_labeled)

# Verify the changes
cat("\n=== After correction ===\n")
cat("NAs in singleR_labels_pruned:", sum(is.na(a@meta.data$singleR_labels_pruned)), "\n")
cat("NAs in singleR_labels_pruned_labeled:", sum(is.na(a@meta.data$singleR_labels_pruned_labeled)), "\n")
cat("\n=== singleR_labels_pruned_labeled counts ===\n")
print(table(a@meta.data$singleR_labels_pruned_labeled))

# Check specifically for Low_Confidence
cat("\nLow_Confidence count:", sum(a@meta.data$singleR_labels_pruned_labeled == "Low_Confidence"), "\n")

cat("\nSingleR labels processing complete!\n")

# ==============================================================================
# SECTION 2: UMAP VISUALIZATIONS - HARMONY
# ==============================================================================

cat("\n==================================================\n")
cat("SECTION 2: CREATING UMAP VISUALIZATIONS (HARMONY)\n")
cat("==================================================\n\n")

# Define all resolutions for Harmony (already defined above, but keeping for clarity)
# resolutions <- c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0")

# Helper function to get the correct cluster column name
# Prioritize harmony_sct_clusters_res* (for 0.1-0.5 and 1.0), fall back to harmony_clusters_res* (for 0.6-0.9)
# Checks if columns exist in meta.data first
get_cluster_col <- function(res, meta_cols) {
  cluster_col_sct <- paste0("harmony_sct_clusters_res", res)
  cluster_col_std <- paste0("harmony_clusters_res", res)
  
  # Check if harmony_sct_clusters_res* exists first (preferred for 0.1-0.5 and 1.0)
  if (cluster_col_sct %in% meta_cols) {
    return(cluster_col_sct)
  } 
  # Fall back to harmony_clusters_res* (for 0.6-0.9)
  else if (cluster_col_std %in% meta_cols) {
    return(cluster_col_std)
  } 
  # Neither exists
  else {
    return(NULL)
  }
}

# Loop through each resolution
for (res in resolutions) {
  
  # Get the correct cluster column name for this resolution
  cluster_col <- get_cluster_col(res, colnames(a@meta.data))
  
  # Check if column exists
  if (is.null(cluster_col)) {
    print(paste0("Skipping resolution ", res, " - column not found"))
    next
  }
  
  tryCatch({
    # Create all 4 plots for this resolution with labels AND legends
    p1 <- DimPlot(a, reduction = "umap.harmony", 
                  group.by = "orig.ident",
                  raster = FALSE) + 
          ggtitle(paste0("Harmony Res ", res, " - Samples"))
    
    # Clusters with labels AND legend
    p2 <- DimPlot(a, reduction = "umap.harmony", 
                  group.by = cluster_col, 
                  label = TRUE,
                  label.size = 6,
                  raster = FALSE) + 
          ggtitle(paste0("Harmony Res ", res, " - Clusters"))
    
    # SingleR labels with labels AND legend
    p3 <- DimPlot(a, reduction = "umap.harmony", 
                  group.by = "singleR_labels",
                  label = TRUE,
                  label.size = 4,
                  repel = TRUE,
                  raster = FALSE) + 
          ggtitle(paste0("Harmony Res ", res, " - SingleR Labels"))
    
    # SingleR labels pruned with labels AND legend
    p4 <- DimPlot(a, reduction = "umap.harmony", 
                  group.by = "singleR_labels_pruned_labeled",
                  label = TRUE,
                  label.size = 4,
                  repel = TRUE,
                  raster = FALSE) + 
          ggtitle(paste0("Harmony Res ", res, " - SingleR Labels (Pruned)"))
    
    # Save combined plot (2x2 grid) - PNG only
    ggsave(file.path(plots_dir, "by_resolution", paste0("harmony_res", res, "_all_comparisons.png")),
           plot = (p1 | p2) / (p3 | p4),
           width = 28, height = 16, dpi = 300)
    
    # Save individual plots - PNG only
    ggsave(file.path(plots_dir, "by_resolution", paste0("harmony_res", res, "_samples.png")),
           plot = p1, width = 12, height = 8, dpi = 300)
    
    ggsave(file.path(plots_dir, "by_resolution", paste0("harmony_res", res, "_clusters.png")),
           plot = p2, width = 12, height = 8, dpi = 300)
    
    ggsave(file.path(plots_dir, "by_resolution", paste0("harmony_res", res, "_singleR.png")),
           plot = p3, width = 12, height = 8, dpi = 300)
    
    ggsave(file.path(plots_dir, "by_resolution", paste0("harmony_res", res, "_singleR_pruned.png")),
           plot = p4, width = 12, height = 8, dpi = 300)
    
    print(paste0("Saved plots for resolution ", res))
    
  }, error = function(e) {
    print(paste0("Error with resolution ", res, ": ", e$message))
  })
}

print(paste0("All resolution plots saved to '", plots_dir, "/by_resolution/' directory!"))

# Create mega comparison plots with labels AND legends - CLUSTERS
plot_list_legend <- list()

for (res in resolutions) {
  cluster_col <- get_cluster_col(res, colnames(a@meta.data))
  
  if (is.null(cluster_col)) {
    next
  }
  
  tryCatch({
    plot_list_legend[[res]] <- DimPlot(a, 
                                       reduction = "umap.harmony", 
                                       group.by = cluster_col, 
                                       label = TRUE,
                                       label.size = 4,
                                       raster = FALSE) + 
                               ggtitle(paste0("Res ", res))
  }, error = function(e) {
    print(paste0("Error creating plot for resolution ", res))
  })
}

# Combine all resolution plots with legends
if (length(plot_list_legend) > 0) {
  ggsave(file.path(plots_dir, "all_resolutions_comparison_with_legend.png"),
         plot = wrap_plots(plot_list_legend, ncol = 3),
         width = 36, height = 32, dpi = 300)
}

# Samples across all resolutions
sample_plots <- list()
for (res in resolutions) {
  cluster_col <- get_cluster_col(res, colnames(a@meta.data))
  
  if (is.null(cluster_col)) {
    next
  }
  
  tryCatch({
    sample_plots[[res]] <- DimPlot(a, 
                                   reduction = "umap.harmony", 
                                   group.by = "orig.ident",
                                   raster = FALSE) + 
                           ggtitle(paste0("Res ", res))
  }, error = function(e) {
    print(paste0("Error creating sample plot for resolution ", res))
  })
}

if (length(sample_plots) > 0) {
  ggsave(file.path(plots_dir, "samples_all_resolutions.png"),
         plot = wrap_plots(sample_plots, ncol = 4),
         width = 32, height = 24, dpi = 300)
}

# Cell types (PRUNED LABELED) across all resolutions - WITH LEGEND
celltype_plots_legend <- list()
for (res in resolutions) {
  cluster_col <- get_cluster_col(res, colnames(a@meta.data))
  
  if (is.null(cluster_col)) {
    next
  }
  
  tryCatch({
    celltype_plots_legend[[res]] <- DimPlot(a, 
                                            reduction = "umap.harmony", 
                                            group.by = "singleR_labels_pruned_labeled",
                                            label = TRUE, 
                                            label.size = 3,
                                            repel = TRUE,
                                            raster = FALSE) +
                                    ggtitle(paste0("Res ", res))
  }, error = function(e) {
    print(paste0("Error creating celltype plot for resolution ", res))
  })
}

if (length(celltype_plots_legend) > 0) {
  ggsave(file.path(plots_dir, "celltypes_pruned_all_resolutions.png"),
         plot = wrap_plots(celltype_plots_legend, ncol = 3),
         width = 36, height = 32, dpi = 300)
}

# SingleR labels (UNPRUNED) across all resolutions - WITH LEGEND
singleR_plots_legend <- list()
for (res in resolutions) {
  cluster_col <- get_cluster_col(res, colnames(a@meta.data))
  
  if (is.null(cluster_col)) {
    next
  }
  
  tryCatch({
    singleR_plots_legend[[res]] <- DimPlot(a, 
                                           reduction = "umap.harmony", 
                                           group.by = "singleR_labels",
                                           label = TRUE,
                                           label.size = 3,
                                           repel = TRUE,
                                           raster = FALSE) + 
                                   ggtitle(paste0("Res ", res))
  }, error = function(e) {
    print(paste0("Error creating singleR plot for resolution ", res))
  })
}

if (length(singleR_plots_legend) > 0) {
  ggsave(file.path(plots_dir, "singleR_unpruned_all_resolutions.png"),
         plot = wrap_plots(singleR_plots_legend, ncol = 3),
         width = 36, height = 32, dpi = 300)
}

print("All mega comparison plots saved!")

# BONUS: Create comprehensive comparisons
p_unpruned <- DimPlot(a, reduction = "umap.harmony", 
                      group.by = "singleR_labels",
                      label = TRUE,
                      label.size = 5,
                      repel = TRUE,
                      raster = FALSE) + 
              ggtitle("SingleR Labels (Unpruned)")

p_pruned <- DimPlot(a, reduction = "umap.harmony", 
                    group.by = "singleR_labels_pruned_labeled",
                    label = TRUE,
                    label.size = 5,
                    repel = TRUE,
                    raster = FALSE) + 
            ggtitle("SingleR Labels (Pruned - Low_Confidence shown)")

ggsave(file.path(plots_dir, "singleR_unpruned_vs_pruned_comparison.png"),
       plot = p_unpruned | p_pruned,
       width = 24, height = 8, dpi = 300)

print("Bonus comparison plots saved!")

cat("\nUMAP visualizations complete!\n")

# ==============================================================================
# SECTION 3: COMPREHENSIVE TABLES GENERATION
# ==============================================================================

cat("\n==================================================\n")
cat("SECTION 3: GENERATING COMPREHENSIVE TABLES\n")
cat("==================================================\n\n")

# Loop through each resolution to create comprehensive tables
for (res in resolutions) {
  
  cluster_col <- get_cluster_col(res, colnames(a@meta.data))
  
  # Check if column exists
  if (is.null(cluster_col)) {
    print(paste0("Skipping resolution ", res, " - column not found"))
    next
  }
  
  cat(paste0("\n=== Processing Resolution ", res, " ===\n"))
  
  # Create a data frame with all relevant columns to ensure alignment
  df_temp <- data.frame(
    cluster = a@meta.data[[cluster_col]],
    condition = a@meta.data$condition,
    singleR_labels = a@meta.data$singleR_labels,
    singleR_labels_pruned_labeled = a@meta.data$singleR_labels_pruned_labeled,
    stringsAsFactors = FALSE
  )
  
  # Remove any rows with NA in cluster (shouldn't happen, but safety check)
  df_temp <- df_temp[!is.na(df_temp$cluster), ]
  
  # Get cluster counts
  cluster_counts <- as.data.frame(table(df_temp$cluster, useNA = "no"))
  colnames(cluster_counts) <- c("Cluster", "Total_Count")
  
  # Get cluster by condition
  cluster_by_condition <- as.data.frame.matrix(table(df_temp$cluster, 
                                                      df_temp$condition, 
                                                      useNA = "no"))
  cluster_by_condition$Cluster <- rownames(cluster_by_condition)
  
  # Get cluster vs singleR
  cluster_vs_singleR <- as.data.frame.matrix(table(df_temp$cluster, 
                                                    df_temp$singleR_labels, 
                                                    useNA = "no"))
  cluster_vs_singleR$Cluster <- rownames(cluster_vs_singleR)
  
  # Get cluster vs singleR pruned
  cluster_vs_singleR_pruned <- as.data.frame.matrix(table(df_temp$cluster, 
                                                           df_temp$singleR_labels_pruned_labeled, 
                                                           useNA = "no"))
  cluster_vs_singleR_pruned$Cluster <- rownames(cluster_vs_singleR_pruned)
  
  # Merge all cluster tables
  table1_cluster <- cluster_counts %>%
    left_join(cluster_by_condition, by = "Cluster") %>%
    left_join(cluster_vs_singleR, by = "Cluster", suffix = c("", "_singleR")) %>%
    left_join(cluster_vs_singleR_pruned, by = "Cluster", suffix = c("", "_pruned"))
  
  # Reorder columns
  condition_cols <- colnames(cluster_by_condition)[colnames(cluster_by_condition) != "Cluster"]
  singleR_cols <- colnames(cluster_vs_singleR)[colnames(cluster_vs_singleR) != "Cluster"]
  singleR_pruned_cols <- colnames(cluster_vs_singleR_pruned)[colnames(cluster_vs_singleR_pruned) != "Cluster"]
  
  col_order <- c("Cluster", "Total_Count", condition_cols, singleR_cols, singleR_pruned_cols)
  col_order <- col_order[col_order %in% colnames(table1_cluster)]
  table1_cluster <- table1_cluster[, col_order]
  
  # Save Table 1
  write.csv(table1_cluster, 
            file.path(tables_dir, "comprehensive", paste0("resolution_", res, "_TABLE1_clusters_comprehensive.csv")), 
            row.names = FALSE)
  
  cat(paste0("  Saved Table 1: Clusters comprehensive (", nrow(table1_cluster), " clusters)\n"))
  
  # Get singleR counts (using aligned data)
  singleR_counts <- as.data.frame(table(df_temp$singleR_labels, useNA = "no"))
  colnames(singleR_counts) <- c("Cell_Type", "Total_Count")
  
  # Get singleR by condition
  singleR_by_condition <- as.data.frame.matrix(table(df_temp$singleR_labels, 
                                                      df_temp$condition, 
                                                      useNA = "no"))
  singleR_by_condition$Cell_Type <- rownames(singleR_by_condition)
  
  # Merge singleR tables
  table2_singleR <- singleR_counts %>%
    left_join(singleR_by_condition, by = "Cell_Type")
  
  col_order2 <- c("Cell_Type", "Total_Count", condition_cols)
  col_order2 <- col_order2[col_order2 %in% colnames(table2_singleR)]
  table2_singleR <- table2_singleR[, col_order2]
  
  # Save Table 2
  write.csv(table2_singleR, 
            file.path(tables_dir, "comprehensive", paste0("resolution_", res, "_TABLE2_singleR_comprehensive.csv")), 
            row.names = FALSE)
  
  cat(paste0("  Saved Table 2: SingleR comprehensive (", nrow(table2_singleR), " cell types)\n"))
  
  # Get singleR pruned counts (using aligned data)
  singleR_pruned_counts <- as.data.frame(table(df_temp$singleR_labels_pruned_labeled, useNA = "no"))
  colnames(singleR_pruned_counts) <- c("Cell_Type", "Total_Count")
  
  # Get singleR pruned by condition
  singleR_pruned_by_condition <- as.data.frame.matrix(table(df_temp$singleR_labels_pruned_labeled, 
                                                             df_temp$condition, 
                                                             useNA = "no"))
  singleR_pruned_by_condition$Cell_Type <- rownames(singleR_pruned_by_condition)
  
  # Merge singleR pruned tables
  table3_singleR_pruned <- singleR_pruned_counts %>%
    left_join(singleR_pruned_by_condition, by = "Cell_Type")
  
  col_order3 <- c("Cell_Type", "Total_Count", condition_cols)
  col_order3 <- col_order3[col_order3 %in% colnames(table3_singleR_pruned)]
  table3_singleR_pruned <- table3_singleR_pruned[, col_order3]
  
  # Save Table 3
  write.csv(table3_singleR_pruned, 
            file.path(tables_dir, "comprehensive", paste0("resolution_", res, "_TABLE3_singleR_pruned_comprehensive.csv")), 
            row.names = FALSE)
  
  cat(paste0("  Saved Table 3: SingleR pruned comprehensive (", nrow(table3_singleR_pruned), " cell types)\n"))
  
  # Save all 3 tables in one Excel file
  comprehensive_tables <- list(
    "1_Clusters_Comprehensive" = table1_cluster,
    "2_SingleR_Comprehensive" = table2_singleR,
    "3_SingleR_Pruned_Comprehensive" = table3_singleR_pruned
  )
  
  write_xlsx(comprehensive_tables, 
             path = file.path(tables_dir, "comprehensive", paste0("resolution_", res, "_ALL_COMPREHENSIVE_TABLES.xlsx")))
  
  cat(paste0("  Saved Excel file with all 3 comprehensive tables\n"))
  
  # Clean up temporary data frame
  rm(df_temp)
  gc(verbose = FALSE)
}

cat("\nComprehensive tables generation complete!\n")

# ==============================================================================
# SECTION 4: QC VIOLIN PLOTS - CUSTOM ggplot2 STYLE
# ==============================================================================

cat("\n==================================================\n")
cat("SECTION 4: CREATING CUSTOM QC VIOLIN PLOTS\n")
cat("==================================================\n\n")

# Define QC metrics
qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hb")

# Function to create ggplot2 violin plots
create_qc_violin_plot <- function(data, group_by, qc_metric, title) {
  p <- ggplot(data, aes_string(x = group_by, y = qc_metric, fill = group_by)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "none",
          plot.title = element_text(size = 10, face = "bold")) +
    labs(title = title, x = "", y = qc_metric) +
    scale_y_continuous(labels = scales::comma)
  
  return(p)
}

# Loop through each resolution for custom QC plots
for (res in resolutions) {
  
  cluster_col <- get_cluster_col(res, colnames(a@meta.data))
  
  if (is.null(cluster_col)) {
    print(paste0("Skipping resolution ", res, " - column not found"))
    next
  }
  
  cat(paste0("\n=== Creating custom QC plots for Resolution ", res, " ===\n"))
  
  # QC plots by HARMONY CLUSTERS
  qc_plots_cluster <- list()
  
  for (qc in qc_metrics) {
    qc_plots_cluster[[qc]] <- create_qc_violin_plot(
      a@meta.data, 
      cluster_col, 
      qc, 
      paste0("Res ", res, " - ", qc)
    )
  }
  
  combined_cluster <- wrap_plots(qc_plots_cluster, ncol = 2)
  
  ggsave(file.path(qc_dir, "by_resolution", paste0("resolution_", res, "_QC_by_harmony_clusters_custom.png")),
         plot = combined_cluster,
         width = 14, height = 12, dpi = 300)
  
  cat(paste0("  Saved: Custom QC by harmony clusters\n"))
  
  # QC plots by singleR_labels
  qc_plots_singleR <- list()
  for (qc in qc_metrics) {
    qc_plots_singleR[[qc]] <- create_qc_violin_plot(
      a@meta.data, 
      "singleR_labels", 
      qc, 
      paste0("Res ", res, " - ", qc)
    )
  }
  combined_singleR <- wrap_plots(qc_plots_singleR, ncol = 2)
  ggsave(file.path(qc_dir, "by_resolution", paste0("resolution_", res, "_QC_by_singleR_custom.png")),
         plot = combined_singleR,
         width = 14, height = 12, dpi = 300)
  
  # QC plots by singleR_labels_pruned_labeled
  qc_plots_singleR_pruned <- list()
  for (qc in qc_metrics) {
    qc_plots_singleR_pruned[[qc]] <- create_qc_violin_plot(
      a@meta.data, 
      "singleR_labels_pruned_labeled", 
      qc, 
      paste0("Res ", res, " - ", qc)
    )
  }
  combined_singleR_pruned <- wrap_plots(qc_plots_singleR_pruned, ncol = 2)
  ggsave(file.path(qc_dir, "by_resolution", paste0("resolution_", res, "_QC_by_singleR_pruned_custom.png")),
         plot = combined_singleR_pruned,
         width = 14, height = 12, dpi = 300)
  
  # QC plots by condition
  qc_plots_condition <- list()
  for (qc in qc_metrics) {
    qc_plots_condition[[qc]] <- create_qc_violin_plot(
      a@meta.data, 
      "condition", 
      qc, 
      paste0("Res ", res, " - ", qc)
    )
  }
  combined_condition <- wrap_plots(qc_plots_condition, ncol = 2)
  ggsave(file.path(qc_dir, "by_resolution", paste0("resolution_", res, "_QC_by_condition_custom.png")),
         plot = combined_condition,
         width = 14, height = 12, dpi = 300)
}

# Create overall comparison plots
for (qc in qc_metrics) {
  
  cat(paste0("\nProcessing custom ", qc, " across all resolutions...\n"))
  
  # Across ALL resolutions for harmony clusters
  plot_list_all_res <- list()
  
  for (res in resolutions) {
    cluster_col <- get_cluster_col(res, colnames(a@meta.data))
    
    if (!is.null(cluster_col)) {
      plot_list_all_res[[res]] <- create_qc_violin_plot(
        a@meta.data, 
        cluster_col, 
        qc, 
        paste0("Res ", res)
      )
    }
  }
  
  if (length(plot_list_all_res) > 0) {
    combined_all_res <- wrap_plots(plot_list_all_res, ncol = 4)
    
    ggsave(file.path(qc_dir, paste0(qc, "_by_harmony_clusters_all_resolutions_custom.png")),
           plot = combined_all_res,
           width = 32, height = 24, dpi = 300)
  }
  
  # By condition
  p_condition <- create_qc_violin_plot(a@meta.data, "condition", qc, 
                                       paste0(qc, " by Condition"))
  
  ggsave(file.path(qc_dir, paste0(qc, "_by_condition_overall_custom.png")),
         plot = p_condition,
         width = 10, height = 8, dpi = 300)
  
  # By singleR labels
  p_singleR <- create_qc_violin_plot(a@meta.data, "singleR_labels", qc, 
                                     paste0(qc, " by SingleR Labels"))
  
  ggsave(file.path(qc_dir, paste0(qc, "_by_singleR_overall_custom.png")),
         plot = p_singleR,
         width = 12, height = 8, dpi = 300)
  
  # By singleR pruned
  p_singleR_pruned <- create_qc_violin_plot(a@meta.data, "singleR_labels_pruned_labeled", qc, 
                                            paste0(qc, " by SingleR Pruned"))
  
  ggsave(file.path(qc_dir, paste0(qc, "_by_singleR_pruned_overall_custom.png")),
         plot = p_singleR_pruned,
         width = 12, height = 8, dpi = 300)
}

cat("\nCustom QC violin plots complete!\n")

# ==============================================================================
# SECTION 5: QC VIOLIN PLOTS - SEURAT STYLE
# ==============================================================================

cat("\n==================================================\n")
cat("SECTION 5: CREATING SEURAT-STYLE QC VIOLIN PLOTS\n")
cat("==================================================\n\n")

# Loop through each resolution for Seurat-style plots
for (res in resolutions) {
  
  cluster_col <- get_cluster_col(res, colnames(a@meta.data))
  
  if (is.null(cluster_col)) {
    print(paste0("Skipping resolution ", res, " - column not found"))
    next
  }
  
  cat(paste0("\n=== Creating Seurat-style QC plots for Resolution ", res, " ===\n"))
  
  # SEURAT-STYLE VlnPlot for HARMONY CLUSTERS
  Idents(a) <- cluster_col
  
  p_vln_cluster <- VlnPlot(a, 
                           features = qc_metrics, 
                           ncol = 3,
                           pt.size = 0.1)
  
  ggsave(file.path(qc_dir, "seurat_style", paste0("resolution_", res, "_VlnPlot_harmony_clusters.png")),
         plot = p_vln_cluster,
         width = 18, height = 12, dpi = 300)
  
  cat(paste0("  Saved: Seurat VlnPlot for harmony clusters\n"))
  
  # Create individual Seurat VlnPlots for each QC metric
  for (qc in qc_metrics) {
    p_vln_individual <- VlnPlot(a, 
                                features = qc,
                                pt.size = 0.1) +
      ggtitle(paste0("Harmony Clusters Res ", res, " - ", qc))
    
    ggsave(file.path(qc_dir, "seurat_style", paste0("resolution_", res, "_VlnPlot_", qc, "_harmony_clusters.png")),
           plot = p_vln_individual,
           width = 12, height = 8, dpi = 300)
  }
  
  # SEURAT-STYLE VlnPlot for singleR_labels
  Idents(a) <- "singleR_labels"
  
  p_vln_singleR <- VlnPlot(a, 
                           features = qc_metrics, 
                           ncol = 3,
                           pt.size = 0.1)
  
  ggsave(file.path(qc_dir, "seurat_style", paste0("resolution_", res, "_VlnPlot_singleR.png")),
         plot = p_vln_singleR,
         width = 18, height = 12, dpi = 300)
  
  # SEURAT-STYLE VlnPlot for singleR_labels_pruned_labeled
  Idents(a) <- "singleR_labels_pruned_labeled"
  
  p_vln_singleR_pruned <- VlnPlot(a, 
                                  features = qc_metrics, 
                                  ncol = 3,
                                  pt.size = 0.1)
  
  ggsave(file.path(qc_dir, "seurat_style", paste0("resolution_", res, "_VlnPlot_singleR_pruned.png")),
         plot = p_vln_singleR_pruned,
         width = 18, height = 12, dpi = 300)
  
  # SEURAT-STYLE VlnPlot for condition
  Idents(a) <- "condition"
  
  p_vln_condition <- VlnPlot(a, 
                             features = qc_metrics, 
                             ncol = 3,
                             pt.size = 0.1)
  
  ggsave(file.path(qc_dir, "seurat_style", paste0("resolution_", res, "_VlnPlot_condition.png")),
         plot = p_vln_condition,
         width = 18, height = 12, dpi = 300)
}

# Overall Seurat VlnPlots
cat("\n=== Creating overall Seurat-style VlnPlots ===\n")

Idents(a) <- "singleR_labels"
p_vln_singleR_overall <- VlnPlot(a, features = qc_metrics, ncol = 3, pt.size = 0.1)
ggsave(file.path(qc_dir, "seurat_style", "VlnPlot_singleR_overall.png"),
       plot = p_vln_singleR_overall,
       width = 18, height = 12, dpi = 300)

Idents(a) <- "singleR_labels_pruned_labeled"
p_vln_singleR_pruned_overall <- VlnPlot(a, features = qc_metrics, ncol = 3, pt.size = 0.1)
ggsave(file.path(qc_dir, "seurat_style", "VlnPlot_singleR_pruned_overall.png"),
       plot = p_vln_singleR_pruned_overall,
       width = 18, height = 12, dpi = 300)

Idents(a) <- "condition"
p_vln_condition_overall <- VlnPlot(a, features = qc_metrics, ncol = 3, pt.size = 0.1)
ggsave(file.path(qc_dir, "seurat_style", "VlnPlot_condition_overall.png"),
       plot = p_vln_condition_overall,
       width = 18, height = 12, dpi = 300)

# Comparison across selected resolutions
for (qc in qc_metrics) {
  plot_list_seurat <- list()
  
  for (res in c("0.1", "0.5", "1.0")) {
    cluster_col <- get_cluster_col(res, colnames(a@meta.data))
    
    if (!is.null(cluster_col)) {
      Idents(a) <- cluster_col
      plot_list_seurat[[res]] <- VlnPlot(a, features = qc, pt.size = 0.1) +
        ggtitle(paste0("Res ", res))
    }
  }
  
  if (length(plot_list_seurat) > 0) {
    combined_seurat <- wrap_plots(plot_list_seurat, ncol = 3)
    ggsave(file.path(qc_dir, "seurat_style", paste0("VlnPlot_", qc, "_comparison_resolutions.png")),
           plot = combined_seurat,
           width = 24, height = 8, dpi = 300)
  }
}

cat("\nSeurat-style QC violin plots complete!\n")

# ==============================================================================
# SECTION 6: QC SUMMARY STATISTICS TABLES
# ==============================================================================

cat("\n==================================================\n")
cat("SECTION 6: GENERATING QC SUMMARY STATISTICS\n")
cat("==================================================\n\n")

# Function to calculate summary statistics
calc_qc_summary <- function(data, group_var, qc_metrics) {
  data %>%
    group_by(across(all_of(group_var))) %>%
    summarise(
      n_cells = n(),
      across(all_of(qc_metrics), 
             list(mean = ~mean(., na.rm = TRUE),
                  median = ~median(., na.rm = TRUE),
                  sd = ~sd(., na.rm = TRUE),
                  min = ~min(., na.rm = TRUE),
                  max = ~max(., na.rm = TRUE)),
             .names = "{.col}_{.fn}")
    ) %>%
    ungroup()
}

# QC summary by condition
qc_summary_condition <- calc_qc_summary(a@meta.data, "condition", qc_metrics)
write.csv(qc_summary_condition, 
          file.path(tables_dir, "QC", "QC_summary_by_condition.csv"), 
          row.names = FALSE)

# QC summary by singleR labels
qc_summary_singleR <- calc_qc_summary(a@meta.data, "singleR_labels", qc_metrics)
write.csv(qc_summary_singleR, 
          file.path(tables_dir, "QC", "QC_summary_by_singleR.csv"), 
          row.names = FALSE)

# QC summary by singleR pruned
qc_summary_singleR_pruned <- calc_qc_summary(a@meta.data, "singleR_labels_pruned_labeled", qc_metrics)
write.csv(qc_summary_singleR_pruned, 
          file.path(tables_dir, "QC", "QC_summary_by_singleR_pruned.csv"), 
          row.names = FALSE)

# QC summary by harmony clusters for each resolution
for (res in resolutions) {
  cluster_col <- get_cluster_col(res, colnames(a@meta.data))
  
  if (!is.null(cluster_col)) {
    qc_summary_cluster <- calc_qc_summary(a@meta.data, cluster_col, qc_metrics)
    write.csv(qc_summary_cluster, 
              file.path(tables_dir, "QC", paste0("QC_summary_by_harmony_sct_clusters_res", res, ".csv")), 
              row.names = FALSE)
    
    cat(paste0("  Saved QC summary for harmony clusters resolution ", res, "\n"))
  }
}

cat("\nQC summary statistics complete!\n")

# ==============================================================================
# SECTION 7: SAVE UPDATED SEURAT OBJECT
# ==============================================================================

cat("\n==================================================\n")
cat("SAVING UPDATED SEURAT OBJECT\n")
cat("==================================================\n\n")

output_filename <- "jk.18.samples.integrate.Seurat5Standard.multiple.cluster.res.normalized_RNA.singleR+.harmony.rds"

cat("Saving Seurat object to:", output_filename, "\n")
cat("Object contains the new metadata column: singleR_labels_pruned_labeled\n")

# saveRDS(a, file = output_filename)

cat("Seurat object saved successfully!\n")

# ==============================================================================
# SECTION 8: FINAL SUMMARY AND COMPLETION
# ==============================================================================

cat("\n\n")
cat("################################################################################\n")
cat("#                                                                              #\n")
cat("#                    ALL ANALYSES COMPLETED SUCCESSFULLY!                      #\n")
cat("#                         HARMONY VERSION                                      #\n")
cat("#                                                                              #\n")
cat("################################################################################\n\n")

cat("==================================================\n")
cat("OUTPUT SUMMARY\n")
cat("==================================================\n\n")

cat("INPUT FILE:\n")
cat("  - jk.18.samples.integrate.SeuratIntegrate.28nov.multiple.cluster.res.normalized_RNA.singleR.rds\n\n")

cat("OUTPUT FILE:\n")
cat("  - jk.18.samples.integrate.SeuratIntegrate.28nov.multiple.cluster.res.normalized_RNA.singleR+.harmony.rds\n")
cat("    (Contains new column: singleR_labels_pruned_labeled)\n\n")

cat("OUTPUT DIRECTORY STRUCTURE:\n")
cat(paste0("  ", output_dir, "/\n"))
cat("    ├── plots/\n")
cat("    │   ├── by_resolution/ (UMAP plots for each resolution)\n")
cat("    │   └── [comparison plots] - ALL WITH LEGENDS\n")
cat("    ├── qc/\n")
cat("    │   ├── by_resolution/ (custom QC violin plots per resolution)\n")
cat("    │   ├── seurat_style/ (Seurat-style QC violin plots)\n")
cat("    │   └── [overall QC comparison plots]\n")
cat("    └── tables/\n")
cat("        ├── comprehensive/ (3 comprehensive tables per resolution)\n")
cat("        └── QC/ (QC summary statistics)\n\n")

cat("INTEGRATION METHOD: Harmony\n")
cat("UMAP REDUCTION: umap.harmony\n")
cat("CLUSTER COLUMNS: harmony_sct_clusters_res* (preferred for 0.1-0.5 and 1.0), harmony_clusters_res* (fallback for 0.6-0.9)\n\n")

cat("QC METRICS ANALYZED:\n")
cat("  - nCount_RNA\n")
cat("  - nFeature_RNA\n")
cat("  - percent.mt\n")
cat("  - percent.ribo\n")
cat("  - percent.hb\n\n")

cat("RESOLUTIONS PROCESSED:\n")
cat("  - 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0\n\n")

cat("GROUPING VARIABLES:\n")
cat("  - Harmony clusters (per resolution)\n")
cat("  - singleR_labels\n")
cat("  - singleR_labels_pruned_labeled (NEW)\n")
cat("  - condition\n\n")

cat("==================================================\n")
cat("SCRIPT EXECUTION COMPLETE\n")
cat("==================================================\n")
cat("\nAll outputs saved successfully!\n")
cat("ALL DIMPLOTS NOW INCLUDE LEGENDS!\n")
cat("Thank you for using this comprehensive analysis pipeline.\n\n")

################################################################################
#                           END OF SCRIPT                                      #
################################################################################

