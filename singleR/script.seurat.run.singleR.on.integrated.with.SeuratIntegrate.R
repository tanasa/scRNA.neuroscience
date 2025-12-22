library(Seurat)
library(Matrix)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(data.table)
library(cowplot)
library(patchwork)
library(fs)
library(grid)  # Required for unit() function

# Load SingleCellExperiment packages EARLY
# Required for colData(), reducedDimNames(), etc.
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SingleR)
  library(scuttle)
  library(scater)
  library(scran)
  library(viridis)
  library(BiocParallel)  # For SerialParam in normalization
})

packageVersion("Seurat")
set.seed(1337)

# ========================================
# RAM MONITORING CONFIGURATION
# ========================================

RAM_THRESHOLD_GB <- 650  # Stop script if RAM usage exceeds this (in GB)
RAM_CHECK_ENABLED <- TRUE  # Set to FALSE to disable RAM monitoring

# Check if /proc/meminfo exists (Linux only)
if (!file.exists("/proc/meminfo")) {
  RAM_CHECK_ENABLED <- FALSE
  warning("⚠ RAM monitoring disabled: /proc/meminfo not found (not on Linux).")
  cat("RAM monitoring is only available on Linux systems.\n")
  cat("Script will continue without RAM checks.\n\n")
}

# Function to check current RAM usage
check_ram_usage <- function() {
  if (!RAM_CHECK_ENABLED) return(NULL)
  
  # Get memory info from /proc/meminfo (Linux)
  tryCatch({
    mem_info <- readLines("/proc/meminfo")
    
    # Extract total and available memory
    mem_total_line <- grep("^MemTotal:", mem_info, value = TRUE)
    mem_avail_line <- grep("^MemAvailable:", mem_info, value = TRUE)
    
    # Parse values (in kB)
    mem_total_kb <- as.numeric(sub(".*:\\s+(\\d+)\\s+kB", "\\1", mem_total_line))
    mem_avail_kb <- as.numeric(sub(".*:\\s+(\\d+)\\s+kB", "\\1", mem_avail_line))
    
    # Convert to GB
    mem_total_gb <- mem_total_kb / (1024^2)
    mem_avail_gb <- mem_avail_kb / (1024^2)
    mem_used_gb <- mem_total_gb - mem_avail_gb
    
    return(list(
      total_gb = mem_total_gb,
      used_gb = mem_used_gb,
      available_gb = mem_avail_gb,
      pct_used = round(100 * mem_used_gb / mem_total_gb, 1)
    ))
  }, error = function(e) {
    warning("Could not read RAM info: ", e$message)
    return(NULL)
  })
}

# Function to check and stop if RAM threshold exceeded
check_ram_and_stop <- function(context = "Processing", force_gc = TRUE) {
  if (!RAM_CHECK_ENABLED) return(invisible(NULL))
  
  # Force garbage collection to free up memory
  if (force_gc) gc(verbose = FALSE)
  
  ram <- check_ram_usage()
  
  # Handle case where RAM check failed
  if (is.null(ram)) {
    cat(sprintf("[RAM CHECK] %s: Unable to check RAM\n", context))
    return(invisible(NULL))
  }
  
  cat(sprintf("[RAM CHECK] %s: %.1f GB used (%.1f%%), %.1f GB available\n", 
              context, ram$used_gb, ram$pct_used, ram$available_gb))
  
  if (ram$used_gb > RAM_THRESHOLD_GB) {
    cat("\n")
    cat("╔════════════════════════════════════════════════════════════╗\n")
    cat("║  🛑 CRITICAL: RAM USAGE EXCEEDED THRESHOLD                ║\n")
    cat("╚════════════════════════════════════════════════════════════╝\n")
    cat("\n")
    cat(sprintf("Current RAM usage: %.1f GB (%.1f%%)\n", ram$used_gb, ram$pct_used))
    cat(sprintf("RAM threshold:     %.1f GB\n", RAM_THRESHOLD_GB))
    cat(sprintf("Total RAM:         %.1f GB\n", ram$total_gb))
    cat("\n")
    cat("STOPPING SCRIPT FOR SAFETY\n")
    cat("\n")
    cat("To resume:\n")
    cat("1. Wait for RAM to be released\n")
    cat("2. Rerun this script (checkpoints will be used to skip completed steps)\n")
    cat("3. Or increase RAM_THRESHOLD_GB in the script\n")
    cat("4. Or delete checkpoint files to restart from scratch\n")
    cat("\n")
    
    # Try to save current state before stopping
    cat("Attempting to save current progress...\n")
    tryCatch({
      if (exists("seurat_obj") && exists("output_dir")) {
        emergency_file <- file.path(output_dir, "jk.18.samples.integrate.SeuratStandard.SCTransform.multiple.cluster.res.normalized_RNA.singleR.EMERGENCY.rds")
        saveRDS(seurat_obj, emergency_file)
        cat("✓ Saved seurat_obj to emergency backup\n")
        cat("  File:", emergency_file, "\n")
      }
    }, error = function(e) {
      cat("Warning: Could not save emergency backup\n")
    })
    
    cat("\n")
    stop("RAM threshold exceeded. Script stopped for safety.")
  }
  
  return(invisible(ram))
}

if (RAM_CHECK_ENABLED) {
  cat("=== RAM MONITORING ENABLED ===\n")
  cat(sprintf("RAM threshold: %d GB\n", RAM_THRESHOLD_GB))
  cat(sprintf("Current RAM status:\n"))
  initial_ram <- check_ram_usage()
  
  if (!is.null(initial_ram)) {
    cat(sprintf("  Total:     %.1f GB\n", initial_ram$total_gb))
    cat(sprintf("  Used:      %.1f GB (%.1f%%)\n", initial_ram$used_gb, initial_ram$pct_used))
    cat(sprintf("  Available: %.1f GB\n", initial_ram$available_gb))
    cat("\n")
    
    if (initial_ram$used_gb > RAM_THRESHOLD_GB * 0.9) {
      warning("⚠️  Warning: RAM usage is already high (>90% of threshold) before starting!")
      cat("Consider closing other applications before continuing.\n\n")
    }
  }
} else {
  cat("=== RAM MONITORING DISABLED ===\n\n")
}

# --- Jupyter inline rendering presets ---
options(repr.plot.width = 8, repr.plot.height = 6, repr.plot.res = 160)
options(bitmapType = "cairo")
theme_set(theme_classic(base_size = 14))
update_geom_defaults("point", list(size = 1.2, alpha = 0.8))
update_geom_defaults("text",  list(size = 4))

getwd()

# ========================================
# CONFIGURATION
# ========================================

# Sample ID
sample_id <- "jk.18.samples.integrate.SeuratIntegrate.28nov.multiple.cluster.res.normalized_RNA"

cat("=== SAMPLE ID ===\n")
cat("Sample ID:", sample_id, "\n\n")

# Base directory
base_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_seurat.BT.analysis.sep30.nov21.18samples"

# Define output directory EARLY (before any processing)
# This ensures emergency saves work even if script crashes early
output_dir <- file.path(base_dir, paste0("singleR_annotation_", sample_id))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Created output folder:\n", output_dir, "\n\n")
} else {
  cat("Output folder already exists:\n", output_dir, "\n\n")
}

# Helper function for output filenames
outfn <- function(name) file.path(output_dir, paste0(sample_id, "_", name))

cat("Output directory:", output_dir, "\n")
cat("Files will be saved with prefix:", sample_id, "\n\n")

# ========================================
# Load MTG-SEA Reference
# ========================================

mtgsea <- "/mnt/nfs/CX000008_DS1/projects/btanasa/brain_refs_MTG_SEA/brain_refs_MTG_SEA.from.matrix.metadata.sce.rds"
mtg_sea <- readRDS(mtgsea)

cat("=== MTG-SEA REFERENCE ===\n")
cat("Class:", class(mtg_sea), "\n")
cat("Dimensions:", nrow(mtg_sea), "genes x", ncol(mtg_sea), "cells\n\n")

cat("Available annotation columns:\n")
print(colnames(colData(mtg_sea)))
cat("\n")

unique_clusters    <- unique(mtg_sea$cluster_label)
unique_classes     <- unique(mtg_sea$class_label)
unique_subclasses  <- unique(mtg_sea$subclass_label)

cat("cluster_label:", length(unique_clusters), "unique\n")
cat("subclass_label:", length(unique_subclasses), "unique\n")
cat("class_label:", length(unique_classes), "unique\n\n")

cat("Dimensionality reductions in reference:\n")
print(reducedDimNames(mtg_sea))
cat("\n")

# ========================================
# Load 18-Sample Seurat Object
# ========================================

input_filename <- paste0(sample_id, ".rds")
input_file <- file.path(base_dir, input_filename)

cat("=== INPUT FILE ===\n")
cat("Input file path:\n", input_file, "\n\n")

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

setwd(base_dir)
cat("Working directory set to:\n", getwd(), "\n\n")

cat("Reading Seurat object...\n")
seurat_obj <- readRDS(input_file)
cat("Seurat object loaded successfully!\n\n")

gc(verbose = FALSE)
check_ram_and_stop("After loading Seurat object")

cat("=== SEURAT OBJECT INFO ===\n")
cat("Class:", class(seurat_obj), "\n")
cat("Dimensions:", nrow(seurat_obj), "genes x", ncol(seurat_obj), "cells\n")
cat("Assays:", paste(names(seurat_obj@assays), collapse = ", "), "\n")
cat("Reductions:", paste(names(seurat_obj@reductions), collapse = ", "), "\n")
cat("Metadata columns:", ncol(seurat_obj@meta.data), "\n\n")

# ========================================
# Convert Seurat to SingleCellExperiment
# ========================================

cat("\n=== Converting Seurat to SingleCellExperiment ===\n")

sce_obj <- as.SingleCellExperiment(seurat_obj, assay = "RNA")

cat("Conversion successful!\n")
cat("SCE dimensions:", nrow(sce_obj), "genes x", ncol(sce_obj), "cells\n")

check_ram_and_stop("After Seurat → SCE conversion")

cat("\nAssays in SCE:\n")
print(assayNames(sce_obj))
cat("\n")

if (!"counts" %in% assayNames(sce_obj)) {
  stop("ERROR: 'counts' assay not found in SCE object. Cannot proceed with normalization.")
}

# Verify counts are raw
sample_counts <- assay(sce_obj, "counts")[1:100, 1:10]
if (any(sample_counts < 0, na.rm = TRUE)) {
  warning("⚠ Negative values found in counts - may not be raw counts!")
}

max_count <- max(sample_counts, na.rm = TRUE)
if (max_count < 10) {
  warning("⚠ Maximum count value is very low (<10) - counts may already be normalized!")
  cat("  Max count in sample:", max_count, "\n")
  cat("  If counts are already log-normalized, SingleR results may be poor.\n")
}

cat("✓ Counts assay verified (max value in sample:", max_count, ")\n")
cat("  First few count values:\n")
print(sample_counts[1:5, 1:3])
cat("\n")

# Free memory
rm(sample_counts)
gc(verbose = FALSE)

if ("logcounts" %in% assayNames(sce_obj)) {
  cat("✓ Logcounts assay already present\n")
} else {
  cat("⚠ Logcounts not found, will be computed\n")
}

cat("\nMetadata columns transferred:", ncol(colData(sce_obj)), "\n")
cat("First few metadata columns:\n")
print(head(colnames(colData(sce_obj)), 10))
cat("\n")

cat("Dimensionality reductions:\n")
print(reducedDimNames(sce_obj))
cat("\n")

# ========================================
# Normalize Both Objects
# ========================================

cat("=== Ensuring Log-Normalization ===\n")

if (!"logcounts" %in% assayNames(sce_obj)) {
  cat("Computing logcounts for sce_obj...\n")
  sce_obj <- scuttle::logNormCounts(sce_obj, BPPARAM = BiocParallel::SerialParam())
  cat("✓ logcounts computed for sce_obj\n")
} else {
  cat("✓ logcounts already present in sce_obj\n")
}

if (!"logcounts" %in% assayNames(mtg_sea)) {
  cat("Computing logcounts for mtg_sea...\n")
  mtg_sea <- scuttle::logNormCounts(mtg_sea, BPPARAM = BiocParallel::SerialParam())
  cat("✓ logcounts computed for mtg_sea\n")
} else {
  cat("✓ logcounts already present in mtg_sea\n")
}

cat("\n")

# ========================================
# Find Common Genes and Subset
# ========================================

cat("=== Finding Common Genes ===\n")

sce_ids <- sub("\\.\\d+$", "", rownames(sce_obj))
mtg_ids <- sub("\\.\\d+$", "", rownames(mtg_sea))

if (anyDuplicated(sce_ids)) {
  n_dup <- length(unique(sce_ids[duplicated(sce_ids)]))
  warning(sprintf("⚠ %d duplicate genes in sce_obj after stripping versions", n_dup))
  cat("Note: If duplicates exist, only the first occurrence will be kept.\n")
  cat("Consider collapsing duplicates (sum counts) for better results.\n")
}

if (anyDuplicated(mtg_ids)) {
  n_dup <- length(unique(mtg_ids[duplicated(mtg_ids)]))
  warning(sprintf("⚠ %d duplicate genes in mtg_sea after stripping versions", n_dup))
  cat("Note: If duplicates exist, only the first occurrence will be kept.\n")
}

common_genes <- intersect(sce_ids, mtg_ids)

cat("sce_obj genes:", length(sce_ids), "\n")
cat("mtg_sea genes:", length(mtg_ids), "\n")
cat("Common genes:", length(common_genes), "\n\n")

if (length(common_genes) < 1000) {
  warning("⚠ Very few common genes found (<1000). Results may be poor.")
}

# Safer subsetting pattern
sce_sub <- sce_obj[sce_ids %in% common_genes, ]
mtg_sub <- mtg_sea[mtg_ids %in% common_genes, ]

ord <- common_genes
sce_sub <- sce_sub[match(ord, sub("\\.\\d+$", "", rownames(sce_sub))), ]
mtg_sub <- mtg_sub[match(ord, sub("\\.\\d+$", "", rownames(mtg_sub))), ]

rownames(sce_sub) <- ord
rownames(mtg_sub) <- ord

if (!identical(rownames(sce_sub), rownames(mtg_sub))) {
  stop("ERROR: Failed to align gene names between sce_sub and mtg_sub")
}

if (any(is.na(rownames(sce_sub))) || any(is.na(rownames(mtg_sub)))) {
  stop("ERROR: NA values in rownames after subsetting")
}

cat("✓ Objects subsetted and aligned\n")
cat("Subset dimensions:\n")
cat("  sce_sub:", nrow(sce_sub), "genes x", ncol(sce_sub), "cells\n")
cat("  mtg_sub:", nrow(mtg_sub), "genes x", ncol(mtg_sub), "cells\n\n")

gc(verbose = FALSE)

# ========================================
# Run SingleR with Checkpointing
# ========================================

cat("=== Running SingleR ===\n")
cat("Using subclass_label from MTG-SEA reference\n")
cat("Method: wilcox\n\n")

check_ram_and_stop("Before SingleR annotation")

pred_file <- outfn("singleR_predictions.rds")

if (file.exists(pred_file)) {
  cat("✓ Found existing SingleR predictions, loading from checkpoint...\n")
  cat("  File:", pred_file, "\n")
  pred2 <- readRDS(pred_file)
  cat("✓ Loaded", length(pred2$labels), "predictions\n\n")
  cat("To rerun SingleR from scratch, delete this file and rerun the script.\n\n")
} else {
  cat("No checkpoint found, running SingleR (this may take a while)...\n\n")
  
  start_time <- Sys.time()
  
  pred2 <- SingleR(
    test   = sce_sub,
    ref    = mtg_sub,
    labels = colData(mtg_sub)$subclass_label,
    de.method = 'wilcox'
  )
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  
  cat("✓ SingleR completed in", round(elapsed, 2), "minutes\n\n")
  
  cat("Saving checkpoint...\n")
  saveRDS(pred2, pred_file)
  cat("✓ Checkpoint saved to:", pred_file, "\n\n")
}

# Free memory after SingleR
rm(mtg_sub)
gc(verbose = FALSE)

check_ram_and_stop("After SingleR annotation")

cat("=== SingleR Results Summary ===\n")
cat("Total cells annotated:", ncol(sce_sub), "\n")
cat("Unique labels:", length(unique(pred2$labels)), "\n\n")

cat("Top 10 most abundant cell types:\n")
print(head(sort(table(pred2$labels), decreasing = TRUE), 10))
cat("\n")

n_pruned <- sum(is.na(pred2$pruned.labels))
pct_pruned <- round(100 * n_pruned / length(pred2$pruned.labels), 2)
cat("Pruned labels (low confidence):", n_pruned, "(", pct_pruned, "%)\n\n")

cat("Saving SingleR predictions...\n")
pred2_df <- data.frame(
  cell_id = rownames(pred2),
  labels = pred2$labels,
  pruned_labels = pred2$pruned.labels,
  delta_next = pred2$delta.next,
  max_score = apply(pred2$scores, 1, max),
  stringsAsFactors = FALSE
)

write.csv(pred2_df, outfn("singleR_predictions.csv"), row.names = FALSE)
cat("✓ Predictions saved\n\n")

# ========================================
# SingleR Annotation Diagnostics
# ========================================

cat("\n=== SingleR Annotation Diagnostics ===\n\n")

current_pred <- pred2

cat("--- 1. Score Heatmap ---\n")
pdf(outfn("singleR_score_heatmap.pdf"), width = 12, height = 10)
plotScoreHeatmap(current_pred)
dev.off()
cat("Saved:", basename(outfn("singleR_score_heatmap.pdf")), "\n\n")

cat("--- 2. Delta Distribution ---\n")
pdf(outfn("singleR_delta_distribution.pdf"), width = 14, height = 40)
plotDeltaDistribution(current_pred, ncol = 4, dots.on.top = FALSE)
dev.off()
cat("Saved:", basename(outfn("singleR_delta_distribution.pdf")), "\n\n")

cat("--- 3. Pruned Labels Analysis ---\n")
n_total <- length(current_pred$pruned.labels)
n_pruned <- sum(is.na(current_pred$pruned.labels))
n_kept <- sum(!is.na(current_pred$pruned.labels))
pct_pruned <- round(100 * n_pruned / n_total, 2)

cat(sprintf("Total cells: %d\n", n_total))
cat(sprintf("High-quality assignments: %d (%.2f%%)\n", n_kept, 100 - pct_pruned))
cat(sprintf("Low-quality/ambiguous (pruned): %d (%.2f%%)\n\n", n_pruned, pct_pruned))

cat("--- 4. Assignment Quality Summary ---\n")
deltas <- current_pred$delta.next
labels <- current_pred$labels

delta_summary <- data.frame(
  label = unique(labels),
  n_cells = as.numeric(table(labels)[unique(labels)]),
  median_delta = sapply(unique(labels), function(l) median(deltas[labels == l], na.rm = TRUE)),
  mean_delta = sapply(unique(labels), function(l) mean(deltas[labels == l], na.rm = TRUE)),
  pct_pruned = sapply(unique(labels), function(l) {
    idx <- which(labels == l)
    100 * sum(is.na(current_pred$pruned.labels[idx])) / length(idx)
  })
)

delta_summary <- delta_summary[order(-delta_summary$n_cells), ]
rownames(delta_summary) <- NULL

cat("\nAssignment quality by cell type (top 20):\n")
print(head(delta_summary, 20))

write.csv(delta_summary, file = outfn("singleR_assignment_quality_summary.csv"), row.names = FALSE)
cat("\nSaved:", basename(outfn("singleR_assignment_quality_summary.csv")), "\n\n")

cat("--- 5. Marker Gene Heatmap (Top Cell Types) ---\n")
top_labels <- head(names(sort(table(current_pred$labels), decreasing = TRUE)), 10)

for (label in top_labels) {
  safe_label <- gsub("[/\\\\:*?\"<>| ]", "_", label)
  pdf_name <- paste0("singleR_markers_", safe_label, ".pdf")
  
  tryCatch({
    pdf(outfn(pdf_name), width = 10, height = 8)
    plotMarkerHeatmap(current_pred, label = label, test = sce_sub)
    dev.off()
    cat(sprintf("  Saved: %s\n", pdf_name))
  }, error = function(e) {
    cat(sprintf("  Warning: Could not generate marker heatmap for '%s'\n", label))
  })
}

cat("\n--- Diagnostics Complete ---\n\n")

gc(verbose = FALSE)

# ========================================
# Transfer Annotations to Full SCE Object
# ========================================

cat("=== Transferring Annotations to Full Object ===\n")

all_labels <- rep(NA_character_, ncol(sce_obj))
all_labels_pruned <- rep(NA_character_, ncol(sce_obj))
names(all_labels) <- colnames(sce_obj)
names(all_labels_pruned) <- colnames(sce_obj)

all_labels[colnames(sce_sub)] <- pred2$labels
all_labels_pruned[colnames(sce_sub)] <- pred2$pruned.labels

sce_obj$singleR_labels <- all_labels
sce_obj$singleR_labels_pruned <- all_labels_pruned

cat("Total cells in sce_obj:", ncol(sce_obj), "\n")
cat("Cells with annotations:", sum(!is.na(all_labels)), "\n")
cat("Cells with high-confidence annotations:", sum(!is.na(all_labels_pruned)), "\n")
cat("Cells removed by pruning:", sum(!is.na(all_labels)) - sum(!is.na(all_labels_pruned)), "\n\n")

# ========================================
# Generate UMAP Visualizations
# ========================================

cat("=== Generating UMAP Visualizations ===\n")

check_ram_and_stop("Before UMAP generation")

if (!"UMAP" %in% reducedDimNames(sce_obj)) {
  cat("Computing UMAP...\n")
  if (!"PCA" %in% reducedDimNames(sce_obj)) {
    cat("  Computing PCA first...\n")
    sce_obj <- runPCA(sce_obj, ncomponents = 50)
  }
  sce_obj <- runUMAP(sce_obj, dimred = "PCA")
  cat("✓ UMAP computed\n\n")
}

umap_coords <- reducedDim(sce_obj, "UMAP")
umap_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  labels = colData(sce_obj)$singleR_labels,
  labels_pruned = colData(sce_obj)$singleR_labels_pruned
)

centroids <- umap_df %>%
  filter(!is.na(labels)) %>%
  group_by(labels) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

centroids_pruned <- umap_df %>%
  filter(!is.na(labels_pruned)) %>%
  group_by(labels_pruned) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

n_types <- length(unique(na.omit(colData(sce_obj)$singleR_labels)))
n_types_pruned <- length(unique(na.omit(colData(sce_obj)$singleR_labels_pruned)))
mycols <- viridis(n_types, option = "turbo")
mycols_pruned <- viridis(n_types_pruned, option = "turbo")

cat("Creating UMAP plots...\n")

# Plot 1: All labels
pdf(outfn("UMAP_singleR_labels.pdf"), width = 12, height = 12)
p1 <- plotReducedDim(sce_obj, dimred = "UMAP", colour_by = "singleR_labels") +
  scale_color_manual(values = mycols, na.value = "grey80") +
  geom_text(data = centroids, aes(x = UMAP1, y = UMAP2, label = labels),
            color = "black", size = 3, fontface = "bold", check_overlap = TRUE, inherit.aes = FALSE) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.box = "horizontal", legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "lines"), plot.margin = margin(10, 10, 10, 10)) +
  guides(colour = guide_legend(nrow = 4, byrow = TRUE, override.aes = list(size = 3), title.position = "top")) +
  labs(title = "SingleR Subclass Annotation (All Predictions)", colour = "Cell Type")
print(p1)
dev.off()
cat("  Saved: UMAP_singleR_labels.pdf\n")

# Plot 2: Pruned labels (high confidence)
pdf(outfn("UMAP_singleR_labels_pruned.pdf"), width = 12, height = 12)
p2 <- plotReducedDim(sce_obj, dimred = "UMAP", colour_by = "singleR_labels_pruned") +
  scale_color_manual(values = mycols_pruned, na.value = "grey80") +
  geom_text(data = centroids_pruned, aes(x = UMAP1, y = UMAP2, label = labels_pruned),
            color = "black", size = 3, fontface = "bold", check_overlap = TRUE, inherit.aes = FALSE) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.box = "horizontal", legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "lines"), plot.margin = margin(10, 10, 10, 10)) +
  guides(colour = guide_legend(nrow = 4, byrow = TRUE, override.aes = list(size = 3), title.position = "top")) +
  labs(title = "SingleR Subclass Annotation (High Confidence Only)", colour = "Cell Type")
print(p2)
dev.off()
cat("  Saved: UMAP_singleR_labels_pruned.pdf\n")

cat("\n✓ UMAP visualizations saved\n\n")

label_counts <- table(sce_obj$singleR_labels)
label_counts_df <- data.frame(cell_type = names(label_counts), count = as.numeric(label_counts))
label_counts_df <- label_counts_df[order(label_counts_df$count, decreasing = TRUE), ]

cat("Top 30 cell types by count:\n")
print(head(label_counts_df, 30))
cat("\n")

write.csv(label_counts_df, file = outfn("label_counts.csv"), row.names = FALSE)

# ========================================
# Save Annotated Seurat Object with Checkpointing
# ========================================

cat("\n=== Saving Annotated Seurat Object ===\n")

output_filename <- "jk.18.samples.integrate.SeuratStandard.SCTransform.multiple.cluster.res.normalized_RNA.singleR.rds"
output_seurat <- file.path(output_dir, output_filename)

# Checkpoint: only save if file doesn't exist
if (file.exists(output_seurat)) {
  cat("✓ Annotated Seurat object already exists, skipping save\n")
  cat("  File:", output_seurat, "\n")
  cat("  To regenerate, delete this file and rerun the script.\n\n")
} else {
  cat("Transferring SingleR annotations to Seurat object...\n")
  singleR_cols <- grep("singleR", colnames(colData(sce_obj)), value = TRUE)
  
  for (col in singleR_cols) {
    seurat_obj@meta.data[colnames(seurat_obj), col] <- colData(sce_obj)[colnames(seurat_obj), col]
  }
  cat("Transferred", length(singleR_cols), "SingleR annotation columns\n")
  cat("Columns:", paste(singleR_cols, collapse = ", "), "\n\n")
  
  saveRDS(seurat_obj, file = output_seurat)
  cat("Seurat object saved to:\n", output_seurat, "\n")
  cat("File size:", format(file.size(output_seurat), units = "auto"), "\n\n")
  
  check_ram_and_stop("After saving Seurat object")
  
  cat("Verifying save...\n")
  test_seurat <- readRDS(output_seurat)
  cat("✓ Seurat object can be loaded successfully\n")
  cat("  Dimensions:", nrow(test_seurat), "genes x", ncol(test_seurat), "cells\n")
  cat("  SingleR columns:", length(grep("singleR", colnames(test_seurat@meta.data), value = TRUE)), "\n")
  
  rm(test_seurat)
  gc(verbose = FALSE)
  cat("\n")
}

# ========================================
# Final Summary
# ========================================

cat("\n=== ANNOTATION COMPLETE ===\n\n")

cat("Input file:\n  ", input_file, "\n")
cat("Output directory:\n  ", output_dir, "\n\n")

cat("Output file:\n")
cat("  ", output_seurat, "\n\n")

cat("To load the annotated object in R:\n")
cat("  seurat <- readRDS('", output_seurat, "')\n", sep = "")
cat("  table(seurat@meta.data$singleR_labels)\n\n")

cat("✓ ALL DONE!\n\n")

# Final cleanup
gc_result <- gc(verbose = FALSE)
cat("Memory cleaned up\n\n")

# Final RAM check with NULL handling
cat("=== FINAL RAM STATUS ===\n")
final_ram <- check_ram_usage()

if (!is.null(final_ram)) {
  cat(sprintf("Total RAM:     %.1f GB\n", final_ram$total_gb))
  cat(sprintf("Used RAM:      %.1f GB (%.1f%%)\n", final_ram$used_gb, final_ram$pct_used))
  cat(sprintf("Available RAM: %.1f GB\n", final_ram$available_gb))
  cat(sprintf("Peak usage stayed below threshold: %s\n\n", 
              ifelse(final_ram$used_gb <= RAM_THRESHOLD_GB, "✓ YES", "✗ NO")))
} else {
  cat("RAM monitoring disabled or unavailable.\n\n")
}

cat("=== GARBAGE COLLECTION SUMMARY ===\n")
cat(sprintf("Ncells:  %.1f MB used, %.1f MB max\n", gc_result[1,2], gc_result[1,6]))
cat(sprintf("Vcells:  %.1f MB used, %.1f MB max\n\n", gc_result[2,2], gc_result[2,6]))
