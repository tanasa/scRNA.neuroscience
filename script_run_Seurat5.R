# sample : AD_AG_2658
# sample : AD_AG_2658.soupx_corrected_counts.rds
set.seed(12345)
# setwd("/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/JK4488_02_combined/pooled/JK4488_02_cellranger_unfilteredMatrix")

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

sample_id <- "AD_AG_2658"

base_dir <- "/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/fastq_file_Dyslexia_r2_cellRanger_local"
path <- file.path(base_dir, sample_id)

# Show the target path
cat("Target path:\n", path, "\n")

# Ensure the parent path exists before setwd (fail fast if not)
if (!dir.exists(path)) {
  stop("Path does not exist: ", path)
}

# Move into the directory and list files
setwd(path)
cat("Working directory set to:\n", getwd(), "\n\n")
print(list.files())

# Create <sample_id>.soupX inside this path
soupX_dir <- file.path(path, paste0(sample_id, ".soupX"))

if (!dir.exists(soupX_dir)) {
  dir.create(soupX_dir, showWarnings = FALSE)
  cat("\nCreated folder:\n", soupX_dir, "\n")
} else {
  cat("\nFolder already exists:\n", soupX_dir, "\n")
}

output_dir = soupX_dir
outfn <- function(name) file.path(output_dir, paste0(sample_id, "_", name))

species    <- "human"      # "human" or "mouse" (affects mito gene pattern)
use_harmony <- FALSE       # set TRUE to run Harmony on orig.ident
npcs       <- 30
resolution <- 0.5 # it was 1
min_features <- 200
max_mt_pct  <- 25

getwd()

rds_file <- paste0(sample_id, ".soupx_corrected_counts.rds")

if (!file.exists(rds_file)) {
  stop("File not found: ", rds_file, "\nCurrent dir: ", getwd())
}

obj <- readRDS(rds_file)
str(obj)

if (!inherits(obj, "Seurat")) {
  if (inherits(obj, "dgCMatrix") || is.matrix(obj)) {
    message("RDS contains a matrix; creating a Seurat object...")
    obj <- CreateSeuratObject(counts = obj, project = "soupx")
  } else {
    stop("RDS is neither a Seurat object nor a matrix. Class: ",
         paste(class(obj), collapse = ", "))
  }
}
str(obj)

# Checking the soupX values
counts_matrix <- Seurat::GetAssayData(object = obj, assay = "RNA", slot = "counts")
print(counts_matrix[1:5][1:5])

cat("Object class:", class(obj), "\n")
cat("Object structure:\n")
str(obj, max.level = 1)

# print(obj@assays$RNA@layers) 

cat("Cells:", ncol(obj), "| Genes:", nrow(obj), "\n\n")

# Peek at metadata and features
print(head(obj@meta.data, 3))
if ("percent.mt" %in% colnames(obj@meta.data)) {
  cat("\npercent.mt summary:\n"); print(summary(obj@meta.data$percent.mt))
}

mt_pat <- if (tolower(species) == "mouse") "^mt-" else "^MT-"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pat)

print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3))

# Print summary stats (overall)
features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
overall_stats <- sapply(
  features,
  function(f) {
    v <- FetchData(obj, vars = f)[[1]]
    c(
      mean   = mean(v, na.rm = TRUE),
      median = median(v, na.rm = TRUE),
      min    = min(v, na.rm = TRUE),
      max    = max(v, na.rm = TRUE)
    )
  }
)
overall_stats <- data.frame(
  feature = features,
  round(t(overall_stats), 3),
  row.names = NULL
)
cat("\n== Overall summary ==\n"); print(overall_stats)


# --- QC metrics: mito, ribo, hemoglobin (human gene name conventions) ---
# obj[["percent.mt"]]  <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RPL|^RPS")
obj[["percent.hb"]]   <- PercentageFeatureSet(obj, pattern = "^HB[ABDEMGZ]")

VlnPlot(obj, features = c("percent.ribo", "percent.hb", "percent.mt"), pt.size = 0.1, ncol = 3)

# Show sharper plots inline
options(repr.plot.width = 18,  # inches
        repr.plot.height = 6,
        repr.plot.res = 120)   # dpi

plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(obj, feature1 = "nFeature_RNA", feature2 = "percent.mt")

p <- (plot1 + plot2 + plot3) + plot_layout(ncol = 3)
p

ggsave(filename = outfn("qc_scatter_triple.png"), plot = p, width = 18, height = 6, dpi = 300)
cat("Saved to:", outfn("qc_scatter_triple.png"), "\n")

# Filter cells (defaults)
obj <- subset(
  obj,
  subset = nFeature_RNA >= min_features & percent.mt <= max_mt_pct
)

# ---- Normalize → HVG → Scale → PCA ----
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj, features = rownames(obj))
# obj <- ScaleData(obj, features = rownames(obj), vars.to.regress = "percent.mt")

obj <- RunPCA(obj, npcs = max(50, npcs), verbose = FALSE)

png(outfn("pca_elbow.png"), width = 800, height = 600)
print(ElbowPlot(obj, ndims = max(50, npcs)))
dev.off()

options(repr.plot.width = 6,   # inches
        repr.plot.height = 6,
        repr.plot.res = 160)   # dpi
print(ElbowPlot(obj, ndims = max(50, npcs)))

# Take top 10 variable genes
top10 <- head(VariableFeatures(obj), 10)

# Base variable feature plot
plot1 <- VariableFeaturePlot(obj)

# Add labels for the top10 (returns a ggplot)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Show inline
options(repr.plot.width = 8,   # inches
        repr.plot.height = 8,
        repr.plot.res = 160)   # dpi
print(plot2)

ggsave(filename = outfn("top10_variable_genes.png"),
       plot = plot2, width = 7, height = 5, dpi = 300)

# ---- (Optional) Harmony, then UMAP/Neighbors/Clusters ----

red.use <- "pca"

# set use_harmony <- TRUE only if we have >1 sample/batch and you see (or expect) batch effects.
# if (use_harmony) {
#  if (!"orig.ident" %in% colnames(obj@meta.data)) {
#    warning("orig.ident not found; Harmony will be skipped.")
#  } else if (length(unique(obj$orig.ident)) <= 1) {
#    message("Only one sample/batch detected; Harmony skipped.")
#  } else {
#    obj <- RunHarmony(
#      obj,
#      group.by.vars = "orig.ident",
#      reduction     = "pca",
#      dims.use      = 1:npcs,
#      assay.use     = DefaultAssay(obj)
#    )
#    red.use <- "harmony"
#  }
# }

print(red.use)

# Show inline
options(repr.plot.width = 14,   # inches
        repr.plot.height = 8,
        repr.plot.res = 160)   # dpi

# Generate the first two plots
plot_loadings <- VizDimLoadings(obj, dims = 1:2, reduction = red.use)
plot_dim      <- DimPlot(obj, reduction = red.use) + NoLegend()

# Show them side by side
plot_loadings | plot_dim

# DimHeatmap(obj, dims = 1, cells = 500, balanced = TRUE)
plot_heat_pca = DimHeatmap(obj, dims = 1:9, cells = 500, balanced = TRUE)

ggsave(
  filename = outfn("PCA_loadings.png"),
  plot     = plot_loadings,
  width    = 7,
  height   = 5,
  dpi      = 300
)

ggsave(
  filename = outfn("PCA_heatmap.png"),
  plot     = plot_heat_pca,
  width    = 7,
  height   = 5,
  dpi      = 300
)

# obj <- RunUMAP(obj, reduction = red.use, dims = 1:npcs, min.dist = 0.3)

# options(repr.plot.width = 8,   # inches
#        repr.plot.height = 6,
#        repr.plot.res = 160)   # dpi

# DimPlot(obj, reduction = "umap", label = TRUE)
# plot_umap = DimPlot(obj, reduction = "umap", label = TRUE)

# ggsave(
#  filename = outfn("UMAP_pre_cluster.png"),
#  plot     = plot_umap,
#  width    = 7,
#  height   = 5,
#  dpi      = 300
# )

obj <- FindNeighbors(obj, reduction = red.use, dims = 1:npcs, k.param = 15)
obj <- FindClusters(obj, resolution = resolution)

# Size in inches; bump res for crispness
options(repr.plot.width = 7, repr.plot.height = 5, repr.plot.res = 160)

# use t-SNE as a temporary solution
obj <- RunTSNE(obj, 
               reduction = red.use, 
               dims = 1:npcs, 
               perplexity = min(30, ncol(obj)/4),
               verbose = TRUE)

# Plot t-SNE instead
DimPlot(obj, reduction = "tsne", label = TRUE)

plot_tsne = DimPlot(obj, reduction = "tsne", label = TRUE)

ggsave(
  filename = outfn("TSNE.png"),
  plot     = plot_tsne,
  width    = 7,
  height   = 5,
  dpi      = 300
)

obj <- RunUMAP(obj, reduction = red.use, dims = 1:npcs, min.dist = 0.3)

options(repr.plot.width = 8,   # inches
        repr.plot.height = 6,
        repr.plot.res = 160)   # dpi

DimPlot(obj, reduction = "umap", label = TRUE)
plot_umap = DimPlot(obj, reduction = "umap", label = TRUE)

ggsave(
  filename = outfn("UMAP_post_cluster.png"),
  plot     = plot_umap,
  width    = 7,
  height   = 5,
  dpi      = 300
)

print(npcs)
# str(obj)

library(DoubletFinder)

# npcs_use = npcs

# the code from Bioturing : it is OK, but slow
# https://colabdev.bioturing.com/notebook/doubletfinder-doublet-detection-in-singlecell-rna--aa4b6bca01

# ---- Step 1: pK parameter sweep ----
# sweep.res.list_obj <- paramSweep(seu = obj, PCs = 1:15, sct = FALSE)
# sweep.stats_obj <- summarizeSweep(sweep.res.list_obj, GT = FALSE)
# bcmvn_obj <- find.pK(sweep.stats_obj)

# Select optimal pK
# bcmvn.max <- bcmvn_obj[which.max(bcmvn_obj$BCmetric), ]
# optimal.pk <- bcmvn.max$pK
# optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
# print(paste0("Optimal pK selected is ", optimal.pk))

# ---- Step 2: Estimate homotypic proportion ----
# annotations <- Idents(obj)   # use cell type / cluster annotations
# homotypic.prop <- modelHomotypic(annotations)

# ---- Step 3: Expected number of doublets ----
# Use 10x Genomics doublet rate estimation (example: 8%)
# nExp_poi <- round(0.08 * nrow(obj@meta.data))
# nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# ---- Step 4: Run DoubletFinder ----
# obj <- doubletFinder(
#  seu = obj,
#  PCs = 1:15,
#  pK = optimal.pk,
#  nExp = nExp_poi.adj
# )

# ---- Inspect results ----
# head(obj@meta.data, 2)

# Find DoubletFinder classification columns
# doubletfinder_cols <- grep("^DF\\.classifications_", colnames(obj@meta.data), value = TRUE)
# if (!length(doubletfinder_cols)) stop("No DoubletFinder classification columns found.")

# message("DF columns: ", paste(doubletfinder_cols, collapse = ", "))

# If multiple runs exist, pick the one with the largest nExp in the suffix
# nExp_from_name <- function(x) suppressWarnings(as.integer(sub(".*_(\\d+)$", "\\1", x)))
# df_col <- if (length(doubletfinder_cols) > 1) {
#  doubletfinder_cols[which.max(nExp_from_name(doubletfinder_cols))]
# } else {
#  doubletfinder_cols[1]
# }
# message("Using DF column: ", df_col)

# Count singlets vs doublets (correct dynamic indexing)
# print(table(obj@meta.data[[df_col]]))

# Optional: check pANN columns
# pann_cols <- grep("^pANN_", colnames(obj@meta.data), value = TRUE)
# if (length(pann_cols)) message("pANN columns: ", paste(pann_cols, collapse = ", "))

# Plot using DoubletFinder classifications
# doublet_finder <- DimPlot(
#  obj,
#  reduction = "umap",
#  group.by = df_col[1],
#  cols = c("Singlet" = "lightblue", "Doublet" = "red")
# ) + ggtitle(paste0("DoubletFinder: ", df_col))

# ggsave(outfn("DoubletFinder_umap_colored.png"),
#        plot = doublet_finder, width = 7, height = 5, dpi = 300)

# print(doublet_finder)

# Use a different doublet detection method temporarily
library(scDblFinder)

# Convert to SingleCellExperiment for scDblFinder
library(SingleCellExperiment)
sce <- as.SingleCellExperiment(obj)

# Run scDblFinder
sce <- scDblFinder(sce, clusters = obj$seurat_clusters)

# Add results back to Seurat object
obj$doublet_scores <- sce$scDblFinder.score
obj$doublet_class <- sce$scDblFinder.class

# Visualize results
scdblfinder <- DimPlot(
  obj,
  # reduction = "tsne",
  reduction = "umap",
  group.by = "doublet_class",
  cols = c("singlet" = "lightblue", "doublet" = "red")
) + ggtitle("scDblFinder: Singlet vs Doublet")

ggsave(outfn("scDblFinder_tsne.png"),
       plot = scdblfinder,
       width = 7, height = 5, dpi = 300)

plot(scdblfinder)

library(cowplot)
# combined <- plot_grid(doublet_finder , scdblfinder, ncol = 2, labels = c("DoubletFinder", "scdblfinder"))
# ggsave(outfn("doublets_umap.png"), scdblfinder, width = 7, height = 5, dpi = 300)
# plot(combined)

head(obj@meta.data)
table(obj@meta.data$doublet_class)
dim(obj)

obj_clean <- subset(obj, subset = doublet_class == "singlet")
dim(obj_clean)

cat("number of singlet cells:")
table(obj_clean$doublet_class)
cat("number of seurat clusters:")
table(obj_clean$seurat_clusters)

# Recompute embeddings & clusters on clean set

obj_clean <- NormalizeData(obj_clean)
obj_clean <- FindVariableFeatures(obj_clean, nfeatures = 3000)
obj_clean <- ScaleData(obj_clean, features = rownames(obj_clean))

obj_clean <- RunPCA(obj_clean, npcs = 50, verbose = FALSE)
obj_clean <- FindNeighbors(obj_clean, dims = 1:30)
obj_clean <- FindClusters(obj_clean, resolution = resolution)

# obj_clean <- RunUMAP(obj_clean, dims = 1:30)

# ---- t-SNE (Seurat -> Rtsne backend) ----
# Perplexity must be < (n_cells - 1) / 3; this helper picks a safe value.
n_cells <- ncol(obj_clean)
perp <- max(5, floor(min(30, (n_cells - 1) / 3)))

obj_clean <- RunTSNE(
  obj_clean,
  dims = 1:30,                 # run t-SNE on the first 30 PCs
  reduction = "pca",
  check_duplicates = FALSE,    # typical for single-cell matrices
  perplexity = perp,
  seed.use = 1337
)

# Plot & save
p_tsne <- DimPlot(obj_clean, reduction = "tsne", label = TRUE) + ggtitle(sprintf("t-SNE (perplexity=%d)", perp))
print(p_tsne)
ggsave(outfn("TSNE_post_doublets.png"), p_tsne, width = 7, height = 5, dpi = 300)

obj_clean <- RunUMAP(obj_clean, reduction = red.use, dims = 1:npcs, min.dist = 0.3)
# options(repr.plot.width = 8,   # inches
#        repr.plot.height = 6,
#        repr.plot.res = 160)   # dpi

DimPlot(obj_clean, reduction = "umap", label = TRUE)
plot_umap = DimPlot(obj_clean, reduction = "umap", label = TRUE)

ggsave(
  filename = outfn("UMAP_post_doublets.png"),
  plot     = plot_umap,
  width    = 7,
  height   = 5,
  dpi      = 300
)

head(obj_clean@meta.data)
table(obj_clean@meta.data$seurat_clusters)
obj = obj_clean

# Cell cycle:

# https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# s.genes <- cc.genes$s.genes
# print(s.genes)
# g2m.genes <- cc.genes$g2m.genes
# print(g2m.genes)

obj <- CellCycleScoring(
  obj,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)

# Continuous score maps on UMAP

# pS   <- FeaturePlot(obj, features = "S.Score",  reduction = "umap",
#                     min.cutoff = "q05", max.cutoff = "q95") + ggtitle("S.Score")

# pG2M <- FeaturePlot(obj, features = "G2M.Score", reduction = "umap",
#                    min.cutoff = "q05", max.cutoff = "q95") + ggtitle("G2M.Score")

# (pS | pG2M)  # patchwork
# ggsave(outfn("UMAP_S_and_G2M_scores.png"), pS | pG2M, width = 12, height = 6, dpi = 300)

# Continuous score maps on UMAP # same scale

# lims <- range(c(obj$S.Score, obj$G2M.Score), na.rm = TRUE)
#
# pS <- FeaturePlot(obj, features = "S.Score", reduction = "umap",
#                  combine = FALSE)[[1]] +
#      scale_color_gradient(low = "lightgrey", high = "blue", limits = lims) +
#      ggtitle("S.Score")

# pG2M <- FeaturePlot(obj, features = "G2M.Score", reduction = "umap",
#                    combine = FALSE)[[1]] +
#        scale_color_gradient(low = "lightgrey", high = "blue", limits = lims) +
#        ggtitle("G2M.Score")

# (pS | pG2M)  # patchwork side by side

# RunPCA(object = obj, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)

# DimPlot(obj, reduction = "umap", group.by = "Phase", label = TRUE) + ggtitle("Cyclone Phase")
# cell_cycle = DimPlot(obj, reduction = "tsne", group.by = "Phase", label = TRUE) + ggtitle("Cyclone Phase")
# print(cell_cycle)
# ggsave(outfn("cell_cycle_TSNE.png"), p_tsne, width = 7, height = 5, dpi = 300)

cell_cycle = DimPlot(obj, reduction = "umap", group.by = "Phase", label = TRUE) + ggtitle("Cyclone Phase")
print(cell_cycle)
ggsave(outfn("cell_cycle_UMAP.png"), p_tsne, width = 7, height = 5, dpi = 300)

head(obj@meta.data)
table(obj@meta.data$Phase)

# Leveraging functions in sce

library(SingleCellExperiment)

# From counts (genes x cells)
# sce <- SingleCellExperiment(assays = list(counts = counts_mat))  # dgCMatrix ok
sce_obj <- as.SingleCellExperiment(obj)

# Quick peek
sce_obj
assayNames(sce_obj)            # "counts", "logcounts", ...
dim(sce_obj)     
# genes, cells

cat("assays in sce:")
assays(sce_obj)$counts[1:3,1:3]
rowData(sce_obj)               # gene-level metadata (DataFrame)
colData(sce_obj)               # cell-level metadata (DataFrame)
# metadata(sce1)              # list for global info


# Standard accessors:
counts(sce_obj)[1:3, 1:3]      # raw counts
logcounts(sce_obj)[1:3, 1:3]   # log-normalized (if present)

inherits(counts(sce),   "dgCMatrix")    # TRUE if sparse
inherits(logcounts(sce),"dgCMatrix")
inherits(assay(sce,"scaledata"), "dgCMatrix")

dim(sce)                                # genes x cells
dim(counts(sce)); dim(logcounts(sce)); dim(assay(sce,"scaledata"))

head(as.data.frame(colData(sce)), 3)    # cell-level metadata
head(as.data.frame(rowData(sce)), 3)    # gene-level metadata (empty in your printout)

reducedDimNames(sce)                    # "PCA", "TSNE"
reducedDim(sce, "PCA")[1:2, 1:2]
reducedDim(sce, "TSNE")[1:2, 1:2]

# from Seurat to SCE: converstion step by step (OK)

# library(SingleCellExperiment)
# library(Matrix)

# counts <- tryCatch(
#  GetAssayData(obj, assay = DefaultAssay(obj), layer = "counts"),
#  error = function(e) GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts")
# )
# logdata <- tryCatch(
#  GetAssayData(obj, assay = DefaultAssay(obj), layer = "data"),
#  error = function(e) GetAssayData(obj, assay = DefaultAssay(obj), slot = "data")
# )

# sce_new <- SingleCellExperiment(
#  assays = list(
#    counts    = as(counts,  "dgCMatrix"),
#    logcounts = as(logdata, "dgCMatrix")
#  ),
#  colData = DataFrame(obj@meta.data)  # attach all Seurat cell metadata here
# )

# ensure row order matches cells
# rownames(colData(sce_new)) <- colnames(sce_new)

# check
# dim(colData(sce_new))
# head(as.data.frame(colData(sce_new)), 3)
# dim(rowData(sce_new))
# head(as.data.frame(rowData(sce_new)), 3)

# https://www.bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html
library("scater")

options(repr.plot.width = 8, repr.plot.height = 10, repr.plot.res = 160)

p <- scater::plotHighestExprs(sce_obj, exprs_values = "counts", n = 50)
print(p) 
ggsave(outfn("highest_exprs_counts.png"), plot = p, width = 8, height = 6, dpi = 300)

library(scuttle)
library(dplyr)

# 1) QC stats
gene_stats <- perFeatureQCMetrics(sce_obj)

# 2) Averages (vectors)
avg_counts <- calculateAverage(sce_obj, exprs_values = "counts")      # numeric
avg_log    <- if ("logcounts" %in% assayNames(sce_obj))
               calculateAverage(sce_obj, exprs_values = "logcounts") else NULL

# 3) Combine & save
gene_df <- as.data.frame(gene_stats)
gene_df$gene         <- rownames(sce_obj)
gene_df$mean_counts  <- as.numeric(avg_counts)
if (!is.null(avg_log)) gene_df$mean_logcounts <- as.numeric(avg_log)

top50 <- gene_df %>% arrange(desc(mean_counts)) %>% slice_head(n = 50)
head(top50)

write.csv(gene_df, outfn("gene_stats_with_stats_means.csv"), row.names = FALSE)

# ?scuttle::calculateAverage
# ?scuttle::aggregateAcrossCells

library(scater)
library(ggplot2)
library(patchwork)

plot_gene_umap_seurat <- function(sce, gene,
                                  assay_name = "logcounts",
                                  viridis_option = "plasma",   # "viridis","magma","inferno","cividis"
                                  cap_quantiles = c(0.01, 0.99),
                                  add_violin = TRUE,
                                  add_box = TRUE) {
  # --- require seurat_clusters in colData ---
  stopifnot("seurat_clusters" %in% colnames(colData(sce)))
  colData(sce)$seurat_clusters <- as.factor(colData(sce)$seurat_clusters)

  # --- ensure UMAP exists (from PCA if available) ---
  if (!"UMAP" %in% reducedDimNames(sce)) {
    set.seed(12345)
    if (!"PCA" %in% reducedDimNames(sce)) {
      sce <- runPCA(sce, ncomponents = 30, exprs_values = assay_name)
    }
    sce <- runUMAP(sce, dimred = "PCA", name = "UMAP", ncomponents = 2)
  }

  # --- resolve gene rowname (optional symbol mapping) ---
  g <- gene
  if (!g %in% rownames(sce) && "symbol" %in% colnames(rowData(sce))) {
    hit <- which(rowData(sce)$symbol == g)[1]
    if (length(hit)) g <- rownames(sce)[hit]
  }
  stopifnot(g %in% rownames(sce))

  # --- expression & limits (cap outliers for nicer colors) ---
  expr <- assay(sce, assay_name)[g, , drop = TRUE]
  lim  <- if (is.null(cap_quantiles)) range(expr, na.rm = TRUE)
          else quantile(expr, probs = cap_quantiles, na.rm = TRUE)

  # expose expression for ggplot via colData
  sce$.expr <- as.numeric(expr)

  # --- UMAP colored by expression ---
  p1 <- plotReducedDim(
          sce, dimred = "UMAP", colour_by = g, by_exprs_values = assay_name,
          point_size = 0.6
        ) +
        ggtitle(paste0("UMAP: ", gene)) +
        scale_colour_viridis_c(
          option = viridis_option, limits = lim,
          name = paste0(gene, " (", assay_name, ")"),
          na.value = "grey85"
        )

  # --- violin/box by Seurat clusters ---
  p2 <- ggcells(sce, mapping = aes(x = seurat_clusters, y = .expr, fill = seurat_clusters)) +
        theme_classic(base_size = 13) +
        theme(axis.text.x = element_text(angle = 40, hjust = 1),
              legend.position = "none") +
        labs(x = "seurat_clusters", y = paste0(gene, " (", assay_name, ")"), title = gene)

  if (add_violin) {
    p2 <- p2 + geom_violin(scale = "width", trim = TRUE)
  }
  if (add_box) {
    p2 <- p2 + geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6, color = "black")
  }

  p1 | p2
}

# example usage
options(repr.plot.width = 12, repr.plot.height = 5, repr.plot.res = 160)
plot_gene_umap_seurat(sce_obj, "APOE")
plot_gene_umap_seurat(sce_obj, "APP")
plot_gene_umap_seurat(sce_obj, "PSEN1")
plot_gene_umap_seurat(sce_obj, "PSEN2")


# Either form works:
class(sce_obj)
methods::slotNames(sce_obj)

# assayNames(sce_obj); 
# counts(sce_obj); 
# logcounts(sce_obj); 
# assay(sce_obj, "scaledata")
# colData(sce_obj); 
# rowData(sce_obj)
# reducedDimNames(sce_obj); 
# reducedDim(sce_obj, "PCA")
# altExpNames(sce_obj); 
# metadata(sce_obj)

# Other functions from scran : https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html



# Find all cluster markers (positive markers only)
obj.markers.wl <- FindAllMarkers(
  obj,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"  
)

# Filter: keep only strong markers (log2FC > 1)
obj.markers.wl.filtered <- obj.markers.wl %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Save full and filtered marker tables
print(obj.markers.wl.filtered)
cat("number of filtered markers:")
dim(obj.markers.wl.filtered)

write.csv(obj.markers.wl, outfn("all_markers.wilcox.csv"), row.names = FALSE)
write.csv(obj.markers.wl.filtered, outfn("all_markers_log2FC_gt1.wilcox.csv"), row.names = FALSE)

# FindMarkers : differential expression tests : test.use :

# "wilcox" : Identifies differentially expressed genes between two groups of cells using a Wilcoxon Rank Sum test (default); will use a fast implementation by Presto if installed
# "wilcox_limma" : Identifies differentially expressed genes between two groups of cells using the limma implementation of the Wilcoxon Rank Sum test; set this option to reproduce results from Seurat v4
# "bimod" : Likelihood-ratio test for single cell gene expression, (McDavid et al., Bioinformatics, 2013)
# "roc" : Identifies 'markers' of gene expression using ROC analysis. For each gene, evaluates (using AUC) a classifier built on that gene alone, to classify between two groups of cells. An AUC value of 1 means that expression values for this gene alone can perfectly classify the two groupings (i.e. Each of the cells in cells.1 exhibit a higher level than each of the cells in cells.2). An AUC value of 0 also means there is perfect classification, but in the other direction. A value of 0.5 implies that the gene has no predictive power to classify the two groups. Returns a 'predictive power' (abs(AUC-0.5) * 2) ranked matrix of putative differentially expressed genes.
# "t" : Identify differentially expressed genes between two groups of cells using the Student's t-test.
# "negbinom" : Identifies differentially expressed genes between two groups of cells using a negative binomial generalized linear model. Use only for UMI-based datasets
# "poisson" : Identifies differentially expressed genes between two groups of cells using a poisson generalized linear model. Use only for UMI-based datasets
# "LR" : Uses a logistic regression framework to determine differentially expressed genes. Constructs a logistic regression model predicting group membership based on each feature individually and compares this to a null model with a likelihood ratio test.
# "MAST" : Identifies differentially expressed genes between two groups of cells using a hurdle model tailored to scRNA-seq data. Utilizes the MAST package to run the DE testing.
# "DESeq2" : Identifies differentially expressed genes between two groups of cells based on a model using DESeq2 which uses a negative binomial distribution (Love et al, Genome Biology, 2014).

# Find all cluster markers (positive markers only) : MAST test
# obj.markers.mast <- FindAllMarkers(
#  obj,
#  only.pos        = TRUE,
#  min.pct         = 0.25,
#  logfc.threshold = 0.25,
#  test.use = "MAST"  
# )

# Filter: keep only strong markers (log2FC > 1)
# obj.markers.mast.filtered <- obj.markers.mast %>%
#  group_by(cluster) %>%
#  dplyr::filter(avg_log2FC > 1)

# Save full and filtered marker tables
# print(obj.markers.mast.filtered)
# cat("number of filtered MAST markers:")
# dim(obj.markers.mast.filtered)

# write.csv(obj.markers.mast, outfn("all_markers.mast.csv"), row.names = FALSE)
# write.csv(obj.markers.mast.filtered, outfn("all_markers_log2FC_gt1.mast.csv"), row.names = FALSE)

# https://bioconductor.org/packages/devel/bioc/vignettes/glmGamPoi/inst/doc/pseudobulk.html

# Pseudobulk : seurat_clusters

agg <- AggregateExpression(
  obj,
  assays = "RNA",
  group.by = "seurat_clusters",
  slot = "counts",
  return.seurat = TRUE,
  verbose = TRUE
)

# list layers present
Layers(agg[["RNA"]])

# pull the matrix stored in the 'counts' layer
mat_counts <- LayerData(agg[["RNA"]], layer = "counts")
mat_counts[1:10, 1:10]

# Export for DE analysis
# head(agg@assays$RNA@layers$counts, 10) 
write.csv(as.data.frame(mat_counts), outfn("pseudobulk_Seurat_Cluster_computed_by_Seurat_AggeregateExpression.csv"))



library(scuttle)
library(SingleCellExperiment)


# Pseudobulk by seurat_clusters only
pseudobulk_sce <- aggregateAcrossCells(
  sce_obj,
  ids = colData(sce_obj)[, "seurat_clusters"],  # just cluster
  use.assay.type = "counts"
)

colData(pseudobulk_sce)  # shows seurat_clusters for each pseudobulk

# Extract counts
counts_mat <- counts(pseudobulk_sce)
dim(counts_mat)  # genes × clusters
counts_mat[1:10, 1:10]

# Export for DE analysis
# head(agg@assays$RNA@layers$counts, 10) 
write.csv(as.data.frame(counts_mat), outfn("pseudobulk_Seurat_Cluster_computed_by_Sce_aggregateAcrossCell.csv"))

# SAVE SEURAT OBJECT

# outfn <- function(name) file.path(output_dir, paste0(sample_id, "_", name))
saveRDS(obj, outfn("seurat_obj.filtered.soupX.rds"))

# SAVE SCE OBJECT

# outfn <- function(name) file.path(output_dir, paste0(sample_id, "_", name))
saveRDS(sce_obj, outfn("sce_obj.filtered.soupX.rds"))



# ENSEMBL_IDs in sce_obj

# Assume sce_obj2 is your SingleCellExperiment
# sce_obj2 <- sce_obj   # as you wrote
# str(sce_obj) 

# Strip version suffix if present (.12 etc.)
ids <- sub("\\.\\d+$", "", rownames(sce_obj))

# Check which look like Ensembl IDs
is_ensg <- grepl("^ENSG", ids)

# Counts
n_total <- length(ids)
n_ensg  <- sum(is_ensg)
n_symbol <- n_total - n_ensg

cat("Total genes:", n_total, "\n")
cat("Ensembl (ENSG) IDs:", n_ensg, "\n")
cat("Other (likely gene symbols):", n_symbol, "\n")

# Optionally show examples
head(ids[is_ensg])
head(ids[!is_ensg])



# Other functions from scran : https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html

# Detecting correlated genes 
# https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html

# library(scran)
# library(scuttle)

# cat("Detecting correlated genes")

## ensure logcounts exist
# if (!"logcounts" %in% assayNames(sce_obj)) {
#  sce_obj <- logNormCounts(sce_obj)
# }

## HVGs
# dec <- modelGeneVar(sce_obj, assay.type = "logcounts")
# top.hvgs1 <- getTopHVGs(dec, prop = 0.1)
# top.hvgs2 <- getTopHVGs(dec, n = 2000)

# Use HVGs for correlation analysis
# cor.pairs.hvsg1 <- correlatePairs(
#  sce_obj,
#  subset.row = top.hvgs1,  # or top.hvgs2
#  assay.type = "logcounts"
# )
# cor.genes.hvsg1 <- correlateGenes(cor.pairs.hvsg1)
# dim(cor.genes.hvsg1)
# head(cor.genes.hvsg1, 2)

# Use HVGs for correlation analysis
# cor.pairs.hvsg2 <- correlatePairs(
#  sce_obj,
#  subset.row = top.hvgs2,  # or top.hvgs2
#  assay.type = "logcounts"
# )
# cor.genes.hvsg2 <- correlateGenes(cor.pairs.hvsg2)
# dim(cor.genes.hvsg2)
# head(cor.genes.hvsg2, 2)

# HVG lists
# write.table(cor.pairs.hvsg1, outfn("topHVGs_prop0.1.txt"),
#            row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(cor.pairs.hvsg2, outfn("correlatePairs_sce_hvg2000.txt"),
#            row.names = FALSE, col.names = FALSE, quote = FALSE)

# correlatePairs / correlateGenes
# write.csv(as.data.frame(cor.pairs.hvsg1),
#          outfn(sprintf("correlatedGenes_HVSG_prop0.1.csv", k)),
#          row.names = FALSE)
# write.csv(as.data.frame(cor.genes.hvsg2),
#          outfn(sprintf("correlateGenes_sce_hvg2000.csv", k)),
#          row.names = FALSE)
