
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

# To increase the limit (e.g., to 70 GiB)
options(future.globals.maxSize = 700 * 1024^3)
setwd("/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_seurat.BT.analysis.sep29")

# Summary function
show_summary <- function(obj, name) {
  cat(sprintf("\n%s:\n", name))
  cat(sprintf("  Genes: %d\n", nrow(obj)))
  cat(sprintf("  Cells: %d\n", ncol(obj)))
  cat(sprintf("  UMAP: %s\n", ifelse("umap" %in% names(obj@reductions), "Yes", "No")))
}

# https://academic.oup.com/bioinformatics/article/41/6/btaf358/8171806
# a new package SEURAT INTEGRATE
# https://github.com/cbib/Seurat-Integrate/blob/main/vignettes/SeuratIntegrate.Rmd

# Point methods to your *path-based* env
UpdateEnvCache("bbknn",
  conda.env = "/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/SeuratIntegrate",
  conda.env.is.path = TRUE
)
UpdateEnvCache("scvi",
  conda.env = "/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/SeuratIntegrate",
  conda.env.is.path = TRUE
)
UpdateEnvCache("scanorama",
  conda.env = "/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/SeuratIntegrate",
  conda.env.is.path = TRUE
)
UpdateEnvCache("trvae",
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

# if (! reticulate::condaenv_exists('umap_0.5.4')) {
#  reticulate::conda_create('umap_0.5.4', packages = 'umap-learn=0.5.4')
# }

# The following integration method functions are available:
#
#        * 'SeuratIntegrate::CombatIntegration'

#        * 'SeuratIntegrate::HarmonyIntegration'

#        * 'SeuratIntegrate::HarmonyIntegration.fix'

#        * 'SeuratIntegrate::MNNIntegration'

#        * 'SeuratIntegrate::ScanoramaIntegration'

#        * 'SeuratIntegrate::bbknnIntegration'

#        * 'SeuratIntegrate::scANVIIntegration'

#        * 'SeuratIntegrate::scVIIntegration'

#        * 'SeuratIntegrate::scVIIntegration.fix'

#        * 'SeuratIntegrate::trVAEIntegration'

#        * 'CCAIntegration'

#        * 'HarmonyIntegration'

#        * 'JointPCAIntegration'

#        * 'RPCAIntegration'

# ---- Load Individual Seurat Objects ----

# To increase the limit (e.g., to 70 GiB)
options(future.globals.maxSize = 700 * 1024^3)
setwd("/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_seurat.BT.analysis.sep29")

# ---- Load Individual Seurat Objects ----
JK4488_01_seurat <- readRDS("JK4488_01_seurat_obj.filtered.soupX.rds")
JK4488_02_seurat <- readRDS("JK4488_02_seurat_obj.filtered.soupX.rds")
JK4488_03_seurat <- readRDS("JK4488_03_seurat_obj.filtered.soupX.rds")
JK4488_04_seurat <- readRDS("JK4488_04_seurat_obj.filtered.soupX.rds")
JK4488_05_seurat <- readRDS("JK4488_05_seurat_obj.filtered.soupX.rds")
JK4488_06_seurat <- readRDS("JK4488_06_seurat_obj.filtered.soupX.rds")

# ---- Setup conda environments ----
library(SeuratIntegrate)
library(Seurat)
library(reticulate)

# Point methods to your path-based env
UpdateEnvCache("bbknn",
               conda.env = "/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/SeuratIntegrate",
               conda.env.is.path = TRUE
)
UpdateEnvCache("scvi",
               conda.env = "/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/SeuratIntegrate",
               conda.env.is.path = TRUE
)
UpdateEnvCache("scanorama",
               conda.env = "/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/SeuratIntegrate",
               conda.env.is.path = TRUE
)
UpdateEnvCache("trvae",
               conda.env = "/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/SeuratIntegrate",
               conda.env.is.path = TRUE
)

# Verify
getCache()

# Use your specific conda environment
use_condaenv("/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/SeuratIntegrate", required = TRUE)
py_config()

# ---- Add sample IDs and merge ----
JK4488_01_seurat$sample_id <- "JK4488_01"
JK4488_02_seurat$sample_id <- "JK4488_02"
JK4488_03_seurat$sample_id <- "JK4488_03"
JK4488_04_seurat$sample_id <- "JK4488_04"
JK4488_05_seurat$sample_id <- "JK4488_05"
JK4488_06_seurat$sample_id <- "JK4488_06"

# Merge all samples
JK4488_merged <- merge(JK4488_01_seurat,
                       y = c(JK4488_02_seurat, JK4488_03_seurat,
                             JK4488_04_seurat, JK4488_05_seurat,
                             JK4488_06_seurat),
                       add.cell.ids = c("JK4488_01", "JK4488_02", "JK4488_03",
                                        "JK4488_04", "JK4488_05", "JK4488_06"))

# Join layers if already split (to avoid error on re-runs)
if (length(Layers(JK4488_merged)) > 1) {
  JK4488_merged <- JoinLayers(JK4488_merged)
}

# Define batch variable
batch.var <- "sample_id"

# ---- Preprocessing with Seurat ----
# Split by sample for independent normalization
JK4488_merged <- split(JK4488_merged, f = JK4488_merged$sample_id)

# SCTransform normalization
JK4488_merged <- SCTransform(JK4488_merged)

# PCA
JK4488_merged <- RunPCA(JK4488_merged, npcs = 50)

# Find neighbors and run UMAP on unintegrated data
JK4488_merged <- FindNeighbors(JK4488_merged, dims = 1:30, k.param = 20L)
JK4488_merged <- RunUMAP(JK4488_merged, dims = 1:30, n.neighbors = 20)

# Create initial clustering on unintegrated data to use as reference
JK4488_merged <- FindClusters(JK4488_merged, resolution = 0.5)
# Use these clusters as a pseudo cell.var for FindOptimalClusters
cell.var <- "seurat_clusters"

# Visualize unintegrated data
DimPlot(JK4488_merged, group.by = "sample_id", label = FALSE)
DimPlot(JK4488_merged, group.by = cell.var, label = TRUE)

# ---- Integration with SeuratIntegrate ----
JK4488_merged <- DoIntegrate(
  object = JK4488_merged,
  
  # integrations
  SeuratIntegrate::CombatIntegration(),
  SeuratIntegrate::bbknnIntegration(orig = "pca", ndims.use = 30),
  SeuratIntegrate::HarmonyIntegration(orig = "pca", dims = 1:30),
  SeuratIntegrate::ScanoramaIntegration(),
  SeuratIntegrate::scVIIntegration(),
  SeuratIntegrate::trVAEIntegration(),
  
  # additional parameters
  use.hvg = TRUE,
  use.future = TRUE
)


# Save the object immediately after integration
saveRDS(JK4488_merged, file = "JK4488_merged_after_DoIntegrate.rds")
cat("Object saved after DoIntegrate: JK4488_merged_after_DoIntegrate.rds\n")


# ---- Post-processing: Combat (corrected counts) ----
DefaultAssay(JK4488_merged) <- "combat.reconstructed"
VariableFeatures(JK4488_merged) <- VariableFeatures(JK4488_merged[["SCT"]])
JK4488_merged <- ScaleData(JK4488_merged)
JK4488_merged <- RunPCA(JK4488_merged, npcs = 50L, reduction.name = "pca.combat")
JK4488_merged <- FindNeighbors(JK4488_merged, reduction = "pca.combat", dims = 1:30,
                               return.neighbor = TRUE, graph.name = "knn.combat")
JK4488_merged <- SymmetrizeKnn(JK4488_merged, graph.name = "knn.combat")
JK4488_merged <- FindOptimalClusters(JK4488_merged, graph.name = "knn.combat_symmetric",
                                     cluster.name = "combat_clusters_{metric}",
                                     cell.var = cell.var,
                                     optimisation.metric = c("nmi", "ari"))

# ---- Post-processing: Unintegrated (dimension reduction) ----
DefaultAssay(JK4488_merged) <- "SCT"
JK4488_merged <- FindNeighbors(JK4488_merged, reduction = "pca", dims = 1:30, k.param = 20L,
                               return.neighbor = TRUE, graph.name = "knn.unintegrated")
JK4488_merged <- SymmetrizeKnn(JK4488_merged, graph.name = "knn.unintegrated")
JK4488_merged <- FindOptimalClusters(JK4488_merged, graph.name = "knn.unintegrated_symmetric",
                                     cluster.name = "unintegrated_clusters_{metric}",
                                     cell.var = cell.var,
                                     optimisation.metric = c("nmi", "ari"))

# ---- Post-processing: Harmony (dimension reduction) ----
JK4488_merged <- FindNeighbors(JK4488_merged, reduction = "harmony", dims = 1:30,
                               return.neighbor = TRUE, graph.name = "knn.harmony")
JK4488_merged <- SymmetrizeKnn(JK4488_merged, graph.name = "knn.harmony")
JK4488_merged <- FindOptimalClusters(JK4488_merged, graph.name = "knn.harmony_symmetric",
                                     cluster.name = "harmony_clusters_{metric}",
                                     cell.var = cell.var,
                                     optimisation.metric = c("nmi", "ari"))

# ---- Post-processing: Scanorama (dimension reduction) ----
JK4488_merged <- FindNeighbors(JK4488_merged, reduction = "integrated.scanorama", dims = 1:30,
                               return.neighbor = TRUE, graph.name = "knn.scanorama")
JK4488_merged <- SymmetrizeKnn(JK4488_merged, graph.name = "knn.scanorama")
JK4488_merged <- FindOptimalClusters(JK4488_merged, graph.name = "knn.scanorama_symmetric",
                                     cluster.name = "scanorama_clusters_{metric}",
                                     cell.var = cell.var,
                                     optimisation.metric = c("nmi", "ari"))

# ---- Post-processing: scVI (dimension reduction) ----
JK4488_merged <- FindNeighbors(JK4488_merged, reduction = "integrated.scVI", dims = 1:30,
                               return.neighbor = TRUE, graph.name = "knn.scvi")
JK4488_merged <- SymmetrizeKnn(JK4488_merged, graph.name = "knn.scvi")
JK4488_merged <- FindOptimalClusters(JK4488_merged, graph.name = "knn.scvi_symmetric",
                                     cluster.name = "scvi_clusters_{metric}",
                                     cell.var = cell.var,
                                     optimisation.metric = c("nmi", "ari"))

# ---- Post-processing: trVAE (dimension reduction) ----
JK4488_merged <- FindNeighbors(JK4488_merged, reduction = "integrated.trVAE", dims = 1:30,
                               return.neighbor = TRUE, graph.name = "knn.trvae")
JK4488_merged <- SymmetrizeKnn(JK4488_merged, graph.name = "knn.trvae")
JK4488_merged <- FindOptimalClusters(JK4488_merged, graph.name = "knn.trvae_symmetric",
                                     cluster.name = "trvae_clusters_{metric}",
                                     cell.var = cell.var,
                                     optimisation.metric = c("nmi", "ari"))

# ---- Post-processing: BBKNN (graph) ----
JK4488_merged <- SymmetrizeKnn(JK4488_merged, graph.name = "bbknn_ridge.residuals_distances")
JK4488_merged <- FindOptimalClusters(JK4488_merged,
                                     graph.name = "bbknn_ridge.residuals_distances_symmetric",
                                     cluster.name = "bbknn_clusters_{metric}",
                                     cell.var = cell.var,
                                     optimisation.metric = c("nmi", "ari"))

# ---- Scoring: BATCH CORRECTION METRICS ONLY ----

# 1. PCA regression scores
JK4488_merged <- AddScoreRegressPC(JK4488_merged, integration = "unintegrated",
                                   batch.var = batch.var, reduction = "pca")
JK4488_merged <- AddScoreRegressPC(JK4488_merged, integration = "combat",
                                   batch.var = batch.var, reduction = "pca.combat")
JK4488_merged <- AddScoreRegressPC(JK4488_merged, integration = "harmony",
                                   batch.var = batch.var, reduction = "harmony")
JK4488_merged <- AddScoreRegressPC(JK4488_merged, integration = "scanorama",
                                   batch.var = batch.var, reduction = "integrated.scanorama")
JK4488_merged <- AddScoreRegressPC(JK4488_merged, integration = "scvi",
                                   batch.var = batch.var, reduction = "integrated.scVI")
JK4488_merged <- AddScoreRegressPC(JK4488_merged, integration = "trvae",
                                   batch.var = batch.var, reduction = "integrated.trVAE")

# 2. PCA density scores
JK4488_merged <- AddScoreDensityPC(JK4488_merged, integration = "unintegrated",
                                   batch.var = batch.var, reduction = "pca")
JK4488_merged <- AddScoreDensityPC(JK4488_merged, integration = "combat",
                                   batch.var = batch.var, reduction = "pca.combat")
JK4488_merged <- AddScoreDensityPC(JK4488_merged, integration = "harmony",
                                   batch.var = batch.var, reduction = "harmony")
JK4488_merged <- AddScoreDensityPC(JK4488_merged, integration = "scanorama",
                                   batch.var = batch.var, reduction = "integrated.scanorama")
JK4488_merged <- AddScoreDensityPC(JK4488_merged, integration = "scvi",
                                   batch.var = batch.var, reduction = "integrated.scVI")
JK4488_merged <- AddScoreDensityPC(JK4488_merged, integration = "trvae",
                                   batch.var = batch.var, reduction = "integrated.trVAE")

# 3. Cell cycle scores
JK4488_merged <- CellCycleScoringPerBatch(JK4488_merged, batch.var = batch.var,
                                          s.features = cc.genes$s.genes,
                                          g2m.features = cc.genes$g2m.genes)
JK4488_merged <- AddScoreRegressPC.CellCycle(JK4488_merged, integration = "unintegrated",
                                             batch.var = batch.var, what = "pca",
                                             compute.cc = FALSE, dims.use = 1:30)
JK4488_merged <- AddScoreRegressPC.CellCycle(JK4488_merged, integration = "combat",
                                             batch.var = batch.var, what = "pca.combat",
                                             compute.cc = FALSE, dims.use = 1:30)
JK4488_merged <- AddScoreRegressPC.CellCycle(JK4488_merged, integration = "harmony",
                                             batch.var = batch.var, what = "harmony",
                                             compute.cc = FALSE, dims.use = 1:30)
JK4488_merged <- AddScoreRegressPC.CellCycle(JK4488_merged, integration = "scanorama",
                                             batch.var = batch.var, what = "integrated.scanorama",
                                             compute.cc = FALSE, dims.use = 1:30)
JK4488_merged <- AddScoreRegressPC.CellCycle(JK4488_merged, integration = "scvi",
                                             batch.var = batch.var, what = "integrated.scVI",
                                             compute.cc = FALSE, dims.use = 1:30)
JK4488_merged <- AddScoreRegressPC.CellCycle(JK4488_merged, integration = "trvae",
                                             batch.var = batch.var, what = "integrated.trVAE",
                                             compute.cc = FALSE, dims.use = 1:30)

# 4. ASW batch scores (without cell.var)
JK4488_merged <- AddScoreASWBatch(JK4488_merged, integration = "unintegrated",
                                  batch.var = batch.var, 
                                  per.cell.var = FALSE,
                                  what = "pca")
JK4488_merged <- AddScoreASWBatch(JK4488_merged, integration = "combat",
                                  batch.var = batch.var,
                                  per.cell.var = FALSE,
                                  what = "pca.combat")
JK4488_merged <- AddScoreASWBatch(JK4488_merged, integration = "harmony",
                                  batch.var = batch.var,
                                  per.cell.var = FALSE,
                                  what = "harmony")
JK4488_merged <- AddScoreASWBatch(JK4488_merged, integration = "scanorama",
                                  batch.var = batch.var,
                                  per.cell.var = FALSE,
                                  what = "integrated.scanorama")
JK4488_merged <- AddScoreASWBatch(JK4488_merged, integration = "scvi",
                                  batch.var = batch.var,
                                  per.cell.var = FALSE,
                                  what = "integrated.scVI")
JK4488_merged <- AddScoreASWBatch(JK4488_merged, integration = "trvae",
                                  batch.var = batch.var,
                                  per.cell.var = FALSE,
                                  what = "integrated.trVAE")

# ---- Scale and visualize scores ----
JK4488_merged <- ScaleScores(JK4488_merged)

# Plot and save score comparison
score_plot <- PlotScores(JK4488_merged)
ggsave("integration_scores_comparison.pdf", plot = score_plot, width = 10, height = 8)
ggsave("integration_scores_comparison.png", plot = score_plot, width = 10, height = 8, dpi = 300)

# ---- Create UMAPs for visual inspection ----
JK4488_merged <- RunUMAP(JK4488_merged, dims = 1:30, reduction = "pca.combat",
                         reduction.name = "umap.combat")
JK4488_merged <- RunUMAP(JK4488_merged, dims = 1:30, reduction = "harmony",
                         reduction.name = "umap.harmony")
JK4488_merged <- RunUMAP(JK4488_merged, dims = 1:30, reduction = "integrated.scanorama",
                         reduction.name = "umap.scanorama")
JK4488_merged <- RunUMAP(JK4488_merged, dims = 1:30, reduction = "integrated.scVI",
                         reduction.name = "umap.scvi")
JK4488_merged <- RunUMAP(JK4488_merged, dims = 1:30, reduction = "integrated.trVAE",
                         reduction.name = "umap.trvae")

# For BBKNN (graph-based)
library(future)
plan(multisession)
JK4488_merged <- future({
  reticulate::use_condaenv("/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/SeuratIntegrate", required = TRUE)
  RunUMAP(JK4488_merged, graph = "bbknn_ridge.residuals_connectivities",
          umap.method = "umap-learn", n.epochs = 200L,
          reduction.name = "umap.bbknn")
})
plan(sequential)

# ---- Save individual UMAP plots ----
# Create output directory for plots
dir.create("integration_plots", showWarnings = FALSE)

# Unintegrated
p1 <- DimPlot(JK4488_merged, group.by = batch.var, reduction = "umap")
ggsave("integration_plots/umap_unintegrated.pdf", plot = p1, width = 8, height = 6)
ggsave("integration_plots/umap_unintegrated.png", plot = p1, width = 8, height = 6, dpi = 300)

# Combat
p2 <- DimPlot(JK4488_merged, group.by = batch.var, reduction = "umap.combat")
ggsave("integration_plots/umap_combat.pdf", plot = p2, width = 8, height = 6)
ggsave("integration_plots/umap_combat.png", plot = p2, width = 8, height = 6, dpi = 300)

# Harmony
p3 <- DimPlot(JK4488_merged, group.by = batch.var, reduction = "umap.harmony")
ggsave("integration_plots/umap_harmony.pdf", plot = p3, width = 8, height = 6)
ggsave("integration_plots/umap_harmony.png", plot = p3, width = 8, height = 6, dpi = 300)

# Scanorama
p4 <- DimPlot(JK4488_merged, group.by = batch.var, reduction = "umap.scanorama")
ggsave("integration_plots/umap_scanorama.pdf", plot = p4, width = 8, height = 6)
ggsave("integration_plots/umap_scanorama.png", plot = p4, width = 8, height = 6, dpi = 300)

# scVI
p5 <- DimPlot(JK4488_merged, group.by = batch.var, reduction = "umap.scvi")
ggsave("integration_plots/umap_scvi.pdf", plot = p5, width = 8, height = 6)
ggsave("integration_plots/umap_scvi.png", plot = p5, width = 8, height = 6, dpi = 300)

# trVAE
p6 <- DimPlot(JK4488_merged, group.by = batch.var, reduction = "umap.trvae")
ggsave("integration_plots/umap_trvae.pdf", plot = p6, width = 8, height = 6)
ggsave("integration_plots/umap_trvae.png", plot = p6, width = 8, height = 6, dpi = 300)

# BBKNN
p7 <- DimPlot(JK4488_merged, group.by = batch.var, reduction = "umap.bbknn")
ggsave("integration_plots/umap_bbknn.pdf", plot = p7, width = 8, height = 6)
ggsave("integration_plots/umap_bbknn.png", plot = p7, width = 8, height = 6, dpi = 300)

# Create a combined plot with all methods
library(patchwork)
combined_plot <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + 
  plot_layout(ncol = 3, guides = 'collect')

ggsave("integration_plots/umap_all_methods_combined.pdf", 
       plot = combined_plot, width = 18, height = 12)
ggsave("integration_plots/umap_all_methods_combined.png", 
       plot = combined_plot, width = 18, height = 12, dpi = 300)

cat("All plots saved to 'integration_plots/' directory\n")

# ---- Save the object ----
saveRDS(JK4488_merged, file = "JK4488_merged_SeuratIntegrate_final.rds")

cat("\n=== Analysis Complete ===\n")
cat("Seurat object saved: JK4488_merged_SeuratIntegrate_final.rds\n")
cat("Score comparison plot: integration_scores_comparison.pdf/png\n")
cat("Individual UMAP plots: integration_plots/ directory\n")
