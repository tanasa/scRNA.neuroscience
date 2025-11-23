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

setwd("/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_seurat.BT.analysis.sep30.nov22.16samples")

list.files()

base_dir <- "." 

# sample_ids <- c("JK4488_01", "JK4488_02", "JK4488_03", 
#                "JK4488_04", "JK4488_05", "JK4488_06")

# ---- Read Seurat objects ----
# seurat_list <- list()
# for (sample in sample_ids) {
#  file_path <- file.path(base_dir, paste0(sample, "_seurat_obj.filtered.soupX.rds"))
#  if (file.exists(file_path)) {
#    cat("Reading Seurat object:", sample, "\n")
#    seurat_list[[sample]] <- readRDS(file_path)
#  } else {
#    warning("File not found: ", file_path)
#  }
# }

# ---- Check what was loaded ----
# cat("Loaded Seurat objects:", length(seurat_list), "\n")
# cat("Loaded SCE objects:", length(sce_list), "\n")

# Summary of each object
# for (sample in names(seurat_list)) {
#  cat("\n", sample, "- Seurat:\n")
#  print(dim(seurat_list[[sample]]))
# }

# for (sample in names(sce_list)) {
#  cat("\n", sample, "- SCE:\n")
#  print(dim(sce_list[[sample]]))
# }

# print(seurat_list)

# ---- Merge samples ----
# cat("Merging samples...\n")
# merged <- merge(seurat_list[[1]], y = seurat_list[-1], 
#                add.cell.ids = names(seurat_list))

# seurat_list[["JK4488_01"]]
# seurat_list[["JK4488_02"]]
# seurat_list[["JK4488_03"]]
# seurat_list[["JK4488_04"]]
# seurat_list[["JK4488_05"]]
# seurat_list[["JK4488_06"]]

# ---- Or read individually with descriptive names ----
# AD_AG_2658_seurat <- readRDS("AD_AG_2658_seurat_obj.filtered.soupX.rds")  # EXCLUDED
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
# Dy_AG_3090_seurat <- readRDS("Dy_AG_3090_seurat_obj.filtered.soupX.rds")  # EXCLUDED
Dy_AG_3296_seurat <- readRDS("Dy_AG_3296_seurat_obj.filtered.soupX.rds")
Healthy_AG_2394_seurat <- readRDS("Healthy_AG_2394_seurat_obj.filtered.soupX.rds")
Healthy_AG_2594_seurat <- readRDS("Healthy_AG_2594_seurat_obj.filtered.soupX.rds")
Healthy_AG_3041_seurat <- readRDS("Healthy_AG_3041_seurat_obj.filtered.soupX.rds")
Healthy_AG_3267_seurat <- readRDS("Healthy_AG_3267_seurat_obj.filtered.soupX.rds")

# Summary function
show_summary <- function(obj, name) {
  cat(sprintf("\n%s:\n", name))
  cat(sprintf("  Genes: %d\n", nrow(obj)))
  cat(sprintf("  Cells: %d\n", ncol(obj)))
  cat(sprintf("  UMAP: %s\n", ifelse("umap" %in% names(obj@reductions), "Yes", "No")))
}

# Show summaries
# show_summary(AD_AG_2658_seurat, "AD_AG_2658")  # EXCLUDED
show_summary(AD_AG_2822_seurat, "AD_AG_2822")
show_summary(AD_AG_3184_seurat, "AD_AG_3184")
show_summary(AD_AG_3188_seurat, "AD_AG_3188")
show_summary(AD_CTL_2273_S2_seurat, "AD_CTL_2273_S2")
show_summary(AD_CTL_2354_S3_seurat, "AD_CTL_2354_S3")
show_summary(AD_CTL_2812_S1_seurat, "AD_CTL_2812_S1")
show_summary(Dy_AD_2328_S5_seurat, "Dy_AD_2328_S5")
show_summary(Dy_AD_2368_S6_seurat, "Dy_AD_2368_S6")
show_summary(Dy_AD_2660_S4_seurat, "Dy_AD_2660_S4")
show_summary(Dy_AG_2659_seurat, "Dy_AG_2659")
show_summary(Dy_AG_2922_seurat, "Dy_AG_2922")
# show_summary(Dy_AG_3090_seurat, "Dy_AG_3090")  # EXCLUDED
show_summary(Dy_AG_3296_seurat, "Dy_AG_3296")
show_summary(Healthy_AG_2394_seurat, "Healthy_AG_2394")
show_summary(Healthy_AG_2594_seurat, "Healthy_AG_2594")
show_summary(Healthy_AG_3041_seurat, "Healthy_AG_3041")
show_summary(Healthy_AG_3267_seurat, "Healthy_AG_3267")

# Create a summary table
summary_df <- data.frame(
  Sample = c("AD_AG_2822", "AD_AG_3184", "AD_AG_3188", 
             "AD_CTL_2273_S2", "AD_CTL_2354_S3", "AD_CTL_2812_S1",
             "Dy_AD_2328_S5", "Dy_AD_2368_S6", "Dy_AD_2660_S4",
             "Dy_AG_2659", "Dy_AG_2922", "Dy_AG_3296",
             "Healthy_AG_2394", "Healthy_AG_2594", "Healthy_AG_3041", "Healthy_AG_3267"),
  Genes = c(nrow(AD_AG_2822_seurat), nrow(AD_AG_3184_seurat), nrow(AD_AG_3188_seurat),
            nrow(AD_CTL_2273_S2_seurat), nrow(AD_CTL_2354_S3_seurat), nrow(AD_CTL_2812_S1_seurat),
            nrow(Dy_AD_2328_S5_seurat), nrow(Dy_AD_2368_S6_seurat), nrow(Dy_AD_2660_S4_seurat),
            nrow(Dy_AG_2659_seurat), nrow(Dy_AG_2922_seurat), nrow(Dy_AG_3296_seurat),
            nrow(Healthy_AG_2394_seurat), nrow(Healthy_AG_2594_seurat), nrow(Healthy_AG_3041_seurat), nrow(Healthy_AG_3267_seurat)),
  Cells = c(ncol(AD_AG_2822_seurat), ncol(AD_AG_3184_seurat), ncol(AD_AG_3188_seurat),
            ncol(AD_CTL_2273_S2_seurat), ncol(AD_CTL_2354_S3_seurat), ncol(AD_CTL_2812_S1_seurat),
            ncol(Dy_AD_2328_S5_seurat), ncol(Dy_AD_2368_S6_seurat), ncol(Dy_AD_2660_S4_seurat),
            ncol(Dy_AG_2659_seurat), ncol(Dy_AG_2922_seurat), ncol(Dy_AG_3296_seurat),
            ncol(Healthy_AG_2394_seurat), ncol(Healthy_AG_2594_seurat), ncol(Healthy_AG_3041_seurat), ncol(Healthy_AG_3267_seurat))
)

print(summary_df)

options(repr.plot.width = 10, repr.plot.height = 8)

# Set orig.ident for each sample
# AD_AG_2658_seurat$orig.ident = "AD_AG_2658"  # EXCLUDED
AD_AG_2822_seurat$orig.ident = "AD_AG_2822"
AD_AG_3184_seurat$orig.ident = "AD_AG_3184"
AD_AG_3188_seurat$orig.ident = "AD_AG_3188"
AD_CTL_2273_S2_seurat$orig.ident = "AD_CTL_2273_S2"
AD_CTL_2354_S3_seurat$orig.ident = "AD_CTL_2354_S3"
AD_CTL_2812_S1_seurat$orig.ident = "AD_CTL_2812_S1"
Dy_AD_2328_S5_seurat$orig.ident = "Dy_AD_2328_S5"
Dy_AD_2368_S6_seurat$orig.ident = "Dy_AD_2368_S6"
Dy_AD_2660_S4_seurat$orig.ident = "Dy_AD_2660_S4"
Dy_AG_2659_seurat$orig.ident = "Dy_AG_2659"
Dy_AG_2922_seurat$orig.ident = "Dy_AG_2922"
# Dy_AG_3090_seurat$orig.ident = "Dy_AG_3090"  # EXCLUDED
Dy_AG_3296_seurat$orig.ident = "Dy_AG_3296"
Healthy_AG_2394_seurat$orig.ident = "Healthy_AG_2394"
Healthy_AG_2594_seurat$orig.ident = "Healthy_AG_2594"
Healthy_AG_3041_seurat$orig.ident = "Healthy_AG_3041"
Healthy_AG_3267_seurat$orig.ident = "Healthy_AG_3267"

# Create list of samples
seurat_list <- list(
  # AD_AG_2658 = AD_AG_2658_seurat,  # EXCLUDED
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
  # Dy_AG_3090 = Dy_AG_3090_seurat,  # EXCLUDED
  Dy_AG_3296 = Dy_AG_3296_seurat,
  Healthy_AG_2394 = Healthy_AG_2394_seurat,
  Healthy_AG_2594 = Healthy_AG_2594_seurat,
  Healthy_AG_3041 = Healthy_AG_3041_seurat,
  Healthy_AG_3267 = Healthy_AG_3267_seurat
)

# Add sample ID to metadata
for (i in seq_along(seurat_list)) {
  seurat_list[[i]]$sample <- names(seurat_list)[i]
}

print(seurat_list)

# ---- Merge samples ----
cat("Merging samples...\n")
rm(merged)
merged <- merge(seurat_list[[1]], y = seurat_list[-1], 
                add.cell.ids = names(seurat_list))

str(merged)

head(merged@meta.data, 2)
tail(merged@meta.data, 2)

# merged[["RNA"]] <- split(merged[["RNA"]], f = merged$sample)

# https://satijalab.org/seurat/articles/seurat5_integration
# merged[["RNA"]] <- split(merged[["RNA"]], f = merged$sample)

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)

merged <- FindNeighbors(merged, 
                        dims = 1:30, 
                        reduction = "pca", 
                        graph.name = "graph_unintegrated")
merged <- FindClusters(merged, 
                       resolution = 0.5, 
                       graph.name   = "graph_unintegrated",
                       cluster.name = "clusters_unintegrated_res0.5") # RESOLUTION 0.5

merged <- FindClusters(merged, 
                       resolution = 0.1, 
                       graph.name   = "graph_unintegrated",
                       cluster.name = "clusters_unintegrated_res0.1") # RESOLUTION 0.1

merged <- FindClusters(merged, 
                       resolution = 1, 
                       graph.name   = "graph_unintegrated",
                       cluster.name = "clusters_unintegrated_res1.0") # RESOLUTION 1

# since the data is split into layers, normalization and variable feature identification is performed for each batch independently 
# (a consensus set of variable features is automatically identified

head(merged@meta.data, 2)
tail(merged@meta.data, 2)

merged <- RunUMAP(merged, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth

# Create output directory for plots
dir.create("plots", showWarnings = FALSE)

# Unintegrated plots
p1 <- DimPlot(merged, reduction = "umap.unintegrated", group.by = "orig.ident")
ggsave("plots/umap_unintegrated_by_sample.pdf", plot = p1, width = 10, height = 8)
ggsave("plots/umap_unintegrated_by_sample.png", plot = p1, width = 10, height = 8, dpi = 300)

p2 <- DimPlot(
  merged,
  reduction = "umap.unintegrated",
  group.by  = "clusters_unintegrated_res0.1",
  label = TRUE, repel = TRUE, raster = TRUE, pt.size = 0.2
)
ggsave("plots/umap_unintegrated_clusters_res0.1.pdf", plot = p2, width = 10, height = 8)
ggsave("plots/umap_unintegrated_clusters_res0.1.png", plot = p2, width = 10, height = 8, dpi = 300)

p3 <- DimPlot(
  merged,
  reduction = "umap.unintegrated",
  group.by  = "clusters_unintegrated_res0.5",
  label = TRUE, repel = TRUE, raster = TRUE, pt.size = 0.2
)
ggsave("plots/umap_unintegrated_clusters_res0.5.pdf", plot = p3, width = 10, height = 8)
ggsave("plots/umap_unintegrated_clusters_res0.5.png", plot = p3, width = 10, height = 8, dpi = 300)

p4 <- DimPlot(
  merged,
  reduction = "umap.unintegrated",
  group.by  = "clusters_unintegrated_res1.0",
  label = TRUE, repel = TRUE, raster = TRUE, pt.size = 0.2
)
ggsave("plots/umap_unintegrated_clusters_res1.0.pdf", plot = p4, width = 10, height = 8)
ggsave("plots/umap_unintegrated_clusters_res1.0.png", plot = p4, width = 10, height = 8, dpi = 300)




# Perform streamlined (one-line) integrative analysis
# https://satijalab.org/seurat/articles/seurat5_integration

# Anchor-based CCA integration (method=CCAIntegration)
# Anchor-based RPCA integration (method=RPCAIntegration)
# Harmony (method=HarmonyIntegration)
# FastMNN (method= FastMNNIntegration)
# scVI (method=scVIIntegration)

obj = merged

# CCA Integration - COMMENTED OUT
# obj <- IntegrateLayers(
#   object = obj, 
#   method = CCAIntegration,
#   orig.reduction = "pca", 
#   new.reduction = "integrated.cca",
#   verbose = FALSE
# )

# slotNames(obj)
(obj@reductions)

library(future)
options(future.globals.maxSize = 256 * 1024^3)  # 128 GB

# plan(sequential)
obj <- IntegrateLayers(
  object = obj, 
  method = RPCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.rpca",
  verbose = FALSE
)
# plan(multicore)  

# slotNames(obj)
(obj@reductions)

library(future)
options(future.globals.maxSize = 256 * 1024^3)  # 128 GB

# plan(sequential)
obj <- IntegrateLayers(
  object = obj, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony",
  verbose = FALSE
)
# plan(multicore)  

# slotNames(obj)
(obj@reductions)

# FastMNN Integration - COMMENTED OUT
# library(future)
# library("batchelor")
# library("scuttle")
# library(SeuratWrappers)

# options(future.globals.maxSize = 256 * 1024^3)  # 128 GB

# plan(sequential)
# obj <- IntegrateLayers(
#   object = obj, 
#   method = FastMNNIntegration,
#   new.reduction = "integrated.mnn",
#   verbose = FALSE
# )
# plan(multicore)  

# slotNames(obj)
(obj@reductions)

getwd()



# CCA downstream analysis - COMMENTED OUT
# obj <- FindNeighbors(obj, 
#                      reduction = "integrated.cca", 
#                      graph.name = "graph_cca",
#                      dims = 1:30)

# obj <- FindClusters(obj, 
#                     resolution = 0.1, 
#                     graph.name   = "graph_cca",
#                     cluster.name = "cca_clusters_res0.1")

# obj <- FindClusters(obj, 
#                     resolution = 0.5, 
#                     graph.name   = "graph_cca",
#                     cluster.name = "cca_clusters_res0.5")

# obj <- FindClusters(obj, 
#                     resolution = 1, 
#                     graph.name   = "graph_cca",
#                     cluster.name = "cca_clusters_res1.0")


# table(obj$cca_clusters_res0.1)
# table(obj$cca_clusters_res0.5)
# table(obj$cca_clusters_res1.0)

# head(obj@meta.data,3)
# Reductions(obj)

# make a UMAP from the integrated CCA space and name it explicitly
# obj <- RunUMAP(
#   obj,
#   reduction      = "integrated.cca",  # your integrated space
#   dims           = 1:30,
#   reduction.name = "umap.cca",        # choose a name
#   verbose        = FALSE
# )

# now plot using that reduction and your cluster columns
# DimPlot(obj, reduction = "umap.cca", group.by = "cca_clusters_res0.1", label = TRUE)
# DimPlot(obj, reduction = "umap.cca", group.by = "cca_clusters_res0.5", label = TRUE)
# DimPlot(obj, reduction = "umap.cca", group.by = "cca_clusters_res1.0",   label = TRUE)



obj <- FindNeighbors(obj, 
                     reduction = "integrated.rpca", 
                     graph.name = "graph_rpca",
                     dims = 1:30)

obj <- FindClusters(obj, 
                    resolution = 0.1, 
                    graph.name   = "graph_rpca",
                    cluster.name = "rpca_clusters_res0.1")

obj <- FindClusters(obj, 
                    resolution = 0.5, 
                    graph.name   = "graph_rpca",
                    cluster.name = "rpca_clusters_res0.5")

obj <- FindClusters(obj, 
                    resolution = 1, 
                    graph.name   = "graph_rpca",
                    cluster.name = "rpca_clusters_res1.0")


obj <- FindClusters(obj, 
                    resolution = 0.2, 
                    graph.name   = "graph_rpca",
                    cluster.name = "rpca_clusters_res0.2")

obj <- FindClusters(obj, 
                    resolution = 0.3, 
                    graph.name   = "graph_rpca",
                    cluster.name = "rpca_clusters_res0.3")

obj <- FindClusters(obj, 
                    resolution = 0.4, 
                    graph.name   = "graph_rpca",
                    cluster.name = "rpca_clusters_res0.4")


table(obj$rpca_clusters_res0.1)
table(obj$rpca_clusters_res0.5)
table(obj$rpca_clusters_res1.0)

head(obj@meta.data,3)
Reductions(obj)

# make a UMAP from the integrated CCA space and name it explicitly
obj <- RunUMAP(
  obj,
  reduction      = "integrated.rpca",  # your integrated space
  dims           = 1:30,
  reduction.name = "umap.rpca",        # choose a name
  verbose        = FALSE
)

# now plot using that reduction and your cluster columns
p5 <- DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters_res0.1", label = TRUE)
ggsave("plots/umap_rpca_clusters_res0.1.pdf", plot = p5, width = 10, height = 8)
ggsave("plots/umap_rpca_clusters_res0.1.png", plot = p5, width = 10, height = 8, dpi = 300)

p6 <- DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters_res0.5", label = TRUE)
ggsave("plots/umap_rpca_clusters_res0.5.pdf", plot = p6, width = 10, height = 8)
ggsave("plots/umap_rpca_clusters_res0.5.png", plot = p6, width = 10, height = 8, dpi = 300)

p7 <- DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters_res1.0",   label = TRUE)
ggsave("plots/umap_rpca_clusters_res1.0.pdf", plot = p7, width = 10, height = 8)
ggsave("plots/umap_rpca_clusters_res1.0.png", plot = p7, width = 10, height = 8, dpi = 300)


# now plot using that reduction and your cluster columns
p8 <- DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters_res0.2", label = TRUE)
ggsave("plots/umap_rpca_clusters_res0.2.pdf", plot = p8, width = 10, height = 8)
ggsave("plots/umap_rpca_clusters_res0.2.png", plot = p8, width = 10, height = 8, dpi = 300)

p9 <- DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters_res0.3", label = TRUE)
ggsave("plots/umap_rpca_clusters_res0.3.pdf", plot = p9, width = 10, height = 8)
ggsave("plots/umap_rpca_clusters_res0.3.png", plot = p9, width = 10, height = 8, dpi = 300)

p10 <- DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters_res0.4",   label = TRUE)
ggsave("plots/umap_rpca_clusters_res0.4.pdf", plot = p10, width = 10, height = 8)
ggsave("plots/umap_rpca_clusters_res0.4.png", plot = p10, width = 10, height = 8, dpi = 300)



obj <- FindNeighbors(obj, 
                     reduction = "harmony", 
                     graph.name = "graph_harmony",
                     dims = 1:30)

obj <- FindClusters(obj, 
                    resolution = 0.1, 
                    graph.name = "graph_harmony",
                    cluster.name = "harmony_clusters_res0.1")
obj <- FindClusters(obj, 
                    resolution = 0.5, 
                    graph.name = "graph_harmony",
                    cluster.name = "harmony_clusters_res0.5")
obj <- FindClusters(obj, 
                    resolution = 1, 
                    graph.name = "graph_harmony",
                    cluster.name = "harmony_clusters_res1.0")

obj <- FindClusters(obj, 
                    resolution = 0.2, 
                    graph.name = "graph_harmony",
                    cluster.name = "harmony_clusters_res0.2")
obj <- FindClusters(obj, 
                    resolution = 0.3, 
                    graph.name = "graph_harmony",
                    cluster.name = "harmony_clusters_res0.3")
obj <- FindClusters(obj, 
                    resolution = 0.4, 
                    graph.name = "graph_harmony",
                    cluster.name = "harmony_clusters_res0.4")

table(obj$harmony_clusters_res0.1)
table(obj$harmony_clusters_res0.5)
table(obj$harmony_clusters_res1.0)

table(obj$harmony_clusters_res0.2)
table(obj$harmony_clusters_res0.3)
table(obj$harmony_clusters_res0.4)

head(obj@meta.data,3)
Reductions(obj)

# make a UMAP from the integrated CCA space and name it explicitly
obj <- RunUMAP(
  obj,
  reduction      = "harmony",  # your integrated space
  dims           = 1:30,
  reduction.name = "umap.harmony",        # choose a name
  verbose        = FALSE
)

# now plot using that reduction and your cluster columns
p11 <- DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters_res0.1", label = TRUE)
ggsave("plots/umap_harmony_clusters_res0.1.pdf", plot = p11, width = 10, height = 8)
ggsave("plots/umap_harmony_clusters_res0.1.png", plot = p11, width = 10, height = 8, dpi = 300)

p12 <- DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters_res0.5", label = TRUE)
ggsave("plots/umap_harmony_clusters_res0.5.pdf", plot = p12, width = 10, height = 8)
ggsave("plots/umap_harmony_clusters_res0.5.png", plot = p12, width = 10, height = 8, dpi = 300)

p13 <- DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters_res1.0",   label = TRUE)
ggsave("plots/umap_harmony_clusters_res1.0.pdf", plot = p13, width = 10, height = 8)
ggsave("plots/umap_harmony_clusters_res1.0.png", plot = p13, width = 10, height = 8, dpi = 300)


# now plot using that reduction and your cluster columns
p14 <- DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters_res0.2", label = TRUE)
ggsave("plots/umap_harmony_clusters_res0.2.pdf", plot = p14, width = 10, height = 8)
ggsave("plots/umap_harmony_clusters_res0.2.png", plot = p14, width = 10, height = 8, dpi = 300)

p15 <- DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters_res0.3", label = TRUE)
ggsave("plots/umap_harmony_clusters_res0.3.pdf", plot = p15, width = 10, height = 8)
ggsave("plots/umap_harmony_clusters_res0.3.png", plot = p15, width = 10, height = 8, dpi = 300)

p16 <- DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters_res0.4",   label = TRUE)
ggsave("plots/umap_harmony_clusters_res0.4.pdf", plot = p16, width = 10, height = 8)
ggsave("plots/umap_harmony_clusters_res0.4.png", plot = p16, width = 10, height = 8, dpi = 300)





# FastMNN downstream analysis - COMMENTED OUT
# obj <- FindNeighbors(obj,
#                      reduction = "integrated.mnn",
#                      graph.name = "graph_mnn",
#                      dims = 1:30)

# obj <- FindClusters(obj,
#                     resolution = 0.1,
#                     graph.name   = "graph_mnn",
#                     cluster.name = "mnn_clusters_res0.1")

# obj <- FindClusters(obj,
#                     resolution = 0.5,
#                     graph.name   = "graph_mnn",
#                     cluster.name = "mnn_clusters_res0.5")

# obj <- FindClusters(obj,
#                     resolution = 1,
#                     graph.name   = "graph_mnn",
#                     cluster.name = "mnn_clusters_res1.0")

# table(obj$mnn_clusters_res0.1)
# table(obj$mnn_clusters_res0.5)
# table(obj$mnn_clusters_res1.0)

# head(obj@meta.data,3)
# Reductions(obj)

# make a UMAP from the integrated CCA space and name it explicitly
# obj <- RunUMAP(
#   obj,
#   reduction      = "integrated.mnn",  # your integrated space
#   dims           = 1:30,
#   reduction.name = "umap.mnn",        # choose a name
#   verbose        = FALSE
# )

# now plot using that reduction and your cluster columns
# DimPlot(obj, reduction = "umap.mnn", group.by = "mnn_clusters_res0.1", label = TRUE)
# DimPlot(obj, reduction = "umap.mnn", group.by = "mnn_clusters_res0.5", label = TRUE)
# DimPlot(obj, reduction = "umap.mnn", group.by = "mnn_clusters_res1.0",   label = TRUE)

getwd()

cat("using scVI")

obj <- IntegrateLayers(
  object = obj, 
  method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/mnt/nfs/CX000008_DS1/projects/btanasa/virtual_env/SeuratIntegrate", 
 verbose = TRUE
)

(obj@reductions)

obj <- FindNeighbors(obj,
                     reduction = "integrated.scvi",
                     graph.name = "graph_scvi",
                     dims = 1:10)

obj <- FindClusters(obj,
                    resolution = 0.1,
                    graph.name   = "graph_scvi",
                    cluster.name = "scvi_clusters_res0.1")

obj <- FindClusters(obj,
                    resolution = 0.5,
                    graph.name   = "graph_scvi",
                    cluster.name = "scvi_clusters_res0.5")

obj <- FindClusters(obj,
                    resolution = 1,
                    graph.name   = "graph_scvi",
                    cluster.name = "scvi_clusters_res1.0")



obj <- FindClusters(obj,
                    resolution = 0.2,
                    graph.name   = "graph_scvi",
                    cluster.name = "scvi_clusters_res0.2")

obj <- FindClusters(obj,
                    resolution = 0.3,
                    graph.name   = "graph_scvi",
                    cluster.name = "scvi_clusters_res0.3")

obj <- FindClusters(obj,
                    resolution = 0.4,
                    graph.name   = "graph_scvi",
                    cluster.name = "scvi_clusters_res0.4")

table(obj$scvi_clusters_res0.1)
table(obj$scvi_clusters_res0.5)
table(obj$scvi_clusters_res1.0)

table(obj$scvi_clusters_res0.2)
table(obj$scvi_clusters_res0.3)
table(obj$scvi_clusters_res0.4)

head(obj@meta.data,3)
Reductions(obj)

# make a UMAP from the integrated CCA space and name it explicitly
obj <- RunUMAP(
  obj,
  reduction      = "integrated.scvi",  # your integrated space
  dims           = 1:10,
  reduction.name = "umap.scvi",        # choose a name
  verbose        = FALSE
)

# now plot using that reduction and your cluster columns
p17 <- DimPlot(obj, reduction = "umap.scvi", group.by = "scvi_clusters_res0.1", label = TRUE)
ggsave("plots/umap_scvi_clusters_res0.1.pdf", plot = p17, width = 10, height = 8)
ggsave("plots/umap_scvi_clusters_res0.1.png", plot = p17, width = 10, height = 8, dpi = 300)

p18 <- DimPlot(obj, reduction = "umap.scvi", group.by = "scvi_clusters_res0.5", label = TRUE)
ggsave("plots/umap_scvi_clusters_res0.5.pdf", plot = p18, width = 10, height = 8)
ggsave("plots/umap_scvi_clusters_res0.5.png", plot = p18, width = 10, height = 8, dpi = 300)

p19 <- DimPlot(obj, reduction = "umap.scvi", group.by = "scvi_clusters_res1.0",   label = TRUE)
ggsave("plots/umap_scvi_clusters_res1.0.pdf", plot = p19, width = 10, height = 8)
ggsave("plots/umap_scvi_clusters_res1.0.png", plot = p19, width = 10, height = 8, dpi = 300)

# now plot using that reduction and your cluster columns
p20 <- DimPlot(obj, reduction = "umap.scvi", group.by = "scvi_clusters_res0.2", label = TRUE)
ggsave("plots/umap_scvi_clusters_res0.2.pdf", plot = p20, width = 10, height = 8)
ggsave("plots/umap_scvi_clusters_res0.2.png", plot = p20, width = 10, height = 8, dpi = 300)

p21 <- DimPlot(obj, reduction = "umap.scvi", group.by = "scvi_clusters_res0.3", label = TRUE)
ggsave("plots/umap_scvi_clusters_res0.3.pdf", plot = p21, width = 10, height = 8)
ggsave("plots/umap_scvi_clusters_res0.3.png", plot = p21, width = 10, height = 8, dpi = 300)

p22 <- DimPlot(obj, reduction = "umap.scvi", group.by = "scvi_clusters_res0.4",   label = TRUE)
ggsave("plots/umap_scvi_clusters_res0.4.pdf", plot = p22, width = 10, height = 8)
ggsave("plots/umap_scvi_clusters_res0.4.png", plot = p22, width = 10, height = 8, dpi = 300)

head(obj@meta.data)

# Save any R object (e.g., a Seurat or SCE object)
saveRDS(obj, file = "zfiles.integrate.seurat.rds")

cat("\n=== Analysis Complete ===\n")
cat("Seurat object saved: zfiles.integrate.seurat.rds\n")
cat("All plots saved to 'plots/' directory\n")

# to continue :
# https://github.com/cbib/Seurat-Integrate/blob/main/vignettes/SeuratIntegrate.Rmd


