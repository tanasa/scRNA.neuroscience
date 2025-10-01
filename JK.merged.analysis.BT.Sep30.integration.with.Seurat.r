setwd("/mnt/nfs/CX000008_DS1/projects/jaeyeon/Dyslexia_Analysis/RDS")

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

setwd("/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/samples.merged_seurat.BT.analysis.sep29")

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
JK4488_01_seurat <- readRDS("JK4488_01_seurat_obj.filtered.soupX.rds")
# JK4488_01_sce <- readRDS("JK4488_01_sce_obj.filtered.soupX.rds")

JK4488_02_seurat <- readRDS("JK4488_02_seurat_obj.filtered.soupX.rds")
# JK4488_02_sce <- readRDS("JK4488_02_sce_obj.filtered.soupX.rds")

JK4488_03_seurat <- readRDS("JK4488_03_seurat_obj.filtered.soupX.rds")
# JK4488_03_sce <- readRDS("JK4488_03_sce_obj.filtered.soupX.rds")

JK4488_04_seurat <- readRDS("JK4488_04_seurat_obj.filtered.soupX.rds")
# JK4488_04_sce <- readRDS("JK4488_04_sce_obj.filtered.soupX.rds")

JK4488_05_seurat <- readRDS("JK4488_05_seurat_obj.filtered.soupX.rds")
# JK4488_05_sce <- readRDS("JK4488_05_sce_obj.filtered.soupX.rds")

JK4488_06_seurat <- readRDS("JK4488_06_seurat_obj.filtered.soupX.rds")
# JK4488_06_sce <- readRDS("JK4488_06_sce_obj.filtered.soupX.rds")

# Summary function
show_summary <- function(obj, name) {
  cat(sprintf("\n%s:\n", name))
  cat(sprintf("  Genes: %d\n", nrow(obj)))
  cat(sprintf("  Cells: %d\n", ncol(obj)))
  cat(sprintf("  UMAP: %s\n", ifelse("umap" %in% names(obj@reductions), "Yes", "No")))
}

# Show summaries
show_summary(JK4488_01_seurat, "JK4488_01")
show_summary(JK4488_02_seurat, "JK4488_02")
show_summary(JK4488_03_seurat, "JK4488_03")
show_summary(JK4488_04_seurat, "JK4488_04")
show_summary(JK4488_05_seurat, "JK4488_05")
show_summary(JK4488_06_seurat, "JK4488_06")

# Create a summary table
summary_df <- data.frame(
  Sample = c("JK4488_01", "JK4488_02", "JK4488_03", "JK4488_04", "JK4488_05", "JK4488_06"),
  Genes = c(nrow(JK4488_01_seurat), nrow(JK4488_02_seurat), nrow(JK4488_03_seurat),
            nrow(JK4488_04_seurat), nrow(JK4488_05_seurat), nrow(JK4488_06_seurat)),
  Cells = c(ncol(JK4488_01_seurat), ncol(JK4488_02_seurat), ncol(JK4488_03_seurat),
            ncol(JK4488_04_seurat), ncol(JK4488_05_seurat), ncol(JK4488_06_seurat))
)

print(summary_df)

options(repr.plot.width = 10, repr.plot.height = 8)

# (1) JK4488_01
cat("\nJK4488_01:\n")
cat(sprintf("  Genes: %d\n", nrow(JK4488_01_seurat)))
cat(sprintf("  Cells: %d\n", ncol(JK4488_01_seurat)))
cat(sprintf("  UMAP reductions: %s\n", paste(Reductions(JK4488_01_seurat), collapse = ", ")))
DimPlot(JK4488_01_seurat, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("JK4488_01 – UMAP by clusters")

# (2) JK4488_02
cat("\nJK4488_02:\n")
cat(sprintf("  Genes: %d\n", nrow(JK4488_02_seurat)))
cat(sprintf("  Cells: %d\n", ncol(JK4488_02_seurat)))
cat(sprintf("  UMAP reductions: %s\n", paste(Reductions(JK4488_02_seurat), collapse = ", ")))
DimPlot(JK4488_02_seurat, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("JK4488_02 – UMAP by clusters")

# (3) JK4488_03
cat("\nJK4488_03:\n")
cat(sprintf("  Genes: %d\n", nrow(JK4488_03_seurat)))
cat(sprintf("  Cells: %d\n", ncol(JK4488_03_seurat)))
cat(sprintf("  UMAP reductions: %s\n", paste(Reductions(JK4488_03_seurat), collapse = ", ")))
DimPlot(JK4488_03_seurat, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("JK4488_03 – UMAP by clusters")

# (4) JK4488_04
cat("\nJK4488_04:\n")
cat(sprintf("  Genes: %d\n", nrow(JK4488_04_seurat)))
cat(sprintf("  Cells: %d\n", ncol(JK4488_04_seurat)))
cat(sprintf("  UMAP reductions: %s\n", paste(Reductions(JK4488_04_seurat), collapse = ", ")))
DimPlot(JK4488_04_seurat, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("JK4488_04 – UMAP by clusters")

# (5) JK4488_05
cat("\nJK4488_05:\n")
cat(sprintf("  Genes: %d\n", nrow(JK4488_05_seurat)))
cat(sprintf("  Cells: %d\n", ncol(JK4488_05_seurat)))
cat(sprintf("  UMAP reductions: %s\n", paste(Reductions(JK4488_05_seurat), collapse = ", ")))
DimPlot(JK4488_05_seurat, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("JK4488_05 – UMAP by clusters")

# (6) JK4488_06
cat("\nJK4488_06:\n")
cat(sprintf("  Genes: %d\n", nrow(JK4488_06_seurat)))
cat(sprintf("  Cells: %d\n", ncol(JK4488_06_seurat)))
cat(sprintf("  UMAP reductions: %s\n", paste(Reductions(JK4488_06_seurat), collapse = ", ")))
DimPlot(JK4488_06_seurat, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("JK4488_06 – UMAP by clusters")


# Show summaries
show_summary(JK4488_01_seurat, "JK4488_01")
show_summary(JK4488_02_seurat, "JK4488_02")
show_summary(JK4488_03_seurat, "JK4488_03")
show_summary(JK4488_04_seurat, "JK4488_04")
show_summary(JK4488_05_seurat, "JK4488_05")
show_summary(JK4488_06_seurat, "JK4488_06")

JK4488_01_seurat$orig.ident = "JK4488_01"
JK4488_02_seurat$orig.ident = "JK4488_02"
JK4488_03_seurat$orig.ident = "JK4488_03"
JK4488_04_seurat$orig.ident = "JK4488_04"
JK4488_05_seurat$orig.ident = "JK4488_05"
JK4488_06_seurat$orig.ident = "JK4488_06"

# Create list of samples
seurat_list <- list(
  JK4488_01 = JK4488_01_seurat,
  JK4488_02 = JK4488_02_seurat,
  JK4488_03 = JK4488_03_seurat,
  JK4488_04 = JK4488_04_seurat,
  JK4488_05 = JK4488_05_seurat,
  JK4488_06 = JK4488_06_seurat
)

# Add sample ID to metadata
for (i in seq_along(seurat_list)) {
  seurat_list[[i]]$sample <- names(seurat_list)[i]
}

print(seurat_list)
print(seurat_list["JK4488_01"])
print(seurat_list["JK4488_02"])
print(seurat_list["JK4488_03"])
print(seurat_list["JK4488_04"])
print(seurat_list["JK4488_05"])
print(seurat_list["JK4488_06"])

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
DimPlot(merged, reduction = "umap.unintegrated", group.by = "orig.ident")

DimPlot(
  merged,
  reduction = "umap.unintegrated",
  group.by  = "clusters_unintegrated_res0.1",
  label = TRUE, repel = TRUE, raster = TRUE, pt.size = 0.2
)

DimPlot(
  merged,
  reduction = "umap.unintegrated",
  group.by  = "clusters_unintegrated_res0.5",
  label = TRUE, repel = TRUE, raster = TRUE, pt.size = 0.2
)

DimPlot(
  merged,
  reduction = "umap.unintegrated",
  group.by  = "clusters_unintegrated_res1.0",
  label = TRUE, repel = TRUE, raster = TRUE, pt.size = 0.2
)




# Perform streamlined (one-line) integrative analysis
# https://satijalab.org/seurat/articles/seurat5_integration

# Anchor-based CCA integration (method=CCAIntegration)
# Anchor-based RPCA integration (method=RPCAIntegration)
# Harmony (method=HarmonyIntegration)
# FastMNN (method= FastMNNIntegration)
# scVI (method=scVIIntegration)

obj = merged

obj <- IntegrateLayers(
  object = obj, 
  method = CCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.cca",
  verbose = FALSE
)

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

library(future)
library("batchelor")
library("scuttle")
library(SeuratWrappers)

options(future.globals.maxSize = 256 * 1024^3)  # 128 GB

# plan(sequential)
obj <- IntegrateLayers(
  object = obj, 
  method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)
# plan(multicore)  

# slotNames(obj)
(obj@reductions)

getwd()



obj <- FindNeighbors(obj, 
                     reduction = "integrated.cca", 
                     graph.name = "graph_cca",
                     dims = 1:30)

obj <- FindClusters(obj, 
                    resolution = 0.1, 
                    graph.name   = "graph_cca",
                    cluster.name = "cca_clusters_res0.1")

obj <- FindClusters(obj, 
                    resolution = 0.5, 
                    graph.name   = "graph_cca",
                    cluster.name = "cca_clusters_res0.5")

obj <- FindClusters(obj, 
                    resolution = 1, 
                    graph.name   = "graph_cca",
                    cluster.name = "cca_clusters_res1.0")


table(obj$cca_clusters_res0.1)
table(obj$cca_clusters_res0.5)
table(obj$cca_clusters_res1.0)

head(obj@meta.data,3)
Reductions(obj)

# make a UMAP from the integrated CCA space and name it explicitly
obj <- RunUMAP(
  obj,
  reduction      = "integrated.cca",  # your integrated space
  dims           = 1:30,
  reduction.name = "umap.cca",        # choose a name
  verbose        = FALSE
)

# now plot using that reduction and your cluster columns
DimPlot(obj, reduction = "umap.cca", group.by = "cca_clusters_res0.1", label = TRUE)
DimPlot(obj, reduction = "umap.cca", group.by = "cca_clusters_res0.5", label = TRUE)
DimPlot(obj, reduction = "umap.cca", group.by = "cca_clusters_res1.0",   label = TRUE)



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
DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters_res0.1", label = TRUE)
DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters_res0.5", label = TRUE)
DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters_res1.0",   label = TRUE)


# now plot using that reduction and your cluster columns
DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters_res0.2", label = TRUE)
DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters_res0.3", label = TRUE)
DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters_res0.4",   label = TRUE)



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
DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters_res0.1", label = TRUE)
DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters_res0.5", label = TRUE)
DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters_res1.0",   label = TRUE)


# now plot using that reduction and your cluster columns
DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters_res0.2", label = TRUE)
DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters_res0.3", label = TRUE)
DimPlot(obj, reduction = "umap.harmony", group.by = "harmony_clusters_res0.4",   label = TRUE)





obj <- FindNeighbors(obj,
                     reduction = "integrated.mnn",
                     graph.name = "graph_mnn",
                     dims = 1:30)

obj <- FindClusters(obj,
                    resolution = 0.1,
                    graph.name   = "graph_mnn",
                    cluster.name = "mnn_clusters_res0.1")

obj <- FindClusters(obj,
                    resolution = 0.5,
                    graph.name   = "graph_mnn",
                    cluster.name = "mnn_clusters_res0.5")

obj <- FindClusters(obj,
                    resolution = 1,
                    graph.name   = "graph_mnn",
                    cluster.name = "mnn_clusters_res1.0")

table(obj$mnn_clusters_res0.1)
table(obj$mnn_clusters_res0.5)
table(obj$mnn_clusters_res1.0)

head(obj@meta.data,3)
Reductions(obj)

# make a UMAP from the integrated CCA space and name it explicitly
obj <- RunUMAP(
  obj,
  reduction      = "integrated.mnn",  # your integrated space
  dims           = 1:30,
  reduction.name = "umap.mnn",        # choose a name
  verbose        = FALSE
)

# now plot using that reduction and your cluster columns
DimPlot(obj, reduction = "umap.mnn", group.by = "mnn_clusters_res0.1", label = TRUE)
DimPlot(obj, reduction = "umap.mnn", group.by = "mnn_clusters_res0.5", label = TRUE)
DimPlot(obj, reduction = "umap.mnn", group.by = "mnn_clusters_res1.0",   label = TRUE)

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
DimPlot(obj, reduction = "umap.scvi", group.by = "scvi_clusters_res0.1", label = TRUE)
DimPlot(obj, reduction = "umap.scvi", group.by = "scvi_clusters_res0.5", label = TRUE)
DimPlot(obj, reduction = "umap.scvi", group.by = "scvi_clusters_res1.0",   label = TRUE)

# now plot using that reduction and your cluster columns
DimPlot(obj, reduction = "umap.scvi", group.by = "scvi_clusters_res0.2", label = TRUE)
DimPlot(obj, reduction = "umap.scvi", group.by = "scvi_clusters_res0.3", label = TRUE)
DimPlot(obj, reduction = "umap.scvi", group.by = "scvi_clusters_res0.4",   label = TRUE)



head(obj@meta.data)

# Save any R object (e.g., a Seurat or SCE object)
saveRDS(obj, file = "JK.samples.merged_seurat.BT.analysis.sep30.v3.rds")

# Load it back later
# obj <- readRDS("obj.rds")





# https://github.com/cbib/Seurat-Integrate/blob/main/vignettes/SeuratIntegrate.Rmd
