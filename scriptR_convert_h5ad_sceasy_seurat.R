############################################
############################################

# convert h5ad into sce or seurat object

############################################ method 1
############################################ sceasy
# https://github.com/cellgeni/sceasy

library(sceasy)
library(reticulate)
use_condaenv('single-cell')
loompy <- reticulate::import('loompy')

library("Seurat")

list.files()
h5ad_file = "hypoxia_anndata_annotated.h5ad"

sceasy::convertFormat(h5ad_file, 
                      from="anndata", to="seurat",
                      outFile='hypoxia_anndata_annotated.sceasy.rds')

# Loading required namespace: Seurat
# X -> counts
# An object of class Seurat
# 34098 features across 104438 samples within 1 assay
# Active assay: RNA (34098 features, 0 variable features)
# 2 layers present: counts, data
# 3 dimensional reductions calculated: pca, scVI, umap

so = readRDS("hypoxia_anndata_annotated.sceasy.rds")

# An object of class Seurat
# 34098 features across 104438 samples within 1 assay
# Active assay: RNA (34098 features, 0 variable features)
# 2 layers present: counts, data
# 3 dimensional reductions calculated: pca, scVI, umap

############################################ method 2 :
############################################ schard

library("schard")

hypoxia.sce = schard::h5ad2sce("hypoxia_anndata_annotated.h5ad")
hypoxia.seurat = schard::h5ad2seurat("hypoxia_anndata_annotated.h5ad")

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
options(encoding = "UTF-8")

saveRDS(hypoxia.sce, file = "hypoxia_anndata_annotated.scahrd.sce.rds")
saveRDS(hypoxia.seurat, file = "hypoxia_anndata_annotated.scahrd.seurat.rds")

Layers(hypoxia.seurat)
Reductions(hypoxia.seurat)

############################################ method 3 :
############################################ zellkonverter
############################################ creates a new environment

#!/usr/bin/env Rscript

## ===== 0) (Optional) Fresh basilisk env for zellkonverter =====
# Comment these two lines out if you don't need a clean rebuild.
Sys.setenv(BASILISK_USE_SYSTEM_DIR = "FALSE", BASILISK_VERBOSE = "1")
try(unlink("~/.cache/R/basilisk", recursive = TRUE, force = TRUE), silent = TRUE)

## ===== 1) Libraries =====

suppressPackageStartupMessages({
  library(zellkonverter)         # readH5AD
  library(SingleCellExperiment)  # SCE container
  library(Seurat)                # Seurat conversion
})

## ===== 2) Read .h5ad into SCE =====
infile <- "hypoxia_anndata_annotated.h5ad"
message("Reading AnnData: ", infile)
sce <- readH5AD(infile)

## ===== 3) Sanitize colData (avoid factor headaches on conversion) =====
if (ncol(colData(sce)) > 0) {
  for (nm in colnames(colData(sce))) {
    x <- colData(sce)[[nm]]
    if (is.factor(x)) colData(sce)[[nm]] <- as.character(x)
  }
}

## ===== 4) Convert SCE -> Seurat =====
# Prefer counts layer if present; fall back to defaults otherwise.
if ("counts" %in% assayNames(sce)) {
  seurat_obj <- as.Seurat(sce, counts = "counts", data = NULL)
} else {
  seurat_obj <- as.Seurat(sce)
}

## ===== 5) Replace/rename assay: originalexp -> RNA (safely) =====
safe_replace_assay <- function(obj, from = "originalexp", to = "RNA") {
  assays_present <- names(obj@assays)

  if (!from %in% assays_present) {
    message(sprintf("Assay '%s' not found. Available assays: %s",
                    from, paste(assays_present, collapse = ", ")))
    # Still set a sensible default assay if RNA exists
    if ("RNA" %in% assays_present) {
      DefaultAssay(obj) <- "RNA"
    }
    return(obj)
  }

  # If target exists already, drop it to avoid duplication/conflict
  if (to %in% assays_present) {
    message(sprintf("Assay '%s' already exists; removing it before replacement.", to))
    obj@assays[[to]] <- NULL
  }

  # Copy and then remove the original
  obj[[to]] <- obj[[from]]
  obj@assays[[from]] <- NULL

  # Set default assay
  SeuratObject::DefaultAssay(obj) <- to
  message(sprintf("Replaced assay '%s' -> '%s' and set DefaultAssay('%s').", from, to, to))
  obj
}

seurat_obj <- safe_replace_assay(seurat_obj, from = "originalexp", to = "RNA")

## ===== 6) Light sanity print =====
message("Assays now in object: ", paste(names(seurat_obj@assays), collapse = ", "))
message("Default assay: ", SeuratObject::DefaultAssay(seurat_obj))
message("Cells: ", ncol(seurat_obj), " | Features (RNA): ",
        if ("RNA" %in% names(seurat_obj@assays)) nrow(seurat_obj@assays$RNA) else NA)

## ===== 7) Save result =====
saveRDS(seurat_obj, file = "hypoxia_anndata_annotated.zellkonverter.seurat.rds")

############################################ method 4
############################################ SeuratDisk
############################################ 
# ============================
# AnnData (.h5ad) → Seurat    |
# Robust loader (SeuratDisk)  |
# ============================

# ---- 0) Setup
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

# ---- 1) Paths
in_h5ad  <- "hypoxia_anndata_annotated.h5ad"
tmp_h5s  <- sub("\\.h5ad$", ".h5seurat", in_h5ad)
out_rds  <- "hypoxia_anndata_annotated.seuratdisk.seurat.rds"

# ---- 2) Convert .h5ad → .h5seurat
# If your AnnData has raw counts in adata.raw or a layer named "counts",
# setting X_name="counts" helps SeuratDisk pick the right matrix as counts.
Convert(in_h5ad,
        dest      = "h5seurat",
        overwrite = TRUE,
        assay     = "RNA",
        X_name    = "counts")

# ---- 3) Load .h5seurat into Seurat with robust metadata handling
# Try normal load first; if 'levels/values' metadata error occurs,
# fall back to meta.data = FALSE to skip problematic obs categoricals.
load_ok <- TRUE
seurat_obj <- tryCatch({
  LoadH5Seurat(tmp_h5s, assays = "RNA", meta.data = TRUE)
}, error = function(e) {
  message("⚠️  Loading with metadata failed (likely obs categoricals). Retrying without metadata...")
  load_ok <<- FALSE
  LoadH5Seurat(tmp_h5s, assays = "RNA", meta.data = FALSE)
})

# ---- 4) Fix assay name if needed (your working approach)
# Some conversions name the sole assay 'originalexp'. Make sure we have an assay called 'RNA'.
assays_now <- Seurat::Assays(seurat_obj)
if (!"RNA" %in% assays_now) {
  if ("originalexp" %in% assays_now) {
    # Your working solution
    seurat_obj[["RNA"]] <- seurat_obj[["originalexp"]]  # copy
    seurat_obj@assays$originalexp <- NULL               # drop old
  } else {
    stop("No 'RNA' or 'originalexp' assay found in the loaded object. Assays present: ",
         paste(assays_now, collapse = ", "))
  }
}

# Set default assay
SeuratObject::DefaultAssay(seurat_obj) <- "RNA"

# ---- 5) Standardize reduction keys (optional, silences those key warnings)
if ("pca"  %in% Seurat::Reductions(seurat_obj)) Seurat::Key(seurat_obj[["pca"]])  <- "pca_"
if ("scVI" %in% Seurat::Reductions(seurat_obj)) Seurat::Key(seurat_obj[["scVI"]]) <- "scvi_"
if ("umap" %in% Seurat::Reductions(seurat_obj)) Seurat::Key(seurat_obj[["umap"]]) <- "umap_"

# ---- 6) Sanity checks (layers, dims)
message("Assays: ", paste(Seurat::Assays(seurat_obj), collapse = ", "))
message("Layers in RNA: ", paste0(tryCatch(SeuratObject::Layers(seurat_obj[["RNA"]]),
                                           error = function(e) "unknown"), collapse = ", "))
message("Cells: ", ncol(seurat_obj), " | Features: ", nrow(seurat_obj))

# ---- 7) (Optional) Recreate nFeature/nCount if missing
if (!all(c("nFeature_RNA","nCount_RNA") %in% colnames(seurat_obj@meta.data))) {
  mat <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  seurat_obj$`nFeature_RNA` <- Matrix::colSums(mat > 0)
  seurat_obj$`nCount_RNA`   <- Matrix::colSums(mat)
}

# ---- 8) Save RDS (avoid copy/paste Unicode issues by building the string)
outfile <- paste0("hypoxia_anndata_annotated", ".SeuratDisk.seurat", ".rds")
saveRDS(seurat_obj, file = outfile)
message("✅ Saved: ", outfile)

# ---- 9) If metadata was skipped, optional: attach a cleaned obs later
if (!load_ok) {
  message("ℹ️  Note: cell-level metadata (obs) was skipped. ",
          "If you need it, export obs from Python (categoricals → strings) ",
          "and add via AddMetaData() matched by cell barcodes.")
}

############################################ 
############################################ 
