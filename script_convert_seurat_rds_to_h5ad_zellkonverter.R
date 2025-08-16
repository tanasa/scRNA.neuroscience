# install.packages("BiocManager")
# BiocManager::install(c("zellkonverter","SingleCellExperiment"))
# install.packages("Seurat")  # if not already

library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

files <- c(
  "hypoxia_anndata_annotated.schard.seurat.AP4765_01.rds",
  "hypoxia_anndata_annotated.schard.seurat.AP4765_02.rds",
  "hypoxia_anndata_annotated.schard.seurat.AP4765_03.rds",
  "hypoxia_anndata_annotated.schard.seurat.AP4765_04.rds",
  "hypoxia_anndata_annotated.schard.seurat.AP4765_05.rds",
  "hypoxia_anndata_annotated.schard.seurat.AP4765_06.rds",
  "hypoxia_anndata_annotated.schard.seurat.AP4765_07.rds",
  "hypoxia_anndata_annotated.schard.seurat.AP4765_08.rds"
)

convert_one <- function(in_rds) {
  message("Converting: ", in_rds)
  out_h5ad <- sub("\\.rds$", ".zellk.h5ad", in_rds)

  obj <- readRDS(in_rds)

  # Convert Seurat -> SingleCellExperiment
  sce <- as.SingleCellExperiment(obj)

  # Write H5AD (AnnData >= 0.10). 
  # By default, zellkonverter writes a sensible X; you can also set X_name="counts" if desired:

  writeH5AD(sce, file = out_h5ad, X_name = "counts")

  # writeH5AD(sce, file = out_h5ad)


  message("Saved: ", out_h5ad)
}

invisible(lapply(files, convert_one))
