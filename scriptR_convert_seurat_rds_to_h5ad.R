# script_convert_sceasy_seurat_rds_to_h5ad.R 
# conda activate single-cell

library(Seurat)
library(sceasy)
library(reticulate)

use_condaenv('single-cell')
loompy <- reticulate::import('loompy')

# List of input files
files <- c(
  "hypoxia_anndata_annotated.sceasy.seurat.AP4765_01.rds",
  "hypoxia_anndata_annotated.sceasy.seurat.AP4765_02.rds",
  "hypoxia_anndata_annotated.sceasy.seurat.AP4765_03.rds",
  "hypoxia_anndata_annotated.sceasy.seurat.AP4765_04.rds",
  "hypoxia_anndata_annotated.sceasy.seurat.AP4765_05.rds",
  "hypoxia_anndata_annotated.sceasy.seurat.AP4765_06.rds",
  "hypoxia_anndata_annotated.sceasy.seurat.AP4765_07.rds",
  "hypoxia_anndata_annotated.sceasy.seurat.AP4765_08.rds"
)

# Loop through files and convert each one
for (f in files) {
  message("Converting: ", f)

  # Read Seurat object
  seurat_obj <- readRDS(f)

  # Define output file (replace .rds with .h5ad)
  outFile <- sub("\\.rds$", ".h5ad", f)

  # Convert
  sceasy::convertFormat(
    seurat_obj,
    from = "seurat",
    to = "anndata",
    outFile = outFile
  )

  message("Saved: ", outFile)
}

