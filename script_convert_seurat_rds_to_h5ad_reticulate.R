# script_convert_seurat_to_h5ad.R

library(reticulate)
use_virtualenv('/mnt/nfs/CX000008_DS1/projects/btanasa/scrna_env', required = TRUE)

library(Seurat)
library(sceasy)
library(reticulate)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if input file is provided
if (length(args) == 0) {
  stop("Usage: Rscript script_convert_seurat_to_h5ad.R <input_file.rds>")
}

input_file <- args[1]

# Extract sample name (everything before "_seurat_obj")
sample_name <- sub("_seurat_obj.*", "", basename(input_file))

message("Sample name: ", sample_name)
message("Input file: ", input_file)

# Read Seurat object
seurat_obj <- readRDS(input_file)

# Convert Assay5 to Assay (v3) format
seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay")

# Define output file
output_file <- paste0(sample_name, "_seurat_obj.filtered.soupX.h5ad")

# Convert
sceasy::convertFormat(
  seurat_obj,
  from = "seurat",
  to = "anndata",
  outFile = output_file
)

message("Saved: ", output_file)
