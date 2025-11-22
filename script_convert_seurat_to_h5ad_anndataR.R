# script_convert_seurat_to_h5ad_anndataR.R

library(Seurat)
library(anndataR)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript script_convert_seurat_to_h5ad_anndataR.R <input_file.rds>")
}

input_file <- args[1]

# Extract sample name
sample_name <- sub("_seurat_obj.*", "", basename(input_file))

message("========================================")
message("Sample name: ", sample_name)
message("Input file: ", input_file)
message("========================================")

# Read Seurat object
message("Reading Seurat object...")
seurat_obj <- readRDS(input_file)

# Convert to AnnData
message("Converting Seurat to AnnData...")
adata <- as_AnnData(seurat_obj)

# Define output file with anndataR suffix
output_file <- paste0(sample_name, "_seurat_obj.filtered.soupX.anndataR.h5ad")

# Write to h5ad
message("Writing to h5ad...")
adata$write_h5ad(output_file)

message("========================================")
message("Saved: ", output_file)
message("========================================")
