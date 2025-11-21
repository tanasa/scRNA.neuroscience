library(SoupX)
library(Matrix)

# Input path
input_path <- './AD_AG_2658/outs/'

# Extract sample name: "AD_AG_2658"
sample_name <- basename(dirname(input_path))

cat("Sample detected:", sample_name, "\n")

# Load data
sc <- load10X(input_path)

# Estimate contamination
sc <- autoEstCont(sc)

# Adjust counts
out <- adjustCounts(sc)

# Build filenames: SAMPLE.soupx_XXXX.ext
rds_file      <- paste0(sample_name, ".soupx_corrected_counts.rds")
features_file <- paste0(sample_name, ".soupx_features.tsv")
barcodes_file <- paste0(sample_name, ".soupx_barcodes.tsv")

# Save corrected matrix
saveRDS(out, rds_file)

# Save features
features <- data.frame(
    gene_id     = rownames(out),
    gene_symbol = rownames(out),
    gene_type   = "Gene Expression"
)
write.table(features, features_file,
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

# Save barcodes
write.table(colnames(out), barcodes_file,
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

cat("SoupX correction complete!\nFiles saved:\n")
cat("- ", rds_file, "\n")
cat("- ", features_file, "\n")
cat("- ", barcodes_file, "\n")
