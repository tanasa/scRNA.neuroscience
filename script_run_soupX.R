library(SoupX)
library(Matrix)

# Load the data
sc <- load10X('../results_cellranger/outs/')

# Automatically estimate contamination
sc <- autoEstCont(sc)

# Get the corrected counts
out <- adjustCounts(sc)

# Save as RDS
saveRDS(out, "soupx_corrected_counts.rds")

# Save features with your specified name
features <- data.frame(
    gene_id = rownames(out),
    gene_symbol = rownames(out),
    gene_type = "Gene Expression"
)
write.table(features, "soupx_features.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Save barcodes with your specified name
write.table(colnames(out), "soupx_barcodes.tsv", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

print("SoupX correction complete!")
print("Files saved:")
print("- soupx_corrected_counts.rds")
print("- soupx_features.tsv")
print("- soupx_barcodes.tsv")
