# library(hdf5r)
library(SoupX)
library(Matrix)
library(Seurat) 
library(DropletUtils)

# Define paths to H5 files
setwd("/mnt/nfs/CX000008_DS1/projects/jaeyeon/fastq_file_Dyslexia_r1/zanalysis_bogdan/fastq_file_Dyslexia_r2_cellRanger/")
list.files()

sample_name <- "AD_AG_2658"
sample_id = sample_name

# Navigate to sample folder
sample_folder <- file.path(getwd(), sample_name)

# Check if folder exists
if (!dir.exists(sample_folder)) {
  stop("Sample folder does not exist: ", sample_folder)
}

# Change to sample folder
setwd(sample_folder)

# Confirm current directory
cat("Current working directory:\n", getwd(), "\n")

# List files in the directory
cat("\nFiles in sample folder:\n")
print(list.files())

# Create output directory for analysis results
output_dir <- file.path(getwd(), paste0(sample_id, ".soupX.gs.analysis")) # gs stands for "gene symbols", no "girls scouts"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE)
  cat("\nCreated output folder:\n", output_dir, "\n")
} else {
  cat("\nOutput folder already exists:\n", output_dir, "\n")
}

# Helper function for output filenames (alternative approach)
outfn <- function(name) {
  filename <- paste(sample_id, name, sep = "_")
  file.path(output_dir, filename)
}

cat("Output files will be saved with prefix:", sample_id, "\n\n")

print(output_dir)

# Now in: /mnt/nfs/.../fastq_file_Dyslexia_r2_cellRanger/AD_AG_2658/

# Define the soupx_corrected folder
soupx_folder <- paste0(sample_name, "_soupx_corrected")

# Check if folder exists
if (!dir.exists(soupx_folder)) {
  stop("SoupX folder does not exist: ", file.path(getwd(), soupx_folder))
}

cat("SoupX folder:", file.path(getwd(), soupx_folder), "\n")

# Define the RDS file path
rds_file <- file.path(soupx_folder, paste0(sample_name, "_soupx_corrected.rds"))

# Check if file exists
if (!file.exists(rds_file)) {
  stop("RDS file not found: ", rds_file)
}

cat("Reading RDS file:", rds_file, "\n\n")

# Read the RDS file
obj <- readRDS(rds_file)

# Check what was loaded
cat("Object class:", class(obj), "\n")
str(obj, max.level = 1)

# Handle different object types
if (!inherits(obj, "Seurat")) {
  if (inherits(obj, "dgCMatrix") || is.matrix(obj)) {
    message("RDS contains a matrix; creating a Seurat object...")
    obj <- CreateSeuratObject(counts = obj, project = sample_name)
  } else {
    stop("RDS is neither a Seurat object nor a matrix. Class: ",
         paste(class(obj), collapse = ", "))
  }
}

cat("\nSeurat object successfully created!\n")
cat("Cells:", ncol(obj), "| Genes:", nrow(obj), "\n")

# Preview counts
counts_matrix <- GetAssayData(object = obj, assay = "RNA", slot = "counts")
print(counts_matrix[1:5, 1:5])

head(obj@meta.data)

cat("\nFirst 20 gene names:\n")
print(head(rownames(obj), 20))

cat("\nLast 20 gene names:\n")
print(tail(rownames(obj), 20))

head(symbols)

library(org.Hs.eg.db)
library(AnnotationDbi)

# Get Ensembl IDs from Seurat object
ensembl_ids <- rownames(obj)

# Clean version numbers (ENSG00000290825.1 -> ENSG00000290825)
ensembl_ids_clean <- sub("\\.\\d+$", "", ensembl_ids)

# Map to gene symbols
symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Create mapping dataframe
gene_mapping <- data.frame(
  ensembl_id = ensembl_ids,
  ensembl_id_clean = ensembl_ids_clean,
  gene_symbol = symbols,
  stringsAsFactors = FALSE
)

head(gene_mapping, 5)
head(gene_mapping, 5)

# Replace NA symbols with original Ensembl ID
gene_mapping$gene_symbol[is.na(gene_mapping$gene_symbol)] <- gene_mapping$ensembl_id[is.na(gene_mapping$gene_symbol)]

# View the mapping
head(gene_mapping, 3)

# View the mapping
head(gene_mapping, 3)

# Define output file name
output_file_gene_mapping <- file.path(output_dir, "gene_mapping_table.csv")

# Write the dataframe to CSV
write.csv(gene_mapping, file = output_file_gene_mapping, row.names = FALSE)

cat("✅ Gene mapping table saved to:", output_file_gene_mapping, "\n")
getwd()

# Summary
cat("\nMapping summary:\n")
cat("Total genes:", nrow(gene_mapping), "\n")
cat("Mapped to symbols:", sum(!is.na(symbols)), "\n")
cat("Kept as Ensembl:", sum(is.na(symbols)), "\n")

# Find duplicate gene symbols
duplicate_symbols <- gene_mapping$gene_symbol[duplicated(gene_mapping$gene_symbol) | 
                                               duplicated(gene_mapping$gene_symbol, fromLast = TRUE)]

cat("\n=== DUPLICATE CHECK ===\n")
cat("Total genes:", nrow(gene_mapping), "\n")
cat("Unique gene symbols:", length(unique(gene_mapping$gene_symbol)), "\n")
cat("Number of duplicates:", length(unique(duplicate_symbols)), "\n\n")


# Use gene symbols as new rownames (make unique if duplicates exist)
# new_rownames <- make.unique(gene_mapping$gene_symbol)
# rownames(obj) <- new_rownames

# Verify
# cat("\nExample conversions:\n")
# print(head(gene_mapping[!is.na(symbols), ], 10))

# code to keep ENSEMBL ID for ALL duplicate gene symbols (not just subsequent occurrences):

if (length(duplicate_symbols) > 0) {
  cat("Duplicate gene symbols found:\n")
  
  # Show all duplicates with their Ensembl IDs
  duplicates_df <- gene_mapping[gene_mapping$gene_symbol %in% unique(duplicate_symbols), ]
  duplicates_df <- duplicates_df[order(duplicates_df$gene_symbol), ]
  
  print(duplicates_df)
  
  # Summary by duplicate
  cat("\n=== Duplicates Summary ===\n")
  dup_summary <- table(duplicates_df$gene_symbol)
  print(dup_summary[dup_summary > 1])
  
  # Handle duplicates: Keep ENSEMBL ID for ALL duplicates
  cat("\n=== Handling Duplicates ===\n")
  cat("Keeping Ensembl ID for all duplicate gene symbols...\n")
  
  # Initialize final_name with gene_symbol
  gene_mapping$final_name <- gene_mapping$gene_symbol
  
  # For ALL duplicates (not just subsequent ones), use Ensembl ID
  is_duplicate <- gene_mapping$gene_symbol %in% unique(duplicate_symbols)
  gene_mapping$final_name[is_duplicate] <- gene_mapping$ensembl_id[is_duplicate]
  
  cat("Replaced", sum(is_duplicate), "duplicate symbols with their Ensembl IDs\n")
  
} else {
  cat("✓ No duplicate gene symbols found!\n")
  gene_mapping$final_name <- gene_mapping$gene_symbol
}

# After handling duplicates, show the results

cat("\n=== DUPLICATED GENES - FINAL MAPPING ===\n\n")

if (length(duplicate_symbols) > 0) {
  # Get all rows with duplicate symbols
  duplicates_final <- gene_mapping[gene_mapping$gene_symbol %in% unique(duplicate_symbols), ]
  duplicates_final <- duplicates_final[order(duplicates_final$gene_symbol), ]
  
  # Display the dataframe
  print(duplicates_final)
  
  cat("\n=== Side-by-Side Comparison ===\n")
  comparison <- duplicates_final[, c("ensembl_id", "gene_symbol", "final_name")]
  print(comparison)
  
  cat("\n=== Summary by Original Symbol ===\n")
  for (sym in unique(duplicate_symbols)) {
    subset_df <- duplicates_final[duplicates_final$gene_symbol == sym, 
                                   c("ensembl_id", "gene_symbol", "final_name")]
    cat("\nOriginal symbol:", sym, "\n")
    print(subset_df, row.names = FALSE)
  }
}

# MITOCHINDRIAL GENES
# ENSG00000198888	MT-ND1	Gene Expression
# ENSG00000198763	MT-ND2	Gene Expression
# ENSG00000198804	MT-CO1	Gene Expression
# ENSG00000198712	MT-CO2	Gene Expression
# ENSG00000228253	MT-ATP8	Gene Expression
# ENSG00000198899	MT-ATP6	Gene Expression
# ENSG00000198938	MT-CO3	Gene Expression
# ENSG00000198840	MT-ND3	Gene Expression
# ENSG00000212907	MT-ND4L	Gene Expression
# ENSG00000198886	MT-ND4	Gene Expression
# ENSG00000198786	MT-ND5	Gene Expression
# ENSG00000198695	MT-ND6	Gene Expression
# ENSG00000198727	MT-CYB	Gene Expression

mt_genes <- c(
  "MT-ND1","MT-ND2","MT-CO1","MT-CO2","MT-ATP8","MT-ATP6",
  "MT-CO3","MT-ND3","MT-ND4L","MT-ND4","MT-ND5","MT-ND6","MT-CYB"
)

# Check presence
found <- mt_genes %in% gene_mapping$final_name

# Summary
cat("Total mitochondrial genes:", length(mt_genes), "\n")
cat("Found in gene_mapping$final_name:", sum(found), "\n\n")

# List missing ones (if any)
if (any(!found)) {
  cat("Missing mitochondrial genes:\n")
  print(mt_genes[!found])
} else {
  cat("✅ All mitochondrial genes are present in gene_mapping$final_name!\n")
}


mt_ensembl_ids <- c(
  "ENSG00000198888",
  "ENSG00000198763",
  "ENSG00000198804",
  "ENSG00000198712",
  "ENSG00000228253",
  "ENSG00000198899",
  "ENSG00000198938",
  "ENSG00000198840",
  "ENSG00000212907",
  "ENSG00000198886",
  "ENSG00000198786",
  "ENSG00000198695",
  "ENSG00000198727"
)

# Boolean presence vector
found <- mt_ensembl_ids %in% gene_mapping$ensembl_id

# Summary
cat("Total mitochondrial genes:", length(mt_ensembl_ids), "\n")
cat("Found in ensembl_id:", sum(found), "\n\n")

# Missing ones (if any)
if (any(!found)) {
  cat("❌ Missing mitochondrial genes:\n")
  print(mt_ensembl_ids[!found])
} else {
  cat("✅ All mitochondrial Ensembl genes are present in gene_mapping$ensembl_id!\n")
}

# Optionally print the matches
if (any(found)) {
  cat("\nMatched mitochondrial genes:\n")
  print(gene_mapping[gene_mapping$final_name %in% mt_ensembl_ids, ])
}


# Create dataframe of mitochondrial genes
mt_genes_df <- data.frame(
  ensembl_id = c(
    "ENSG00000198888",
    "ENSG00000198763",
    "ENSG00000198804",
    "ENSG00000198712",
    "ENSG00000228253",
    "ENSG00000198899",
    "ENSG00000198938",
    "ENSG00000198840",
    "ENSG00000212907",
    "ENSG00000198886",
    "ENSG00000198786",
    "ENSG00000198695",
    "ENSG00000198727"
  ),
  gene_symbol_MT = c(
    "MT-ND1",
    "MT-ND2",
    "MT-CO1",
    "MT-CO2",
    "MT-ATP8",
    "MT-ATP6",
    "MT-CO3",
    "MT-ND3",
    "MT-ND4L",
    "MT-ND4",
    "MT-ND5",
    "MT-ND6",
    "MT-CYB"
  ),
  category = "Gene Expression",
  stringsAsFactors = FALSE
)

# Preview
print(mt_genes_df)

# Perfect — you want to update your gene_mapping$final_name for the 13 mitochondrial genes based on your mt_genes_df,
# i.e. whenever gene_mapping$ensembl_id_clean matches mt_genes_df$ensembl_id, set final_name = gene_symbol_MT.

# Make sure mt_genes_df and gene_mapping exist
stopifnot(exists("gene_mapping"), exists("mt_genes_df"))

# Create a lookup named vector: Ensembl → MT gene symbol
mt_lookup <- setNames(mt_genes_df$gene_symbol_MT, mt_genes_df$ensembl_id)

# Find matching indices
matches <- gene_mapping$ensembl_id_clean %in% names(mt_lookup)

cat("Found", sum(matches), "mitochondrial genes in gene_mapping\n")

# Update the final_name for matched genes
gene_mapping$final_name[matches] <- mt_lookup[gene_mapping$ensembl_id_clean[matches]]

# Verify the updates
cat("\n✅ Updated mitochondrial gene names:\n")
print(gene_mapping[matches, c("ensembl_id_clean", "gene_symbol", "final_name")])


# Define key hemoglobin genes (adult + fetal + embryonic)
hb_genes <- c(
  "HBA1", "HBA2", "HBB",  # adult
  "HBD",                  # delta (minor adult HbA2)
  "HBG1", "HBG2",         # fetal
  "HBZ", "HBE1"           # embryonic
)

# Check if they are in gene_mapping$final_name
present <- hb_genes %in% gene_mapping$final_name

# Create a results dataframe
hb_check <- data.frame(
  gene = hb_genes,
  present = present
)

# Print summary
cat("✅ Hemoglobin gene presence summary:\n")
print(hb_check)

cat("\nFound", sum(present), "of", length(hb_genes), "hemoglobin genes in gene_mapping$final_name\n")

# Optionally, list missing ones
if (any(!present)) {
  cat("\n❌ Missing hemoglobin genes:\n")
  print(hb_genes[!present])
}


# Check for genes starting with RPS or RPL (case-sensitive)
ribosomal_matches <- grep("^RPS|^RPL", gene_mapping$final_name, value = TRUE)

# Count them
cat("✅ Total ribosomal genes (RPS/RPL) found:", length(ribosomal_matches), "\n")

# Preview first few
head(ribosomal_matches)



# Display the complete gene mapping dataframe

cat("\n=== COMPLETE GENE METADATA ===\n\n")

# Show first 50 rows
cat("First 50 genes:\n")
print(head(gene_mapping, 50))

cat("\n--- --- ---\n\n")

# Show last 20 rows
cat("Last 20 genes:\n")
print(tail(gene_mapping, 20))

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Total genes:", nrow(gene_mapping), "\n")
cat("Genes with symbols:", sum(gene_mapping$gene_symbol != gene_mapping$ensembl_id), "\n")
cat("Genes kept as Ensembl ID:", sum(gene_mapping$gene_symbol == gene_mapping$ensembl_id), "\n")
cat("Duplicate symbols replaced:", sum(gene_mapping$final_name != gene_mapping$gene_symbol & 
                                        grepl("^ENSG", gene_mapping$final_name)), "\n")

# Check specific genes
cat("\n=== EXAMPLE GENE LOOKUPS ===\n")
cat("\nAPOE:\n")
print(gene_mapping[grep("APOE", gene_mapping$gene_symbol), ])

cat("\nAPP:\n")
print(gene_mapping[grep("^APP$", gene_mapping$gene_symbol), ])

cat("\nPSEN1:\n")
print(gene_mapping[grep("PSEN1", gene_mapping$gene_symbol), ])

cat("\nMitochondrial genes (MT-):\n")
print(head(gene_mapping[grep("^MT-", gene_mapping$gene_symbol), ], 20))

# Save complete mapping to file
# Define proper output path
output_file_gene_mapping <- file.path(output_dir, paste0(sample_id, "_complete_gene_metadata.csv"))

# Save dataframe
write.csv(gene_mapping, file = output_file_gene_mapping, row.names = FALSE)

# Confirmation message
cat("✅ Complete gene metadata saved to:\n", output_file_gene_mapping, "\n")

# Check for duplicates in final_name column

cat("\n=== CHECKING DUPLICATES IN FINAL_NAME ===\n\n")

# Find duplicate final_names
duplicate_final_names <- gene_mapping$final_name[duplicated(gene_mapping$final_name) | 
                                                  duplicated(gene_mapping$final_name, fromLast = TRUE)]

# Check if any duplicates exist
if (length(duplicate_final_names) > 0) {
  cat("⚠️  WARNING: Duplicates found in final_name!\n\n")
  
  # Show all rows with duplicate final_names
  duplicates_final_df <- gene_mapping[gene_mapping$final_name %in% unique(duplicate_final_names), ]
  duplicates_final_df <- duplicates_final_df[order(duplicates_final_df$final_name), ]
  
  cat("Duplicate final_names found:\n")
  print(duplicates_final_df)
  
  # Summary
  cat("\n=== Duplicate Summary ===\n")
  cat("Total unique final_names:", length(unique(gene_mapping$final_name)), "\n")
  cat("Expected unique names:", nrow(gene_mapping), "\n")
  cat("Number of duplicated names:", length(unique(duplicate_final_names)), "\n")
  
  # Count by final_name
  dup_counts <- table(duplicates_final_df$final_name)
  cat("\nDuplication counts:\n")
  print(dup_counts)
  
} else {
  cat("✅ SUCCESS: No duplicates found in final_name!\n\n")
  cat("Total genes:", nrow(gene_mapping), "\n")
  cat("Unique final_names:", length(unique(gene_mapping$final_name)), "\n")
  cat("All gene names are unique! ✓\n")
}

# Verify uniqueness
cat("\n=== VERIFICATION ===\n")
cat("Are all final_names unique?", !any(duplicated(gene_mapping$final_name)), "\n")
cat("Number of rows:", nrow(gene_mapping), "\n")
cat("Number of unique final_names:", length(unique(gene_mapping$final_name)), "\n")

# Check for mitochondrial genes by pattern "MT-" (case-sensitive)
mt_pattern_matches <- grep("^MT-", gene_mapping$final_name, value = TRUE)

# Summary
cat("✅ Mitochondrial gene pattern check (pattern: ^MT-)\n")
cat("Total mitochondrial genes found:", length(mt_pattern_matches), "\n\n")

# Preview
if (length(mt_pattern_matches) > 0) {
  print(mt_pattern_matches)
} else {
  cat("⚠️ No mitochondrial genes with 'MT-' prefix found in gene_mapping$final_name\n")
}

# Optional: case-insensitive check too
mt_pattern_matches_ic <- grep("^MT-", gene_mapping$final_name, value = TRUE, ignore.case = TRUE)
cat("\nIncluding case-insensitive matches:", length(mt_pattern_matches_ic), "\n")

# Compare strict vs relaxed
if (length(mt_pattern_matches_ic) > length(mt_pattern_matches)) {
  cat("⚠️ Some mitochondrial genes appear in lowercase or mixed case.\n")
}


# rownames(obj)

# Assign final_name as rownames in Seurat object

cat("\n=== UPDATING SEURAT OBJECT ROWNAMES ===\n\n")

# Verify final_names are unique
if (any(duplicated(gene_mapping$final_name))) {
  stop("ERROR: final_name contains duplicates! Cannot assign as rownames.")
}

# Verify number of rows match
if (nrow(gene_mapping) != nrow(obj)) {
  stop("ERROR: Number of genes in gene_mapping (", nrow(gene_mapping), 
       ") does not match Seurat object (", nrow(obj), ")")
}

# Store old rownames for reference
# old_rownames <- rownames(obj)

# Assign new rownames
# rownames(obj) <- gene_mapping$final_name

# cat("✅ Rownames successfully updated!\n\n")

# Verification
# cat("=== VERIFICATION ===\n")
# cat("Total genes:", nrow(obj), "\n")
# cat("Unique rownames:", length(unique(rownames(obj))), "\n")
# cat("All rownames unique:", !any(duplicated(rownames(obj))), "\n\n")

# Show examples
# cat("=== EXAMPLES OF NEW ROWNAMES ===\n")
# cat("\nFirst 20 genes:\n")
# print(head(rownames(obj), 20))

# cat("\nLast 20 genes:\n")
# print(tail(rownames(obj), 20))

# Check specific genes
# cat("\n=== KEY GENE CHECKS ===\n")
# cat("APOE present:", "APOE" %in% rownames(obj), "\n")
# cat("APP present:", "APP" %in% rownames(obj), "\n")
# cat("PSEN1 present:", "PSEN1" %in% rownames(obj), "\n")
# cat("PSEN2 present:", "PSEN2" %in% rownames(obj), "\n")

# Show some mitochondrial genes
# cat("\nMitochondrial genes:\n")
# mt_genes <- grep("^MT-", rownames(obj), value = TRUE)
# print(head(mt_genes, 13))

# Show conversion summary
# cat("\n=== CONVERSION SUMMARY ===\n")
# converted <- sum(old_rownames != rownames(obj))
# cat("Genes converted from Ensembl to Symbol:", converted, "\n")
# cat("Genes kept as Ensembl ID:", nrow(obj) - converted, "\n")

# cat("\n✓ Seurat object updated successfully!\n")

# Feature-level metadata lives on the assay:
head(obj[["RNA"]]@meta.data)        # data.frame with one row per feature (gene)
nrow(obj[["RNA"]]@meta.data)        # should equal nrow(obj)

gene_meta <- data.frame(
  gene_id = rownames(obj),        # keep the gene names themselves
  stringsAsFactors = FALSE
)
head(gene_meta, 2)
dim(gene_meta)
dim(gene_mapping)
head(gene_mapping, 2)

# Confirm both dataframes have matching rows
nrow(gene_meta) == nrow(gene_mapping)   # should be TRUE

# Merge gene_meta and gene_mapping

library(dplyr)

# Left join: keep all rows from gene_meta, match on gene_id = final_name
merged_meta <- gene_meta %>%
  left_join(gene_mapping, by = c("gene_id" = "ensembl_id"))

# Verify dimensions
cat("Merged metadata dimensions:", dim(merged_meta), "\n")

# Preview first few rows
head(merged_meta, 5)

dim(merged_meta)

cat("Set Seurat rownames from merged_meta$final_name")

# Safety checks first
stopifnot(nrow(obj) == nrow(merged_meta))

# Confirm that gene order matches
all(rownames(obj) == merged_meta$gene_id)  # should return TRUE

# Replace rownames with the final_name column
rownames(obj) <- merged_meta$final_name

# Verification
cat("✅ Rownames successfully updated!\n")
cat("Total genes:", nrow(obj), "\n")
cat("Unique rownames:", length(unique(rownames(obj))), "\n")
cat("All rownames unique:", !any(duplicated(rownames(obj))), "\n\n")

# Preview
cat("First 10 new rownames:\n")
print(head(rownames(obj), 10))

cat("\nLast 10 new rownames:\n")
print(tail(rownames(obj), 10))

# Confirm mitochondrial and ribosomal features are now readable
grep("^MT-", rownames(obj), value = TRUE)
grep("^RPS|^RPL", rownames(obj), value = TRUE)

# Display the counts matrix with new gene symbols

cat("\n=== COUNTS MATRIX WITH GENE SYMBOLS ===\n\n")

# Get counts matrix
counts_matrix <- GetAssayData(object = obj, assay = "RNA", slot = "counts")

cat("Matrix dimensions:", nrow(counts_matrix), "genes x", ncol(counts_matrix), "cells\n")
cat("Matrix class:", class(counts_matrix), "\n\n")

# Show first 10 genes, first 10 cells
cat("First 10 genes x 10 cells:\n")
print(counts_matrix[1:10, 1:10])

cat("\n--- --- ---\n\n")

# Show specific genes of interest
cat("Key genes (first 10 cells):\n")
genes_of_interest <- c("APOE", "APP", "PSEN1", "PSEN2", "CD3D", "CD79A", "LYZ")
genes_present <- genes_of_interest[genes_of_interest %in% rownames(counts_matrix)]

if (length(genes_present) > 0) {
  print(counts_matrix[genes_present, 1:10])
} else {
  cat("None of the key genes found in the matrix\n")
}

cat("\n--- --- ---\n\n")

# Show mitochondrial genes
cat("Mitochondrial genes (first 10 cells):\n")
mt_genes <- grep("^MT-", rownames(counts_matrix), value = TRUE)
if (length(mt_genes) > 0) {
  print(counts_matrix[head(mt_genes, 13), 1:10])
} else {
  cat("No MT- genes found\n")
}

cat("\n--- --- ---\n\n")

# Show some highly expressed genes
cat("Top 10 highly expressed genes (first 5 cells):\n")
gene_means <- Matrix::rowMeans(counts_matrix)
top_genes <- names(sort(gene_means, decreasing = TRUE)[1:10])
print(counts_matrix[top_genes, 1:5])

cat("\n--- --- ---\n\n")

# Summary statistics
cat("=== MATRIX SUMMARY ===\n")
cat("Total counts:", sum(counts_matrix), "\n")
cat("Mean counts per cell:", mean(Matrix::colSums(counts_matrix)), "\n")
cat("Mean counts per gene:", mean(Matrix::rowSums(counts_matrix)), "\n")
cat("Sparsity:", round(100 * (1 - nnzero(counts_matrix) / length(counts_matrix)), 2), "%\n")

# Show cell barcodes
cat("\nFirst 10 cell barcodes:\n")
print(head(colnames(counts_matrix), 10))

cat("\nLast 10 cell barcodes:\n")
print(tail(colnames(counts_matrix), 10))

# Distribution of counts for all genes (Claude)

cat("\n=== GENE COUNTS DISTRIBUTION ANALYSIS ===\n\n")

# Get counts matrix
counts_matrix <- GetAssayData(object = obj, assay = "RNA", slot = "counts")

# Calculate gene-level statistics
gene_total_counts <- Matrix::rowSums(counts_matrix)
gene_mean_counts <- Matrix::rowMeans(counts_matrix)
gene_detection_rate <- Matrix::rowMeans(counts_matrix > 0) * 100  # % of cells expressing
gene_max_counts <- apply(counts_matrix, 1, max)

# Create summary dataframe
gene_stats <- data.frame(
  gene = rownames(counts_matrix),
  total_counts = gene_total_counts,
  mean_counts = gene_mean_counts,
  detection_rate = gene_detection_rate,
  max_counts = gene_max_counts,
  stringsAsFactors = FALSE
)

# Summary statistics
cat("=== SUMMARY STATISTICS ===\n")
summary(gene_stats[, -1])  # Exclude gene name column

cat("\n=== DISTRIBUTION BREAKDOWN ===\n")
cat("Genes with 0 counts:", sum(gene_stats$total_counts == 0), "\n")
cat("Genes with >0 counts:", sum(gene_stats$total_counts > 0), "\n")
cat("Genes detected in >50% cells:", sum(gene_stats$detection_rate > 50), "\n")
cat("Genes detected in >10% cells:", sum(gene_stats$detection_rate > 10), "\n")
cat("Genes detected in <1% cells:", sum(gene_stats$detection_rate < 1), "\n")

# Top 20 most expressed genes
cat("\n=== TOP 20 MOST EXPRESSED GENES ===\n")
top20 <- gene_stats[order(-gene_stats$total_counts), ]
print(head(top20, 20))

# Top 20 most detected genes
cat("\n=== TOP 20 MOST DETECTED GENES (% cells) ===\n")
top20_detected <- gene_stats[order(-gene_stats$detection_rate), ]
print(head(top20_detected, 20))

# Additional functions from scater packages

library(Seurat)
library(SingleCellExperiment)
library(scater)

sce <- as.SingleCellExperiment(obj)

sce
assayNames(sce)
colnames(colData(sce))
head(rowData(sce))

# Add QC metrics if not yet done
library(scater)
sce <- addPerCellQC(sce)

# Plot two histograms side-by-side
par(mfrow = c(1, 2))

# Histogram of total counts per cell (library size)
hist(
  sce$sum / 1e6,
  xlab = "Library sizes (millions)",
  ylab = "Number of cells",
  main = "",
  breaks = 20,
  col = "grey80",
  border = "white"
)

# Histogram of number of detected features per cell
hist(
  sce$detected,
  xlab = "Number of expressed genes",
  ylab = "Number of cells",
  main = "",
  breaks = 20,
  col = "grey80",
  border = "white"
)


getwd()


library(ggplot2)

# 1️⃣ Total counts per cell
p_counts <- ggplot(as.data.frame(colData(sce)), aes(x = sum / 1e6)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(x = "Library sizes (millions)", y = "Number of cells") +
  theme_minimal()

# Save in the same folder as analysis outputs
ggsave(
  filename = file.path(output_dir, paste0(sample_id, "_sce_hist_library_sizes.png")),
  plot = p_counts,
  width = 7, height = 6, units = "in", dpi = 300
)

# 2️⃣ Number of detected features per cell
p_genes <- ggplot(as.data.frame(colData(sce)), aes(x = detected)) +
  geom_histogram(bins = 30, fill = "grey60", color = "white") +
  labs(x = "Number of expressed genes", y = "Number of cells") +
  theme_minimal()

# Save in the same folder as analysis outputs
ggsave(
  filename = file.path(output_dir, paste0(sample_id, "_sce_hist_expressed_genes.png")),
  plot = p_genes,
  width = 7, height = 6, units = "in", dpi = 300
)

cat("✅ Saved histograms to:", output_dir, "\n")

library(Matrix)
library(ggplot2)

# 1) Per-gene average (from counts)
cnts <- assay(sce, "counts")                    # genes x cells (sparse ok)
avg   <- Matrix::rowMeans(cnts)
log10_avg <- log10(pmax(avg, 1e-8))             # avoid -Inf

# 2) Choose a threshold (pick one)
# thr <- -0.7                                     # A) set manually
# thr <- quantile(log10_avg, 0.10)              # B) 10th percentile
thr <- median(log10_avg) - 3*mad(log10_avg)   # C) robust MAD rule

# 3) Plot
df <- data.frame(log10_avg = log10_avg)
p <- ggplot(df, aes(x = log10_avg)) +
  geom_histogram(bins = 80, fill = "grey70", color = "grey20") +
  geom_vline(xintercept = thr, linetype = "dashed", color = "dodgerblue", linewidth = 1) +
  labs(x = expression(Log[10]~average~count), y = "Frequency") +
  theme_classic()

print(p)

# 4) Save
# Save histogram of log10 average counts in the same folder
ggsave(
  filename = file.path(output_dir, paste0(sample_id, "_log10_average_counts_hist.png")),
  plot = p,
  width = 7, height = 6, units = "in", dpi = 300
)

cat("✅ Saved log10 average counts histogram to:", output_dir, "\n")


# code GPT5

library(SingleCellExperiment)
library(scran)
library(scater)
library(ggplot2)

# Assume `sce` is your SingleCellExperiment
# Normalize if needed
sce <- logNormCounts(sce)

# Model the variance for each gene
dec <- modelGeneVar(sce)

# === Save the plot in the same folder as other analysis outputs ===
png(
  file.path(output_dir, paste0(sample_id, "_sce_variance_mean_relationship.png")),
  width = 7, height = 6, units = "in", res = 300
)

plot(
  dec$mean, dec$total,
  pch = 16, cex = 0.5,
  xlab = "Mean log-expression",
  ylab = "Variance of log-expression"
)
curve(metadata(dec)$trend(x), col = "red", add = TRUE, lwd = 2)

dev.off()

cat("✅ Saved variance–mean relationship plot to:",
    file.path(output_dir, paste0(sample_id, "_sce_variance_mean_relationship.png")), "\n")

plot(
  dec$mean, dec$total,
  pch = 16, cex = 0.5,
  xlab = "Mean log-expression",
  ylab = "Variance of log-expression"
)

# https://pmc.ncbi.nlm.nih.gov/articles/PMC5112579/pdf/f1000research-5-10712.pdf

# Define filename
seurat_file <- file.path(output_dir, paste0(sample_id, "_soupX.gs.rds"))

# Save Seurat object
saveRDS(obj, seurat_file)

cat("✅ Seurat object saved successfully!\n")
cat("File path:", seurat_file, "\n")

# ==============================================================
# ✅ Read back the saved Seurat object and verify gene naming
# ==============================================================
cat("\n--- Verifying saved Seurat object ---\n")

# Read the saved Seurat object
obj_check <- readRDS(seurat_file)

cat("Loaded Seurat object from:", seurat_file, "\n")
cat("Object class:", class(obj_check), "\n")
cat("Genes:", nrow(obj_check), "| Cells:", ncol(obj_check), "\n\n")

# ==============================================================
# Check for ENSEMBL-style gene names
# ==============================================================
ensembl_pattern <- "^ENSG"
ensembl_genes <- grep(ensembl_pattern, rownames(obj_check), value = TRUE)

cat("Number of ENSEMBL-style gene names:", length(ensembl_genes), "\n")
if (length(ensembl_genes) > 0) {
  cat("Example ENSEMBL gene names:\n")
  print(head(ensembl_genes, 10))
} else {
  cat("✅ No ENSEMBL-style gene names found — all converted to symbols.\n")
}

# ==============================================================
# Check for mitochondrial genes (MT-)
# ==============================================================
mt_genes_present <- grep("^MT-", rownames(obj_check), value = TRUE)

cat("\nNumber of mitochondrial genes (MT-):", length(mt_genes_present), "\n")
if (length(mt_genes_present) > 0) {
  cat("✅ Mitochondrial genes found:\n")
  print(mt_genes_present)
} else {
  cat("⚠️ No MT- genes found — check mapping or capitalization.\n")
}

cat("\n--- Verification complete ---\n")



