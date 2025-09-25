library(scRNAseq)

# Brain resources

## 1) Search: only "brain", human datasets
# brain <- searchDatasets(defineTextQuery("brain", partial = TRUE))

## 2) Tidy summar
# unlist_first <- function(x) paste0(unlist(x), collapse = ", ")
# summ <- data.frame(
#  name    = brain$name,
#  version = brain$version,
#  title   = brain$title,
#  taxonomy_id = sapply(brain$taxonomy_id, unlist_first),  # 9606 = human; 10090 = mouse
#  genome      = sapply(brain$genome,      unlist_first),  # GRCh* human; GRCm* mouse
#  n_genes = brain$rows,
#  n_cells = brain$columns,
#  assays  = sapply(brain$assays,  unlist_first),
#  coldata = sapply(brain$column_annotations, unlist_first),
#  stringsAsFactors = FALSE
# )
# summ <- summ[order(summ$name, summ$version), ]
# summ
# If you want only human brain:
# subset(summ, grepl("9606", taxonomy_id) | grepl("^GRCh", genome))

## 3) Pretty “details per study”
# for (i in seq_len(nrow(brain))) {
#  src <- tryCatch(as.data.frame(brain$sources[[i]]), error = function(e) NULL)
#  src_str <- if (!is.null(src) && nrow(src)) {
#    paste(sprintf("%s:%s", src$provider, src$id), collapse = " | ")
#  } else NA_character_

#  cat(sprintf("\n[%d] %s (version %s)\n", i, brain$name[i], brain$version[i]))
#  cat(sprintf("Title:   %s\n", brain$title[i]))
#  cat(sprintf("Genome:  %s   |   Taxonomy: %s\n",
#              paste(unlist(brain$genome[[i]]), collapse=", "),
#              paste(unlist(brain$taxonomy_id[[i]]), collapse=", ")))
#  cat(sprintf("Cells:   %s   |   Genes: %s\n", brain$columns[i], brain$rows[i]))
#  cat(sprintf("Assays:  %s\n",
#              paste(unlist(brain$assays[[i]]), collapse=", ")))
#  cat(sprintf("colData: %s\n",
#              paste(unlist(brain$column_annotations[[i]]), collapse=", ")))
#  cat(sprintf("Sources: %s\n", src_str))
# }

## 4) (Optional)Save the table
# write.csv(summ, "scRNAseq_brain_datasets.csv", row.names = FALSE)

# library(scRNAseq)
# name <- "zhong-prefrontal-2018"
# ver  <- "2023-12-22"

# File-backed by default; set realize.assays=TRUE to load into RAM (dgCMatrix)
# sce <- fetchDataset(name, ver, realize.assays = TRUE)
# sce

# assayNames(sce)              # counts
# dim(sce)                     # 24153 genes x 2394 cells (per your summary)
# names(colData(sce))          # developmental_stage, gender, sample, cell_types, week, ...

# Peek at labels
# head(unique(colData(sce)$cell_types))
# table(colData(sce)$developmental_stage, useNA = "ifany")

# lab <- colData(sce)[["cell_types"]]
# unique(lab)
# [1] "Neurons"           "GABAergic neurons" "Microglia"        
# [4] "Stem cells"        NA                  "Astrocytes"       
# [7] "OPC"       

# SeuratData humancortexref.SeuratData	

# library(SeuratData)
# AvailableData() 
# humancortexref.SeuratData	humancortexref	1.0.0	Azimuth Reference: humancortex	human	motor cortex	76533	cells

# CellXGene 

# Seattle AZ Alzheimer Disease Atlas (SEA-AD) : 
# https://cellxgene.cziscience.com/collections/1ca90a2d-2943-483d-b678-b809bf464c30
# https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443

# Human Brain Atlas 
# https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443

library("SingleR")
library(SingleCellExperiment)

# sce <- as.SingleCellExperiment(obj)
# ref <- celldex::HumanPrimaryCellAtlasData()

# pred_clust <- SingleR(test = sce, 
#                        ref = ref, 
#                        labels = ref$label.main,
#                        clusters = sce$seurat_clusters)

base_dir <- "/mnt/nfs/CX000008_DS1/projects/btanasa/brain_refs_scRNAseq"

# Human Motor Cortex : 
# https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x
# https://celltypes.brain-map.org/rnaseq/human_m1_10x?selectedVisualization=Heatmap&colorByFeature=Cell+Type&colorByFeatureValue=GAD1
# https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad
# https://cellxgene.cziscience.com/collections/1ca90a2d-2943-483d-b678-b809bf464c30

# Mouse annotations
zeisel   <- readRDS(file.path(base_dir, "ZeiselBrainData.rds"))
tasic    <- readRDS(file.path(base_dir, "TasicBrainData.rds"))
romanov  <- readRDS(file.path(base_dir, "RomanovBrainData.rds"))

# sanity checks
stopifnot(inherits(zeisel, "SingleCellExperiment"))
stopifnot(inherits(tasic, "SingleCellExperiment"))
stopifnot(inherits(romanov, "SingleCellExperiment"))

show_sce <- function(sce, name) {
  cat("\n================ ", name, " ================\n", sep = "")
  cat("assays: ", paste(assayNames(sce), collapse = ", "), "\n", sep = "")

## 1) logcounts
  if ("logcounts" %in% assayNames(sce)) {
    logc <- assay(sce, "logcounts")
    cat("logcounts  dim: ", paste(dim(logc), collapse = " x "),
        "  class: ", class(logc), "\n", sep = "")
    print(as.matrix(logc[1:min(5, nrow(logc)), 1:min(5, ncol(logc))]))
  } else {
    cat("No 'logcounts' assay found.\n")
  }

## 2) raw counts (if available)
  raw_name <- if ("counts" %in% assayNames(sce)) "counts" else NA
  if (is.na(raw_name)) {
    cat("No raw 'counts' assay in this reference.\n")
  } else {
    cnt <- assay(sce, raw_name)
    cat(raw_name, " dim: ", paste(dim(cnt), collapse = " x "),
        "  class: ", class(cnt), "\n", sep = "")
    print(as.matrix(cnt[1:min(5, nrow(cnt)), 1:min(5, ncol(cnt))]))
  }

## 3) colData
  cd <- colData(sce)
  cat("colData columns: ", paste(colnames(cd), collapse = ", "), "\n", sep = "")
  print(head(as.data.frame(cd), 10))
 }

## run one by one
show_sce(zeisel,   "ZeiselBrainData")
show_sce(tasic,    "TasicBrainData")
show_sce(romanov,  "RomanovBrainData")

## run one by one
colnames(colData(zeisel))
show_sce(zeisel,   "ZeiselBrainData")
unique(colData(zeisel)[, 8])

## run one by one
colnames(colData(tasic))
show_sce(tasic,   "TasicBrainData")
unique(colData(tasic)[, 1])

## run one by one
colnames(colData(romanov))
show_sce(romanov,   "RomanovBrainData")
unique(colData(tasic)[, 1])

# Show unique values for categorical columns in colData
df <- as.data.frame(SummarizedExperiment::colData(romanov))
cat_cols <- names(Filter(function(x) is.factor(x) || is.character(x) || is.logical(x), df))

invisible(lapply(cat_cols, function(nm) {
  cat("\n==", nm, "==\n", sep = "")
  print(unique(df[[nm]]))
}))

# Show unique values for categorical columns in colData
df <- as.data.frame(SummarizedExperiment::colData(tasic))
cat_cols <- names(Filter(function(x) is.factor(x) || is.character(x) || is.logical(x), df))

invisible(lapply(cat_cols, function(nm) {
  cat("\n==", nm, "==\n", sep = "")
  print(unique(df[[nm]]))
}))

# Show unique values for categorical columns in colData
df <- as.data.frame(SummarizedExperiment::colData(zeisel))
cat_cols <- names(Filter(function(x) is.factor(x) || is.character(x) || is.logical(x), df))

invisible(lapply(cat_cols, function(nm) {
  cat("\n==", nm, "==\n", sep = "")
  print(unique(df[[nm]]))
}))

# make rds objects from MTG_SEA resources

# /mnt/nfs/CX000008_DS1/projects/btanasa/brain_refs_MTG_SE

# working with a file from : /mnt/nfs/CX000008_DS1/projects/btanasa/brain_refs_MTG_SE

# setwd("/mnt/nfs/CX000008_DS1/projects/btanasa/brain_refs_MTG_SEA")
# dir="./"

# suppressPackageStartupMessages({
#  library(data.table)
#  library(Matrix)
#  library(SingleCellExperiment)
#  library(SummarizedExperiment)
# })

# mfile <- file.path(dir, "matrix.csv")
# meta  <- file.path(dir, "metadata.csv")

# dt <- fread(mfile)   
# head(dt[1:10,1:10])

# md <- fread(meta)
# colnames(md)
# head(md,2)

# unique(md$cluster_label)
# unique(md$class_label)
# unique(md$subclass_label)

# cell_ids <- as.character(dt[["sample_name"]])
# gene_ids <- as.character(names(dt)[-1])

# fast path (dense -> transpose -> sparse)
# counts <- try({
#  M <- as.matrix(dt[, -1, with = FALSE])   # cells x genes
#  mode(M) <- "numeric"
#  M <- t(M)                                # genes x cells
#  colnames(M) <- cell_ids
#  rownames(M) <- gene_ids
#  Matrix(M, sparse = TRUE)
# }, silent = TRUE)

# o <- match(colnames(counts), md$sample_name)

# dim(counts)
# head(o)

# o <- match(colnames(counts), md$sample_name)
# if (anyNA(o))
#  stop("Some cells in matrix.csv are missing from metadata.csv via 'sample_name'.")
# md <- as.data.frame(md[o, , drop = FALSE])

#  Build SCE and add logcounts
# sce <- SingleCellExperiment(
#  assays  = list(counts = counts),
#  colData = S4Vectors::DataFrame(md)
# )

# sce <- logNormCounts(sce)

# sce                             # class, dims, assays, reduced dims, altExps
# class(sce)
# dim(sce)                        # genes x cells
# assayNames(sce)                 # e.g., "counts", "logcounts"
# head(rownames(sce), 5)          # first genes
# head(colnames(sce), 5)          # first cells (if present)

# names(rowData(sce))             # feature metadata fields
# names(colData(sce))             # cell metadata fields
# head(as.data.frame(colData(sce))[ , 1:min(10, ncol(colData(sce))), drop=FALSE])

# reducedDimNames(sce)            # embeddings stored (e.g., PCA/TSNE/UMAP)
# altExpNames(sce)                # alternative experiments (e.g., ERCC, repeats)
# metadata(sce)                   # list of extra metadata (taxonomy, dendrogram, etc.)

# --- Basic counts/logcounts checks ---
# has_counts   <- "counts"    %in% assayNames(sce)
# has_logcounts<- "logcounts" %in% assayNames(sce)
# cat("Has counts:", has_counts, " | Has logcounts:", has_logcounts, "\n")

# saveRDS(sce, file.path(dir, "brain_refs_MTG_SEA.from.matrix.metadata.sce.rds"))



mtgsea <- "/mnt/nfs/CX000008_DS1/projects/btanasa/brain_refs_MTG_SEA/brain_refs_MTG_SEA.from.matrix.metadata.sce.rds"
mtg_sea <- readRDS(mtgsea)
# Inspect its structure
class(mtg_sea)
# str(mtg_sea)

# .$ cluster_label              : chr [1:76533] "Inh L1-2 SST CCNJL" "Exc L5-6 FEZF2 IFNG-AS1" "Exc L3-5 RORB LINC01202" "Exc L2 LINC00507 GLRA3" ...
# .$ class_label                : chr [1:76533] "GABAergic" "Glutamatergic" "Glutamatergic" "Glutamatergic" ...
# .$ subclass_label             : chr [1:76533] "Sst" "L5/6 NP" "L5 IT" "L2/3 IT" ...

# Check available columns in colData
colnames(colData(mtg_sea))

# Get unique values for each label
unique_clusters    <- unique(mtg_sea$cluster_label)
unique_classes     <- unique(mtg_sea$class_label)
unique_subclasses  <- unique(mtg_sea$subclass_label)

# Print counts + first few
cat("cluster_label:", length(unique_clusters), "unique\n")
print(unique_clusters)

cat("class_label:", length(unique_classes), "unique\n")
print(unique_classes)

cat("subclass_label:", length(unique_subclasses), "unique\n")
print(unique_subclasses)

# Inspect the object type
class(mtg_sea)

# Extract rownames (gene IDs)
ids <- sub("\\.\\d+$", "", rownames(mtg_sea))   # strip version suffixes like .12

# Flag Ensembl IDs
is_ensg <- grepl("^ENSG", ids)

# Count totals
n_total  <- length(ids)
n_ensg   <- sum(is_ensg)
n_symbol <- n_total - n_ensg

cat("Object class:", class(mtg_sea), "\n")
cat("Total genes:", n_total, "\n")
cat("Ensembl (ENSG) IDs:", n_ensg, "\n")
cat("Other (likely gene symbols):", n_symbol, "\n")

# Peek at examples
cat("\nExample Ensembl IDs:\n")
print(head(ids[is_ensg]))

cat("\nExample non-Ensembl IDs:\n")
print(head(ids[!is_ensg]))

mtgseaad <- "/mnt/nfs/CX000008_DS1/projects/btanasa/brain_refs_MTG_SEA_AD/Reference_MTG_RNAseq_all-nuclei.2022-06-07.sce.rds"
mtg_sea_ad <- readRDS(mtgseaad)

# Inspect its structure
class(mtg_sea_ad)

#.$ cluster_label            : chr [1:166868] "Pax6_1" NA "L5/6 NP_1" "L5 IT_7" ...
#.$ subclass_label           : chr [1:166868] "Pax6" NA "L5/6 NP" "L5 IT" ...
#.$ class_label              : chr [1:166868] "Neuronal: GABAergic" "NA" "Neuronal: Glutamatergic" "Neuronal: Glutamatergic" ...
#.$ GA_cluster_label         : chr [1:166868] "Pax6_1" NA "L5/6 NP_1" "L5 IT_7" ...
#.$ GA_subclass_label        : chr [1:166868] "Pax6" NA "L5/6 NP" "L5 IT" ...
#.$ GA_neighborhood_label    : chr [1:166868] "lamp5_sncg_vip" NA "l5et_l56np_l6ct_l6b" "it_types" ...
#.$ CA_cluster_label         : chr [1:166868] "Pax6_1" NA "L5/6 NP_3" "L5 IT_1" ...
#.$ CA_subclass_label        : chr [1:166868] "Pax6" NA "L5/6 NP" "L5 IT" ...
#.$ CA_neighborhood_label    : chr [1:166868] "CGE Inh" NA "Deep Exc" "IT types" ..

# Check available columns in colData
colnames(colData(mtg_sea_ad))

# Get unique values for each label
unique_clusters    <- unique(mtg_sea_ad$cluster_label)
unique_subclasses  <- unique(mtg_sea_ad$subclass_label)
unique_classes     <- unique(mtg_sea_ad$class_label)

# Print counts + first few
cat("cluster_label:", length(unique_clusters), "unique\n")
print(unique_clusters)

cat("class_label:", length(unique_classes), "unique\n")
print(unique_classes)

cat("subclass_label:", length(unique_subclasses), "unique\n")
print(unique_subclasses)

# Inspect class
class(mtg_sea_ad)

# Extract feature IDs, strip version suffixes
ids <- sub("\\.\\d+$", "", rownames(mtg_sea_ad))

# Flag Ensembl IDs
is_ensg <- grepl("^ENSG", ids)

# Counts
n_total  <- length(ids)
n_ensg   <- sum(is_ensg)
n_symbol <- n_total - n_ensg

cat("Object class:", class(mtg_sea_ad), "\n")
cat("Total features:", n_total, "\n")
cat("Ensembl (ENSG) IDs:", n_ensg, "\n")
cat("Other (likely gene symbols):", n_symbol, "\n")

# Peek at examples
cat("\nExample Ensembl IDs:\n")
print(head(ids[is_ensg]))

cat("\nExample non-Ensembl IDs:\n")
print(head(ids[!is_ensg]))

# To work with : 

# cluster_label
# subclass_label


