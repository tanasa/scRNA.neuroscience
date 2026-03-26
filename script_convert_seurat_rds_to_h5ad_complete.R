library(Seurat)
library(anndataR)

# ── Read RDS ──────────────────────────────────────────────────
cat("Reading RDS...\n")
a <- readRDS("18_merged_AG_031826.rds")
message("✔ Loaded: ", ncol(a), " cells, ", nrow(a), " genes")

# ── Extract all components ────────────────────────────────────
cat("Extracting components...\n")

# 1. Raw counts
raw_counts <- GetAssayData(a, assay = "RNA", layer = "counts")

# 2. Normalized counts
norm_counts <- GetAssayData(a, assay = "RNA", layer = "data")

# 3. Metadata
metadata <- a@meta.data

# 4. Reductions
pca_coords          <- as.matrix(a@reductions$pca@cell.embeddings)
umap_coords         <- as.matrix(a@reductions$umap@cell.embeddings)
harmony_coords      <- as.matrix(a@reductions$harmony@cell.embeddings)
umap_harmony_coords <- as.matrix(a@reductions$umap.harmony@cell.embeddings)

cat("Raw counts dim:        ", dim(raw_counts), "\n")
cat("Normalized counts dim: ", dim(norm_counts), "\n")
cat("PCA dim:               ", dim(pca_coords), "\n")
cat("UMAP dim:              ", dim(umap_coords), "\n")
cat("Harmony dim:           ", dim(harmony_coords), "\n")
cat("UMAP harmony dim:      ", dim(umap_harmony_coords), "\n")

# ── Build AnnData object ──────────────────────────────────────
cat("Building AnnData...\n")

adata <- AnnData(
  X      = t(norm_counts),
  obs    = metadata,
  var    = data.frame(gene = rownames(raw_counts),
                      row.names = rownames(raw_counts)),
  layers = list(
    counts     = t(raw_counts),
    normalized = t(norm_counts)
  ),
  obsm   = list(
    X_pca          = pca_coords,
    X_umap         = umap_coords,
    X_harmony      = harmony_coords,
    X_umap_harmony = umap_harmony_coords
  )
)

print(adata)

# ── Save ──────────────────────────────────────────────────────
cat("Writing h5ad...\n")
adata$write_h5ad("18_merged_AG_031826_full.h5ad")
message("✔ Saved: 18_merged_AG_031826_full.h5ad")
