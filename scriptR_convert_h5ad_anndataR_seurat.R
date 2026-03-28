library(anndataR)
library(Seurat)

h5ad_file <- "18_merged_AG_031826.26celltypes.min50cells.18K.proportional.h5ad"

adata <- read_h5ad(h5ad_file)

# Check what layers and reductions exist
print(adata$layers_keys())
print(names(adata$obsm))

# Convert to Seurat
seurat_obj <- adata$as_Seurat(
  layers_mapping = TRUE,
  reduction_mapping = list(
    pca          = c(key = "PC_",      embeddings = "X_pca"),
    harmony      = c(key = "harmony_", embeddings = "X_harmony"),
    umap         = c(key = "UMAP_",    embeddings = "X_umap"),
    umap_harmony = c(key = "UMAPH_",   embeddings = "X_umap_harmony")
  )
)

# Rename normalized layer to data (Seurat convention)
seurat_obj <- SetAssayData(seurat_obj, 
                           layer    = "data", 
                           new.data = LayerData(seurat_obj, layer = "normalized"))

# Verify
print(seurat_obj)
print(Layers(seurat_obj))

# Save
saveRDS(seurat_obj, "18_merged_AG_031826.26celltypes.min50cells.18K.proportional.anndataR.rds")
