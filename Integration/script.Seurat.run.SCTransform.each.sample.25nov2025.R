# Load required libraries
library(patchwork)
library(Seurat)
library(ggplot2)

# Define all sample files
sample_files <- c(
  "AD_AG_2658_seurat_obj.filtered.soupX.rds",
  "AD_AG_2822_seurat_obj.filtered.soupX.rds",
  "AD_AG_3184_seurat_obj.filtered.soupX.rds",
  "AD_AG_3188_seurat_obj.filtered.soupX.rds",
  "AD_CTL_2273_S2_seurat_obj.filtered.soupX.rds",
  "AD_CTL_2354_S3_seurat_obj.filtered.soupX.rds",
  "AD_CTL_2812_S1_seurat_obj.filtered.soupX.rds",
  "Dy_AD_2328_S5_seurat_obj.filtered.soupX.rds",
  "Dy_AD_2368_S6_seurat_obj.filtered.soupX.rds",
  "Dy_AD_2660_S4_seurat_obj.filtered.soupX.rds",
  "Dy_AG_2659_seurat_obj.filtered.soupX.rds",
  "Dy_AG_2922_seurat_obj.filtered.soupX.rds",
  "Dy_AG_3090_seurat_obj.filtered.soupX.rds",
  "Dy_AG_3296_seurat_obj.filtered.soupX.rds",
  "Healthy_AG_2394_seurat_obj.filtered.soupX.rds",
  "Healthy_AG_2594_seurat_obj.filtered.soupX.rds",
  "Healthy_AG_3041_seurat_obj.filtered.soupX.rds",
  "Healthy_AG_3267_seurat_obj.filtered.soupX.rds"
)

# Define AD-related genes
ad_genes <- c("APP", "PSEN1", "PSEN2", "APOE", "TREM2", "BIN1")

# Define dyslexia-related genes
dyslexia_genes <- c("DCDC2", "KIAA0319", "DYX1C1", "ROBO1")

# Directory for QC plots
qc_dir <- "QC_plots"
if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)

# Loop through each sample
for (rds_file in sample_files) {
  
  # Check if file exists
  if (!file.exists(rds_file)) {
    print(paste("WARNING: File not found:", rds_file))
    next
  }
  
  print("========================================")
  print(paste("PROCESSING:", rds_file))
  print("========================================")
  
  # Extract sample name (everything before "_seurat_obj")
  sample_name <- sub("_seurat_obj.*", "", basename(rds_file))
  print(paste("Sample name:", sample_name))
  
  # Read the object
  tryCatch({
    obj <- readRDS(rds_file)
    
    # Ensure the object has sample_id metadata
    if (!"sample_id" %in% colnames(obj@meta.data)) {
      obj$sample_id <- sample_name
      print(paste("Added sample_id metadata:", sample_name))
    } else {
      # Update sample_id if it exists but is different
      if (unique(obj$sample_id)[1] != sample_name) {
        obj$sample_id <- sample_name
        print(paste("Updated sample_id metadata to:", sample_name))
      } else {
        print(paste("sample_id metadata already exists:", unique(obj$sample_id)[1]))
      }
    }
    
    # Compute percent.mt if not present
    if (!"percent.mt" %in% colnames(obj@meta.data)) {
      print("Computing percent.mt...")
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
      print("percent.mt computed and added to metadata")
    } else {
      print("percent.mt already exists in metadata")
    }
    
    # ============================================================
    # QC violin plot (RNA assay) BEFORE SCTransform
    # ============================================================
    print("Creating QC violin plot...")
    DefaultAssay(obj) <- "RNA"
    qc_plot <- VlnPlot(
      obj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      ncol = 3,
      pt.size = 0
    ) + ggtitle(paste(sample_name, "- QC metrics"))
    
    ggsave(
      file.path(qc_dir, paste0(sample_name, "_QC_metrics.png")),
      qc_plot,
      width = 12,
      height = 4,
      dpi = 300
    )
    rm(qc_plot)
    gc(verbose = FALSE)
    print("QC violin plot saved.")
    
    # ============================================================
    # SCTransform workflow (SCT assay created here)
    # ============================================================
    print("Running SCTransform...")
    obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
    
    # ============================================================
    # RNA-based PCA / Neighbors / Clusters / UMAP
    # ============================================================
    print("----------------------------------------")
    print("RNA-based PCA / Neighbors / Clusters / UMAP")
    print("----------------------------------------")
    
    if (!"RNA" %in% Assays(obj)) {
      stop("RNA assay not found in object: ", rds_file)
    }
    
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    
    obj <- RunPCA(
      obj,
      reduction.name = "pca_RNA",
      npcs = 30,
      verbose = FALSE
    )
    
    obj <- FindNeighbors(
      obj,
      reduction  = "pca_RNA",
      dims       = 1:30,
      graph.name = "RNA_snn"
    )
    
    obj <- FindClusters(
      obj,
      graph.name   = "RNA_snn",
      resolution   = 0.5,
      algorithm    = 1,
      cluster.name = "RNA_snn_res.0.5"
    )
    
    obj <- RunUMAP(
      obj,
      reduction      = "pca_RNA",
      dims           = 1:30,
      reduction.name = "umap_RNA",
      verbose        = FALSE
    )
    
    # ============================================================
    # SCT-based PCA / Neighbors / Clusters / UMAP
    # ============================================================
    print("----------------------------------------")
    print("SCT-based PCA / Neighbors / Clusters / UMAP")
    print("----------------------------------------")
    
    if (!"SCT" %in% Assays(obj)) {
      stop("SCT assay not found in object: ", rds_file)
    }
    
    DefaultAssay(obj) <- "SCT"
    
    obj <- RunPCA(
      obj,
      reduction.name = "pca_SCT",
      npcs = 30,
      verbose = FALSE
    )
    
    obj <- FindNeighbors(
      obj,
      reduction  = "pca_SCT",
      dims       = 1:30,
      graph.name = "SCT_snn"
    )
    
    obj <- FindClusters(
      obj,
      graph.name   = "SCT_snn",
      resolution   = 0.5,
      algorithm    = 1,
      cluster.name = "SCT_snn_res.0.5"
    )
    
    obj <- RunUMAP(
      obj,
      reduction      = "pca_SCT",
      dims           = 1:30,
      reduction.name = "umap_SCT",
      verbose        = FALSE
    )
    
    # Set SCT clusters as default identities for downstream plots
    Idents(obj) <- obj$SCT_snn_res.0.5
    
    # ============================================================================
    # UMAP WITH CLUSTER LABELS (RNA vs SCT side-by-side)
    # ============================================================================
    print("----------------------------------------")
    print("CREATING UMAP WITH CLUSTER LABELS (RNA vs SCT)")
    print("----------------------------------------")
    tryCatch({
      print("Creating UMAP DimPlots with cluster labels...")
      
      p_umap_rna <- DimPlot(
        obj,
        reduction = "umap_RNA",
        group.by  = "RNA_snn_res.0.5",
        label     = TRUE,
        label.size = 4,
        pt.size   = 0.5
      ) + ggtitle(paste(sample_name, "- RNA UMAP (RNA_snn_res.0.5)"))
      
      p_umap_sct <- DimPlot(
        obj,
        reduction = "umap_SCT",
        group.by  = "SCT_snn_res.0.5",
        label     = TRUE,
        label.size = 4,
        pt.size   = 0.5
      ) + ggtitle(paste(sample_name, "- SCT UMAP (SCT_snn_res.0.5)"))
      
      p_umap_combined <- p_umap_rna | p_umap_sct
      
      ggsave(
        paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform_UMAP_clusters_RNA_vs_SCT.png"),
        p_umap_combined,
        width  = 14,
        height = 8,
        dpi    = 300
      )
      
      rm(p_umap_rna, p_umap_sct, p_umap_combined)
      gc(verbose = FALSE)
      print("UMAP with cluster labels (RNA vs SCT) saved successfully")
    }, error = function(e) {
      message("Error creating UMAP DimPlot: ", e$message)
    })
    
    # ============================================================================
    # SECTION 1: ALZHEIMER'S DISEASE GENES
    # ============================================================================
    print("----------------------------------------")
    print("ANALYZING ALZHEIMER'S DISEASE GENES")
    print("----------------------------------------")
    
    # Check which genes are present
    genes_present_ad <- ad_genes[ad_genes %in% rownames(obj)]
    genes_missing_ad <- ad_genes[!ad_genes %in% rownames(obj)]
    
    if (length(genes_missing_ad) > 0) {
      print(paste("Missing AD genes:", paste(genes_missing_ad, collapse = ", ")))
    }
    
    if (length(genes_present_ad) > 0) {
      print(paste("Found AD genes:", paste(genes_present_ad, collapse = ", ")))
      
      # UMAP Feature Plots - default colors - SCT assay (use umap_SCT)
      tryCatch({
        print("Creating AD FeaturePlot (default colors) - SCT assay...")
        DefaultAssay(obj) <- "SCT"
        p1 <- FeaturePlot(
          obj,
          features   = genes_present_ad,
          ncol       = 3,
          pt.size    = 0.5,
          reduction  = "umap_SCT"
        )
        ggsave(
          paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform_AD_genes_FeaturePlot_default_SCT.png"),
          p1, width = 12, height = 8, dpi = 300
        )
      }, error = function(e) {
        message("Error in AD FeaturePlot (default): ", e$message)
      })
      
      # UMAP plots with red gradient - SCT assay (umap_SCT)
      tryCatch({
        print("Creating AD FeaturePlot (red gradient) - SCT assay...")
        DefaultAssay(obj) <- "SCT"
        p2 <- FeaturePlot(
          obj,
          features   = genes_present_ad,
          ncol       = 3,
          pt.size    = 0.5,
          cols       = c("lightgrey", "red"),
          reduction  = "umap_SCT"
        )
        ggsave(
          paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform_AD_genes_FeaturePlot_red_SCT.png"),
          p2, width = 12, height = 8, dpi = 300
        )
      }, error = function(e) {
        message("Error in AD FeaturePlot (red gradient): ", e$message)
      })
      
      # Individual gene comparisons: SCT vs RNA (explicit reductions)
      tryCatch({
        print("Creating individual AD gene comparisons (SCT vs RNA)...")
        for (gene in genes_present_ad) {
          # SCT assay
          DefaultAssay(obj) <- "SCT"
          p_sct <- FeaturePlot(
            obj,
            features  = gene,
            pt.size   = 0.5,
            cols      = c("lightgrey", "red"),
            reduction = "umap_SCT"
          ) + ggtitle(paste(gene, "- SCT assay"))
          
          # RNA assay
          DefaultAssay(obj) <- "RNA"
          p_rna <- FeaturePlot(
            obj,
            features  = gene,
            pt.size   = 0.5,
            cols      = c("lightgrey", "red"),
            reduction = "umap_RNA"
          ) + ggtitle(paste(gene, "- RNA assay"))
          
          # Combine side by side
          p_combined <- p_sct | p_rna
          
          # Save
          ggsave(
            paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform_AD_", gene, "_SCT_vs_RNA.png"),
            p_combined, width = 14, height = 6, dpi = 300
          )
        }
      }, error = function(e) {
        message("Error in individual AD gene comparisons: ", e$message)
      })
      
      # Set back to SCT and SCT-derived clusters for DotPlot/VlnPlot
      DefaultAssay(obj) <- "SCT"
      Idents(obj) <- obj$SCT_snn_res.0.5
      
      # Dot Plot by cluster (unscaled)
      tryCatch({
        print("Creating AD DotPlot (unscaled)...")
        p3 <- DotPlot(
          obj,
          features = genes_present_ad
        ) +
          RotatedAxis() +
          ggtitle(paste(sample_name, "- AD-related genes across SCT clusters"))
        ggsave(
          paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform_AD_genes_DotPlot_unscaled.png"),
          p3, width = 8, height = 6, dpi = 300
        )
      }, error = function(e) {
        message("Error in AD DotPlot (unscaled): ", e$message)
      })
      
      # Dot Plot with scaling
      tryCatch({
        print("Creating AD DotPlot (scaled)...")
        p4 <- DotPlot(
          obj,
          features = genes_present_ad,
          scale    = TRUE
        ) +
          RotatedAxis() +
          ggtitle(paste(sample_name, "- AD-related genes (scaled expression, SCT clusters)"))
        ggsave(
          paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform_AD_genes_DotPlot_scaled.png"),
          p4, width = 8, height = 6, dpi = 300
        )
      }, error = function(e) {
        message("Error in AD DotPlot (scaled): ", e$message)
      })
      
      # Violin plots
      tryCatch({
        print("Creating AD VlnPlot...")
        p5 <- VlnPlot(
          obj,
          features = genes_present_ad,
          ncol     = 3,
          pt.size  = 0
        )
        ggsave(
          paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform_AD_genes_VlnPlot.png"),
          p5, width = 12, height = 8, dpi = 300
        )
      }, error = function(e) {
        message("Error in AD VlnPlot: ", e$message)
      })
      
      print("AD gene plots completed!")
      
    } else {
      print("No AD-related genes found in the dataset!")
    }
    
    # ============================================================================
    # SECTION 2: DYSLEXIA GENES
    # ============================================================================
    print("----------------------------------------")
    print("ANALYZING DYSLEXIA GENES")
    print("----------------------------------------")
    
    # Check which genes are present
    genes_present_dys <- dyslexia_genes[dyslexia_genes %in% rownames(obj)]
    genes_missing_dys <- dyslexia_genes[!dyslexia_genes %in% rownames(obj)]
    
    if (length(genes_missing_dys) > 0) {
      print(paste("Missing dyslexia genes:", paste(genes_missing_dys, collapse = ", ")))
    }
    
    if (length(genes_present_dys) > 0) {
      print(paste("Found dyslexia genes:", paste(genes_present_dys, collapse = ", ")))
      
      # UMAP Feature Plots - default colors - SCT assay (umap_SCT)
      tryCatch({
        print("Creating dyslexia FeaturePlot (default colors) - SCT assay...")
        DefaultAssay(obj) <- "SCT"
        p1 <- FeaturePlot(
          obj,
          features   = genes_present_dys,
          ncol       = 2,
          pt.size    = 0.5,
          reduction  = "umap_SCT"
        )
        ggsave(
          paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform_dyslexia_genes_FeaturePlot_default_SCT.png"),
          p1, width = 10, height = 8, dpi = 300
        )
      }, error = function(e) {
        message("Error in dyslexia FeaturePlot (default): ", e$message)
      })
      
      # UMAP plots with red gradient - SCT assay (umap_SCT)
      tryCatch({
        print("Creating dyslexia FeaturePlot (red gradient) - SCT assay...")
        DefaultAssay(obj) <- "SCT"
        p2 <- FeaturePlot(
          obj,
          features   = genes_present_dys,
          ncol       = 2,
          pt.size    = 0.5,
          cols       = c("lightgrey", "red"),
          reduction  = "umap_SCT"
        )
        ggsave(
          paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform_dyslexia_genes_FeaturePlot_red_SCT.png"),
          p2, width = 10, height = 8, dpi = 300
        )
      }, error = function(e) {
        message("Error in dyslexia FeaturePlot (red gradient): ", e$message)
      })
      
      # Individual gene comparisons: SCT vs RNA (explicit reductions)
      tryCatch({
        print("Creating individual dyslexia gene comparisons (SCT vs RNA)...")
        for (gene in genes_present_dys) {
          # SCT assay
          DefaultAssay(obj) <- "SCT"
          p_sct <- FeaturePlot(
            obj,
            features  = gene,
            pt.size   = 0.5,
            cols      = c("lightgrey", "red"),
            reduction = "umap_SCT"
          ) + ggtitle(paste(gene, "- SCT assay"))
          
          # RNA assay
          DefaultAssay(obj) <- "RNA"
          p_rna <- FeaturePlot(
            obj,
            features  = gene,
            pt.size   = 0.5,
            cols      = c("lightgrey", "red"),
            reduction = "umap_RNA"
          ) + ggtitle(paste(gene, "- RNA assay"))
          
          # Combine side by side
          p_combined <- p_sct | p_rna
          
          # Save
          ggsave(
            paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform_dyslexia_", gene, "_SCT_vs_RNA.png"),
            p_combined, width = 14, height = 6, dpi = 300
          )
        }
      }, error = function(e) {
        message("Error in individual dyslexia gene comparisons: ", e$message)
      })
      
      # Set back to SCT clusters
      DefaultAssay(obj) <- "SCT"
      Idents(obj) <- obj$SCT_snn_res.0.5
      
      # Dot Plot by cluster (unscaled)
      tryCatch({
        print("Creating dyslexia DotPlot (unscaled)...")
        p3 <- DotPlot(
          obj,
          features = genes_present_dys
        ) +
          RotatedAxis() +
          ggtitle(paste(sample_name, "- Dyslexia-related genes across SCT clusters"))
        ggsave(
          paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform_dyslexia_genes_DotPlot_unscaled.png"),
          p3, width = 8, height = 6, dpi = 300
        )
      }, error = function(e) {
        message("Error in dyslexia DotPlot (unscaled): ", e$message)
      })
      
      # Dot Plot with scaling
      tryCatch({
        print("Creating dyslexia DotPlot (scaled)...")
        p4 <- DotPlot(
          obj,
          features = genes_present_dys,
          scale    = TRUE
        ) +
          RotatedAxis() +
          ggtitle(paste(sample_name, "- Dyslexia-related genes (scaled expression, SCT clusters)"))
        ggsave(
          paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform_dyslexia_genes_DotPlot_scaled.png"),
          p4, width = 8, height = 6, dpi = 300
        )
      }, error = function(e) {
        message("Error in dyslexia DotPlot (scaled): ", e$message)
      })
      
      # Violin plots
      tryCatch({
        print("Creating dyslexia VlnPlot...")
        p5 <- VlnPlot(
          obj,
          features = genes_present_dys,
          ncol     = 2,
          pt.size  = 0
        )
        ggsave(
          paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform_dyslexia_genes_VlnPlot.png"),
          p5, width = 10, height = 8, dpi = 300
        )
      }, error = function(e) {
        message("Error in dyslexia VlnPlot: ", e$message)
      })
      
      print("Dyslexia gene plots completed!")
      
    } else {
      print("No dyslexia-related genes found in the dataset!")
    }
    
    print(paste("COMPLETED:", sample_name))
    print("========================================")
    
    # Save the processed object (now includes SCT + RNA & SCT UMAPs/clusters)
    saveRDS(obj, paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform.rds"))
    print(paste("Saved:", paste0(sample_name, "_seurat_obj.filtered.soupX.sctransform.rds")))
    
    # Clean up memory
    rm(obj)
    gc()
    
  }, error = function(e) {
    message(paste("ERROR processing", rds_file, ":", e$message))
  })
}

print("========================================")
print("ALL SAMPLES PROCESSED!")
print("========================================")
