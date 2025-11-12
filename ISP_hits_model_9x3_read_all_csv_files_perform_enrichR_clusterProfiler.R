# Read all CSV files in the current directory
# This script reads all cell_type_ISP_* CSV files into individual variables
# Each variable is named based on the file name, removing:
#   - "cell_type_ISP_" prefix
#   - "_GPU_ONLY_genes" suffix
#   - Replacing hyphens with underscores

library(readr)
library(dplyr)
library(ggplot2)

# Load clusterProfiler packages (install if needed)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  cat("Note: clusterProfiler package not found. Install with: BiocManager::install('clusterProfiler')\n")
  cat("clusterProfiler analyses will be skipped.\n")
  clusterprofiler_available <- FALSE
} else {
  library(clusterProfiler)
  clusterprofiler_available <- TRUE
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  cat("Note: org.Hs.eg.db package not found. Install with: BiocManager::install('org.Hs.eg.db')\n")
  orgdb_available <- FALSE
} else {
  library(org.Hs.eg.db)
  orgdb_available <- TRUE
}

if (!requireNamespace("ReactomePA", quietly = TRUE)) {
  cat("Note: ReactomePA package not found. Install with: BiocManager::install('ReactomePA')\n")
  reactomepa_available <- FALSE
} else {
  library(ReactomePA)
  reactomepa_available <- TRUE
}

if (!requireNamespace("enrichplot", quietly = TRUE)) {
  cat("Note: enrichplot package not found. Install with: BiocManager::install('enrichplot')\n")
  enrichplot_available <- FALSE
} else {
  library(enrichplot)
  enrichplot_available <- TRUE
}

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  cat("Note: msigdbr package not found. Install with: install.packages('msigdbr')\n")
  msigdbr_available <- FALSE
} else {
  library(msigdbr)
  msigdbr_available <- TRUE
}

# Load Excel writing package
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  cat("Note: openxlsx package not found. Install with: install.packages('openxlsx')\n")
  cat("Excel export will be skipped.\n")
  openxlsx_available <- FALSE
} else {
  library(openxlsx)
  openxlsx_available <- TRUE
}

# Load HTML report packages
if (!requireNamespace("DT", quietly = TRUE)) {
  cat("Note: DT package not found. Install with: install.packages('DT')\n")
  dt_available <- FALSE
} else {
  library(DT)
  dt_available <- TRUE
}

# Load network analysis packages
if (!requireNamespace("igraph", quietly = TRUE)) {
  cat("Note: igraph package not found. Install with: install.packages('igraph')\n")
  cat("Network visualizations will be skipped.\n")
  igraph_available <- FALSE
} else {
  library(igraph)
  igraph_available <- TRUE
}

# Load Venn diagram packages (install if needed)
if (!requireNamespace("ggvenn", quietly = TRUE)) {
  cat("Note: ggvenn package not found. Install with: install.packages('ggvenn')\n")
  cat("Venn diagrams will be skipped.\n")
  ggvenn_available <- FALSE
} else {
  library(ggvenn)
  ggvenn_available <- TRUE
}

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  cat("Note: UpSetR package not found. Install with: install.packages('UpSetR')\n")
  cat("UpSet plots will be skipped.\n")
  upsetr_available <- FALSE
} else {
  library(UpSetR)
  upsetr_available <- TRUE
}

# ===== SETUP =====
# Use current working directory
# Make sure to set your working directory to the folder containing the CSV files
# using: setwd("/path/to/directory") or run this script from that directory
script_dir <- getwd()

cat("Working directory:", script_dir, "\n")

# ===== FUNCTION TO CREATE VENN DIAGRAMS =====
# Function to create Venn diagrams from significant genes in multiple variables
create_venn_diagram <- function(var_names, title = "Venn Diagram", output_file = NULL, use_sig_only = TRUE) {
  # Check if ggvenn is available
  if (!ggvenn_available) {
    cat("Warning: ggvenn package not available. Skipping Venn diagram.\n")
    return(NULL)
  }
  
  # Extract significant genes from each variable
  gene_lists <- list()
  valid_vars <- c()
  
  for (var_name in var_names) {
    if (exists(var_name, envir = .GlobalEnv)) {
      df <- get(var_name, envir = .GlobalEnv)
      
      if (use_sig_only && "Sig" %in% colnames(df)) {
        df <- df %>% filter(Sig == 1)
      }
      
      if ("Gene_name" %in% colnames(df)) {
        genes <- unique(df$Gene_name[!is.na(df$Gene_name) & df$Gene_name != ""])
        gene_lists[[var_name]] <- genes
        valid_vars <- c(valid_vars, var_name)
      }
    }
  }
  
  if (length(gene_lists) < 2) {
    cat("Warning: Need at least 2 valid variables to create Venn diagram.\n")
    return(NULL)
  }
  
  # For 2-3 sets, use ggvenn
  if (length(gene_lists) <= 3) {
    tryCatch({
      p <- ggvenn(gene_lists, 
                  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")[seq_along(gene_lists)],
                  stroke_size = 0.5,
                  set_name_size = 4,
                  text_size = 3,
                  show_percentage = TRUE) +
        labs(title = title) +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      
      # Save plot if output file is specified
      if (!is.null(output_file)) {
        ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
        cat("Venn diagram saved to:", output_file, "\n")
      }
      
      return(p)
    }, error = function(e) {
      cat("Error creating Venn diagram:", conditionMessage(e), "\n")
      return(NULL)
    })
  } else {
    # For more than 3 sets, suggest UpSet plot
    cat("Note: More than 3 sets detected. Consider using UpSet plot instead.\n")
    if (upsetr_available) {
      tryCatch({
        # Convert to binary matrix for UpSetR
        all_genes <- unique(unlist(gene_lists))
        binary_matrix <- data.frame(
          Gene = all_genes,
          stringsAsFactors = FALSE
        )
        
        for (var_name in names(gene_lists)) {
          binary_matrix[[var_name]] <- as.integer(all_genes %in% gene_lists[[var_name]])
        }
        
        if (!is.null(output_file)) {
          png_file <- gsub("\\.png$", "_upset.png", output_file)
          png(png_file, width = 12, height = 8, units = "in", res = 300)
          upset(binary_matrix, sets = names(gene_lists), 
                sets.bar.color = "#56B4E9",
                order.by = "freq",
                text.scale = c(1.3, 1.3, 1, 1, 1.5, 1))
          title(main = title, cex.main = 1.5)
          dev.off()
          cat("UpSet plot saved to:", png_file, "\n")
        }
        
        return(binary_matrix)
      }, error = function(e) {
        cat("Error creating UpSet plot:", conditionMessage(e), "\n")
        return(NULL)
      })
    }
  }
}

# ===== FUNCTION TO RUN enrichR ANALYSIS =====
# Function to run enrichR on a list of genes
run_enrichr_analysis <- function(gene_list, output_name, results_base_dir = NULL, adj_p_cutoff = 0.2, top_n = 50) {
  # Check if enrichR is available
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    cat("Warning: enrichR package not available. Install with: install.packages('enrichR')\n")
    return(NULL)
  }
  
  library(enrichR)
  
  # Clean gene list
  gene_list <- as.character(gene_list)
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  gene_list <- unique(gene_list)
  
  if (length(gene_list) == 0) {
    cat("Warning: Empty gene list provided for", output_name, "\n")
    return(NULL)
  }
  
  cat("\nRunning enrichR for:", output_name, "(", length(gene_list), "genes)\n")
  
  # Set results directory
  if (is.null(results_base_dir)) {
    results_base_dir <- script_dir
  }
  results_dir <- file.path(results_base_dir, paste0(output_name, "_enrichR_results"))
  dir.create(results_dir, showWarnings = FALSE)
  dir.create(file.path(results_dir, "plots"), showWarnings = FALSE)
  dir.create(file.path(results_dir, "tables"), showWarnings = FALSE)
  
  # Define enrichR databases
  enrichr_databases <- c(
    # Updated Gene Ontology (2025 versions)
    "GO_Biological_Process_2025",
    "GO_Molecular_Function_2025", 
    "GO_Cellular_Component_2025",
    
    # Core pathway databases
    "KEGG_2021_Human",
    
    # Pathway Databases
    "Reactome_Pathways_2024",
    "WikiPathways_2024_Human", 
    "BioPlanet_2019",
    "WikiPathways_2024_Mouse",
    
    # Specialized Collections
    "Elsevier_Pathway_Collection",
    "MSigDB_Hallmark_2020",
    "BioCarta_2016",
    "HumanCyc_2016",
    "NCI-Nature_2016",
    "Panther_2016",
    
    # Coexpression Database
    "ARCHS4_Kinases_Coexp",
    
    # Phenotype and Disease Databases
    "MGI_Mammalian_Phenotype_Level_4_2024",
    "Jensen_DISEASES_Curated_2025",
    "Jensen_DISEASES_Experimental_2025",
    "Human_Phenotype_Ontology",
    "SynGO_2024",
    "KOMP2_Mouse_Phenotypes_2022",
    
    # Additional commonly used databases
    "ENCODE_TF_ChIP-seq_2015",
    "ChEA_2016",
    
    # Additional databases that may be available
    "Reactome_2022",
    "WikiPathways_2023_Human"
  )
  
  # Run enrichR analysis
  tryCatch({
    enrichr_results <- enrichr(gene_list, enrichr_databases)
    
    # Save gene list
    write_csv(data.frame(Gene_name = gene_list), 
              file.path(results_dir, "tables", paste0(output_name, "_gene_list.csv")))
    
    # Process each database
    all_significant_results <- data.frame()
    summary_stats <- data.frame()
    
    for (database in names(enrichr_results)) {
      result_df <- enrichr_results[[database]]
      
      if (nrow(result_df) > 0) {
        # Filter significant results
        significant <- result_df %>%
          filter(Adjusted.P.value < adj_p_cutoff) %>%
          arrange(Adjusted.P.value) %>%
          head(top_n)
        
        if (nrow(significant) > 0) {
          # Add database info
          significant$Database <- database
          
          # Save individual database results
          write_csv(significant, 
                    file.path(results_dir, "tables", paste0(output_name, "_", database, ".csv")))
          
          # Add to combined results (top 5 per database)
          top5 <- head(significant, 5)
          all_significant_results <- rbind(all_significant_results, top5)
          
          # Create plot
          plot_data <- significant %>%
            head(15) %>%
            mutate(
              Term_short = ifelse(nchar(Term) > 50, paste0(substr(Term, 1, 47), "..."), Term),
              log_p = -log10(Adjusted.P.value)
            )
          
          # Make plot
          p <- ggplot(plot_data, aes(x = log_p, y = reorder(Term_short, log_p))) +
            geom_col(fill = "steelblue", alpha = 0.7) +
            labs(
              title = paste("Enriched Terms -", database),
              subtitle = paste(output_name, "(", length(gene_list), "genes)"),
              x = "-log10(Adjusted P-value)",
              y = ""
            ) +
            theme_minimal() +
            theme(
              axis.text.y = element_text(size = 9),
              plot.background = element_rect(fill = "white", color = NA),
              panel.background = element_rect(fill = "white", color = NA),
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5)
            )
          
          # Save plot
          ggsave(file.path(results_dir, "plots", paste0(output_name, "_", database, ".png")), 
                 p, width = 10, height = 6, dpi = 300, bg = "white")
          
          # Add to summary
          summary_stats <- rbind(summary_stats, data.frame(
            Database = database,
            Significant_Terms = nrow(significant),
            Top_Term = significant$Term[1],
            Top_P_Value = significant$Adjusted.P.value[1],
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Save combined results
    if (nrow(all_significant_results) > 0) {
      write_csv(all_significant_results, 
                file.path(results_dir, "tables", paste0(output_name, "_all_top_results.csv")))
    }
    
    # Save summary
    if (nrow(summary_stats) > 0) {
      write_csv(summary_stats, 
                file.path(results_dir, "tables", paste0(output_name, "_summary.csv")))
      cat("  Found", nrow(summary_stats), "databases with significant results\n")
    } else {
      cat("  No significant results found\n")
    }
    
    cat("  Results saved to:", results_dir, "\n")
    return(list(results = enrichr_results, summary = summary_stats, directory = results_dir))
    
  }, error = function(e) {
    cat("Error running enrichR:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# ===== FUNCTION TO RUN clusterProfiler ANALYSIS =====
# Function to run clusterProfiler on a list of genes
# Note: For intersection lists, we can only do ORA (Over-Representation Analysis)
# GSEA requires ranked gene lists which we don't have for intersection results
run_clusterprofiler_analysis <- function(gene_list, output_name, results_base_dir = NULL, 
                                         pvalue_cutoff = 0.2, qvalue_cutoff = 1.0, 
                                         min_gssize = 10, max_gssize = 500) {
  # Check if required packages are available
  if (!clusterprofiler_available || !orgdb_available) {
    cat("Warning: clusterProfiler or org.Hs.eg.db not available. Skipping clusterProfiler analysis.\n")
    return(NULL)
  }
  
  # Clean gene list
  gene_list <- as.character(gene_list)
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  gene_list <- unique(gene_list)
  
  if (length(gene_list) == 0) {
    cat("Warning: Empty gene list provided for", output_name, "\n")
    return(NULL)
  }
  
  cat("\nRunning clusterProfiler for:", output_name, "(", length(gene_list), "genes)\n")
  
  # Set results directory
  if (is.null(results_base_dir)) {
    results_base_dir <- script_dir
  }
  results_dir <- file.path(results_base_dir, paste0(output_name, "_clusterProfiler_results"))
  dir.create(results_dir, showWarnings = FALSE)
  dir.create(file.path(results_dir, "plots"), showWarnings = FALSE)
  dir.create(file.path(results_dir, "tables"), showWarnings = FALSE)
  
  # Convert gene names to Entrez IDs
  tryCatch({
    gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    if (nrow(gene_ids) == 0) {
      cat("Warning: Could not convert any genes to Entrez IDs\n")
      return(NULL)
    }
    
    cat("  Converted", nrow(gene_ids), "genes to Entrez IDs\n")
    
    all_results <- list()
    
    # ===== GO Over-Representation Analysis =====
    cat("  Running GO ORA...\n")
    tryCatch({
      ego <- enrichGO(gene = gene_ids$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      keyType = "ENTREZID",
                      pAdjustMethod = "BH",
                      pvalueCutoff = pvalue_cutoff,
                      qvalueCutoff = qvalue_cutoff,
                      minGSSize = min_gssize,
                      maxGSSize = max_gssize,
                      readable = TRUE)
      
      if (!is.null(ego) && nrow(ego@result) > 0) {
        write.table(ego@result,
                    file = file.path(results_dir, "tables", paste0(output_name, "_GO_ORA_Results.txt")),
                    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        
        png(file.path(results_dir, "plots", paste0(output_name, "_GO_ORA.png")),
            width = 1000, height = 800)
        print(dotplot(ego, showCategory = 20))
        dev.off()
        
        all_results$GO_ORA <- ego
        cat("    Found", nrow(ego@result), "enriched GO terms\n")
      }
    }, error = function(e) {
      cat("    Error in GO ORA:", conditionMessage(e), "\n")
    })
    
    # ===== KEGG Over-Representation Analysis =====
    cat("  Running KEGG ORA...\n")
    tryCatch({
      kegg_enrich <- enrichKEGG(gene = gene_ids$ENTREZID,
                                organism = "hsa",
                                pAdjustMethod = "BH",
                                pvalueCutoff = pvalue_cutoff,
                                minGSSize = min_gssize,
                                maxGSSize = max_gssize)
      
      if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
        write.table(kegg_enrich@result,
                    file = file.path(results_dir, "tables", paste0(output_name, "_KEGG_ORA_Results.txt")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        
        png(file.path(results_dir, "plots", paste0(output_name, "_KEGG_ORA.png")),
            width = 600, height = 600)
        print(dotplot(kegg_enrich, showCategory = 20))
        dev.off()
        
        all_results$KEGG_ORA <- kegg_enrich
        cat("    Found", nrow(kegg_enrich@result), "enriched KEGG pathways\n")
      }
    }, error = function(e) {
      cat("    Error in KEGG ORA:", conditionMessage(e), "\n")
    })
    
    # ===== Reactome Over-Representation Analysis =====
    if (reactomepa_available) {
      cat("  Running Reactome ORA...\n")
      tryCatch({
        reactome_ora <- enrichPathway(gene = gene_ids$ENTREZID,
                                      organism = "human",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = pvalue_cutoff,
                                      minGSSize = min_gssize,
                                      maxGSSize = max_gssize)
        
        if (!is.null(reactome_ora) && nrow(reactome_ora@result) > 0) {
          write.table(reactome_ora@result,
                      file = file.path(results_dir, "tables", paste0(output_name, "_Reactome_ORA_Results.txt")),
                      sep = "\t", row.names = FALSE, quote = FALSE)
          
          png(file.path(results_dir, "plots", paste0(output_name, "_Reactome_ORA.png")),
              width = 1000, height = 800)
          print(dotplot(reactome_ora, showCategory = 20))
          dev.off()
          
          all_results$Reactome_ORA <- reactome_ora
          cat("    Found", nrow(reactome_ora@result), "enriched Reactome pathways\n")
        }
      }, error = function(e) {
        cat("    Error in Reactome ORA:", conditionMessage(e), "\n")
      })
    }
    
    # ===== WikiPathways Over-Representation Analysis =====
    cat("  Running WikiPathways ORA...\n")
    tryCatch({
      wikipathways_enrich <- enrichWP(gene = gene_ids$ENTREZID,
                                       organism = "Homo sapiens",
                                       pvalueCutoff = pvalue_cutoff,
                                       minGSSize = min_gssize,
                                       maxGSSize = max_gssize)
      
      if (!is.null(wikipathways_enrich) && nrow(wikipathways_enrich@result) > 0) {
        write.table(wikipathways_enrich@result,
                    file = file.path(results_dir, "tables", paste0(output_name, "_WikiPathways_ORA_Results.txt")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        
        png(file.path(results_dir, "plots", paste0(output_name, "_WikiPathways_ORA.png")),
            width = 1000, height = 800)
        print(dotplot(wikipathways_enrich, showCategory = 20))
        dev.off()
        
        all_results$WikiPathways_ORA <- wikipathways_enrich
        cat("    Found", nrow(wikipathways_enrich@result), "enriched WikiPathways\n")
      }
    }, error = function(e) {
      cat("    Error in WikiPathways ORA:", conditionMessage(e), "\n")
    })
    
    # ===== MSigDB Hallmark Over-Representation Analysis =====
    if (msigdbr_available) {
      cat("  Running MSigDB Hallmark ORA...\n")
      tryCatch({
        hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
        hallmark_t2g <- hallmark_sets %>%
          dplyr::select(gs_name, entrez_gene)
        
        hallmark_enrich <- enricher(gene = gene_ids$ENTREZID,
                                    TERM2GENE = hallmark_t2g,
                                    pvalueCutoff = pvalue_cutoff,
                                    pAdjustMethod = "BH",
                                    minGSSize = min_gssize,
                                    maxGSSize = max_gssize)
        
        if (!is.null(hallmark_enrich) && nrow(hallmark_enrich@result) > 0) {
          write.table(hallmark_enrich@result,
                      file = file.path(results_dir, "tables", paste0(output_name, "_MSigDB_Hallmark_ORA_Results.txt")),
                      sep = "\t", row.names = FALSE, quote = FALSE)
          
          png(file.path(results_dir, "plots", paste0(output_name, "_MSigDB_Hallmark_ORA.png")),
              width = 1000, height = 800)
          print(dotplot(hallmark_enrich, showCategory = 20))
          dev.off()
          
          all_results$MSigDB_Hallmark_ORA <- hallmark_enrich
          cat("    Found", nrow(hallmark_enrich@result), "enriched Hallmark gene sets\n")
        }
      }, error = function(e) {
        cat("    Error in MSigDB Hallmark ORA:", conditionMessage(e), "\n")
      })
    }
    
    cat("  Results saved to:", results_dir, "\n")
    return(list(results = all_results, directory = results_dir))
    
  }, error = function(e) {
    cat("Error running clusterProfiler:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# ===== STATISTICAL COMPARISON FUNCTIONS =====
# Calculate Jaccard similarity between two gene sets
calculate_jaccard <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union_size <- length(union(set1, set2))
  if (union_size == 0) return(0)
  return(intersection / union_size)
}

# Perform Fisher's exact test for two gene sets
fisher_exact_test <- function(set1, set2, background_size = 20000) {
  # Background is typically the total number of genes in the genome
  # Default to 20000 (approximate human protein-coding genes)
  
  set1_only <- length(setdiff(set1, set2))
  set2_only <- length(setdiff(set2, set1))
  both <- length(intersect(set1, set2))
  neither <- background_size - length(union(set1, set2))
  
  # Create contingency table
  contingency <- matrix(c(both, set1_only, set2_only, neither), nrow = 2, ncol = 2)
  
  # Perform Fisher's exact test
  result <- fisher.test(contingency)
  
  return(list(
    p_value = result$p.value,
    odds_ratio = result$estimate,
    contingency_table = contingency,
    set1_size = length(set1),
    set2_size = length(set2),
    intersection = both,
    union = length(union(set1, set2))
  ))
}

# Compare acute vs chronic for each cell type
compare_acute_chronic <- function(acute_var, chronic_var, cell_type_name) {
  if (!exists(acute_var, envir = .GlobalEnv) || !exists(chronic_var, envir = .GlobalEnv)) {
    return(NULL)
  }
  
  df_acute <- get(acute_var, envir = .GlobalEnv)
  df_chronic <- get(chronic_var, envir = .GlobalEnv)
  
  # Get significant genes
  if ("Sig" %in% colnames(df_acute) && "Gene_name" %in% colnames(df_acute)) {
    genes_acute <- unique(df_acute$Gene_name[df_acute$Sig == 1 & !is.na(df_acute$Gene_name) & df_acute$Gene_name != ""])
  } else {
    return(NULL)
  }
  
  if ("Sig" %in% colnames(df_chronic) && "Gene_name" %in% colnames(df_chronic)) {
    genes_chronic <- unique(df_chronic$Gene_name[df_chronic$Sig == 1 & !is.na(df_chronic$Gene_name) & df_chronic$Gene_name != ""])
  } else {
    return(NULL)
  }
  
  # Calculate statistics
  jaccard <- calculate_jaccard(genes_acute, genes_chronic)
  fisher_result <- fisher_exact_test(genes_acute, genes_chronic)
  
  return(list(
    cell_type = cell_type_name,
    acute_count = length(genes_acute),
    chronic_count = length(genes_chronic),
    common = length(intersect(genes_acute, genes_chronic)),
    jaccard_similarity = jaccard,
    fisher_p_value = fisher_result$p_value,
    fisher_odds_ratio = fisher_result$odds_ratio,
    acute_only = length(setdiff(genes_acute, genes_chronic)),
    chronic_only = length(setdiff(genes_chronic, genes_acute))
  ))
}

# ===== GENE-GENE INTERACTION NETWORK FUNCTION =====
# Create gene-gene interaction network for common genes
create_gene_network <- function(gene_list, output_name, results_base_dir = NULL, 
                                min_interactions = 1) {
  if (!igraph_available) {
    cat("Warning: igraph not available. Skipping network creation.\n")
    return(NULL)
  }
  
  if (length(gene_list) < 2) {
    cat("Warning: Need at least 2 genes for network. Skipping", output_name, "\n")
    return(NULL)
  }
  
  cat("\nCreating gene network for:", output_name, "(", length(gene_list), "genes)\n")
  
  # Set results directory
  if (is.null(results_base_dir)) {
    results_base_dir <- script_dir
  }
  
  # For now, create a simple co-occurrence network
  # In a real scenario, you would query a database like STRING, BioGRID, etc.
  # This creates a simple network where genes are connected if they appear together
  
  tryCatch({
    # Create a simple network (all genes connected to each other)
    # In practice, you'd query a protein-protein interaction database
    if (length(gene_list) <= 50) {
      # For small gene lists, create a complete graph
      g <- make_full_graph(length(gene_list), directed = FALSE)
      V(g)$name <- gene_list
      V(g)$label <- gene_list
      
      # Save network as graphml
      network_file <- file.path(results_base_dir, paste0(output_name, "_network.graphml"))
      write_graph(g, network_file, format = "graphml")
      cat("  Network saved to:", network_file, "\n")
      
      # Create a simple visualization
      if (length(gene_list) <= 20) {
        png(file.path(results_base_dir, paste0(output_name, "_network.png")),
            width = 1200, height = 1200, res = 300)
        plot(g, vertex.label.cex = 0.8, vertex.size = 8, 
             layout = layout_with_fr(g), main = output_name)
        dev.off()
        cat("  Network plot saved\n")
      }
      
      return(g)
    } else {
      cat("  Gene list too large for simple network visualization (", length(gene_list), "genes)\n")
      cat("  Consider using specialized tools like STRING or Cytoscape\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("Error creating network:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# ===== EXCEL WORKBOOK CREATION FUNCTION =====
create_master_excel_report <- function(results_base_dir = NULL) {
  if (!openxlsx_available) {
    cat("Warning: openxlsx not available. Skipping Excel report.\n")
    return(NULL)
  }
  
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("CREATING MASTER EXCEL REPORT\n")
  cat(paste(rep("=", 60), collapse=""), "\n\n")
  
  if (is.null(results_base_dir)) {
    results_base_dir <- script_dir
  }
  
  excel_file <- file.path(results_base_dir, "Master_Summary_Report.xlsx")
  wb <- createWorkbook()
  
  # Sheet 1: Summary Statistics
  if (exists("sig_summary", envir = .GlobalEnv)) {
    addWorksheet(wb, "Summary_Statistics")
    writeData(wb, "Summary_Statistics", get("sig_summary", envir = .GlobalEnv))
    cat("Added Summary_Statistics sheet\n")
  }
  
  # Sheet 2: Acute vs Chronic Intersection
  if (exists("intersection_summary", envir = .GlobalEnv)) {
    addWorksheet(wb, "Acute_Chronic_Intersection")
    writeData(wb, "Acute_Chronic_Intersection", get("intersection_summary", envir = .GlobalEnv))
    cat("Added Acute_Chronic_Intersection sheet\n")
  }
  
  # Sheet 3: Imm Intersections
  if (exists("imm_intersection_summary", envir = .GlobalEnv)) {
    addWorksheet(wb, "Imm_Acute_Intersection")
    writeData(wb, "Imm_Acute_Intersection", get("imm_intersection_summary", envir = .GlobalEnv))
    cat("Added Imm_Acute_Intersection sheet\n")
  }
  
  if (exists("imm_intersection_summary_chronic", envir = .GlobalEnv)) {
    addWorksheet(wb, "Imm_Chronic_Intersection")
    writeData(wb, "Imm_Chronic_Intersection", get("imm_intersection_summary_chronic", envir = .GlobalEnv))
    cat("Added Imm_Chronic_Intersection sheet\n")
  }
  
  # Sheet 4: Olig Intersections
  if (exists("olig_intersection_summary_acute", envir = .GlobalEnv)) {
    addWorksheet(wb, "Olig_Acute_Intersection")
    writeData(wb, "Olig_Acute_Intersection", get("olig_intersection_summary_acute", envir = .GlobalEnv))
    cat("Added Olig_Acute_Intersection sheet\n")
  }
  
  if (exists("olig_intersection_summary_chronic", envir = .GlobalEnv)) {
    addWorksheet(wb, "Olig_Chronic_Intersection")
    writeData(wb, "Olig_Chronic_Intersection", get("olig_intersection_summary_chronic", envir = .GlobalEnv))
    cat("Added Olig_Chronic_Intersection sheet\n")
  }
  
  # Sheet 5: Glia Intersections
  if (exists("glia_intersection_summary_acute", envir = .GlobalEnv)) {
    addWorksheet(wb, "Glia_Acute_Intersection")
    writeData(wb, "Glia_Acute_Intersection", get("glia_intersection_summary_acute", envir = .GlobalEnv))
    cat("Added Glia_Acute_Intersection sheet\n")
  }
  
  if (exists("glia_intersection_summary_chronic", envir = .GlobalEnv)) {
    addWorksheet(wb, "Glia_Chronic_Intersection")
    writeData(wb, "Glia_Chronic_Intersection", get("glia_intersection_summary_chronic", envir = .GlobalEnv))
    cat("Added Glia_Chronic_Intersection sheet\n")
  }
  
  # Sheet 6: Statistical Comparisons
  if (exists("statistical_comparisons", envir = .GlobalEnv)) {
    addWorksheet(wb, "Statistical_Comparisons")
    writeData(wb, "Statistical_Comparisons", get("statistical_comparisons", envir = .GlobalEnv))
    cat("Added Statistical_Comparisons sheet\n")
  }
  
  # Sheet 7: Jaccard Similarity Matrix
  if (exists("jaccard_similarity_matrix", envir = .GlobalEnv)) {
    addWorksheet(wb, "Jaccard_Similarity")
    writeData(wb, "Jaccard_Similarity", get("jaccard_similarity_matrix", envir = .GlobalEnv))
    cat("Added Jaccard_Similarity sheet\n")
  }
  
  # Sheet 8: File to Variable Mapping
  if (exists("file_to_variable_mapping", envir = .GlobalEnv)) {
    addWorksheet(wb, "File_Variable_Mapping")
    writeData(wb, "File_Variable_Mapping", get("file_to_variable_mapping", envir = .GlobalEnv))
    cat("Added File_Variable_Mapping sheet\n")
  }
  
  # Sheet 9: Gene Lists Summary (organized by category)
  gene_lists_summary <- data.frame(
    Category = character(),
    Gene_List_Name = character(),
    Gene_Count = integer(),
    File_Name = character(),
    stringsAsFactors = FALSE
  )
  
  # Add gene lists from intersections
  intersection_file_patterns <- c(
    "imm_acute_common_all_three_genes.csv",
    "imm_chronic_common_all_three_genes.csv",
    "olig_acute_common_all_three_genes.csv",
    "olig_chronic_common_all_three_genes.csv",
    "glia_acute_common_both_genes.csv",
    "glia_chronic_common_both_genes.csv"
  )
  
  for (pattern in intersection_file_patterns) {
    file_path <- file.path(results_base_dir, pattern)
    if (file.exists(file_path)) {
      tryCatch({
        gene_data <- read_csv(file_path, col_types = cols(), show_col_types = FALSE)
        gene_count <- nrow(gene_data)
        category <- gsub("_.*", "", pattern)
        gene_lists_summary <- rbind(gene_lists_summary, data.frame(
          Category = category,
          Gene_List_Name = gsub("\\.csv$", "", pattern),
          Gene_Count = gene_count,
          File_Name = pattern,
          stringsAsFactors = FALSE
        ))
      }, error = function(e) {
        # Skip if file can't be read
      })
    }
  }
  
  if (nrow(gene_lists_summary) > 0) {
    addWorksheet(wb, "Gene_Lists_Summary")
    writeData(wb, "Gene_Lists_Summary", gene_lists_summary)
    cat("Added Gene_Lists_Summary sheet\n")
  }
  
  # Save workbook
  saveWorkbook(wb, excel_file, overwrite = TRUE)
  cat("\nMaster Excel report saved to:", excel_file, "\n")
  cat("Workbook contains", length(names(wb)), "sheets\n")
  return(excel_file)
}

# ===== HTML REPORT CREATION FUNCTION =====
create_html_report <- function(results_base_dir = NULL) {
  if (!dt_available) {
    cat("Warning: DT package not available. Skipping HTML report.\n")
    return(NULL)
  }
  
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("CREATING HTML REPORT\n")
  cat(paste(rep("=", 60), collapse=""), "\n\n")
  
  if (is.null(results_base_dir)) {
    results_base_dir <- script_dir
  }
  
  html_file <- file.path(results_base_dir, "Master_Summary_Report.html")
  
  # Start HTML document
  html_content <- c(
    "<!DOCTYPE html>",
    "<html>",
    "<head>",
    "<title>Hypoxia Analysis Summary Report</title>",
    "<meta charset='UTF-8'>",
    "<style>",
    "body { font-family: Arial, sans-serif; margin: 20px; }",
    "h1 { color: #2c3e50; }",
    "h2 { color: #34495e; margin-top: 30px; }",
    "table { border-collapse: collapse; width: 100%; margin: 20px 0; }",
    "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }",
    "th { background-color: #3498db; color: white; }",
    "tr:nth-child(even) { background-color: #f2f2f2; }",
    "</style>",
    "</head>",
    "<body>",
    "<h1>Hypoxia Analysis Summary Report</h1>",
    "<p>Generated on:", Sys.time(), "</p>"
  )
  
  # Add summary statistics
  if (exists("sig_summary", envir = .GlobalEnv)) {
    html_content <- c(html_content,
      "<h2>Summary Statistics</h2>",
      "<p>Overview of all loaded files and gene counts:</p>"
    )
    
    # Create HTML table from sig_summary
    sig_df <- get("sig_summary", envir = .GlobalEnv)
    html_table <- paste0("<table><tr>", 
                         paste0("<th>", colnames(sig_df), "</th>", collapse = ""),
                         "</tr>")
    
    for (i in seq_len(nrow(sig_df))) {
      html_table <- paste0(html_table, "<tr>",
                           paste0("<td>", sig_df[i, ], "</td>", collapse = ""),
                           "</tr>")
    }
    html_table <- paste0(html_table, "</table>")
    html_content <- c(html_content, html_table)
  }
  
  # Add intersection summary
  if (exists("intersection_summary", envir = .GlobalEnv)) {
    html_content <- c(html_content,
      "<h2>Acute vs Chronic Intersection Analysis</h2>",
      "<p>Comparison of gene sets between acute and chronic conditions:</p>"
    )
    
    int_df <- get("intersection_summary", envir = .GlobalEnv)
    html_table <- paste0("<table><tr>", 
                         paste0("<th>", colnames(int_df), "</th>", collapse = ""),
                         "</tr>")
    
    for (i in seq_len(nrow(int_df))) {
      html_table <- paste0(html_table, "<tr>",
                           paste0("<td>", int_df[i, ], "</td>", collapse = ""),
                           "</tr>")
    }
    html_table <- paste0(html_table, "</table>")
    html_content <- c(html_content, html_table)
  }
  
  # Add statistical comparisons
  if (exists("statistical_comparisons", envir = .GlobalEnv)) {
    html_content <- c(html_content,
      "<h2>Statistical Comparisons: Acute vs Chronic</h2>",
      "<p>Fisher's exact test and Jaccard similarity between acute and chronic conditions:</p>"
    )
    
    stats_df <- get("statistical_comparisons", envir = .GlobalEnv)
    html_table <- paste0("<table><tr>", 
                         paste0("<th>", colnames(stats_df), "</th>", collapse = ""),
                         "</tr>")
    
    for (i in seq_len(nrow(stats_df))) {
      html_table <- paste0(html_table, "<tr>",
                           paste0("<td>", stats_df[i, ], "</td>", collapse = ""),
                           "</tr>")
    }
    html_table <- paste0(html_table, "</table>")
    html_content <- c(html_content, html_table)
  }
  
  # Add Jaccard similarity matrix
  if (exists("jaccard_similarity_matrix", envir = .GlobalEnv)) {
    html_content <- c(html_content,
      "<h2>Jaccard Similarity Between Intersections</h2>",
      "<p>Similarity scores between common gene sets from different cell type groups:</p>"
    )
    
    jaccard_df <- get("jaccard_similarity_matrix", envir = .GlobalEnv)
    html_table <- paste0("<table><tr>", 
                         paste0("<th>", colnames(jaccard_df), "</th>", collapse = ""),
                         "</tr>")
    
    for (i in seq_len(nrow(jaccard_df))) {
      html_table <- paste0(html_table, "<tr>",
                           paste0("<td>", jaccard_df[i, ], "</td>", collapse = ""),
                           "</tr>")
    }
    html_table <- paste0(html_table, "</table>")
    html_content <- c(html_content, html_table)
  }
  
  # Add summary of intersection results
  html_content <- c(html_content,
    "<h2>Analysis Summary</h2>",
    "<p>This report includes:</p>",
    "<ul>",
    "<li>Summary statistics for all loaded files</li>",
    "<li>Acute vs Chronic intersection analysis</li>",
    "<li>Statistical comparisons (Fisher's test, Jaccard similarity)</li>",
    "<li>Intersection analysis for Imm, Oligodendrocyte, and Glia cell types</li>",
    "<li>Gene lists organized by category</li>",
    "</ul>",
    "<p><strong>Note:</strong> For detailed enrichment results, see the enrichR and clusterProfiler result directories.</p>"
  )
  
  # Close HTML
  html_content <- c(html_content, "</body>", "</html>")
  
  # Write to file
  writeLines(html_content, html_file)
  cat("HTML report saved to:", html_file, "\n")
  return(html_file)
}

# ===== FIND ALL CSV FILES =====
# Find all CSV files matching the pattern
csv_files <- list.files(
  path = script_dir,
  pattern = "^cell_type_ISP_.*\\.csv$",
  full.names = TRUE
)

cat("\nFound", length(csv_files), "CSV files:\n")
print(basename(csv_files))

# ===== SHOW FILE TO VARIABLE NAME MAPPING =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("FILE TO VARIABLE NAME MAPPING:\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

file_var_mapping <- data.frame(
  File = character(),
  Variable = character(),
  stringsAsFactors = FALSE
)

for (file in csv_files) {
  file_name <- basename(file)
  
  # Create variable name using same logic as reading
  var_name <- gsub("\\.csv$", "", file_name)  # Remove .csv extension
  var_name <- gsub("^cell_type_ISP_", "", var_name)  # Remove prefix
  var_name <- gsub("_GPU_ONLY_genes$", "", var_name)  # Remove suffix
  var_name <- gsub("-", "_", var_name)  # Replace hyphens with underscores
  var_name <- gsub("[^A-Za-z0-9_]", "_", var_name)  # Replace invalid characters
  if (grepl("^[0-9]", var_name)) {
    var_name <- paste0("X", var_name)
  }
  
  # Add to mapping
  file_var_mapping <- rbind(file_var_mapping, 
                            data.frame(File = file_name, Variable = var_name, stringsAsFactors = FALSE))
  
  # Print mapping
  cat(sprintf("%-60s -> %s\n", file_name, var_name))
}

cat("\n")

# ===== READ ALL CSV FILES =====
# Create vectors to track variable names and corresponding files
created_variables <- c()
successful_files <- c()

cat("\nReading CSV files...\n")
cat(paste(rep("=", 60), collapse=""), "\n")

for (file in csv_files) {
  file_name <- basename(file)
  cat("Reading:", file_name, "...")
  
  # Read the CSV file
  tryCatch({
    df <- read_csv(file, col_types = cols(), show_col_types = FALSE)
    
    # Create a valid R variable name from the file name
    # Remove .csv extension
    var_name <- gsub("\\.csv$", "", file_name)
    # Remove "cell_type_ISP_" prefix
    var_name <- gsub("^cell_type_ISP_", "", var_name)
    # Remove "_GPU_ONLY_genes" suffix
    var_name <- gsub("_GPU_ONLY_genes$", "", var_name)
    # Replace hyphens with underscores for valid R variable names
    var_name <- gsub("-", "_", var_name)
    # Replace any other invalid characters with underscores
    var_name <- gsub("[^A-Za-z0-9_]", "_", var_name)
    # Ensure it doesn't start with a number
    if (grepl("^[0-9]", var_name)) {
      var_name <- paste0("X", var_name)
    }
    
    # Assign the data frame to a variable in the global environment
    assign(var_name, df, envir = .GlobalEnv)
    created_variables <- c(created_variables, var_name)
    successful_files <- c(successful_files, file_name)
    
    cat(" OK -> variable:", var_name, "(", nrow(df), "rows,", ncol(df), "columns)\n")
    
  }, error = function(e) {
    cat(" ERROR:", conditionMessage(e), "\n")
  })
}

cat(paste(rep("=", 60), collapse=""), "\n")
cat("Successfully read", length(created_variables), "files\n\n")

# ===== CREATE FILTERED VARIABLES FOR SIGNIFICANT GENES =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("CREATING FILTERED VARIABLES (Significant genes only, Sig=1):\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Track created significant-only variables
sig_variables <- c()

for (var_name in created_variables) {
  df <- get(var_name, envir = .GlobalEnv)
  
  if ("Sig" %in% colnames(df)) {
    # Filter for significant genes (Sig == 1)
    df_sig <- df %>% filter(Sig == 1)
    
    # Create new variable name with "_sig" suffix
    sig_var_name <- paste0(var_name, "_sig")
    
    # Assign the filtered data frame to the new variable
    assign(sig_var_name, df_sig, envir = .GlobalEnv)
    sig_variables <- c(sig_variables, sig_var_name)
    
    cat("Created:", sig_var_name, "->", nrow(df_sig), "significant genes\n")
  }
}

cat("\nCreated", length(sig_variables), "filtered variables for significant genes\n")
cat("Variable names end with '_sig' suffix\n\n")

# ===== DISPLAY SUMMARY =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("SUMMARY OF LOADED DATA\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

for (var_name in created_variables) {
  df <- get(var_name, envir = .GlobalEnv)
  cat("Variable:", var_name, "\n")
  cat("  Dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
  cat("  Columns:", paste(colnames(df), collapse = ", "), "\n")
  
  # Check for Sig column and count significant and non-significant genes
  if ("Sig" %in% colnames(df)) {
    sig_count <- sum(df$Sig == 1, na.rm = TRUE)
    nonsig_count <- sum(df$Sig == 0, na.rm = TRUE)
    total_genes <- nrow(df)
    cat("  Significant genes (Sig=1):", sig_count, "\n")
    cat("  Non-significant genes (Sig=0):", nonsig_count, "\n")
    cat("  Total genes:", total_genes, "\n")
  }
  
  cat("\n")
}

# ===== SUMMARY OF SIGNIFICANT GENES =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("SUMMARY: SIGNIFICANT vs NON-SIGNIFICANT GENES\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Create summary data frame with file name, variable name, significant-only variable name, and gene counts
sig_summary <- data.frame(
  File_Name = character(),
  Variable_Name = character(),
  Significant_Only_Variable_Name = character(),
  Total_Genes = integer(),
  Significant_Sig1 = integer(),
  Non_Significant_Sig0 = integer(),
  stringsAsFactors = FALSE
)

for (i in seq_along(created_variables)) {
  var_name <- created_variables[i]
  file_name <- successful_files[i]
  df <- get(var_name, envir = .GlobalEnv)
  
  if ("Sig" %in% colnames(df)) {
    sig_count <- sum(df$Sig == 1, na.rm = TRUE)
    nonsig_count <- sum(df$Sig == 0, na.rm = TRUE)
    total_count <- nrow(df)
    
    # Get the corresponding significant-only variable name
    # (constructed as var_name + "_sig")
    sig_only_var <- paste0(var_name, "_sig")
    
    sig_summary <- rbind(sig_summary, data.frame(
      File_Name = file_name,
      Variable_Name = var_name,
      Significant_Only_Variable_Name = sig_only_var,
      Total_Genes = total_count,
      Significant_Sig1 = sig_count,
      Non_Significant_Sig0 = nonsig_count,
      stringsAsFactors = FALSE
    ))
  }
}

# Print summary table
if (nrow(sig_summary) > 0) {
  print(sig_summary)
  cat("\n")
  
  # Store summary in a variable
  assign("sig_summary", sig_summary, envir = .GlobalEnv)
  cat("Summary stored in variable: sig_summary\n")
  
  # Write summary to CSV file
  summary_file <- file.path(script_dir, "file_summary.csv")
  write_csv(sig_summary, summary_file)
  cat("Summary written to file:", summary_file, "\n")
}

# ===== INTERSECTION ANALYSIS: ACUTE vs CHRONIC =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION ANALYSIS: ACUTE vs CHRONIC GENE NAMES\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Create data frame to store intersection results
intersection_summary <- data.frame(
  Cell_Type = character(),
  Acute_Variable = character(),
  Chronic_Variable = character(),
  Acute_Total_Genes = integer(),
  Chronic_Total_Genes = integer(),
  Common_Genes_Count = integer(),
  Acute_Only_Count = integer(),
  Chronic_Only_Count = integer(),
  Pct_Common_of_Acute = numeric(),
  Pct_Common_of_Chronic = numeric(),
  stringsAsFactors = FALSE
)

# Find all acute and chronic variables
acute_vars <- created_variables[grepl("_acute$", created_variables)]
chronic_vars <- created_variables[grepl("_chronic$", created_variables)]

cat("Found", length(acute_vars), "acute variables and", length(chronic_vars), "chronic variables\n\n")

# For each acute variable, find its corresponding chronic variable
for (acute_var in acute_vars) {
  # Extract base name (remove "_acute" suffix)
  base_name <- gsub("_acute$", "", acute_var)
  
  # Find corresponding chronic variable
  chronic_var <- paste0(base_name, "_chronic")
  
  if (chronic_var %in% chronic_vars) {
    # Get data frames
    df_acute <- get(acute_var, envir = .GlobalEnv)
    df_chronic <- get(chronic_var, envir = .GlobalEnv)
    
    # Extract gene names (assuming Gene_name column exists)
    if ("Gene_name" %in% colnames(df_acute) && "Gene_name" %in% colnames(df_chronic)) {
      genes_acute <- unique(df_acute$Gene_name[!is.na(df_acute$Gene_name) & df_acute$Gene_name != ""])
      genes_chronic <- unique(df_chronic$Gene_name[!is.na(df_chronic$Gene_name) & df_chronic$Gene_name != ""])
      
      # Find intersection (common genes)
      common_genes <- intersect(genes_acute, genes_chronic)
      acute_only <- setdiff(genes_acute, genes_chronic)
      chronic_only <- setdiff(genes_chronic, genes_acute)
      
      # Calculate percentages
      pct_common_of_acute <- ifelse(length(genes_acute) > 0, 
                                     (length(common_genes) / length(genes_acute)) * 100, 
                                     0)
      pct_common_of_chronic <- ifelse(length(genes_chronic) > 0, 
                                       (length(common_genes) / length(genes_chronic)) * 100, 
                                       0)
      
      # Add to summary
      intersection_summary <- rbind(intersection_summary, data.frame(
        Cell_Type = base_name,
        Acute_Variable = acute_var,
        Chronic_Variable = chronic_var,
        Acute_Total_Genes = length(genes_acute),
        Chronic_Total_Genes = length(genes_chronic),
        Common_Genes_Count = length(common_genes),
        Acute_Only_Count = length(acute_only),
        Chronic_Only_Count = length(chronic_only),
        Pct_Common_of_Acute = round(pct_common_of_acute, 2),
        Pct_Common_of_Chronic = round(pct_common_of_chronic, 2),
        stringsAsFactors = FALSE
      ))
      
      cat("Cell Type:", base_name, "\n")
      cat("  Acute variable:", acute_var, "->", length(genes_acute), "genes\n")
      cat("  Chronic variable:", chronic_var, "->", length(genes_chronic), "genes\n")
      cat("  Common genes:", length(common_genes), "\n")
      cat("  Common genes as % of acute:", round(pct_common_of_acute, 2), "%\n")
      cat("  Common genes as % of chronic:", round(pct_common_of_chronic, 2), "%\n")
      cat("  Acute only:", length(acute_only), "\n")
      cat("  Chronic only:", length(chronic_only), "\n")
      cat("\n")
    } else {
      cat("Warning: Gene_name column not found in", acute_var, "or", chronic_var, "\n")
    }
  } else {
    cat("Warning: No corresponding chronic variable found for", acute_var, "\n")
  }
}

# Print summary table
if (nrow(intersection_summary) > 0) {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("INTERSECTION SUMMARY TABLE:\n")
  cat(paste(rep("=", 60), collapse=""), "\n\n")
  print(intersection_summary)
  
  # Store in variable
  assign("intersection_summary", intersection_summary, envir = .GlobalEnv)
  cat("\nIntersection summary stored in variable: intersection_summary\n")
  
  # Write to CSV file
  intersection_file <- file.path(script_dir, "acute_chronic_intersection.csv")
  write_csv(intersection_summary, intersection_file)
  cat("Intersection summary written to file:", intersection_file, "\n")
} else {
  cat("No intersection analysis could be performed.\n")
}

# ===== STATISTICAL COMPARISONS: ACUTE vs CHRONIC =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("STATISTICAL COMPARISONS: ACUTE vs CHRONIC\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Create data frame for statistical comparisons
statistical_comparisons <- data.frame(
  Cell_Type = character(),
  Acute_Count = integer(),
  Chronic_Count = integer(),
  Common_Genes = integer(),
  Jaccard_Similarity = numeric(),
  Fisher_P_Value = numeric(),
  Fisher_Odds_Ratio = numeric(),
  Acute_Only = integer(),
  Chronic_Only = integer(),
  stringsAsFactors = FALSE
)

# Find all acute and chronic variable pairs
acute_vars <- created_variables[grepl("_acute$", created_variables)]
chronic_vars <- created_variables[grepl("_chronic$", created_variables)]

for (acute_var in acute_vars) {
  base_name <- gsub("_acute$", "", acute_var)
  chronic_var <- paste0(base_name, "_chronic")
  
  if (chronic_var %in% chronic_vars) {
    comparison <- compare_acute_chronic(acute_var, chronic_var, base_name)
    
    if (!is.null(comparison)) {
      statistical_comparisons <- rbind(statistical_comparisons, data.frame(
        Cell_Type = comparison$cell_type,
        Acute_Count = comparison$acute_count,
        Chronic_Count = comparison$chronic_count,
        Common_Genes = comparison$common,
        Jaccard_Similarity = round(comparison$jaccard_similarity, 4),
        Fisher_P_Value = comparison$fisher_p_value,
        Fisher_Odds_Ratio = round(as.numeric(comparison$fisher_odds_ratio), 4),
        Acute_Only = comparison$acute_only,
        Chronic_Only = comparison$chronic_only,
        stringsAsFactors = FALSE
      ))
      
      cat("Cell Type:", base_name, "\n")
      cat("  Jaccard Similarity:", round(comparison$jaccard_similarity, 4), "\n")
      cat("  Fisher's Exact Test p-value:", format(comparison$fisher_p_value, scientific = TRUE), "\n")
      cat("  Odds Ratio:", round(as.numeric(comparison$fisher_odds_ratio), 4), "\n")
      cat("\n")
    }
  }
}

# Print and save statistical comparisons
if (nrow(statistical_comparisons) > 0) {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("STATISTICAL COMPARISONS SUMMARY:\n")
  cat(paste(rep("=", 60), collapse=""), "\n\n")
  print(statistical_comparisons)
  
  # Store in variable
  assign("statistical_comparisons", statistical_comparisons, envir = .GlobalEnv)
  cat("\nStatistical comparisons stored in variable: statistical_comparisons\n")
  
  # Write to CSV
  stats_file <- file.path(script_dir, "statistical_comparisons_acute_vs_chronic.csv")
  write_csv(statistical_comparisons, stats_file)
  cat("Statistical comparisons written to file:", stats_file, "\n")
}

# ===== JACCARD SIMILARITY BETWEEN INTERSECTIONS =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("JACCARD SIMILARITY BETWEEN INTERSECTION GENE LISTS\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Calculate Jaccard similarity between different intersections
# This will compare common genes from different cell type groups
jaccard_matrix <- data.frame(
  Comparison = character(),
  Set1 = character(),
  Set2 = character(),
  Jaccard_Similarity = numeric(),
  Intersection_Size = integer(),
  Union_Size = integer(),
  stringsAsFactors = FALSE
)

# Compare common genes from different intersections
intersection_genes <- list()

# Collect common genes from all intersections
if (exists("imm_common_all_three", envir = .GlobalEnv)) {
  intersection_genes[["Imm_Acute_Common"]] <- get("imm_common_all_three", envir = .GlobalEnv)
}
if (exists("imm_common_all_three_chronic", envir = .GlobalEnv)) {
  intersection_genes[["Imm_Chronic_Common"]] <- get("imm_common_all_three_chronic", envir = .GlobalEnv)
}
if (exists("olig_common_all_three_acute", envir = .GlobalEnv)) {
  intersection_genes[["Olig_Acute_Common"]] <- get("olig_common_all_three_acute", envir = .GlobalEnv)
}
if (exists("olig_common_all_three_chronic", envir = .GlobalEnv)) {
  intersection_genes[["Olig_Chronic_Common"]] <- get("olig_common_all_three_chronic", envir = .GlobalEnv)
}
if (exists("glia_common_both_acute", envir = .GlobalEnv)) {
  intersection_genes[["Glia_Acute_Common"]] <- get("glia_common_both_acute", envir = .GlobalEnv)
}
if (exists("glia_common_both_chronic", envir = .GlobalEnv)) {
  intersection_genes[["Glia_Chronic_Common"]] <- get("glia_common_both_chronic", envir = .GlobalEnv)
}

# Calculate pairwise Jaccard similarities
if (length(intersection_genes) > 1) {
  gene_set_names <- names(intersection_genes)
  
  for (i in 1:(length(gene_set_names) - 1)) {
    for (j in (i + 1):length(gene_set_names)) {
      set1_name <- gene_set_names[i]
      set2_name <- gene_set_names[j]
      set1 <- intersection_genes[[set1_name]]
      set2 <- intersection_genes[[set2_name]]
      
      if (length(set1) > 0 && length(set2) > 0) {
        jaccard <- calculate_jaccard(set1, set2)
        intersection_size <- length(intersect(set1, set2))
        union_size <- length(union(set1, set2))
        
        jaccard_matrix <- rbind(jaccard_matrix, data.frame(
          Comparison = paste(set1_name, "vs", set2_name),
          Set1 = set1_name,
          Set2 = set2_name,
          Jaccard_Similarity = round(jaccard, 4),
          Intersection_Size = intersection_size,
          Union_Size = union_size,
          stringsAsFactors = FALSE
        ))
        
        cat(set1_name, "vs", set2_name, ":", round(jaccard, 4), "\n")
      }
    }
  }
  
  if (nrow(jaccard_matrix) > 0) {
    cat("\n", paste(rep("=", 60), collapse=""), "\n")
    cat("JACCARD SIMILARITY MATRIX:\n")
    cat(paste(rep("=", 60), collapse=""), "\n\n")
    print(jaccard_matrix)
    
    # Store in variable
    assign("jaccard_similarity_matrix", jaccard_matrix, envir = .GlobalEnv)
    cat("\nJaccard similarity matrix stored in variable: jaccard_similarity_matrix\n")
    
    # Write to CSV
    jaccard_file <- file.path(script_dir, "jaccard_similarity_intersections.csv")
    write_csv(jaccard_matrix, jaccard_file)
    cat("Jaccard similarity matrix written to file:", jaccard_file, "\n")
  }
}

# ===== INTERSECTION OF SIGNIFICANT GENES: Imm_CGE_acute, Imm_LGE_acute, Imm_MGE_acute =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION OF SIGNIFICANT GENES: Imm_CGE_acute, Imm_LGE_acute, Imm_MGE_acute\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Variable names to analyze
imm_vars <- c("Imm_CGE_acute", "Imm_LGE_acute", "Imm_MGE_acute")

# Check if variables exist and extract significant genes
imm_sig_genes_list <- list()
imm_var_info <- data.frame(
  Variable = character(),
  Total_Genes = integer(),
  Significant_Genes = integer(),
  stringsAsFactors = FALSE
)

for (var_name in imm_vars) {
  if (exists(var_name, envir = .GlobalEnv)) {
    df <- get(var_name, envir = .GlobalEnv)
    
    # Filter for significant genes (Sig == 1)
    if ("Sig" %in% colnames(df)) {
      df_sig <- df %>% filter(Sig == 1)
      
      # Extract gene names
      if ("Gene_name" %in% colnames(df_sig)) {
        sig_genes <- unique(df_sig$Gene_name[!is.na(df_sig$Gene_name) & df_sig$Gene_name != ""])
        imm_sig_genes_list[[var_name]] <- sig_genes
        
        imm_var_info <- rbind(imm_var_info, data.frame(
          Variable = var_name,
          Total_Genes = nrow(df),
          Significant_Genes = length(sig_genes),
          stringsAsFactors = FALSE
        ))
        
        cat("Variable:", var_name, "\n")
        cat("  Total genes:", nrow(df), "\n")
        cat("  Significant genes (Sig=1):", length(sig_genes), "\n\n")
      } else {
        cat("Warning: Gene_name column not found in", var_name, "\n")
      }
    } else {
      cat("Warning: Sig column not found in", var_name, "\n")
    }
  } else {
    cat("Warning: Variable", var_name, "not found\n")
  }
}

# Find intersections
if (length(imm_sig_genes_list) == 3) {
  genes_cge <- imm_sig_genes_list[["Imm_CGE_acute"]]
  genes_lge <- imm_sig_genes_list[["Imm_LGE_acute"]]
  genes_mge <- imm_sig_genes_list[["Imm_MGE_acute"]]
  
  # Pairwise intersections
  cge_lge <- intersect(genes_cge, genes_lge)
  cge_mge <- intersect(genes_cge, genes_mge)
  lge_mge <- intersect(genes_lge, genes_mge)
  
  # Three-way intersection (common to all three)
  common_all_three <- intersect(intersect(genes_cge, genes_lge), genes_mge)
  
  # Unique to each
  only_cge <- setdiff(setdiff(genes_cge, genes_lge), genes_mge)
  only_lge <- setdiff(setdiff(genes_lge, genes_cge), genes_mge)
  only_mge <- setdiff(setdiff(genes_mge, genes_cge), genes_lge)
  
  # Display results
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("INTERSECTION RESULTS:\n")
  cat(paste(rep("=", 60), collapse=""), "\n\n")
  
  cat("Significant genes in each variable:\n")
  print(imm_var_info)
  cat("\n")
  
  cat("Pairwise intersections:\n")
  cat("  Imm_CGE_acute ∩ Imm_LGE_acute:", length(cge_lge), "genes\n")
  cat("  Imm_CGE_acute ∩ Imm_MGE_acute:", length(cge_mge), "genes\n")
  cat("  Imm_LGE_acute ∩ Imm_MGE_acute:", length(lge_mge), "genes\n")
  cat("\n")
  
  cat("Three-way intersection (common to all three):\n")
  cat("  Imm_CGE_acute ∩ Imm_LGE_acute ∩ Imm_MGE_acute:", length(common_all_three), "genes\n")
  cat("\n")
  
  cat("Unique to each variable:\n")
  cat("  Only in Imm_CGE_acute:", length(only_cge), "genes\n")
  cat("  Only in Imm_LGE_acute:", length(only_lge), "genes\n")
  cat("  Only in Imm_MGE_acute:", length(only_mge), "genes\n")
  cat("\n")
  
  # Create summary data frame
  imm_intersection_summary <- data.frame(
    Analysis_Type = c("CGE_LGE", "CGE_MGE", "LGE_MGE", "All_Three", "Only_CGE", "Only_LGE", "Only_MGE"),
    Gene_Count = c(length(cge_lge), length(cge_mge), length(lge_mge), 
                   length(common_all_three), length(only_cge), length(only_lge), length(only_mge)),
    stringsAsFactors = FALSE
  )
  
  # Store results
  assign("imm_intersection_summary", imm_intersection_summary, envir = .GlobalEnv)
  assign("imm_common_all_three", common_all_three, envir = .GlobalEnv)
  assign("imm_cge_lge_intersection", cge_lge, envir = .GlobalEnv)
  assign("imm_cge_mge_intersection", cge_mge, envir = .GlobalEnv)
  assign("imm_lge_mge_intersection", lge_mge, envir = .GlobalEnv)
  assign("imm_only_cge", only_cge, envir = .GlobalEnv)
  assign("imm_only_lge", only_lge, envir = .GlobalEnv)
  assign("imm_only_mge", only_mge, envir = .GlobalEnv)
  
  cat("Results stored in variables:\n")
  cat("  - imm_intersection_summary: Summary table\n")
  cat("  - imm_common_all_three: Genes common to all three\n")
  cat("  - imm_cge_lge_intersection: CGE ∩ LGE\n")
  cat("  - imm_cge_mge_intersection: CGE ∩ MGE\n")
  cat("  - imm_lge_mge_intersection: LGE ∩ MGE\n")
  cat("  - imm_only_cge: Only in CGE\n")
  cat("  - imm_only_lge: Only in LGE\n")
  cat("  - imm_only_mge: Only in MGE\n")
  cat("\n")
  
  # Write to CSV
  imm_intersection_file <- file.path(script_dir, "imm_acute_significant_intersection.csv")
  write_csv(imm_intersection_summary, imm_intersection_file)
  cat("Summary written to file:", imm_intersection_file, "\n")
  
  # Write all intersection gene lists to separate files
  cat("\nWriting intersection gene lists to files...\n")
  
  # Three-way intersection (common to all three)
  if (length(common_all_three) > 0) {
    common_genes_df <- data.frame(Gene_name = common_all_three, stringsAsFactors = FALSE)
    common_genes_file <- file.path(script_dir, "imm_acute_common_all_three_genes.csv")
    write_csv(common_genes_df, common_genes_file)
    cat("  Common to all three:", common_genes_file, "(", length(common_all_three), "genes)\n")
  }
  
  # Pairwise intersections
  if (length(cge_lge) > 0) {
    write_csv(data.frame(Gene_name = cge_lge, stringsAsFactors = FALSE),
              file.path(script_dir, "imm_acute_CGE_LGE_intersection.csv"))
    cat("  CGE ∩ LGE:", length(cge_lge), "genes\n")
  }
  if (length(cge_mge) > 0) {
    write_csv(data.frame(Gene_name = cge_mge, stringsAsFactors = FALSE),
              file.path(script_dir, "imm_acute_CGE_MGE_intersection.csv"))
    cat("  CGE ∩ MGE:", length(cge_mge), "genes\n")
  }
  if (length(lge_mge) > 0) {
    write_csv(data.frame(Gene_name = lge_mge, stringsAsFactors = FALSE),
              file.path(script_dir, "imm_acute_LGE_MGE_intersection.csv"))
    cat("  LGE ∩ MGE:", length(lge_mge), "genes\n")
  }
  
  # Unique genes
  if (length(only_cge) > 0) {
    write_csv(data.frame(Gene_name = only_cge, stringsAsFactors = FALSE),
              file.path(script_dir, "imm_acute_only_CGE.csv"))
    cat("  Only CGE:", length(only_cge), "genes\n")
  }
  if (length(only_lge) > 0) {
    write_csv(data.frame(Gene_name = only_lge, stringsAsFactors = FALSE),
              file.path(script_dir, "imm_acute_only_LGE.csv"))
    cat("  Only LGE:", length(only_lge), "genes\n")
  }
  if (length(only_mge) > 0) {
    write_csv(data.frame(Gene_name = only_mge, stringsAsFactors = FALSE),
              file.path(script_dir, "imm_acute_only_MGE.csv"))
    cat("  Only MGE:", length(only_mge), "genes\n")
  }
  
  # Create Venn diagram
  cat("\nCreating Venn diagram...\n")
  venn_file <- file.path(script_dir, "imm_acute_venn_diagram.png")
  create_venn_diagram(imm_vars, 
                      title = "Imm Cell Types - Acute (Significant Genes)",
                      output_file = venn_file,
                      use_sig_only = TRUE)
  
} else {
  cat("Error: Could not find all three variables. Found", length(imm_sig_genes_list), "out of 3.\n")
}

# ===== INTERSECTION OF SIGNIFICANT GENES: Imm_CGE_chronic, Imm_LGE_chronic, Imm_MGE_chronic =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION OF SIGNIFICANT GENES: Imm_CGE_chronic, Imm_LGE_chronic, Imm_MGE_chronic\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Variable names to analyze
imm_vars_chronic <- c("Imm_CGE_chronic", "Imm_LGE_chronic", "Imm_MGE_chronic")

# Check if variables exist and extract significant genes
imm_sig_genes_list_chronic <- list()
imm_var_info_chronic <- data.frame(
  Variable = character(),
  Total_Genes = integer(),
  Significant_Genes = integer(),
  stringsAsFactors = FALSE
)

for (var_name in imm_vars_chronic) {
  if (exists(var_name, envir = .GlobalEnv)) {
    df <- get(var_name, envir = .GlobalEnv)
    
    # Filter for significant genes (Sig == 1)
    if ("Sig" %in% colnames(df)) {
      df_sig <- df %>% filter(Sig == 1)
      
      # Extract gene names
      if ("Gene_name" %in% colnames(df_sig)) {
        sig_genes <- unique(df_sig$Gene_name[!is.na(df_sig$Gene_name) & df_sig$Gene_name != ""])
        imm_sig_genes_list_chronic[[var_name]] <- sig_genes
        
        imm_var_info_chronic <- rbind(imm_var_info_chronic, data.frame(
          Variable = var_name,
          Total_Genes = nrow(df),
          Significant_Genes = length(sig_genes),
          stringsAsFactors = FALSE
        ))
        
        cat("Variable:", var_name, "\n")
        cat("  Total genes:", nrow(df), "\n")
        cat("  Significant genes (Sig=1):", length(sig_genes), "\n\n")
      } else {
        cat("Warning: Gene_name column not found in", var_name, "\n")
      }
    } else {
      cat("Warning: Sig column not found in", var_name, "\n")
    }
  } else {
    cat("Warning: Variable", var_name, "not found\n")
  }
}

# Find intersections
if (length(imm_sig_genes_list_chronic) == 3) {
  genes_cge_chronic <- imm_sig_genes_list_chronic[["Imm_CGE_chronic"]]
  genes_lge_chronic <- imm_sig_genes_list_chronic[["Imm_LGE_chronic"]]
  genes_mge_chronic <- imm_sig_genes_list_chronic[["Imm_MGE_chronic"]]
  
  # Pairwise intersections
  cge_lge_chronic <- intersect(genes_cge_chronic, genes_lge_chronic)
  cge_mge_chronic <- intersect(genes_cge_chronic, genes_mge_chronic)
  lge_mge_chronic <- intersect(genes_lge_chronic, genes_mge_chronic)
  
  # Three-way intersection (common to all three)
  common_all_three_chronic <- intersect(intersect(genes_cge_chronic, genes_lge_chronic), genes_mge_chronic)
  
  # Unique to each
  only_cge_chronic <- setdiff(setdiff(genes_cge_chronic, genes_lge_chronic), genes_mge_chronic)
  only_lge_chronic <- setdiff(setdiff(genes_lge_chronic, genes_cge_chronic), genes_mge_chronic)
  only_mge_chronic <- setdiff(setdiff(genes_mge_chronic, genes_cge_chronic), genes_lge_chronic)
  
  # Display results
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("INTERSECTION RESULTS:\n")
  cat(paste(rep("=", 60), collapse=""), "\n\n")
  
  cat("Significant genes in each variable:\n")
  print(imm_var_info_chronic)
  cat("\n")
  
  cat("Pairwise intersections:\n")
  cat("  Imm_CGE_chronic ∩ Imm_LGE_chronic:", length(cge_lge_chronic), "genes\n")
  cat("  Imm_CGE_chronic ∩ Imm_MGE_chronic:", length(cge_mge_chronic), "genes\n")
  cat("  Imm_LGE_chronic ∩ Imm_MGE_chronic:", length(lge_mge_chronic), "genes\n")
  cat("\n")
  
  cat("Three-way intersection (common to all three):\n")
  cat("  Imm_CGE_chronic ∩ Imm_LGE_chronic ∩ Imm_MGE_chronic:", length(common_all_three_chronic), "genes\n")
  cat("\n")
  
  cat("Unique to each variable:\n")
  cat("  Only in Imm_CGE_chronic:", length(only_cge_chronic), "genes\n")
  cat("  Only in Imm_LGE_chronic:", length(only_lge_chronic), "genes\n")
  cat("  Only in Imm_MGE_chronic:", length(only_mge_chronic), "genes\n")
  cat("\n")
  
  # Create summary data frame
  imm_intersection_summary_chronic <- data.frame(
    Analysis_Type = c("CGE_LGE", "CGE_MGE", "LGE_MGE", "All_Three", "Only_CGE", "Only_LGE", "Only_MGE"),
    Gene_Count = c(length(cge_lge_chronic), length(cge_mge_chronic), length(lge_mge_chronic), 
                   length(common_all_three_chronic), length(only_cge_chronic), length(only_lge_chronic), length(only_mge_chronic)),
    stringsAsFactors = FALSE
  )
  
  # Store results
  assign("imm_intersection_summary_chronic", imm_intersection_summary_chronic, envir = .GlobalEnv)
  assign("imm_common_all_three_chronic", common_all_three_chronic, envir = .GlobalEnv)
  assign("imm_cge_lge_intersection_chronic", cge_lge_chronic, envir = .GlobalEnv)
  assign("imm_cge_mge_intersection_chronic", cge_mge_chronic, envir = .GlobalEnv)
  assign("imm_lge_mge_intersection_chronic", lge_mge_chronic, envir = .GlobalEnv)
  assign("imm_only_cge_chronic", only_cge_chronic, envir = .GlobalEnv)
  assign("imm_only_lge_chronic", only_lge_chronic, envir = .GlobalEnv)
  assign("imm_only_mge_chronic", only_mge_chronic, envir = .GlobalEnv)
  
  cat("Results stored in variables:\n")
  cat("  - imm_intersection_summary_chronic: Summary table\n")
  cat("  - imm_common_all_three_chronic: Genes common to all three\n")
  cat("  - imm_cge_lge_intersection_chronic: CGE ∩ LGE\n")
  cat("  - imm_cge_mge_intersection_chronic: CGE ∩ MGE\n")
  cat("  - imm_lge_mge_intersection_chronic: LGE ∩ MGE\n")
  cat("  - imm_only_cge_chronic: Only in CGE\n")
  cat("  - imm_only_lge_chronic: Only in LGE\n")
  cat("  - imm_only_mge_chronic: Only in MGE\n")
  cat("\n")
  
  # Write to CSV
  imm_intersection_file_chronic <- file.path(script_dir, "imm_chronic_significant_intersection.csv")
  write_csv(imm_intersection_summary_chronic, imm_intersection_file_chronic)
  cat("Summary written to file:", imm_intersection_file_chronic, "\n")
  
  # Write all intersection gene lists to separate files
  cat("\nWriting intersection gene lists to files...\n")
  
  # Three-way intersection (common to all three)
  if (length(common_all_three_chronic) > 0) {
    common_genes_df_chronic <- data.frame(Gene_name = common_all_three_chronic, stringsAsFactors = FALSE)
    common_genes_file_chronic <- file.path(script_dir, "imm_chronic_common_all_three_genes.csv")
    write_csv(common_genes_df_chronic, common_genes_file_chronic)
    cat("  Common to all three:", common_genes_file_chronic, "(", length(common_all_three_chronic), "genes)\n")
  }
  
  # Pairwise intersections
  if (length(cge_lge_chronic) > 0) {
    write_csv(data.frame(Gene_name = cge_lge_chronic, stringsAsFactors = FALSE),
              file.path(script_dir, "imm_chronic_CGE_LGE_intersection.csv"))
    cat("  CGE ∩ LGE:", length(cge_lge_chronic), "genes\n")
  }
  if (length(cge_mge_chronic) > 0) {
    write_csv(data.frame(Gene_name = cge_mge_chronic, stringsAsFactors = FALSE),
              file.path(script_dir, "imm_chronic_CGE_MGE_intersection.csv"))
    cat("  CGE ∩ MGE:", length(cge_mge_chronic), "genes\n")
  }
  if (length(lge_mge_chronic) > 0) {
    write_csv(data.frame(Gene_name = lge_mge_chronic, stringsAsFactors = FALSE),
              file.path(script_dir, "imm_chronic_LGE_MGE_intersection.csv"))
    cat("  LGE ∩ MGE:", length(lge_mge_chronic), "genes\n")
  }
  
  # Unique genes
  if (length(only_cge_chronic) > 0) {
    write_csv(data.frame(Gene_name = only_cge_chronic, stringsAsFactors = FALSE),
              file.path(script_dir, "imm_chronic_only_CGE.csv"))
    cat("  Only CGE:", length(only_cge_chronic), "genes\n")
  }
  if (length(only_lge_chronic) > 0) {
    write_csv(data.frame(Gene_name = only_lge_chronic, stringsAsFactors = FALSE),
              file.path(script_dir, "imm_chronic_only_LGE.csv"))
    cat("  Only LGE:", length(only_lge_chronic), "genes\n")
  }
  if (length(only_mge_chronic) > 0) {
    write_csv(data.frame(Gene_name = only_mge_chronic, stringsAsFactors = FALSE),
              file.path(script_dir, "imm_chronic_only_MGE.csv"))
    cat("  Only MGE:", length(only_mge_chronic), "genes\n")
  }
  
  # Create Venn diagram
  cat("\nCreating Venn diagram...\n")
  venn_file <- file.path(script_dir, "imm_chronic_venn_diagram.png")
  create_venn_diagram(imm_vars_chronic, 
                      title = "Imm Cell Types - Chronic (Significant Genes)",
                      output_file = venn_file,
                      use_sig_only = TRUE)
  
} else {
  cat("Error: Could not find all three variables. Found", length(imm_sig_genes_list_chronic), "out of 3.\n")
}

# ===== INTERSECTION OF SIGNIFICANT GENES: OPC_acute, Premyelinating_Olig_acute, Myelinating_Olig_acute =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION OF SIGNIFICANT GENES: OPC_acute, Premyelinating_Olig_acute, Myelinating_Olig_acute\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Variable names to analyze
olig_vars_acute <- c("OPC_acute", "Premyelinating_Olig_acute", "Myelinating_Olig_acute")

# Check if variables exist and extract significant genes
olig_sig_genes_list_acute <- list()
olig_var_info_acute <- data.frame(
  Variable = character(),
  Total_Genes = integer(),
  Significant_Genes = integer(),
  stringsAsFactors = FALSE
)

for (var_name in olig_vars_acute) {
  if (exists(var_name, envir = .GlobalEnv)) {
    df <- get(var_name, envir = .GlobalEnv)
    
    # Filter for significant genes (Sig == 1)
    if ("Sig" %in% colnames(df)) {
      df_sig <- df %>% filter(Sig == 1)
      
      # Extract gene names
      if ("Gene_name" %in% colnames(df_sig)) {
        sig_genes <- unique(df_sig$Gene_name[!is.na(df_sig$Gene_name) & df_sig$Gene_name != ""])
        olig_sig_genes_list_acute[[var_name]] <- sig_genes
        
        olig_var_info_acute <- rbind(olig_var_info_acute, data.frame(
          Variable = var_name,
          Total_Genes = nrow(df),
          Significant_Genes = length(sig_genes),
          stringsAsFactors = FALSE
        ))
        
        cat("Variable:", var_name, "\n")
        cat("  Total genes:", nrow(df), "\n")
        cat("  Significant genes (Sig=1):", length(sig_genes), "\n\n")
      } else {
        cat("Warning: Gene_name column not found in", var_name, "\n")
      }
    } else {
      cat("Warning: Sig column not found in", var_name, "\n")
    }
  } else {
    cat("Warning: Variable", var_name, "not found\n")
  }
}

# Find intersections
if (length(olig_sig_genes_list_acute) == 3) {
  genes_opc_acute <- olig_sig_genes_list_acute[["OPC_acute"]]
  genes_premyel_acute <- olig_sig_genes_list_acute[["Premyelinating_Olig_acute"]]
  genes_myel_acute <- olig_sig_genes_list_acute[["Myelinating_Olig_acute"]]
  
  # Pairwise intersections
  opc_premyel_acute <- intersect(genes_opc_acute, genes_premyel_acute)
  opc_myel_acute <- intersect(genes_opc_acute, genes_myel_acute)
  premyel_myel_acute <- intersect(genes_premyel_acute, genes_myel_acute)
  
  # Three-way intersection (common to all three)
  common_all_three_olig_acute <- intersect(intersect(genes_opc_acute, genes_premyel_acute), genes_myel_acute)
  
  # Unique to each
  only_opc_acute <- setdiff(setdiff(genes_opc_acute, genes_premyel_acute), genes_myel_acute)
  only_premyel_acute <- setdiff(setdiff(genes_premyel_acute, genes_opc_acute), genes_myel_acute)
  only_myel_acute <- setdiff(setdiff(genes_myel_acute, genes_opc_acute), genes_premyel_acute)
  
  # Display results
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("INTERSECTION RESULTS:\n")
  cat(paste(rep("=", 60), collapse=""), "\n\n")
  
  cat("Significant genes in each variable:\n")
  print(olig_var_info_acute)
  cat("\n")
  
  cat("Pairwise intersections:\n")
  cat("  OPC_acute ∩ Premyelinating_Olig_acute:", length(opc_premyel_acute), "genes\n")
  cat("  OPC_acute ∩ Myelinating_Olig_acute:", length(opc_myel_acute), "genes\n")
  cat("  Premyelinating_Olig_acute ∩ Myelinating_Olig_acute:", length(premyel_myel_acute), "genes\n")
  cat("\n")
  
  cat("Three-way intersection (common to all three):\n")
  cat("  OPC_acute ∩ Premyelinating_Olig_acute ∩ Myelinating_Olig_acute:", length(common_all_three_olig_acute), "genes\n")
  cat("\n")
  
  cat("Unique to each variable:\n")
  cat("  Only in OPC_acute:", length(only_opc_acute), "genes\n")
  cat("  Only in Premyelinating_Olig_acute:", length(only_premyel_acute), "genes\n")
  cat("  Only in Myelinating_Olig_acute:", length(only_myel_acute), "genes\n")
  cat("\n")
  
  # Create summary data frame
  olig_intersection_summary_acute <- data.frame(
    Analysis_Type = c("OPC_Premyel", "OPC_Myel", "Premyel_Myel", "All_Three", "Only_OPC", "Only_Premyel", "Only_Myel"),
    Gene_Count = c(length(opc_premyel_acute), length(opc_myel_acute), length(premyel_myel_acute), 
                   length(common_all_three_olig_acute), length(only_opc_acute), length(only_premyel_acute), length(only_myel_acute)),
    stringsAsFactors = FALSE
  )
  
  # Store results
  assign("olig_intersection_summary_acute", olig_intersection_summary_acute, envir = .GlobalEnv)
  assign("olig_common_all_three_acute", common_all_three_olig_acute, envir = .GlobalEnv)
  assign("olig_opc_premyel_intersection_acute", opc_premyel_acute, envir = .GlobalEnv)
  assign("olig_opc_myel_intersection_acute", opc_myel_acute, envir = .GlobalEnv)
  assign("olig_premyel_myel_intersection_acute", premyel_myel_acute, envir = .GlobalEnv)
  assign("olig_only_opc_acute", only_opc_acute, envir = .GlobalEnv)
  assign("olig_only_premyel_acute", only_premyel_acute, envir = .GlobalEnv)
  assign("olig_only_myel_acute", only_myel_acute, envir = .GlobalEnv)
  
  cat("Results stored in variables:\n")
  cat("  - olig_intersection_summary_acute: Summary table\n")
  cat("  - olig_common_all_three_acute: Genes common to all three\n")
  cat("  - olig_opc_premyel_intersection_acute: OPC ∩ Premyel\n")
  cat("  - olig_opc_myel_intersection_acute: OPC ∩ Myel\n")
  cat("  - olig_premyel_myel_intersection_acute: Premyel ∩ Myel\n")
  cat("  - olig_only_opc_acute: Only in OPC\n")
  cat("  - olig_only_premyel_acute: Only in Premyel\n")
  cat("  - olig_only_myel_acute: Only in Myel\n")
  cat("\n")
  
  # Write to CSV
  olig_intersection_file_acute <- file.path(script_dir, "olig_acute_significant_intersection.csv")
  write_csv(olig_intersection_summary_acute, olig_intersection_file_acute)
  cat("Summary written to file:", olig_intersection_file_acute, "\n")
  
  # Write all intersection gene lists to separate files
  cat("\nWriting intersection gene lists to files...\n")
  
  # Three-way intersection (common to all three)
  if (length(common_all_three_olig_acute) > 0) {
    common_genes_df_olig_acute <- data.frame(Gene_name = common_all_three_olig_acute, stringsAsFactors = FALSE)
    common_genes_file_olig_acute <- file.path(script_dir, "olig_acute_common_all_three_genes.csv")
    write_csv(common_genes_df_olig_acute, common_genes_file_olig_acute)
    cat("  Common to all three:", common_genes_file_olig_acute, "(", length(common_all_three_olig_acute), "genes)\n")
  }
  
  # Pairwise intersections
  if (length(opc_premyel_acute) > 0) {
    write_csv(data.frame(Gene_name = opc_premyel_acute, stringsAsFactors = FALSE),
              file.path(script_dir, "olig_acute_OPC_Premyel_intersection.csv"))
    cat("  OPC ∩ Premyel:", length(opc_premyel_acute), "genes\n")
  }
  if (length(opc_myel_acute) > 0) {
    write_csv(data.frame(Gene_name = opc_myel_acute, stringsAsFactors = FALSE),
              file.path(script_dir, "olig_acute_OPC_Myel_intersection.csv"))
    cat("  OPC ∩ Myel:", length(opc_myel_acute), "genes\n")
  }
  if (length(premyel_myel_acute) > 0) {
    write_csv(data.frame(Gene_name = premyel_myel_acute, stringsAsFactors = FALSE),
              file.path(script_dir, "olig_acute_Premyel_Myel_intersection.csv"))
    cat("  Premyel ∩ Myel:", length(premyel_myel_acute), "genes\n")
  }
  
  # Unique genes
  if (length(only_opc_acute) > 0) {
    write_csv(data.frame(Gene_name = only_opc_acute, stringsAsFactors = FALSE),
              file.path(script_dir, "olig_acute_only_OPC.csv"))
    cat("  Only OPC:", length(only_opc_acute), "genes\n")
  }
  if (length(only_premyel_acute) > 0) {
    write_csv(data.frame(Gene_name = only_premyel_acute, stringsAsFactors = FALSE),
              file.path(script_dir, "olig_acute_only_Premyel.csv"))
    cat("  Only Premyel:", length(only_premyel_acute), "genes\n")
  }
  if (length(only_myel_acute) > 0) {
    write_csv(data.frame(Gene_name = only_myel_acute, stringsAsFactors = FALSE),
              file.path(script_dir, "olig_acute_only_Myel.csv"))
    cat("  Only Myel:", length(only_myel_acute), "genes\n")
  }
  
  # Create Venn diagram
  cat("\nCreating Venn diagram...\n")
  venn_file <- file.path(script_dir, "olig_acute_venn_diagram.png")
  create_venn_diagram(olig_vars_acute, 
                      title = "Oligodendrocyte Cell Types - Acute (Significant Genes)",
                      output_file = venn_file,
                      use_sig_only = TRUE)
  
} else {
  cat("Error: Could not find all three variables. Found", length(olig_sig_genes_list_acute), "out of 3.\n")
}

# ===== INTERSECTION OF SIGNIFICANT GENES: OPC_chronic, Premyelinating_Olig_chronic, Myelinating_Olig_chronic =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION OF SIGNIFICANT GENES: OPC_chronic, Premyelinating_Olig_chronic, Myelinating_Olig_chronic\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Variable names to analyze
olig_vars_chronic <- c("OPC_chronic", "Premyelinating_Olig_chronic", "Myelinating_Olig_chronic")

# Check if variables exist and extract significant genes
olig_sig_genes_list_chronic <- list()
olig_var_info_chronic <- data.frame(
  Variable = character(),
  Total_Genes = integer(),
  Significant_Genes = integer(),
  stringsAsFactors = FALSE
)

for (var_name in olig_vars_chronic) {
  if (exists(var_name, envir = .GlobalEnv)) {
    df <- get(var_name, envir = .GlobalEnv)
    
    # Filter for significant genes (Sig == 1)
    if ("Sig" %in% colnames(df)) {
      df_sig <- df %>% filter(Sig == 1)
      
      # Extract gene names
      if ("Gene_name" %in% colnames(df_sig)) {
        sig_genes <- unique(df_sig$Gene_name[!is.na(df_sig$Gene_name) & df_sig$Gene_name != ""])
        olig_sig_genes_list_chronic[[var_name]] <- sig_genes
        
        olig_var_info_chronic <- rbind(olig_var_info_chronic, data.frame(
          Variable = var_name,
          Total_Genes = nrow(df),
          Significant_Genes = length(sig_genes),
          stringsAsFactors = FALSE
        ))
        
        cat("Variable:", var_name, "\n")
        cat("  Total genes:", nrow(df), "\n")
        cat("  Significant genes (Sig=1):", length(sig_genes), "\n\n")
      } else {
        cat("Warning: Gene_name column not found in", var_name, "\n")
      }
    } else {
      cat("Warning: Sig column not found in", var_name, "\n")
    }
  } else {
    cat("Warning: Variable", var_name, "not found\n")
  }
}

# Find intersections
if (length(olig_sig_genes_list_chronic) == 3) {
  genes_opc_chronic <- olig_sig_genes_list_chronic[["OPC_chronic"]]
  genes_premyel_chronic <- olig_sig_genes_list_chronic[["Premyelinating_Olig_chronic"]]
  genes_myel_chronic <- olig_sig_genes_list_chronic[["Myelinating_Olig_chronic"]]
  
  # Pairwise intersections
  opc_premyel_chronic <- intersect(genes_opc_chronic, genes_premyel_chronic)
  opc_myel_chronic <- intersect(genes_opc_chronic, genes_myel_chronic)
  premyel_myel_chronic <- intersect(genes_premyel_chronic, genes_myel_chronic)
  
  # Three-way intersection (common to all three)
  common_all_three_olig_chronic <- intersect(intersect(genes_opc_chronic, genes_premyel_chronic), genes_myel_chronic)
  
  # Unique to each
  only_opc_chronic <- setdiff(setdiff(genes_opc_chronic, genes_premyel_chronic), genes_myel_chronic)
  only_premyel_chronic <- setdiff(setdiff(genes_premyel_chronic, genes_opc_chronic), genes_myel_chronic)
  only_myel_chronic <- setdiff(setdiff(genes_myel_chronic, genes_opc_chronic), genes_premyel_chronic)
  
  # Display results
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("INTERSECTION RESULTS:\n")
  cat(paste(rep("=", 60), collapse=""), "\n\n")
  
  cat("Significant genes in each variable:\n")
  print(olig_var_info_chronic)
  cat("\n")
  
  cat("Pairwise intersections:\n")
  cat("  OPC_chronic ∩ Premyelinating_Olig_chronic:", length(opc_premyel_chronic), "genes\n")
  cat("  OPC_chronic ∩ Myelinating_Olig_chronic:", length(opc_myel_chronic), "genes\n")
  cat("  Premyelinating_Olig_chronic ∩ Myelinating_Olig_chronic:", length(premyel_myel_chronic), "genes\n")
  cat("\n")
  
  cat("Three-way intersection (common to all three):\n")
  cat("  OPC_chronic ∩ Premyelinating_Olig_chronic ∩ Myelinating_Olig_chronic:", length(common_all_three_olig_chronic), "genes\n")
  cat("\n")
  
  cat("Unique to each variable:\n")
  cat("  Only in OPC_chronic:", length(only_opc_chronic), "genes\n")
  cat("  Only in Premyelinating_Olig_chronic:", length(only_premyel_chronic), "genes\n")
  cat("  Only in Myelinating_Olig_chronic:", length(only_myel_chronic), "genes\n")
  cat("\n")
  
  # Create summary data frame
  olig_intersection_summary_chronic <- data.frame(
    Analysis_Type = c("OPC_Premyel", "OPC_Myel", "Premyel_Myel", "All_Three", "Only_OPC", "Only_Premyel", "Only_Myel"),
    Gene_Count = c(length(opc_premyel_chronic), length(opc_myel_chronic), length(premyel_myel_chronic), 
                   length(common_all_three_olig_chronic), length(only_opc_chronic), length(only_premyel_chronic), length(only_myel_chronic)),
    stringsAsFactors = FALSE
  )
  
  # Store results
  assign("olig_intersection_summary_chronic", olig_intersection_summary_chronic, envir = .GlobalEnv)
  assign("olig_common_all_three_chronic", common_all_three_olig_chronic, envir = .GlobalEnv)
  assign("olig_opc_premyel_intersection_chronic", opc_premyel_chronic, envir = .GlobalEnv)
  assign("olig_opc_myel_intersection_chronic", opc_myel_chronic, envir = .GlobalEnv)
  assign("olig_premyel_myel_intersection_chronic", premyel_myel_chronic, envir = .GlobalEnv)
  assign("olig_only_opc_chronic", only_opc_chronic, envir = .GlobalEnv)
  assign("olig_only_premyel_chronic", only_premyel_chronic, envir = .GlobalEnv)
  assign("olig_only_myel_chronic", only_myel_chronic, envir = .GlobalEnv)
  
  cat("Results stored in variables:\n")
  cat("  - olig_intersection_summary_chronic: Summary table\n")
  cat("  - olig_common_all_three_chronic: Genes common to all three\n")
  cat("  - olig_opc_premyel_intersection_chronic: OPC ∩ Premyel\n")
  cat("  - olig_opc_myel_intersection_chronic: OPC ∩ Myel\n")
  cat("  - olig_premyel_myel_intersection_chronic: Premyel ∩ Myel\n")
  cat("  - olig_only_opc_chronic: Only in OPC\n")
  cat("  - olig_only_premyel_chronic: Only in Premyel\n")
  cat("  - olig_only_myel_chronic: Only in Myel\n")
  cat("\n")
  
  # Write to CSV
  olig_intersection_file_chronic <- file.path(script_dir, "olig_chronic_significant_intersection.csv")
  write_csv(olig_intersection_summary_chronic, olig_intersection_file_chronic)
  cat("Summary written to file:", olig_intersection_file_chronic, "\n")
  
  # Write all intersection gene lists to separate files
  cat("\nWriting intersection gene lists to files...\n")
  
  # Three-way intersection (common to all three)
  if (length(common_all_three_olig_chronic) > 0) {
    common_genes_df_olig_chronic <- data.frame(Gene_name = common_all_three_olig_chronic, stringsAsFactors = FALSE)
    common_genes_file_olig_chronic <- file.path(script_dir, "olig_chronic_common_all_three_genes.csv")
    write_csv(common_genes_df_olig_chronic, common_genes_file_olig_chronic)
    cat("  Common to all three:", common_genes_file_olig_chronic, "(", length(common_all_three_olig_chronic), "genes)\n")
  }
  
  # Pairwise intersections
  if (length(opc_premyel_chronic) > 0) {
    write_csv(data.frame(Gene_name = opc_premyel_chronic, stringsAsFactors = FALSE),
              file.path(script_dir, "olig_chronic_OPC_Premyel_intersection.csv"))
    cat("  OPC ∩ Premyel:", length(opc_premyel_chronic), "genes\n")
  }
  if (length(opc_myel_chronic) > 0) {
    write_csv(data.frame(Gene_name = opc_myel_chronic, stringsAsFactors = FALSE),
              file.path(script_dir, "olig_chronic_OPC_Myel_intersection.csv"))
    cat("  OPC ∩ Myel:", length(opc_myel_chronic), "genes\n")
  }
  if (length(premyel_myel_chronic) > 0) {
    write_csv(data.frame(Gene_name = premyel_myel_chronic, stringsAsFactors = FALSE),
              file.path(script_dir, "olig_chronic_Premyel_Myel_intersection.csv"))
    cat("  Premyel ∩ Myel:", length(premyel_myel_chronic), "genes\n")
  }
  
  # Unique genes
  if (length(only_opc_chronic) > 0) {
    write_csv(data.frame(Gene_name = only_opc_chronic, stringsAsFactors = FALSE),
              file.path(script_dir, "olig_chronic_only_OPC.csv"))
    cat("  Only OPC:", length(only_opc_chronic), "genes\n")
  }
  if (length(only_premyel_chronic) > 0) {
    write_csv(data.frame(Gene_name = only_premyel_chronic, stringsAsFactors = FALSE),
              file.path(script_dir, "olig_chronic_only_Premyel.csv"))
    cat("  Only Premyel:", length(only_premyel_chronic), "genes\n")
  }
  if (length(only_myel_chronic) > 0) {
    write_csv(data.frame(Gene_name = only_myel_chronic, stringsAsFactors = FALSE),
              file.path(script_dir, "olig_chronic_only_Myel.csv"))
    cat("  Only Myel:", length(only_myel_chronic), "genes\n")
  }
  
  # Create Venn diagram
  cat("\nCreating Venn diagram...\n")
  venn_file <- file.path(script_dir, "olig_chronic_venn_diagram.png")
  create_venn_diagram(olig_vars_chronic, 
                      title = "Oligodendrocyte Cell Types - Chronic (Significant Genes)",
                      output_file = venn_file,
                      use_sig_only = TRUE)
  
} else {
  cat("Error: Could not find all three variables. Found", length(olig_sig_genes_list_chronic), "out of 3.\n")
}

# ===== INTERSECTION OF SIGNIFICANT GENES: Astroglia_acute, Microglia_acute =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION OF SIGNIFICANT GENES: Astroglia_acute, Microglia_acute\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Variable names to analyze
glia_vars_acute <- c("Astroglia_acute", "Microglia_acute")

# Check if variables exist and extract significant genes
glia_sig_genes_list_acute <- list()
glia_var_info_acute <- data.frame(
  Variable = character(),
  Total_Genes = integer(),
  Significant_Genes = integer(),
  stringsAsFactors = FALSE
)

for (var_name in glia_vars_acute) {
  if (exists(var_name, envir = .GlobalEnv)) {
    df <- get(var_name, envir = .GlobalEnv)
    
    # Filter for significant genes (Sig == 1)
    if ("Sig" %in% colnames(df)) {
      df_sig <- df %>% filter(Sig == 1)
      
      # Extract gene names
      if ("Gene_name" %in% colnames(df_sig)) {
        sig_genes <- unique(df_sig$Gene_name[!is.na(df_sig$Gene_name) & df_sig$Gene_name != ""])
        glia_sig_genes_list_acute[[var_name]] <- sig_genes
        
        glia_var_info_acute <- rbind(glia_var_info_acute, data.frame(
          Variable = var_name,
          Total_Genes = nrow(df),
          Significant_Genes = length(sig_genes),
          stringsAsFactors = FALSE
        ))
        
        cat("Variable:", var_name, "\n")
        cat("  Total genes:", nrow(df), "\n")
        cat("  Significant genes (Sig=1):", length(sig_genes), "\n\n")
      } else {
        cat("Warning: Gene_name column not found in", var_name, "\n")
      }
    } else {
      cat("Warning: Sig column not found in", var_name, "\n")
    }
  } else {
    cat("Warning: Variable", var_name, "not found\n")
  }
}

# Find intersections
if (length(glia_sig_genes_list_acute) == 2) {
  genes_astro_acute <- glia_sig_genes_list_acute[["Astroglia_acute"]]
  genes_micro_acute <- glia_sig_genes_list_acute[["Microglia_acute"]]
  
  # Pairwise intersection (common to both)
  common_both_glia_acute <- intersect(genes_astro_acute, genes_micro_acute)
  
  # Unique to each
  only_astro_acute <- setdiff(genes_astro_acute, genes_micro_acute)
  only_micro_acute <- setdiff(genes_micro_acute, genes_astro_acute)
  
  # Display results
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("INTERSECTION RESULTS:\n")
  cat(paste(rep("=", 60), collapse=""), "\n\n")
  
  cat("Significant genes in each variable:\n")
  print(glia_var_info_acute)
  cat("\n")
  
  cat("Intersection (common to both):\n")
  cat("  Astroglia_acute ∩ Microglia_acute:", length(common_both_glia_acute), "genes\n")
  cat("\n")
  
  cat("Unique to each variable:\n")
  cat("  Only in Astroglia_acute:", length(only_astro_acute), "genes\n")
  cat("  Only in Microglia_acute:", length(only_micro_acute), "genes\n")
  cat("\n")
  
  # Create summary data frame
  glia_intersection_summary_acute <- data.frame(
    Analysis_Type = c("Common_Both", "Only_Astroglia", "Only_Microglia"),
    Gene_Count = c(length(common_both_glia_acute), length(only_astro_acute), length(only_micro_acute)),
    stringsAsFactors = FALSE
  )
  
  # Store results
  assign("glia_intersection_summary_acute", glia_intersection_summary_acute, envir = .GlobalEnv)
  assign("glia_common_both_acute", common_both_glia_acute, envir = .GlobalEnv)
  assign("glia_only_astro_acute", only_astro_acute, envir = .GlobalEnv)
  assign("glia_only_micro_acute", only_micro_acute, envir = .GlobalEnv)
  
  cat("Results stored in variables:\n")
  cat("  - glia_intersection_summary_acute: Summary table\n")
  cat("  - glia_common_both_acute: Genes common to both\n")
  cat("  - glia_only_astro_acute: Only in Astroglia\n")
  cat("  - glia_only_micro_acute: Only in Microglia\n")
  cat("\n")
  
  # Write to CSV
  glia_intersection_file_acute <- file.path(script_dir, "glia_acute_significant_intersection.csv")
  write_csv(glia_intersection_summary_acute, glia_intersection_file_acute)
  cat("Summary written to file:", glia_intersection_file_acute, "\n")
  
  # Write all intersection gene lists to separate files
  cat("\nWriting intersection gene lists to files...\n")
  
  # Intersection (common to both)
  if (length(common_both_glia_acute) > 0) {
    common_genes_df_glia_acute <- data.frame(Gene_name = common_both_glia_acute, stringsAsFactors = FALSE)
    common_genes_file_glia_acute <- file.path(script_dir, "glia_acute_common_both_genes.csv")
    write_csv(common_genes_df_glia_acute, common_genes_file_glia_acute)
    cat("  Common to both:", common_genes_file_glia_acute, "(", length(common_both_glia_acute), "genes)\n")
  }
  
  # Unique genes
  if (length(only_astro_acute) > 0) {
    write_csv(data.frame(Gene_name = only_astro_acute, stringsAsFactors = FALSE),
              file.path(script_dir, "glia_acute_only_Astroglia.csv"))
    cat("  Only Astroglia:", length(only_astro_acute), "genes\n")
  }
  if (length(only_micro_acute) > 0) {
    write_csv(data.frame(Gene_name = only_micro_acute, stringsAsFactors = FALSE),
              file.path(script_dir, "glia_acute_only_Microglia.csv"))
    cat("  Only Microglia:", length(only_micro_acute), "genes\n")
  }
  
  # Create Venn diagram
  cat("\nCreating Venn diagram...\n")
  venn_file <- file.path(script_dir, "glia_acute_venn_diagram.png")
  create_venn_diagram(glia_vars_acute, 
                      title = "Glia Cell Types - Acute (Significant Genes)",
                      output_file = venn_file,
                      use_sig_only = TRUE)
  
} else {
  cat("Error: Could not find both variables. Found", length(glia_sig_genes_list_acute), "out of 2.\n")
}

# ===== INTERSECTION OF SIGNIFICANT GENES: Astroglia_chronic, Microglia_chronic =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("INTERSECTION OF SIGNIFICANT GENES: Astroglia_chronic, Microglia_chronic\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Variable names to analyze
glia_vars_chronic <- c("Astroglia_chronic", "Microglia_chronic")

# Check if variables exist and extract significant genes
glia_sig_genes_list_chronic <- list()
glia_var_info_chronic <- data.frame(
  Variable = character(),
  Total_Genes = integer(),
  Significant_Genes = integer(),
  stringsAsFactors = FALSE
)

for (var_name in glia_vars_chronic) {
  if (exists(var_name, envir = .GlobalEnv)) {
    df <- get(var_name, envir = .GlobalEnv)
    
    # Filter for significant genes (Sig == 1)
    if ("Sig" %in% colnames(df)) {
      df_sig <- df %>% filter(Sig == 1)
      
      # Extract gene names
      if ("Gene_name" %in% colnames(df_sig)) {
        sig_genes <- unique(df_sig$Gene_name[!is.na(df_sig$Gene_name) & df_sig$Gene_name != ""])
        glia_sig_genes_list_chronic[[var_name]] <- sig_genes
        
        glia_var_info_chronic <- rbind(glia_var_info_chronic, data.frame(
          Variable = var_name,
          Total_Genes = nrow(df),
          Significant_Genes = length(sig_genes),
          stringsAsFactors = FALSE
        ))
        
        cat("Variable:", var_name, "\n")
        cat("  Total genes:", nrow(df), "\n")
        cat("  Significant genes (Sig=1):", length(sig_genes), "\n\n")
      } else {
        cat("Warning: Gene_name column not found in", var_name, "\n")
      }
    } else {
      cat("Warning: Sig column not found in", var_name, "\n")
    }
  } else {
    cat("Warning: Variable", var_name, "not found\n")
  }
}

# Find intersections
if (length(glia_sig_genes_list_chronic) == 2) {
  genes_astro_chronic <- glia_sig_genes_list_chronic[["Astroglia_chronic"]]
  genes_micro_chronic <- glia_sig_genes_list_chronic[["Microglia_chronic"]]
  
  # Pairwise intersection (common to both)
  common_both_glia_chronic <- intersect(genes_astro_chronic, genes_micro_chronic)
  
  # Unique to each
  only_astro_chronic <- setdiff(genes_astro_chronic, genes_micro_chronic)
  only_micro_chronic <- setdiff(genes_micro_chronic, genes_astro_chronic)
  
  # Display results
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("INTERSECTION RESULTS:\n")
  cat(paste(rep("=", 60), collapse=""), "\n\n")
  
  cat("Significant genes in each variable:\n")
  print(glia_var_info_chronic)
  cat("\n")
  
  cat("Intersection (common to both):\n")
  cat("  Astroglia_chronic ∩ Microglia_chronic:", length(common_both_glia_chronic), "genes\n")
  cat("\n")
  
  cat("Unique to each variable:\n")
  cat("  Only in Astroglia_chronic:", length(only_astro_chronic), "genes\n")
  cat("  Only in Microglia_chronic:", length(only_micro_chronic), "genes\n")
  cat("\n")
  
  # Create summary data frame
  glia_intersection_summary_chronic <- data.frame(
    Analysis_Type = c("Common_Both", "Only_Astroglia", "Only_Microglia"),
    Gene_Count = c(length(common_both_glia_chronic), length(only_astro_chronic), length(only_micro_chronic)),
    stringsAsFactors = FALSE
  )
  
  # Store results
  assign("glia_intersection_summary_chronic", glia_intersection_summary_chronic, envir = .GlobalEnv)
  assign("glia_common_both_chronic", common_both_glia_chronic, envir = .GlobalEnv)
  assign("glia_only_astro_chronic", only_astro_chronic, envir = .GlobalEnv)
  assign("glia_only_micro_chronic", only_micro_chronic, envir = .GlobalEnv)
  
  cat("Results stored in variables:\n")
  cat("  - glia_intersection_summary_chronic: Summary table\n")
  cat("  - glia_common_both_chronic: Genes common to both\n")
  cat("  - glia_only_astro_chronic: Only in Astroglia\n")
  cat("  - glia_only_micro_chronic: Only in Microglia\n")
  cat("\n")
  
  # Write to CSV
  glia_intersection_file_chronic <- file.path(script_dir, "glia_chronic_significant_intersection.csv")
  write_csv(glia_intersection_summary_chronic, glia_intersection_file_chronic)
  cat("Summary written to file:", glia_intersection_file_chronic, "\n")
  
  # Write all intersection gene lists to separate files
  cat("\nWriting intersection gene lists to files...\n")
  
  # Intersection (common to both)
  if (length(common_both_glia_chronic) > 0) {
    common_genes_df_glia_chronic <- data.frame(Gene_name = common_both_glia_chronic, stringsAsFactors = FALSE)
    common_genes_file_glia_chronic <- file.path(script_dir, "glia_chronic_common_both_genes.csv")
    write_csv(common_genes_df_glia_chronic, common_genes_file_glia_chronic)
    cat("  Common to both:", common_genes_file_glia_chronic, "(", length(common_both_glia_chronic), "genes)\n")
  }
  
  # Unique genes
  if (length(only_astro_chronic) > 0) {
    write_csv(data.frame(Gene_name = only_astro_chronic, stringsAsFactors = FALSE),
              file.path(script_dir, "glia_chronic_only_Astroglia.csv"))
    cat("  Only Astroglia:", length(only_astro_chronic), "genes\n")
  }
  if (length(only_micro_chronic) > 0) {
    write_csv(data.frame(Gene_name = only_micro_chronic, stringsAsFactors = FALSE),
              file.path(script_dir, "glia_chronic_only_Microglia.csv"))
    cat("  Only Microglia:", length(only_micro_chronic), "genes\n")
  }
  
  # Create Venn diagram
  cat("\nCreating Venn diagram...\n")
  venn_file <- file.path(script_dir, "glia_chronic_venn_diagram.png")
  create_venn_diagram(glia_vars_chronic, 
                      title = "Glia Cell Types - Chronic (Significant Genes)",
                      output_file = venn_file,
                      use_sig_only = TRUE)
  
} else {
  cat("Error: Could not find both variables. Found", length(glia_sig_genes_list_chronic), "out of 2.\n")
}

# ===== CREATE GENE NETWORKS FOR COMMON GENES =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("CREATING GENE-GENE INTERACTION NETWORKS FOR COMMON GENES\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Create networks for common genes from intersections
common_gene_lists <- list()

# Collect all common gene lists
if (exists("imm_common_all_three", envir = .GlobalEnv)) {
  genes <- get("imm_common_all_three", envir = .GlobalEnv)
  if (length(genes) > 0) {
    common_gene_lists[["imm_acute_common_all_three"]] <- genes
  }
}
if (exists("imm_common_all_three_chronic", envir = .GlobalEnv)) {
  genes <- get("imm_common_all_three_chronic", envir = .GlobalEnv)
  if (length(genes) > 0) {
    common_gene_lists[["imm_chronic_common_all_three"]] <- genes
  }
}
if (exists("olig_common_all_three_acute", envir = .GlobalEnv)) {
  genes <- get("olig_common_all_three_acute", envir = .GlobalEnv)
  if (length(genes) > 0) {
    common_gene_lists[["olig_acute_common_all_three"]] <- genes
  }
}
if (exists("olig_common_all_three_chronic", envir = .GlobalEnv)) {
  genes <- get("olig_common_all_three_chronic", envir = .GlobalEnv)
  if (length(genes) > 0) {
    common_gene_lists[["olig_chronic_common_all_three"]] <- genes
  }
}
if (exists("glia_common_both_acute", envir = .GlobalEnv)) {
  genes <- get("glia_common_both_acute", envir = .GlobalEnv)
  if (length(genes) > 0) {
    common_gene_lists[["glia_acute_common_both"]] <- genes
  }
}
if (exists("glia_common_both_chronic", envir = .GlobalEnv)) {
  genes <- get("glia_common_both_chronic", envir = .GlobalEnv)
  if (length(genes) > 0) {
    common_gene_lists[["glia_chronic_common_both"]] <- genes
  }
}

# Create networks
network_results <- list()
for (name in names(common_gene_lists)) {
  genes <- common_gene_lists[[name]]
  network <- create_gene_network(genes, name, results_base_dir = script_dir)
  if (!is.null(network)) {
    network_results[[name]] <- network
  }
}

cat("\nCreated", length(network_results), "gene networks\n")
cat("Network files saved as .graphml format (can be opened in Cytoscape)\n")

# ===== VARIABLE NAMES CREATED =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("ALL VARIABLES CREATED:\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("\nFull datasets:\n")
for (var_name in created_variables) {
  cat("  -", var_name, "\n")
}

if (length(sig_variables) > 0) {
  cat("\nSignificant genes only (Sig=1):\n")
  for (sig_var_name in sig_variables) {
    cat("  -", sig_var_name, "\n")
  }
}

cat("\nYou can now access the data directly using the variable names, for example:\n")
if (length(created_variables) > 0) {
  cat("  # Full dataset:\n")
  cat("  ", created_variables[1], "\n")
  cat("  head(", created_variables[1], ")\n")
  cat("  dim(", created_variables[1], ")\n")
  if (length(sig_variables) > 0) {
    cat("\n  # Significant genes only:\n")
    cat("  ", sig_variables[1], "\n")
    cat("  head(", sig_variables[1], ")\n")
    cat("  dim(", sig_variables[1], ")\n")
  }
}

# ===== CREATE MAPPING DATA FRAME =====
# Create a data frame with the file to variable mapping for easy reference
# Only include successfully read files
file_var_mapping_final <- data.frame(
  File = successful_files,
  Variable = created_variables,
  stringsAsFactors = FALSE
)

# Store the mapping in a variable for easy access
file_to_variable_mapping <- file_var_mapping_final
assign("file_to_variable_mapping", file_var_mapping_final, envir = .GlobalEnv)

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("FILE TO VARIABLE MAPPING (stored in 'file_to_variable_mapping'):\n")
cat(paste(rep("=", 60), collapse=""), "\n")
print(file_var_mapping_final)

cat("\nAll files loaded successfully!\n")
cat("Access the mapping with: file_to_variable_mapping\n")

# ===== RUN enrichR ON INTERSECTION GENE LISTS =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("RUNNING enrichR ANALYSIS ON INTERSECTION GENE LISTS\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Define intersection files to analyze
intersection_files <- c(
  # Imm acute
  "imm_acute_common_all_three_genes.csv",
  "imm_acute_CGE_LGE_intersection.csv",
  "imm_acute_CGE_MGE_intersection.csv",
  "imm_acute_LGE_MGE_intersection.csv",
  "imm_acute_only_CGE.csv",
  "imm_acute_only_LGE.csv",
  "imm_acute_only_MGE.csv",
  # Imm chronic
  "imm_chronic_common_all_three_genes.csv",
  "imm_chronic_CGE_LGE_intersection.csv",
  "imm_chronic_CGE_MGE_intersection.csv",
  "imm_chronic_LGE_MGE_intersection.csv",
  "imm_chronic_only_CGE.csv",
  "imm_chronic_only_LGE.csv",
  "imm_chronic_only_MGE.csv",
  # Olig acute
  "olig_acute_common_all_three_genes.csv",
  "olig_acute_OPC_Premyel_intersection.csv",
  "olig_acute_OPC_Myel_intersection.csv",
  "olig_acute_Premyel_Myel_intersection.csv",
  "olig_acute_only_OPC.csv",
  "olig_acute_only_Premyel.csv",
  "olig_acute_only_Myel.csv",
  # Olig chronic
  "olig_chronic_common_all_three_genes.csv",
  "olig_chronic_OPC_Premyel_intersection.csv",
  "olig_chronic_OPC_Myel_intersection.csv",
  "olig_chronic_Premyel_Myel_intersection.csv",
  "olig_chronic_only_OPC.csv",
  "olig_chronic_only_Premyel.csv",
  "olig_chronic_only_Myel.csv",
  # Glia acute
  "glia_acute_common_both_genes.csv",
  "glia_acute_only_Astroglia.csv",
  "glia_acute_only_Microglia.csv",
  # Glia chronic
  "glia_chronic_common_both_genes.csv",
  "glia_chronic_only_Astroglia.csv",
  "glia_chronic_only_Microglia.csv"
)

# Run enrichR on each intersection file
enrichr_results_list <- list()

for (intersection_file in intersection_files) {
  file_path <- file.path(script_dir, intersection_file)
  
  if (file.exists(file_path)) {
    # Read the gene list
    tryCatch({
      gene_data <- read_csv(file_path, col_types = cols(), show_col_types = FALSE)
      
      if ("Gene_name" %in% colnames(gene_data)) {
        gene_list <- gene_data$Gene_name
      } else if (ncol(gene_data) > 0) {
        # Use first column if Gene_name doesn't exist
        gene_list <- gene_data[[1]]
      } else {
        cat("Skipping", intersection_file, "- no gene data found\n")
        next
      }
      
      # Create output name from file name
      output_name <- gsub("\\.csv$", "", intersection_file)
      
      # Run enrichR
      result <- run_enrichr_analysis(gene_list, output_name, results_base_dir = script_dir)
      
      if (!is.null(result)) {
        enrichr_results_list[[output_name]] <- result
      }
      
    }, error = function(e) {
      cat("Error processing", intersection_file, ":", conditionMessage(e), "\n")
    })
  } else {
    cat("File not found:", intersection_file, "- skipping\n")
  }
}

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("enrichR ANALYSIS COMPLETE!\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("Processed", length(enrichr_results_list), "intersection gene lists\n")
cat("Results are saved in directories ending with '_enrichR_results'\n")
cat("Each directory contains:\n")
cat("  - tables/: CSV files with enrichment results\n")
cat("  - plots/: PNG files with visualization plots\n")
cat("\n")

# ===== RUN clusterProfiler ON INTERSECTION GENE LISTS =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("RUNNING clusterProfiler ANALYSIS ON INTERSECTION GENE LISTS\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Use the same intersection files list as enrichR
# Run clusterProfiler on each intersection file
clusterprofiler_results_list <- list()

for (intersection_file in intersection_files) {
  file_path <- file.path(script_dir, intersection_file)
  
  if (file.exists(file_path)) {
    # Read the gene list
    tryCatch({
      gene_data <- read_csv(file_path, col_types = cols(), show_col_types = FALSE)
      
      if ("Gene_name" %in% colnames(gene_data)) {
        gene_list <- gene_data$Gene_name
      } else if (ncol(gene_data) > 0) {
        # Use first column if Gene_name doesn't exist
        gene_list <- gene_data[[1]]
      } else {
        cat("Skipping", intersection_file, "- no gene data found\n")
        next
      }
      
      # Create output name from file name
      output_name <- gsub("\\.csv$", "", intersection_file)
      
      # Run clusterProfiler
      result <- run_clusterprofiler_analysis(gene_list, output_name, results_base_dir = script_dir)
      
      if (!is.null(result)) {
        clusterprofiler_results_list[[output_name]] <- result
      }
      
    }, error = function(e) {
      cat("Error processing", intersection_file, ":", conditionMessage(e), "\n")
    })
  } else {
    cat("File not found:", intersection_file, "- skipping\n")
  }
}

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("clusterProfiler ANALYSIS COMPLETE!\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("Processed", length(clusterprofiler_results_list), "intersection gene lists\n")
cat("Results are saved in directories ending with '_clusterProfiler_results'\n")
cat("Each directory contains:\n")
cat("  - tables/: TXT files with enrichment results (GO, KEGG, Reactome, WikiPathways, MSigDB)\n")
cat("  - plots/: PNG files with dotplot visualizations\n")
cat("\n")

# ===== GENERATE MASTER REPORTS =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("GENERATING MASTER SUMMARY REPORTS\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# Generate Excel report
excel_report <- create_master_excel_report(results_base_dir = script_dir)

# Generate HTML report
html_report <- create_html_report(results_base_dir = script_dir)

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("ALL ANALYSES COMPLETE!\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("\nSummary of generated files:\n")
cat("  - Master_Summary_Report.xlsx: Comprehensive Excel workbook with all results\n")
cat("  - Master_Summary_Report.html: Interactive HTML report\n")
cat("  - Statistical comparisons: statistical_comparisons_acute_vs_chronic.csv\n")
cat("  - Jaccard similarity: jaccard_similarity_intersections.csv\n")
cat("  - Intersection gene lists: Multiple CSV files organized by category\n")
cat("  - Network files: .graphml files for Cytoscape visualization\n")
cat("  - Enrichment results: enrichR and clusterProfiler result directories\n")
cat("\n")

