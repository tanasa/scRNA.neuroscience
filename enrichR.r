# Read file, filter Sig=1 genes, run enrichR, save results

# Load required libraries
library(enrichR)
library(readr)
library(dplyr)
library(ggplot2)

# ===== SETUP FILES AND DIRECTORIES =====
input_file <- "hypoxia9000proportional_ISP_Imm-MGE_chronic.txt"

# Create results directory
base_name <- gsub("\\.txt$", "", basename(input_file))
results_dir <- paste0(base_name, "_enrichR_results")
dir.create(results_dir, showWarnings = FALSE)
dir.create(file.path(results_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(results_dir, "tables"), showWarnings = FALSE)

cat("Results will be saved to:", results_dir, "\n")

# ===== READ THE DATA FILE =====
data <- read_csv(input_file, col_types = cols())
cat("Data dimensions:", dim(data), "\n")
cat("Column names:", paste(colnames(data), collapse = ", "), "\n")
cat("First few rows:\n")
print(head(data, 3))

head(data)

# ===== FILTER FOR SIGNIFICANT GENES (Sig=1) =====
cat("\nFiltering for significant genes...\n")

# Check Sig column
if ("Sig" %in% colnames(data)) {
  sig_counts <- table(data$Sig, useNA = "always")
  cat("Sig column distribution:\n")
  print(sig_counts)
  
  # Filter for Sig=1
  sig_data <- data %>% filter(Sig == 1)
  cat("Genes with Sig=1:", nrow(sig_data), "out of", nrow(data), "total\n")
} else {
  cat("No 'Sig' column found! Available columns:", paste(colnames(data), collapse = ", "), "\n")
  stop("Cannot find Sig column for filtering")
}

head(sig_data)
tail(sig_data)

# ===== EXTRACT GENE NAMES =====
if ("Gene_name" %in% colnames(sig_data)) {
  gene_list <- sig_data$Gene_name
  gene_column <- "Gene_name"
} else if ("Gene" %in% colnames(sig_data)) {
  gene_list <- sig_data$Gene
  gene_column <- "Gene"
} else {
  # Use second column if first is row numbers
  gene_list <- sig_data[[2]]
  gene_column <- colnames(sig_data)[2]
}

# Clean gene names
gene_list <- as.character(gene_list)
gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
gene_list <- unique(gene_list)

cat("Using gene column:", gene_column, "\n")
cat("Final gene list length:", length(gene_list), "\n")
cat("The gene list:\n")
print(gene_list)



# ===== DEFINE enrichR DATABASES =====

enrichr_databases <- c(
  # Updated Gene Ontology (2025 versions)
  "GO_Biological_Process_2025",
  "GO_Molecular_Function_2025", 
  "GO_Cellular_Component_2025",
  
  # Core pathway databases
  "KEGG_2021_Human",
  
  # Pathway Databases (from your images)
  "Reactome_Pathways_2024",
  "WikiPathways_2024_Human", 
  "BioPlanet_2019",
  "WikiPathways_2024_Mouse",
  
  # Specialized Collections (from your images)
  "Elsevier_Pathway_Collection",
  "MSigDB_Hallmark_2020",
  "BioCarta_2016",
  "HumanCyc_2016",
  "NCI-Nature_2016",
  "Panther_2016",
  
  # Coexpression Database (from your images)
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
  
  # Additional databases that may be available with different names
  "Reactome_2022",
  "WikiPathways_2023_Human"
)

cat("Databases to query:", length(enrichr_databases), "\n")


# ===== RUN enrichR ANALYSIS =====
cat("\nRunning enrichR analysis...\n")
enrichr_results <- enrichr(gene_list, enrichr_databases)

# ===== PROCESS AND SAVE RESULTS =====
# adj_p_cutoff <- 0.05
# top_n <- 20

# ===== PROCESS AND SAVE RESULTS =====
adj_p_cutoff <- 0.2
top_n <- 50

# Save significant genes list
write_csv(data.frame(Gene_name = gene_list), 
          file.path(results_dir, "tables", paste0(base_name, "_significant_genes.csv")))

# Process each database
all_significant_results <- data.frame()
summary_stats <- data.frame()

for (database in names(enrichr_results)) {
  cat("Processing", database, "...\n")
  
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
                file.path(results_dir, "tables", paste0(base_name, "_", database, ".csv")))
      
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
          x = "-log10(Adjusted P-value)",
          y = ""
        ) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 9),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)
        )
      
      # Save plot with white background
      ggsave(file.path(results_dir, "plots", paste0(base_name, "_", database, ".png")), 
             p, width = 10, height = 6, dpi = 300, bg = "white")
      
      # Add to summary
      summary_stats <- rbind(summary_stats, data.frame(
        Database = database,
        Significant_Terms = nrow(significant),
        Top_Term = significant$Term[1],
        Top_P_Value = significant$Adjusted.P.value[1]
      ))
    }
  }
}

# ===== SAVE COMBINED RESULTS =====
# Save combined top results
if (nrow(all_significant_results) > 0) {
  write_csv(all_significant_results, 
            file.path(results_dir, "tables", paste0(base_name, "_all_top_results.csv")))
}

# Save summary
if (nrow(summary_stats) > 0) {
  write_csv(summary_stats, 
            file.path(results_dir, "tables", paste0(base_name, "_summary.csv")))
}

# ===== PRINT FINAL SUMMARY =====
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("ANALYSIS COMPLETE!\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("Input file:", input_file, "\n")
cat("Total genes in file:", nrow(data), "\n")
cat("Genes with Sig=1:", length(gene_list), "\n")
cat("Results directory:", results_dir, "\n")
cat("Databases with significant results:", nrow(summary_stats), "\n")

if (nrow(summary_stats) > 0) {
  cat("\nTop results by database:\n")
  print(summary_stats)
}

cat("\nFiles created:\n")
cat("- Significant genes: tables/", base_name, "_significant_genes.csv\n", sep="")
cat("- Combined results: tables/", base_name, "_all_top_results.csv\n", sep="")
cat("- Summary: tables/", base_name, "_summary.csv\n", sep="")
cat("- Individual database files: tables/[database_name].csv\n")
cat("- Plots: plots/[database_name].png\n")
cat("Done!\n")




