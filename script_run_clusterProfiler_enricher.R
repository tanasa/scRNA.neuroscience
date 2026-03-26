#!/usr/bin/env Rscript
# Custom pathway enrichment with clusterProfiler (Geneformer ISP hits).
#
# Analysis folder (default: "enrichment 2_only_pos" under this script's directory):
#   - *.csv  = genes of interest (column Gene_name)
#   - Pathways = fixed list PATHWAY_FILENAMES (FerrDb gene_symbols, GOBP, Hallmark) — see below
#   - LIST.GENES.in.all.ISP.hits.for.background.with.clusterProfiler.txt = full background universe
#
# Outputs (under <ANALYSIS_DIR>/clusterProfiler_results/):
#   - <basename>_enrichment_all_pathways.csv (all tested pathways)
#   - <basename>_enrichment_fdr_filtered.csv (BH p.adjust < PVALUE_CUTOFF; enricher filter)
#   - <basename>_dotplot_all_categories.png (all tested pathways, no p.adjust filter)
#   - <basename>_dotplot_BH_padj_lt_<PVALUE_CUTOFF>.png (BH p.adjust < PVALUE_CUTOFF only)
#   - <basename>_dotplot_nominal_p_lt_<NOMINAL_P_CUTOFF>.png (raw hypergeometric p-value < NOMINAL_P_CUTOFF only)
#   - <basename>_enrichment_all_pathways_neglog10_pvalue.png (horizontal −log10(pvalue); order by BH p.adjust)
#   - <basename>_enrichment_all_pathways_neglog10_p_adjust.png (−log10(p.adjust))
#   - <basename>_enrichment_all_pathways_neglog10_qvalue.png (−log10(qvalue); pathways with NA q-value omitted)
#   - <basename>_enrichment_all_pathways_neglog10_*_refs_p005_p02.png — same as above + dashed vertical lines at p=0.05 and p=0.2 (−log10 scale)
#
# Optional: pass analysis directory as first argument:
#   Rscript run_clusterProfiler_custom_pathways.R "/path/to/enrichment 2_only_pos"

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
})

UNIVERSE_FILENAME <- "LIST.GENES.in.all.ISP.hits.for.background.with.clusterProfiler.txt"

# BH threshold for filtered table + second dotplot. enricher() applies pvalueCutoff to raw & adjusted p (?enricher).
PVALUE_CUTOFF <- 0.2
QVALUE_CUTOFF <- 1
# Third dotplot: filter by nominal (raw) hypergeometric p-value only — ignores BH p.adjust and q-value.
NOMINAL_P_CUTOFF <- 0.2
# Show every tested pathway on the “all categories” dotplot (enricher must not pre-filter).
PVALUE_CUTOFF_ALL <- 1

# Include all user-defined pathways: clusterProfiler defaults (minGSSize=10, maxGSSize=500) drop small/large sets.
MIN_GS_SIZE <- 1L
MAX_GS_SIZE <- 10000L

# Vertical reference lines on −log10 barplots (x = −log10(p) for each nominal p).
NEGLOG_P_REF_LINES <- c(0.05, 0.2)

# Exact pathway gene-list files (must exist inside ANALYSIS_DIR)
PATHWAY_FILENAMES <- c(
  "ferrdb_driver.APOPTOSIS_gene_symbols.txt",
  "ferrdb_driver.FERROPTOSIS_gene_symbols.txt",
  "ferrdb_marker.APOPTOSIS_gene_symbols.txt",
  "ferrdb_marker.FERROPTOSIS_gene_symbols.txt",
  "ferrdb_suppressor.APOPTOSIS_gene_symbols.txt",
  "ferrdb_suppressor.FERROPTOSIS_gene_symbols.txt",
  "ferrdb_unclassified.APOPTOSIS_gene_symbols.txt",
  "ferrdb_unclassified.FERROPTOSIS_gene_symbols.txt",
  "GOBP_NEGATIVE_REGULATION_OF_AUTOPHAGY.txt",
  "GOBP_NEGATIVE_REGULATION_OF_CYTOSKELETON_ORGANIZATION.txt",
  "GOBP_POSITIVE_REGULATION_OF_AUTOPHAGY.txt",
  "GOBP_POSITIVE_REGULATION_OF_CYTOSKELETON_ORGANIZATION.txt",
  "GOBP_REGULATION_OF_AUTOPHAGY.txt",
  "GOBP_REGULATION_OF_CYTOSKELETON_ORGANIZATION.txt",
  "GOBP_REGULATION_OF_MICROTUBULE_CYTOSKELETON_ORGANIZATION.txt",
  "GOBP_REGULATION_OF_MODIFICATION_OF_POSTSYNAPTIC_ACTIN_CYTOSKELETON.txt",
  "HALLMARK_APOPTOSIS.txt",
  "HALLMARK_GLYCOLYSIS.txt",
  "HALLMARK_HEDGEHOG_SIGNALING.txt",
  "HALLMARK_HYPOXIA.txt",
  "HALLMARK_IL2_STAT5_SIGNALING.txt",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING.txt",
  "HALLMARK_INFLAMMATORY_RESPONSE.txt",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE.txt",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE.txt",
  "HALLMARK_MTORC1_SIGNALING.txt",
  "HALLMARK_NOTCH_SIGNALING.txt",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION.txt",
  "HALLMARK_PEROXISOME.txt",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING.txt",
  "HALLMARK_PROTEIN_SECRETION.txt",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY.txt",
  "HALLMARK_TGF_BETA_SIGNALING.txt",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING.txt"
)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
BASE <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()

args_trailing <- commandArgs(trailingOnly = TRUE)
if (length(args_trailing) >= 1L) {
  ANALYSIS_DIR <- normalizePath(args_trailing[1])
} else {
  ANALYSIS_DIR <- file.path(BASE, "enrichment 2_only_pos")
}

universe_file <- file.path(ANALYSIS_DIR, UNIVERSE_FILENAME)
out_root <- file.path(ANALYSIS_DIR, "clusterProfiler_results")

if (!dir.exists(ANALYSIS_DIR)) {
  stop("Analysis directory not found: ", ANALYSIS_DIR, call. = FALSE)
}
if (!file.exists(universe_file)) {
  stop("Universe file not found: ", universe_file, call. = FALSE)
}

message("Analysis directory: ", ANALYSIS_DIR)
message(
  "enricher gene-set size filter: minGSSize=", MIN_GS_SIZE, ", maxGSSize=", MAX_GS_SIZE,
  " (clusterProfiler defaults 10/500 would omit small or large pathways)"
)

# Horizontal barplot of −log10(metric) for each pathway; y order = increasing BH p.adjust (best at top).
# If p_ref_lines is non-empty, draw vertical lines at x = −log10(p) for each p (e.g. 0.05 and 0.2).
plot_enrichment_all_neglog10 <- function(df, value_col, out_png, p_ref_lines = NULL) {
  req <- c("ID", "p.adjust", value_col)
  if (!all(req %in% names(df))) {
    message("  skip ", basename(out_png), " (missing column(s): ", value_col, ")")
    return(invisible(NULL))
  }
  d <- df[!is.na(df[[value_col]]) & as.numeric(df[[value_col]]) > 0, , drop = FALSE]
  if (nrow(d) < 1L) {
    message("  skip ", basename(out_png), " (no finite ", value_col, ")")
    return(invisible(NULL))
  }
  d$neglog10 <- -log10(as.numeric(d[[value_col]]))
  ord <- order(d$p.adjust, na.last = TRUE)
  d$ID <- factor(d$ID, levels = rev(unique(as.character(d$ID)[ord])))
  xl <- sprintf("-log10(%s)", value_col)
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$neglog10, y = .data$ID)) +
    ggplot2::geom_col(width = 0.85, fill = "#2E5F8C") +
    ggplot2::labs(x = xl, y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

  pr <- as.numeric(p_ref_lines)
  pr <- unique(sort(pr[!is.na(pr) & pr > 0 & pr < 1]))
  if (length(pr) > 0L) {
    ref_cols <- c("#C44E52", "#4C72B0", "#55A868", "#8172B3")
    for (i in seq_along(pr)) {
      g <- g + ggplot2::geom_vline(
        xintercept = -log10(pr[i]),
        linetype = "dashed",
        color = ref_cols[((i - 1L) %% length(ref_cols)) + 1L],
        linewidth = 0.55,
        alpha = 0.95
      )
    }
    cap <- paste(
      "Dashed lines: p = ",
      paste(sprintf("%g (x = −log10(p) = %.3f)", pr, -log10(pr)), collapse = "; "),
      sep = ""
    )
    g <- g + ggplot2::labs(caption = cap) +
      ggplot2::theme(plot.caption = ggplot2::element_text(size = 9, hjust = 0))
  }

  nr <- nrow(d)
  grDevices::png(out_png, width = 2800, height = max(2000L, 600L + 45L * nr), res = 200)
  print(g)
  grDevices::dev.off()
  message("  wrote ", out_png, " (", nr, " pathways, ", value_col, ")")
  invisible(NULL)
}

read_genes_from_isp_csv <- function(csv_path) {
  df <- utils::read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"Gene_name" %in% names(df)) {
    stop("Missing column Gene_name in ", csv_path, call. = FALSE)
  }
  g <- unique(trimws(as.character(df[["Gene_name"]])))
  g[!is.na(g) & nzchar(g)]
}

# --- 1. Universe ---
universe_genes <- trimws(readLines(universe_file, warn = FALSE))
universe_genes <- universe_genes[nzchar(universe_genes)]
universe_genes <- unique(universe_genes)
message("Universe (ISP background): ", length(universe_genes), " genes")

# --- 2. TERM2GENE: fixed pathway list (PATHWAY_FILENAMES) ---
pathway_files <- file.path(ANALYSIS_DIR, PATHWAY_FILENAMES)
miss <- pathway_files[!file.exists(pathway_files)]
if (length(miss)) {
  stop(
    "Missing pathway file(s) in analysis directory:\n",
    paste(basename(miss), collapse = "\n"),
    call. = FALSE
  )
}
message("Pathway files (fixed list): ", length(pathway_files))

pathway_parts <- lapply(pathway_files, function(f) {
  genes <- trimws(readLines(f, warn = FALSE))
  genes <- genes[nzchar(genes)]
  genes <- unique(genes)
  if (length(genes) < 1L) {
    warning("Skipping empty pathway file: ", basename(f), call. = FALSE)
    return(NULL)
  }
  data.frame(
    term = tools::file_path_sans_ext(basename(f)),
    gene = genes,
    stringsAsFactors = FALSE
  )
})
pathway_parts <- Filter(Negate(is.null), pathway_parts)
if (length(pathway_parts) < 1L) {
  stop("No non-empty pathway files in PATHWAY_FILENAMES (check empty .txt).", call. = FALSE)
}
term2gene <- do.call(rbind, pathway_parts)
message(
  "TERM2GENE rows: ", nrow(term2gene), " (",
  length(unique(term2gene$term)), " pathways after dropping empty lists)"
)

run_one <- function(csv_path) {
  base_name <- tools::file_path_sans_ext(basename(csv_path))
  out_dir <- out_root
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  query_genes <- read_genes_from_isp_csv(csv_path)
  message("[", base_name, "] query genes (raw): ", length(query_genes))

  drop_q <- setdiff(query_genes, universe_genes)
  if (length(drop_q) > 0) {
    message(
      "  dropped ", length(drop_q), " gene(s) not in universe (showing up to 5): ",
      paste(head(drop_q, 5), collapse = ", ")
    )
  }
  query_genes <- intersect(query_genes, universe_genes)
  message("  query genes (in universe): ", length(query_genes))

  out_csv_all <- file.path(out_dir, paste0(base_name, "_enrichment_all_pathways.csv"))
  out_csv_sig <- file.path(out_dir, paste0(base_name, "_enrichment_fdr_filtered.csv"))
  out_png_all <- file.path(out_dir, paste0(base_name, "_dotplot_all_categories.png"))
  out_png_sig <- file.path(
    out_dir,
    paste0(base_name, "_dotplot_BH_padj_lt_", PVALUE_CUTOFF, ".png")
  )
  out_png_nominal <- file.path(
    out_dir,
    paste0(base_name, "_dotplot_nominal_p_lt_", NOMINAL_P_CUTOFF, ".png")
  )

  if (length(query_genes) < 1L) {
    message("  skip enrichment (no genes in universe).")
    write.csv(
      data.frame(note = "no query genes after intersecting with universe"),
      out_csv_all,
      row.names = FALSE
    )
    return(invisible(NULL))
  }

  result_all <- enricher(
    gene          = query_genes,
    TERM2GENE     = term2gene,
    universe      = universe_genes,
    pAdjustMethod = "BH",
    pvalueCutoff  = PVALUE_CUTOFF_ALL,
    qvalueCutoff  = QVALUE_CUTOFF,
    minGSSize     = MIN_GS_SIZE,
    maxGSSize     = MAX_GS_SIZE
  )

  if (is.null(result_all)) {
    message("  enricher returned NULL (no testable pathways).")
    write.csv(data.frame(note = "enricher returned NULL"), out_csv_all, row.names = FALSE)
    return(invisible(NULL))
  }

  df_all <- as.data.frame(result_all@result)
  write.csv(df_all, out_csv_all, row.names = FALSE)
  message("  wrote ", out_csv_all, " (", nrow(df_all), " pathways, all tested)")

  plot_enrichment_all_neglog10(
    df_all,
    "pvalue",
    file.path(out_dir, paste0(base_name, "_enrichment_all_pathways_neglog10_pvalue.png"))
  )
  plot_enrichment_all_neglog10(
    df_all,
    "p.adjust",
    file.path(out_dir, paste0(base_name, "_enrichment_all_pathways_neglog10_p_adjust.png"))
  )
  plot_enrichment_all_neglog10(
    df_all,
    "qvalue",
    file.path(out_dir, paste0(base_name, "_enrichment_all_pathways_neglog10_qvalue.png"))
  )

  plot_enrichment_all_neglog10(
    df_all,
    "pvalue",
    file.path(out_dir, paste0(base_name, "_enrichment_all_pathways_neglog10_pvalue_refs_p005_p02.png")),
    p_ref_lines = NEGLOG_P_REF_LINES
  )
  plot_enrichment_all_neglog10(
    df_all,
    "p.adjust",
    file.path(out_dir, paste0(base_name, "_enrichment_all_pathways_neglog10_p_adjust_refs_p005_p02.png")),
    p_ref_lines = NEGLOG_P_REF_LINES
  )
  plot_enrichment_all_neglog10(
    df_all,
    "qvalue",
    file.path(out_dir, paste0(base_name, "_enrichment_all_pathways_neglog10_qvalue_refs_p005_p02.png")),
    p_ref_lines = NEGLOG_P_REF_LINES
  )

  result_sig <- enricher(
    gene          = query_genes,
    TERM2GENE     = term2gene,
    universe      = universe_genes,
    pAdjustMethod = "BH",
    pvalueCutoff  = PVALUE_CUTOFF,
    qvalueCutoff  = QVALUE_CUTOFF,
    minGSSize     = MIN_GS_SIZE,
    maxGSSize     = MAX_GS_SIZE
  )

  if (is.null(result_sig)) {
    message("  enricher (BH < ", PVALUE_CUTOFF, ") returned NULL; writing empty filtered CSV.")
    write.csv(data.frame(note = "enricher returned NULL for BH cutoff"), out_csv_sig, row.names = FALSE)
  } else {
    df_sig <- as.data.frame(result_sig)
    write.csv(df_sig, out_csv_sig, row.names = FALSE)
    message(
      "  wrote ", out_csv_sig, " (", nrow(df_sig),
      " pathways after enricher filter: BH p.adjust < ", PVALUE_CUTOFF,
      " (and raw p < ", PVALUE_CUTOFF, " per ?enricher))"
    )
  }

  n_all <- nrow(df_all)
  png(out_png_all, width = 2800, height = max(2000L, 800L + 65L * n_all), res = 200)
  print(dotplot(result_all, showCategory = n_all))
  dev.off()
  message("  wrote ", out_png_all, " (", n_all, " categories)")

  if (!is.null(result_sig) && nrow(as.data.frame(result_sig)) > 0L) {
    n_sig <- nrow(as.data.frame(result_sig))
    png(out_png_sig, width = 2800, height = max(2000L, 800L + 65L * n_sig), res = 200)
    print(dotplot(result_sig, showCategory = n_sig))
    dev.off()
    message("  wrote ", out_png_sig, " (", n_sig, " categories, BH p.adjust < ", PVALUE_CUTOFF, ")")
  } else {
    message("  no terms with BH p.adjust < ", PVALUE_CUTOFF, "; skipped ", basename(out_png_sig))
  }

  res_tbl <- result_all@result
  keep_nom <- res_tbl$pvalue < NOMINAL_P_CUTOFF & !is.na(res_tbl$pvalue)
  n_nom <- sum(keep_nom)
  if (n_nom > 0L) {
    result_nominal <- result_all
    result_nominal@result <- res_tbl[keep_nom, , drop = FALSE]
    png(out_png_nominal, width = 2800, height = max(2000L, 800L + 65L * n_nom), res = 200)
    print(dotplot(result_nominal, showCategory = n_nom))
    dev.off()
    message(
      "  wrote ", out_png_nominal, " (", n_nom, " categories, nominal p < ", NOMINAL_P_CUTOFF, ")"
    )
  } else {
    message(
      "  no terms with nominal p < ", NOMINAL_P_CUTOFF, "; skipped ", basename(out_png_nominal)
    )
  }
}

csv_files <- list.files(ANALYSIS_DIR, pattern = "\\.csv$", full.names = TRUE)
if (length(csv_files) == 0L) {
  stop("No .csv files in ", ANALYSIS_DIR, call. = FALSE)
}

message("--- Enrichment (", length(csv_files), " query CSVs) ---")
for (f in csv_files) run_one(f)

message("Done. Results under: ", out_root)
