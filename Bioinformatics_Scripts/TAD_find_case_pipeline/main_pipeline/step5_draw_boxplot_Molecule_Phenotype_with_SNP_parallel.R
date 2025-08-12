# SNP-TAD-Gene Expression-Phenotype Visualization (Refactored)
# Author: Refactored for improved maintainability and robustness
# Date: 2025
# Version: 2.7.0 (CLI-controlled cores and logging)

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(dplyr)
  library(scales)
  library(foreach)
  library(doParallel)
})

# --- Command-Line Argument Parsing ---
args <- commandArgs(trailingOnly = TRUE)
# --- 修改点②: 参数数量检查从3改为4 ---
if (length(args) != 4) {
  stop("Usage: Rscript your_script_name.R <path_to_triplets_file> <path_to_snp_data_file> <path_to_output_directory> <num_cores>", call. = FALSE)
}
# --- End of Argument Parsing ---


# Configuration constants
CONFIG <- list(
  # --- 修改点②: N_CORES 从硬编码改为从命令行参数读取 ---
  N_CORES = as.integer(args[4]), # <--- 从第4个命令行参数获取核心数
  BASE_SIZE = 8,
  NATURE_COLORS = c("#1f77b4", "#d62728", "#ff7f0e", "#2ca02c", "#9467bd"),
  MIN_SAMPLES = 10,
  MIN_DATA_ROWS = 5,
  EXPECTED_GENOTYPES = c(0, 2),  # Homozygous only
  PLOT_WIDTH = 5,
  PLOT_HEIGHT = 8,
  DPI = 300
)

# File paths
PATHS <- list(
  TRIPLETS = args[1],
  SNP_DATA = args[2],
  OUTPUT_DIR = args[3],
  GENE_EXP = "/share/home/cyzhao24/workplace/01.ML/02.loop_association/coloc/3d_with_gene_new/tad_gene_vaild/raw_data/formatted_rna_data_for_qtl.bed",
  TAD_DATA = "/share/home/cyzhao24/workplace/01.ML/02.loop_association/coloc/3d_with_gene_new/tad_gene_vaild/raw_data/pan3Dgenome_TAD_boundary_IS_244sample.csv",
  PHENOTYPE = "/share/home/cyzhao24//workplace/01.ML/02.loop_association/wmj_structure/copy_pipeline_huang/pipeline_04_run_fastlmm/all_data_for_fastlmm/Lsample/phenotype/ezhou_blup_results_L244.txt"
)

# --- Helper Functions (All Included) ---

# Enhanced Nature/Cell theme
create_nature_theme <- function(base_size = CONFIG$BASE_SIZE) {
  theme_classic(base_size = base_size) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(2, "pt"),
      axis.text = element_text(color = "black", size = base_size, family = "sans"),
      axis.title = element_text(color = "black", size = base_size + 1, face = "bold", family = "sans"),
      axis.title.x = element_text(margin = margin(t = 4, b = 0)),
      axis.title.y = element_text(margin = margin(r = 4, l = 0)),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5,
                                family = "sans", margin = margin(b = 8)),
      plot.subtitle = element_blank(),
      legend.position = "none",
      plot.margin = margin(8, 8, 8, 8, "pt"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", family = "sans")
    )
}

# Data loading functions
load_all_data <- function() {
  # --- 修改点⑤: 将所有cat输出重定向到stderr ---
  cat("Loading data files...\n", file = stderr())
  
  data_list <- list(
    triplets = fread(PATHS$TRIPLETS, header = TRUE, sep = '\t'),
    snp_data = fread(PATHS$SNP_DATA, header = TRUE),
    gene_exp = fread(PATHS$GENE_EXP, header = TRUE, sep = '\t'),
    tad_data = fread(PATHS$TAD_DATA, header = TRUE, sep = ','),
    phenotype = fread(PATHS$PHENOTYPE, header = TRUE, sep = '\t')
  )
  
  # Clean column names
  colnames(data_list$gene_exp)[1] <- gsub("^#", "", colnames(data_list$gene_exp)[1])
  
  # --- 修改点⑤: 将所有cat输出重定向到stderr ---
  cat("Data loading complete:", nrow(data_list$triplets), "triplets\n", file = stderr())
  return(data_list)
}

# Utility functions
extract_alleles <- function(snp_id) {
  parts <- strsplit(snp_id, ":")[[1]]
  if (length(parts) >= 3) {
    ref_allele <- parts[length(parts) - 1]
    alt_allele <- parts[length(parts)]
    return(list(ref = ref_allele, alt = alt_allele))
  }
  return(list(ref = "A", alt = "G"))  # Default fallback
}

convert_to_homozygous <- function(genotypes, ref, alt) {
  sapply(genotypes, function(x) {
    if (is.na(x) || x == "NA" || x == "") return(NA)
    x_num <- as.numeric(x)
    if (is.na(x_num)) return(NA)
    if (x_num == 2) return(paste0(ref, "/", ref)) # REF/REF
    if (x_num == 0) return(paste0(alt, "/", alt)) # ALT/ALT
    return(NA)  # Exclude heterozygous
  })
}

perform_wilcox_test <- function(data, value_col, group_col) {
  complete_data <- data[complete.cases(data[, c(value_col, group_col)]), ]
  
  if (nrow(complete_data) < 3 || length(unique(complete_data[[group_col]])) < 2) {
    return(list(p_value = NA, results = NULL))
  }
  
  test_result <- tryCatch({
    wilcox.test(complete_data[[value_col]] ~ complete_data[[group_col]])
  }, error = function(e) NULL)
  
  if (is.null(test_result)) {
    return(list(p_value = NA, results = NULL))
  }
  
  groups <- sort(unique(complete_data[[group_col]]))
  pairwise_df <- data.frame(
    group1 = as.character(groups[1]),
    group2 = as.character(groups[2]),
    p.adj = test_result$p.value,
    stringsAsFactors = FALSE
  )
  
  return(list(p_value = test_result$p.value, results = pairwise_df))
}

calculate_correlation <- function(data, x_col, y_col) {
  complete_data <- data[complete.cases(data[, c(x_col, y_col)]), ]
  
  if (nrow(complete_data) < 3) {
    return(list(cor = NA, p_value = NA))
  }
  
  cor_result <- tryCatch({
    cor.test(complete_data[[x_col]], complete_data[[y_col]], method = "pearson")
  }, error = function(e) NULL)
  
  if (is.null(cor_result)) {
    return(list(cor = NA, p_value = NA))
  }
  
  return(list(cor = cor_result$estimate, p_value = cor_result$p.value))
}

# Data processing functions
find_common_samples <- function(data_list) {
  snp_samples <- setdiff(colnames(data_list$snp_data), colnames(data_list$snp_data)[1])
  gene_samples <- setdiff(colnames(data_list$gene_exp), c("chr", "start", "end", "pid", "gid", "strand"))
  tad_samples <- setdiff(colnames(data_list$tad_data), "GeneId")
  pheno_samples <- data_list$phenotype$id
  
  common_samples <- intersect(
    intersect(intersect(snp_samples, gene_samples), tad_samples), 
    pheno_samples
  )
  
  return(common_samples)
}

extract_sample_data <- function(data_list, sample_cols, snp_row, gene_row, tad_row) {
  snp_values <- as.numeric(unlist(snp_row[, .SD, .SDcols = sample_cols]))
  gene_values <- as.numeric(unlist(gene_row[, .SD, .SDcols = sample_cols]))
  tad_values <- as.numeric(unlist(tad_row[, .SD, .SDcols = sample_cols]))
  pheno_match <- data_list$phenotype[match(sample_cols, data_list$phenotype$id), ]
  
  if (any(is.na(pheno_match$FL)) || any(is.na(pheno_match$FS))) {
    return(NULL)
  }
  
  return(list(
    snp = snp_values, gene = gene_values, tad = tad_values,
    fl = pheno_match$FL, fs = pheno_match$FS
  ))
}

# Plotting functions
create_significance_symbol <- function(p_val) {
  if (is.na(p_val)) return("ns")
  if (p_val < 0.001) return("***")
  if (p_val < 0.01) return("**")
  if (p_val < 0.05) return("*")
  return("ns")
}

add_significance_annotations <- function(plot, pairwise_results, y_range, genotype_levels) {
  if (is.null(pairwise_results) || nrow(pairwise_results) == 0) return(plot)
  
  y_diff <- diff(y_range)
  if (is.na(y_diff) || y_diff == 0) y_diff <- 1
  base_height <- y_diff * 0.15
  line_spacing <- y_diff * 0.12
  
  for (i in 1:nrow(pairwise_results)) {
    comparison <- pairwise_results[i, ]
    group1_num <- as.numeric(as.character(comparison$group1))
    group2_num <- as.numeric(as.character(comparison$group2))
    
    x1_pos <- match(group1_num, genotype_levels)
    x2_pos <- match(group2_num, genotype_levels)
    y_pos <- y_range[2] + base_height + (line_spacing * (i - 1))
    sig_symbol <- create_significance_symbol(comparison$p.adj)
    
    plot <- plot +
      annotate("segment", x = x1_pos, xend = x2_pos, y = y_pos, yend = y_pos, size = 0.6) +
      annotate("segment", x = x1_pos, xend = x1_pos, y = y_pos, yend = y_pos - (line_spacing * 0.25), size = 0.6) +
      annotate("segment", x = x2_pos, xend = x2_pos, y = y_pos, yend = y_pos - (line_spacing * 0.25), size = 0.6) +
      annotate("text", x = (x1_pos + x2_pos) / 2, y = y_pos + (line_spacing * 0.18), label = sig_symbol, size = 3, fontface = "bold")
  }
  return(plot)
}

calculate_plot_limits <- function(data_range, num_comparisons) {
  base_padding <- 0.12
  count_space <- 0.18
  annotation_space <- 0.12
  cor_space <- 0.08
  range_diff <- diff(data_range)
  if (is.na(range_diff) || range_diff == 0) range_diff <- 1
  
  y_min <- data_range[1] - (range_diff * (base_padding + count_space))
  y_max <- data_range[2] + (range_diff * (base_padding + 0.08 + (annotation_space * num_comparisons) + cor_space + 0.05))
  return(c(y_min, y_max))
}

create_boxplot <- function(plot_data, y_var, y_label, title, genotype_counts, cor_stats, data_range, genotype_levels, ylim) {
  range_diff <- diff(data_range)
  if (is.na(range_diff) || range_diff == 0) range_diff <- 1
  base_padding <- 0.12
  count_space <- 0.18
  cor_space <- 0.08
  
  p <- ggplot(plot_data, aes(x = Genotype_factor, y = .data[[y_var]], fill = Genotype_factor)) +
    geom_boxplot(alpha = 0.85, width = 0.7, linewidth = 0.6, outlier.size = 1) +
    scale_fill_manual(values = CONFIG$NATURE_COLORS[1:length(genotype_levels)]) +
    scale_y_continuous(limits = ylim, expand = c(0, 0), labels = function(x) format(x, digits = 2, nsmall = 1)) +
    labs(title = title, x = "Genotype", y = y_label) +
    create_nature_theme() +
    geom_text(data = genotype_counts,
              aes(x = Genotype_factor,
                  y = data_range[1] - (range_diff * (base_padding + count_space * 0.6)),
                  label = paste0("n=", Count)),
              inherit.aes = FALSE, size = 2.8, color = "gray30") +
    annotate("text", x = length(genotype_levels)/2 + 0.5,
             y = data_range[2] + (range_diff * (base_padding + 0.08 + cor_space * 0.5)),
             label = sprintf("r=%.3f, P=%.2e", cor_stats$cor, cor_stats$p_value),
             size = 2.8, hjust = 0.5, fontface = "bold")
  return(p)
}


# --- Main Processing Function (Calculates, Filters, Plots) ---
process_triplet <- function(data_list, triplet_data, common_samples) {
  # Stage 1: Data Extraction and Validation
  gene_id <- triplet_data$Gene_ID
  tad_id <- triplet_data$TAD_ID
  snp_id <- triplet_data$SNP_ID
  
  snp_col_name <- colnames(data_list$snp_data)[1]
  escaped_snp_id <- gsub("([\\[\\]\\(\\)\\{\\}\\^\\$\\.\\*\\+\\?\\|\\\\])", "\\\\\\1", snp_id)
  pattern <- paste0("^", escaped_snp_id, ":")
  snp_row <- data_list$snp_data[grepl(pattern, data_list$snp_data[[snp_col_name]]), ]
  
  if (nrow(snp_row) == 0) return(list(success = FALSE, message = paste("SNP not found:", snp_id)))
  if (nrow(snp_row) > 1) snp_row <- snp_row[1, ]
  actual_snp_id <- snp_row[[snp_col_name]]
  
  gene_row <- data_list$gene_exp[pid == gene_id, ]
  tad_row <- data_list$tad_data[GeneId == tad_id, ]
  if (nrow(gene_row) == 0 || nrow(tad_row) == 0) return(list(success = FALSE, message = paste("Gene or TAD data not found:", gene_id, "/", tad_id)))
  
  sample_data <- extract_sample_data(data_list, common_samples, snp_row, gene_row, tad_row)
  if (is.null(sample_data)) return(list(success = FALSE, message = "Phenotype data matching failed"))
  
  alleles <- extract_alleles(actual_snp_id)
  plot_data <- data.frame(
    Genotype_num = sample_data$snp,
    Genotype = convert_to_homozygous(sample_data$snp, alleles$ref, alleles$alt),
    TAD_IS_Score = sample_data$tad, Gene_Exp = sample_data$gene,
    FL = sample_data$fl, FS = sample_data$fs, stringsAsFactors = FALSE
  )
  plot_data <- plot_data[complete.cases(plot_data), ]
  
  unique_genotypes <- unique(plot_data$Genotype_num)
  if (nrow(plot_data) < CONFIG$MIN_DATA_ROWS || length(unique_genotypes) != 2 || !all(unique_genotypes %in% CONFIG$EXPECTED_GENOTYPES)) {
    return(list(success = FALSE, message = "Insufficient data or incorrect genotype distribution"))
  }
  
  # Stage 2: Always Calculate All Statistics
  stats_list <- list(
    tad = perform_wilcox_test(plot_data, "TAD_IS_Score", "Genotype_num"),
    gene = perform_wilcox_test(plot_data, "Gene_Exp", "Genotype_num"),
    fl = perform_wilcox_test(plot_data, "FL", "Genotype_num"),
    fs = perform_wilcox_test(plot_data, "FS", "Genotype_num")
  )
  cor_list <- list(
    tad = calculate_correlation(plot_data, "Genotype_num", "TAD_IS_Score"),
    gene = calculate_correlation(plot_data, "Genotype_num", "Gene_Exp"),
    fl = calculate_correlation(plot_data, "Genotype_num", "FL"),
    fs = calculate_correlation(plot_data, "Genotype_num", "FS"),
    gene_tad = calculate_correlation(plot_data, "TAD_IS_Score", "Gene_Exp")
  )
  
  all_stats <- list(
    SNP_ID = snp_id, Gene_ID = gene_id, TAD_ID = tad_id, Full_SNP_ID = actual_snp_id,
    Wilcox_P_TAD = stats_list$tad$p_value,
    Wilcox_P_Gene = stats_list$gene$p_value,
    Wilcox_P_FL = stats_list$fl$p_value,
    Wilcox_P_FS = stats_list$fs$p_value,
    Cor_R_Gene_TAD = cor_list$gene_tad$cor,
    Cor_P_Gene_TAD = cor_list$gene_tad$p_value
  )
  
  # Stage 3: Filtering and Plotting Decision
  is_gene_significant <- !is.na(stats_list$gene$p_value) && stats_list$gene$p_value < 0.05
  is_corr_negative <- !is.na(cor_list$gene_tad$cor) && cor_list$gene_tad$cor < 0
  
  combined_plot <- NULL
  if (is_gene_significant && is_corr_negative) {
    genotype_levels <- sort(unique(plot_data$Genotype_num), decreasing = TRUE)
    plot_data$Genotype_factor <- factor(plot_data$Genotype_num, levels = genotype_levels, labels = sapply(genotype_levels, function(gt) unique(plot_data$Genotype[plot_data$Genotype_num == gt])[1]))
    genotype_counts <- plot_data %>% group_by(Genotype_factor) %>% summarise(Count = n(), .groups = 'drop')
    
    ranges <- list(tad = range(plot_data$TAD_IS_Score, na.rm = TRUE), gene = range(plot_data$Gene_Exp, na.rm = TRUE), fl = range(plot_data$FL, na.rm = TRUE), fs = range(plot_data$FS, na.rm = TRUE))
    num_comparisons <- sapply(stats_list, function(x) ifelse(is.null(x$results), 0, nrow(x$results)))
    ylims <- mapply(calculate_plot_limits, ranges, num_comparisons, SIMPLIFY = FALSE)
    
    plots <- list(
      p1 = create_boxplot(plot_data, "TAD_IS_Score", "TAD predicted IS-score", paste(tad_id), genotype_counts, cor_list$tad, ranges$tad, genotype_levels, ylims$tad),
      p2 = create_boxplot(plot_data, "Gene_Exp", "Gene expression", paste("Gene", gene_id), genotype_counts, cor_list$gene, ranges$gene, genotype_levels, ylims$gene),
      p_fl = create_boxplot(plot_data, "FL", "Phenotype FL", "Phenotype FL", genotype_counts, cor_list$fl, ranges$fl, genotype_levels, ylims$fl),
      p_fs = create_boxplot(plot_data, "FS", "Phenotype FS", "Phenotype FS", genotype_counts, cor_list$fs, ranges$fs, genotype_levels, ylims$fs)
    )
    plots$p1 <- add_significance_annotations(plots$p1, stats_list$tad$results, ranges$tad, genotype_levels)
    plots$p2 <- add_significance_annotations(plots$p2, stats_list$gene$results, ranges$gene, genotype_levels)
    plots$p_fl <- add_significance_annotations(plots$p_fl, stats_list$fl$results, ranges$fl, genotype_levels)
    plots$p_fs <- add_significance_annotations(plots$p_fs, stats_list$fs$results, ranges$fs, genotype_levels)
    
    main_title <- sprintf("Gene-TAD correlation: r=%.3f, P=%.2e (SNP: %s)", cor_list$gene_tad$cor, cor_list$gene_tad$p_value, actual_snp_id)
    combined_plot <- arrangeGrob(plots$p1, plots$p_fl, plots$p2, plots$p_fs, ncol = 2, padding = unit(0.5, "line"), top = textGrob(main_title, gp = gpar(fontsize = 10, fontface = "bold")))
  }
  
  # Stage 4: Return Comprehensive Results
  return(list(
    success = TRUE,
    plot = combined_plot, 
    stats = all_stats,     
    ids = list(gene_id = gene_id, snp_id = snp_id, tad_id = tad_id)
  ))
}

# --- Main Execution Function ---
main <- function() {
  # Setup
  data_list <- load_all_data()
  common_samples <- find_common_samples(data_list)
  if (length(common_samples) < CONFIG$MIN_SAMPLES) stop("Insufficient common samples")
  # --- 修改点⑤: 将所有cat输出重定向到stderr ---
  cat("Found", length(common_samples), "common samples\n", file = stderr())
  
  dir.create(PATHS$OUTPUT_DIR, showWarnings = FALSE)
  unique_tads <- unique(data_list$triplets$TAD_ID)
  for (tad in unique_tads) dir.create(file.path(PATHS$OUTPUT_DIR, tad), showWarnings = FALSE, recursive = TRUE)

  # Parallel Processing
  # --- 修改点⑤: 将所有cat输出重定向到stderr ---
  cat("Setting up parallel environment using", CONFIG$N_CORES, "CPU cores...\n", file = stderr())
  cl <- makeCluster(CONFIG$N_CORES)
  registerDoParallel(cl)
  
  # --- 修改点⑤: 将所有cat输出重定向到stderr ---
  cat("Starting parallel processing of all triplets...\n", file = stderr())
  total_count <- nrow(data_list$triplets)
  
  # --- OPTIMIZATION: Explicitly define objects to export ---
  # This is more robust and efficient than exporting the entire global environment.
  export_objects <- c(
    "process_triplet", "data_list", "common_samples", "PATHS", "CONFIG", 
    "create_nature_theme", "extract_alleles", "convert_to_homozygous", 
    "perform_wilcox_test", "calculate_correlation", "extract_sample_data", 
    "create_significance_symbol", "add_significance_annotations", 
    "calculate_plot_limits", "create_boxplot"
  )
  
  results_list <- tryCatch({
    foreach(
      i = 1:total_count,
      .export = export_objects,
      .packages = c("data.table", "ggplot2", "gridExtra", "grid", "dplyr", "scales")
    ) %dopar% {
      tryCatch({
        process_triplet(data_list, data_list$triplets[i, ], common_samples)
      }, error = function(e) {
        list(success = FALSE, message = paste("Error processing triplet", i, ":", e$message))
      })
    }
  }, finally = {
    # --- 修改点⑤: 将所有cat输出重定向到stderr ---
    cat("Shutting down parallel environment...\n", file = stderr())
    stopCluster(cl)
  })
  
  # Post-processing results, saving plots and CSV
  # --- 修改点⑤: 将所有cat输出重定向到stderr ---
  cat("Parallel processing complete. Consolidating results and filtering plots...\n", file = stderr())
  
  success_count <- 0
  plots_saved_count <- 0
  stats_to_export <- list()
  filtered_stats_to_export <- list()

  for (i in 1:length(results_list)) {
    res <- results_list[[i]]
    if (res$success) {
      success_count <- success_count + 1
      stats_to_export[[length(stats_to_export) + 1]] <- as.data.table(res$stats)

      if (!is.null(res$plot)) {
        plots_saved_count <- plots_saved_count + 1
        filtered_stats_to_export[[length(filtered_stats_to_export) + 1]] <- as.data.table(res$stats)
        
        safe_gene <- gsub("[^A-Za-z0-9_.-]", "_", res$ids$gene_id)
        safe_snp <- gsub("[^A-Za-z0-9_.-]", "_", res$ids$snp_id)
        filename <- file.path(PATHS$OUTPUT_DIR, res$ids$tad_id, paste0(safe_gene, "_", safe_snp, ".png"))
        
        ggsave(filename, res$plot, 
               width = CONFIG$PLOT_WIDTH, height = CONFIG$PLOT_HEIGHT, 
               dpi = CONFIG$DPI, units = "in", bg = "white", device = "png")
      }
    } else {
      # --- 修改点⑤: 将所有cat输出重定向到stderr ---
      cat(paste0("  - 警告: ", res$message, "\n"), file = stderr())
    }
  }

  # Export all statistics to a CSV file
  if (length(stats_to_export) > 0) {
    stats_df <- rbindlist(stats_to_export, use.names = TRUE, fill = TRUE)
    csv_filename <- file.path(PATHS$OUTPUT_DIR, "all_triplet_statistics.csv")
    fwrite(stats_df, csv_filename)
    # --- 修改点⑤: 将所有cat输出重定向到stderr ---
    cat("All triplet statistics saved to:", csv_filename, "\n", file = stderr())
  }
  
  # Export filtered statistics to a separate CSV file
  if (length(filtered_stats_to_export) > 0) {
    filtered_stats_df <- rbindlist(filtered_stats_to_export, use.names = TRUE, fill = TRUE)
    filtered_csv_filename <- file.path(PATHS$OUTPUT_DIR, "filtered_plotted_statistics.csv")
    fwrite(filtered_stats_df, filtered_csv_filename)
    # --- 修改点⑤: 将所有cat输出重定向到stderr ---
    cat("Statistics for plotted triplets saved to:", filtered_csv_filename, "\n", file = stderr())
  }

  # Updated final summary message
  # --- 修改点⑤: 将所有cat输出重定向到stderr ---
  cat("Analysis complete:\n", file = stderr())
  cat("  - Successfully analyzed", success_count, "of", total_count, "triplets.\n", file = stderr())
  cat("  - Found and plotted", plots_saved_count, "triplets that met the filtering criteria.\n", file = stderr())
  cat("  - Plots are saved in subdirectories within", PATHS$OUTPUT_DIR, "\n", file = stderr())
}

# Execute main function
if (!interactive()) {
  main()
}