# =============================================================================
# R脚本：根据基因型绘制Loop强度箱形图并进行关联检验
# 版本: 3.0 (修正版: 真正优化并行内存管理，解决unserialize错误)
# =============================================================================

# --- 1. 加载必要的R包 ---
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(foreach)
  library(doParallel)
})

# --- 解析命令行参数 (保持不变) ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("
  错误：提供的参数数量不正确！需要8个参数。
  用法: Rscript analyze_and_plot_loops.R ... <num_cores>
  ", call. = FALSE)
}

# --- 2. 定义文件路径和全局配置 (保持不变) ---
PATHS <- list(
  LOOP_ANALYSIS   = args[1],
  STATS_BRIDGE    = args[2],
  LOOP_STRENGTH   = args[3],
  GENOTYPE        = args[4],
  LOOP_PRESENCE   = args[5]
)
MAIN_OUTPUT_DIR <- args[6]
SUMMARY_CSV_PATH <- args[7]
NUM_CORES <- as.integer(args[8])

CONFIG <- list(
  BASE_SIZE = 8,
  NATURE_COLORS = c("#1f77b4", "#d62728", "#ff7f0e", "#2ca02c", "#9467bd"),
  MIN_SAMPLES = 3,
  EXPECTED_GENOTYPES = c(0, 2)
)

# =============================================================================
# 3. 辅助函数 (保持不变)
# =============================================================================
create_nature_theme <- function(base_size = CONFIG$BASE_SIZE) {
  theme_classic(base_size = base_size) +
    theme(
      panel.background = element_blank(), panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(2, "pt"),
      axis.text = element_text(color = "black", size = base_size),
      axis.title = element_text(color = "black", size = base_size + 1, face = "bold"),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5, margin = margin(b = 8)),
      legend.position = "none"
    )
}

extract_alleles <- function(snp_id) {
  parts <- strsplit(snp_id, ":")[[1]]
  if (length(parts) >= 3) {
    return(list(ref = parts[2], alt = parts[3]))
  }
  return(list(ref = "REF", alt = "ALT"))
}

convert_genotype_to_allele <- function(genotypes, ref, alt) {
  sapply(genotypes, function(g) {
    if (is.na(g)) return(NA)
    if (g == 2) return(paste0(ref, "/", ref))
    if (g == 0) return(paste0(alt, "/", alt))
    return(NA)
  })
}

perform_wilcox_test <- function(data, formula) {
  if (nrow(data) < 3 || length(unique(data$Genotype)) < 2) {
    return(list(p.value = NA))
  }
  tryCatch({
    wilcox.test(formula, data = data)
  }, error = function(e) list(p.value = NA))
}

calculate_correlation <- function(data, x_var, y_var) {
  if (nrow(data) < 3) {
    return(list(estimate = NA, p.value = NA))
  }
  tryCatch({
    cor.test(data[[x_var]], data[[y_var]], method = "pearson")
  }, error = function(e) list(estimate = NA, p.value = NA))
}

get_significance_stars <- function(p_val) {
  if (is.na(p_val)) return("ns")
  if (p_val < 0.001) return("***")
  if (p_val < 0.01) return("**")
  if (p_val < 0.05) return("*")
  return("ns")
}

perform_categorical_test <- function(data) {
  data[, Genotype := factor(Genotype, levels = c(2, 0))]
  data[, Loop_Presence := factor(Loop_Presence, levels = c(0, 1))]
  
  contingency_table <- table(data$Genotype, data$Loop_Presence)
  
  if (sum(contingency_table) < 4 || any(rowSums(contingency_table) == 0) || any(colSums(contingency_table) == 0)) {
    return(list(table = contingency_table, test_name = "Not enough data", result = NULL))
  }
  
  expected_counts <- tryCatch(chisq.test(contingency_table)$expected, warning = function(w) { return(NULL) })
  
  if (is.null(expected_counts) || any(expected_counts < 5)) {
    test_name <- "Fisher's Exact Test"
    test_result <- fisher.test(contingency_table)
  } else {
    test_name <- "Chi-squared Test"
    test_result <- chisq.test(contingency_table, correct = FALSE)
  }
  
  return(list(table = contingency_table, test_name = test_name, result = test_result))
}

# =============================================================================
# 4. 核心绘图函数 (保持不变)
# =============================================================================
create_loop_boxplot <- function(plot_data, title_text, p_val, cor_stats) {
  genotype_levels <- sort(unique(plot_data$Genotype_Allele), decreasing = TRUE)
  plot_data$Genotype_Allele <- factor(plot_data$Genotype_Allele, levels = genotype_levels)
  n_counts <- plot_data[, .N, by=Genotype_Allele]; y_min <- min(plot_data$Loop_Strength, na.rm=T); y_max <- max(plot_data$Loop_Strength, na.rm=T); y_range <- y_max - y_min; if(y_range==0) y_range <- 1
  p <- ggplot(plot_data, aes(x=Genotype_Allele, y=Loop_Strength, fill=Genotype_Allele)) + geom_boxplot(width=0.5, outlier.shape=1, outlier.size=1) + scale_fill_manual(values=CONFIG$NATURE_COLORS[1:2]) + labs(title=title_text, x="Genotype", y="Loop Strength") + create_nature_theme() + geom_text(data=n_counts, aes(x=Genotype_Allele, y=y_min-y_range*0.1, label=paste0("n=", N)), size=2.8, color="gray30") + annotate("text", x=1.5, y=y_max+y_range*0.18, label=sprintf("r=%.3f, P=%.2e", cor_stats$estimate, cor_stats$p.value), size=2.8, fontface="bold")
  if (!is.na(p_val) && length(genotype_levels) == 2) {
    stars <- get_significance_stars(p_val); y_anno <- y_max + y_range * 0.08
    p <- p + geom_segment(aes(x=1, xend=2, y=y_anno, yend=y_anno), color="black", linewidth=0.6) + geom_segment(aes(x=1, xend=1, y=y_anno, yend=y_anno-y_range*0.02), color="black", linewidth=0.6) + geom_segment(aes(x=2, xend=2, y=y_anno, yend=y_anno-y_range*0.02), color="black", linewidth=0.6) + annotate("text", x=1.5, y=y_anno+y_range*0.01, label=stars, size=4, fontface="bold", vjust="bottom")
  }
  return(p)
}

# =============================================================================
# 5. 单个任务处理函数 (保持不变)
# =============================================================================
process_single_loop <- function(current_row, stats_df, loop_strengths_wide, 
                                loop_presence_wide, genotype_df, MAIN_OUTPUT_DIR, CONFIG) {
  
  gene_id <- current_row$Gene_ID
  tad_id <- current_row$TAD_ID
  loop_id <- current_row$Matched_Loop_ID
  
  loop_results_collector <- list()
  
  loop_dir <- file.path(MAIN_OUTPUT_DIR, paste0("Loop_", loop_id))
  dir.create(loop_dir, showWarnings = FALSE, recursive = TRUE)
  
  associated_snps <- unique(stats_df[Gene_ID %in% gene_id & TAD_ID %in% tad_id, Full_SNP_ID])
  if (length(associated_snps) == 0) {
    return(list(data.table(Gene_ID=gene_id, TAD_ID=tad_id, Loop_ID=loop_id, 
                           Full_SNP_ID=NA, Wilcoxon_P_Value_Loop_Strength=NA, 
                           Categorical_P_Value_Loop_Presence=NA, 
                           warning_msg="No associated SNPs found for this Gene-TAD pair")))
  }
  
  if (nrow(loop_strengths_wide) == 0 || nrow(loop_presence_wide) == 0) {
    return(list(data.table(Gene_ID=gene_id, TAD_ID=tad_id, Loop_ID=loop_id, 
                           Full_SNP_ID=NA, Wilcoxon_P_Value_Loop_Strength=NA, 
                           Categorical_P_Value_Loop_Presence=NA, 
                           warning_msg=paste("Loop ID", loop_id, "not found in strength/presence files"))))
  }

  for (snp_full_id in associated_snps) {
    genotypes_wide <- genotype_df[.(snp_full_id)]
    if (nrow(genotypes_wide) == 0) {
      loop_results_collector[[length(loop_results_collector) + 1]] <- 
        data.table(Gene_ID=gene_id, TAD_ID=tad_id, Loop_ID=loop_id, 
                   Full_SNP_ID=snp_full_id, Wilcoxon_P_Value_Loop_Strength=NA, 
                   Categorical_P_Value_Loop_Presence=NA, 
                   warning_msg=paste("SNP ID", snp_full_id, "not found in genotype file"))
      next
    }
    
    safe_snp_id <- gsub(":", "_", snp_full_id)
    base_filename <- paste0(gene_id, "_", tad_id, "_", safe_snp_id)

    p_value_wilcox <- NA
    p_value_categorical <- NA
    
    # Part A: 箱形图绘制
    loop_long_str <- melt(loop_strengths_wide, 
                          id.vars=names(loop_strengths_wide)[1:7], 
                          measure.vars=patterns("^L[0-9]"), 
                          variable.name="Sample", value.name="Loop_Strength")
    geno_long <- melt(genotypes_wide, id.vars="Full_SNP_ID", 
                      measure.vars=patterns("^L[0-9]"), 
                      variable.name="Sample", value.name="Genotype")
    plot_df <- merge(loop_long_str[, .(Sample, Loop_Strength)], 
                     geno_long[, .(Sample, Genotype)], by="Sample")
    plot_df <- plot_df[Genotype %in% CONFIG$EXPECTED_GENOTYPES]
    plot_df <- na.omit(plot_df)
    
    if (nrow(plot_df) >= CONFIG$MIN_SAMPLES && length(unique(plot_df$Genotype)) == 2) {
      cat(sprintf("  - Plotting boxplot for Loop %d & SNP %s\n", loop_id, snp_full_id), file=stderr())
      wilcox_res <- perform_wilcox_test(plot_df, Loop_Strength ~ Genotype)
      cor_res <- calculate_correlation(plot_df, "Genotype", "Loop_Strength")
      p_value_wilcox <- wilcox_res$p.value
    
      alleles <- extract_alleles(snp_full_id)
      plot_df[, Genotype_Allele := convert_genotype_to_allele(Genotype, alleles$ref, alleles$alt)]
      final_plot <- create_loop_boxplot(plot_df, paste("Loop", loop_id, "Strength"), 
                                        wilcox_res$p.value, cor_res)
      suppressMessages(ggsave(file.path(loop_dir, paste0(base_filename, ".png")), 
                              final_plot, width=5, height=8, dpi=300, units="in", bg="white"))
    }
    
    # Part B: 卡方/Fisher检验
    loop_long_pre <- melt(loop_presence_wide, id.vars="Loop_ID", 
                          measure.vars=patterns("^L[0-9]"), 
                          variable.name="Sample", value.name="Loop_Presence")
    presence_df <- merge(loop_long_pre[, .(Sample, Loop_Presence)], 
                         geno_long[, .(Sample, Genotype)], by="Sample")
    presence_df <- presence_df[Genotype %in% CONFIG$EXPECTED_GENOTYPES]
    presence_df <- na.omit(presence_df)
    
    if (nrow(presence_df) > 0 && length(unique(presence_df$Genotype)) == 2) {
      cat(sprintf("  - Performing test for Loop %d & SNP %s\n", loop_id, snp_full_id), file=stderr())
      test_output <- perform_categorical_test(presence_df)
      
      if (!is.null(test_output$result)) {
          p_value_categorical <- test_output$result$p.value
      }
      
      log_path <- file.path(loop_dir, paste0(base_filename, "_chisq_test.log"))
      sink(log_path)
      cat("====================================================\n"); cat(" Loop Presence vs. Genotype Association Test\n"); cat("====================================================\n\n")
      cat("Gene ID:", gene_id, "\n"); cat("TAD ID:", tad_id, "\n"); cat("Loop ID:", loop_id, "\n"); cat("SNP ID:", snp_full_id, "\n\n")
      cat("--- Contingency Table ---\n")
      dimnames(test_output$table) <- list(Genotype = c("REF/REF (2)", "ALT/ALT (0)"), Loop_Status = c("Absent (0)", "Present (1)"))
      print(test_output$table)
      cat("\n--- Test Result ---\n"); cat("Test Used:", test_output$test_name, "\n\n")
      if (!is.null(test_output$result)) { print(test_output$result) }
      sink()
    }
    
    current_result <- data.table(
        Gene_ID = gene_id, TAD_ID = tad_id, Loop_ID = loop_id, Full_SNP_ID = snp_full_id,
        Wilcoxon_P_Value_Loop_Strength = p_value_wilcox,
        Categorical_P_Value_Loop_Presence = p_value_categorical,
        warning_msg = NA_character_
    )
    loop_results_collector[[length(loop_results_collector) + 1]] <- current_result
  }
  return(loop_results_collector)
}

# =============================================================================
# 6. 主逻辑 (【修正后】: 采用任务预处理模式)
# =============================================================================
main <- function() {
  dir.create(MAIN_OUTPUT_DIR, showWarnings = FALSE)
  
  cat("--- 阶段1: 正在加载所有必需的数据文件... ---\n", file=stderr())
  tryCatch({
    loop_strength_df <- fread(PATHS$LOOP_STRENGTH); setkey(loop_strength_df, Loop_ID)
    loop_presence_df <- fread(PATHS$LOOP_PRESENCE); setkey(loop_presence_df, Loop_ID)
    genotype_df <- fread(PATHS$GENOTYPE); setnames(genotype_df, old="SNP_ID", new="Full_SNP_ID"); setkey(genotype_df, Full_SNP_ID)
    stats_df <- fread(PATHS$STATS_BRIDGE)
    loop_analysis_df <- fread(PATHS$LOOP_ANALYSIS)
    valid_loops_df <- loop_analysis_df[!is.na(Matched_Loop_ID)]
    
    cat(sprintf("成功加载 %d 条Loop强度信息。\n", nrow(loop_strength_df)), file=stderr())
    cat(sprintf("成功加载 %d 条Loop存在/不存在信息。\n", nrow(loop_presence_df)), file=stderr())
    cat(sprintf("成功加载 %d 条基因型信息。\n", nrow(genotype_df)), file=stderr())
    cat(sprintf("成功加载 %d 条统计关联信息。\n", nrow(stats_df)), file=stderr())
    cat(sprintf("成功加载主分析文件，筛选出 %d 条有效的TAD-Gene-Loop连接。\n", nrow(valid_loops_df)), file=stderr())
    
  }, error = function(e) {
    stop(paste("错误：加载文件时失败。\n原始错误:", e$message))
  })
  
  # 【修正后】: 阶段2, 在主进程中预先创建任务列表
  cat(sprintf("\n--- 阶段2: 正在准备并行任务列表... ---\n"), file=stderr())
  tasks_list <- list()
  if (nrow(valid_loops_df) > 0) {
    for (i in seq_len(nrow(valid_loops_df))) {
      current_row_data <- valid_loops_df[i, ]
      current_loop_id <- current_row_data$Matched_Loop_ID
      
      # 在主进程中提取数据子集
      strength_subset <- loop_strength_df[.(current_loop_id)]
      presence_subset <- loop_presence_df[.(current_loop_id)]
      
      # 只有当数据子集有效时，才创建任务
      if (nrow(strength_subset) > 0 && nrow(presence_subset) > 0) {
        tasks_list[[length(tasks_list) + 1]] <- list(
          row_data = current_row_data,
          strength = strength_subset,
          presence = presence_subset
        )
      }
    }
  }
  cat(sprintf("已准备 %d 个有效任务进入并行处理。\n", length(tasks_list)), file=stderr())
  
  # 【修正后】: 阶段3, 执行并行计算
  cat(sprintf("\n--- 阶段3: 开始并行处理 (使用 %d CPU核心)... ---\n", NUM_CORES), file=stderr())
  if (length(tasks_list) > 0) {
    cl <- makeCluster(NUM_CORES)
    registerDoParallel(cl)
    
    # 【修正后】: 导出对象列表不再需要包含大型数据框
    export_objects <- c(
        "stats_df", "genotype_df", "MAIN_OUTPUT_DIR", "CONFIG",
        "process_single_loop",
        "create_nature_theme", "extract_alleles", "convert_genotype_to_allele",
        "perform_wilcox_test", "calculate_correlation", "get_significance_stars",
        "perform_categorical_test", "create_loop_boxplot"
    )
    
    results_list <- tryCatch({
      # 【修正后】: foreach 现在迭代的是预处理好的任务列表
      foreach(
        task = tasks_list,
        .packages = c("data.table", "ggplot2", "scales", "dplyr"),
        .export = export_objects,
        .errorhandling = 'pass'
      ) %dopar% {
        # 【修正后】: 调用处理函数，传入任务包中的小数据块
        process_single_loop(
            current_row = task$row_data, 
            stats_df = stats_df,
            loop_strengths_wide = task$strength, 
            loop_presence_wide = task$presence, 
            genotype_df = genotype_df, 
            MAIN_OUTPUT_DIR = MAIN_OUTPUT_DIR, 
            CONFIG = CONFIG
        )
      }
    }, finally = {
      cat("--- 关闭并行环境... ---\n", file=stderr())
      stopCluster(cl)
    })
  } else {
    cat("警告: 没有有效的任务可供处理，跳过并行分析步骤。\n", file=stderr())
    results_list <- list()
  }
  
  cat("\n--- 阶段4: 处理完成，整理结果... ---\n", file=stderr())
  cat(sprintf("所有图表和日志已按Loop ID保存在: %s\n", MAIN_OUTPUT_DIR), file=stderr())
  
  errors <- results_list[sapply(results_list, function(x) "error" %in% class(x))]
  valid_results <- results_list[sapply(results_list, function(x) !"error" %in% class(x))]
  
  if (length(errors) > 0) {
    cat(sprintf("\n--- 检测到 %d 个任务在并行处理中发生严重错误: ---\n", length(errors)), file=stderr())
    for (err in errors) {
      cat(paste0("- ", conditionMessage(err), "\n"), file=stderr())
    }
  }

  results_collector <- unlist(valid_results, recursive = FALSE)
  
  if (length(results_collector) > 0) {
      final_summary_df <- rbindlist(results_collector, fill=TRUE)
      
      warnings_df <- final_summary_df[!is.na(warning_msg), ]
      if(nrow(warnings_df) > 0){
        cat("\n--- 检测到处理警告 (汇总): ---\n", file=stderr())
        warning_summary <- warnings_df[, .N, by=warning_msg]
        for(i in 1:nrow(warning_summary)){
          cat(sprintf("- %s (出现 %d 次)\n", warning_summary$warning_msg[i], warning_summary$N[i]), file=stderr())
        }
      }

      final_summary_df[, warning_msg := NULL]
      final_summary_df <- final_summary_df[!(is.na(Wilcoxon_P_Value_Loop_Strength) & 
                                              is.na(Categorical_P_Value_Loop_Presence))]

      if (nrow(final_summary_df) > 0) {
        tryCatch({
            fwrite(final_summary_df, file = SUMMARY_CSV_PATH, row.names = FALSE)
            cat(sprintf("成功将 %d 条有效统计检验结果汇总到: %s\n", 
                        nrow(final_summary_df), SUMMARY_CSV_PATH), file=stderr())
        }, error = function(e) {
            cat(sprintf("错误：无法写入汇总CSV文件到 %s\n", SUMMARY_CSV_PATH), file=stderr())
            cat(paste("原始错误:", e$message, "\n"), file=stderr())
        })
      } else {
        cat("未生成任何有效的统计结果可供汇总。\n", file=stderr())
      }
  } else {
      cat("未生成任何可供汇总的统计结果。\n", file=stderr())
  }
}

# =============================================================================
# 7. 执行主函数
# =============================================================================
if (!interactive()) {
  main()
}