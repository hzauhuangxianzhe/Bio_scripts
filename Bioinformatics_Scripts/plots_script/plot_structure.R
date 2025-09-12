# Structure群体结构分析结果可视化脚本（命令行版）
# 功能：读取Structure输出文件并绘制堆叠柱状图
# 特点：命令行参数、自定义输出、简化图表、按群体排序
# 修正：确保概率最大的部分在柱子底部，移除所有网格线，零误差容忍

# 加载必要的R包
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales) # 用于y轴格式
library(optparse) # 用于解析命令行参数

# 1. 命令行参数解析
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Structure input file (K*.outfile)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="structure_plot.pdf", 
              help="Output PDF file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# 检查输入文件是否已指定
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Input file must be specified. (--input FILE)", call.=FALSE)
}

file_path <- opt$input
output_file <- opt$output

# 2. 数据读取和解析函数（零误差版本）
parse_structure_file <- function(file_path) {
  cat("正在读取Structure输出文件:", file_path, "\n")
  
  lines <- readLines(file_path)
  cat("文件包含", length(lines), "行数据\n")
  
  result_list <- list()
  
  cat("正在解析数据...\n")
  for (i in 1:length(lines)) {
    line <- trimws(lines[i])
    
    if (line == "") next
    
    colon_pos <- regexpr(":", line)
    if (colon_pos == -1) {
      cat("警告：第", i, "行格式异常，跳过\n")
      next
    }
    
    before_colon <- substr(line, 1, colon_pos - 1)
    after_colon <- substr(line, colon_pos + 1, nchar(line))
    
    before_parts <- strsplit(trimws(before_colon), "\\s+")[[1]]
    if (length(before_parts) == 0) next
    
    individual_id <- as.numeric(before_parts[1])
    if (is.na(individual_id)) next
    
    prob_values <- as.numeric(strsplit(trimws(after_colon), "\\s+")[[1]])
    if (any(is.na(prob_values))) {
      cat("警告：个体", individual_id, "概率数据异常\n")
      next
    }
    
    original_sum <- sum(prob_values)
    if (abs(original_sum - 1.0) > .Machine$double.eps) {
      cat("个体", individual_id, "概率和:", sprintf("%.15f", original_sum), "，进行高精度标准化\n")
    }
    
    # 高精度标准化：使用更严格的数值处理
    standardized_probs <- prob_values / original_sum
    
    # 使用高精度计算确保严格等于1.0
    current_sum <- sum(standardized_probs)
    diff <- 1.0 - current_sum
    
    # 将剩余差值分配给最大概率的成分
    max_prob_index <- which.max(standardized_probs)
    standardized_probs[max_prob_index] <- standardized_probs[max_prob_index] + diff
    
    # 最终验证：必须严格等于1.0，零误差容忍
    new_sum <- sum(standardized_probs)
    if (abs(new_sum - 1.0) > .Machine$double.eps) {
      cat("致命错误：个体", individual_id, "标准化后和不等于1.0，值为", sprintf("%.15f", new_sum), "\n")
      cat("差值:", sprintf("%.15e", new_sum - 1.0), "\n")
      stop("高精度标准化失败，程序终止")
    }
    
    # 强制最终修正：确保绝对精确
    if (new_sum != 1.0) {
      final_diff <- 1.0 - new_sum
      standardized_probs[max_prob_index] <- standardized_probs[max_prob_index] + final_diff
      
      # 终极验证
      final_sum <- sum(standardized_probs)
      if (final_sum != 1.0) {
        stop(paste("无法达到严格的1.0标准，最终和为:", sprintf("%.15f", final_sum)))
      }
    }
    
    cluster_names <- paste0("Cluster", 1:length(standardized_probs))
    new_row <- data.frame(
      Individual = individual_id
    )
    new_row[cluster_names] <- as.list(standardized_probs)
    result_list[[length(result_list) + 1]] <- new_row
  }
  
  if (length(result_list) == 0) {
    stop("没有成功解析任何数据")
  }
  
  parsed_df <- do.call(rbind, result_list)
  cat("成功解析", nrow(parsed_df), "个个体的数据\n")
  
  return(parsed_df)
}

# 3. 数据转换函数（零误差版本）
convert_to_long_format <- function(wide_data) {
  cat("正在转换数据格式...\n")

  long_data <- wide_data %>%
    pivot_longer(
      cols = starts_with("Cluster"),
      names_to = "Cluster",
      values_to = "Probability"
    ) %>%
    arrange(Individual, Cluster)

  prob_check <- long_data %>%
    group_by(Individual) %>%
    summarise(total_prob = sum(Probability), .groups = 'drop')

  # 严格验证：零误差容忍
  tolerance <- .Machine$double.eps
  prob_deviations <- abs(prob_check$total_prob - 1.0)
  max_deviation <- max(prob_deviations)
  
  if (max_deviation > tolerance) {
    cat("数据转换后发现概率和异常，最大偏差:", sprintf("%.15e", max_deviation), "\n")
    problem_individuals <- prob_check$Individual[prob_deviations > tolerance]
    cat("问题个体:", paste(problem_individuals, collapse = ", "), "\n")
    stop("概率和验证失败：存在不等于1.0的个体")
  } else {
    cat("数据转换完成，所有概率和均严格等于1.0（最大偏差:", sprintf("%.2e", max_deviation), ")\n")
  }

  return(long_data)
}

# 4. 群体分配函数
assign_main_clusters <- function(long_data) {
  cat("正在确定主要群体归属...\n")
  
  main_clusters <- long_data %>%
    group_by(Individual) %>%
    slice_max(Probability, n = 1, with_ties = FALSE) %>%
    select(Individual, MainCluster = Cluster) %>%
    ungroup()
  
  cluster_stats <- main_clusters %>%
    count(MainCluster, name = "Size") %>%
    arrange(desc(Size))
  
  return(list(assignments = main_clusters, stats = cluster_stats))
}

# 5. 个体排序函数
order_individuals <- function(long_data, main_assignments, cluster_stats) {
  cat("正在排序个体...\n")
  
  cluster_priority <- cluster_stats$MainCluster
  
  main_probs <- long_data %>%
    inner_join(main_assignments, by = "Individual") %>%
    filter(Cluster == MainCluster) %>%
    select(Individual, MainCluster, MainProbability = Probability)
  
  ordered_individuals <- main_probs %>%
    mutate(MainCluster = factor(MainCluster, levels = cluster_priority)) %>%
    arrange(MainCluster, desc(MainProbability)) %>%
    pull(Individual)
  
  cat("个体排序完成：", length(ordered_individuals), "个个体\n")
  return(ordered_individuals)
}

# 6. 绘图数据准备函数（零误差版本）
# 在个体级别计算堆叠区间（ymin/ymax），并为不同主要群体之间插入空隙
prepare_plot_data <- function(long_data, individual_order) {
  cat("正在准备绘图数据...\n")
  
  # 最终概率验证
  final_check <- long_data %>%
    group_by(Individual) %>%
    summarise(total = sum(Probability), .groups = 'drop')
  
  max_error <- max(abs(final_check$total - 1.0))
  if (max_error > .Machine$double.eps) {
    stop(paste("绘图数据准备失败：概率和不等于1.0，最大误差:", sprintf("%.15e", max_error)))
  }
  
  # 计算每个个体的主要群体（与 assign_main_clusters 同逻辑，局部复用）
  main_assign <- long_data %>%
    group_by(Individual) %>%
    slice_max(Probability, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(Individual, MainCluster = Cluster) %>%
    mutate(Individual = factor(Individual, levels = individual_order))
  
  # 个体顺序（因子），并在"每个个体内部"按概率从大到小排序，计算手工堆叠区间
  base <- long_data %>%
    mutate(Individual = factor(Individual, levels = individual_order)) %>%
    left_join(main_assign, by = "Individual") %>%
    arrange(Individual) %>%
    group_by(Individual) %>%
    arrange(desc(Probability), .by_group = TRUE) %>%
    mutate(
      ymax = cumsum(Probability),
      ymin = ymax - Probability
    ) %>%
    ungroup()
  
  # 为不同主要群体之间添加 x 轴空隙
  gap_size <- 3
  
  ind_main <- base %>%
    distinct(Individual, MainCluster) %>%
    arrange(Individual)
  
  group_change <- c(0, as.integer(ind_main$MainCluster[-1] != ind_main$MainCluster[-nrow(ind_main)]))
  gap_cumsum <- cumsum(group_change)
  x_pos <- seq_len(nrow(ind_main)) + gap_size * gap_cumsum
  
  ind_x_map <- ind_main %>%
    mutate(x_idx = x_pos)
  
  plot_data <- base %>%
    left_join(ind_x_map %>% select(Individual, x_idx), by = "Individual")
  
  cat("绘图数据准备完成（已插入群体间隙），概率和验证通过（零误差容忍）\n")
  return(plot_data)
}

# 7. 绘图函数（完全无网格线版本）
# 使用 geom_rect 根据 ymin/ymax 手工堆叠，保证个体级别最大块在底部
create_structure_plot <- function(plot_data) {
  cat("正在生成堆叠柱状图（无网格线版本）...\n")
  
  num_clusters <- length(unique(plot_data$Cluster))
  cluster_colors <- scales::hue_pal()(num_clusters)
  xmax_val <- max(plot_data$x_idx)
  
  p <- ggplot(plot_data) +
    geom_rect(
      aes(xmin = x_idx - 0.5, xmax = x_idx + 0.5,
          ymin = ymin, ymax = ymax, fill = Cluster),
      color = "white", linewidth = 0.1
    ) +
    scale_fill_manual(values = cluster_colors, name = "Population") +
    labs(
      x = "Individual",
      y = "Ancestry Proportion"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      # 移除所有网格线
      panel.grid.major = element_blank(),    # 移除所有主要网格线（包括水平和垂直）
      panel.grid.minor = element_blank(),    # 移除所有次要网格线
      panel.grid = element_blank(),          # 额外保险：移除所有网格
      plot.title = element_blank(), 
      plot.subtitle = element_blank(), 
      plot.caption = element_blank(), 
      legend.title = element_text(size = 12),
      axis.title = element_text(size = 12)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),  # 保持y轴刻度但不显示网格线
      expand = c(0, 0)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.5, xmax_val + 0.5), ylim = c(0, 1))
    
  return(p)
}

# 8. 主执行流程
main <- function() {
  cat("=== Structure群体结构分析可视化程序（无网格线版本）===\n\n")
  
  # 步骤1：解析数据文件
  wide_data <- parse_structure_file(file_path)
  
  # 2. 数据转换和排序
  long_data <- convert_to_long_format(wide_data)
  cluster_info <- assign_main_clusters(long_data)
  individual_order <- order_individuals(long_data, cluster_info$assignments, cluster_info$stats)
  
  # 3. 准备绘图数据（个体级别手工堆叠）
  plot_data <- prepare_plot_data(long_data, individual_order)
  
  # 4. 创建并保存图形
  p <- create_structure_plot(plot_data)
  
  # 保存图形（只保存PDF）
  ggsave(output_file, plot = p, width = 6, height = 4)
  
  cat("\n=== 处理完成 ===\n")
  cat("图形已保存为:", output_file, "\n")
  
  return(list(
    data = plot_data,
    assignments = cluster_info$assignments,
    stats = cluster_info$stats
  ))
}

# 执行主程序
result <- main()