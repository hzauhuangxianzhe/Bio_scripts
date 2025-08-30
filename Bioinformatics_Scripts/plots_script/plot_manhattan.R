# 棉花多组学曼哈顿图绘制脚本
# 作者: 基因组分析脚本
# 日期: 2025-08-30

# 加载必需的R包
library(data.table)
library(ggplot2)
library(gridExtra)
library(dplyr)

# 设置工作目录和定义文件路径
cat("开始处理棉花多组学数据...\n")

# 定义表型列表，只处理FL和FS
traits <- c("FL", "FS")

# 读取基因组信息文件
cat("读取染色体长度信息...\n")
genome_info <- fread("l128.genome.lst", header = FALSE, col.names = c("chr", "length"))

# 读取位置信息文件
cat("读取分子表型位置信息...\n")
gene_pos <- fread("final_gene_pos_info.tsv")
loop_pos <- fread("final_loop_pos_info.tsv")
tad_pos <- fread("final_tad_pos_info.tsv")

# 数字染色体到HC04格式的转换函数
convert_chr_format <- function(chr_num) {
  chr_num <- as.numeric(chr_num)
  if (chr_num <= 13) {
    return(paste0("HC04_A", sprintf("%02d", chr_num)))
  } else {
    return(paste0("HC04_D", sprintf("%02d", chr_num - 13)))
  }
}

# 处理每个表型的数据
all_data <- list()

for (trait_name in traits) {
  cat(paste("处理表型", trait_name, "的数据...\n"))
  
  # 1. 处理GWAS数据
  cat(paste("  处理", trait_name, "的GWAS数据...\n"))
  gwas_file <- paste0("gwas_", trait_name, ".tsv")
  gwas_data <- fread(gwas_file)
  
  # 提取SNP的染色体和位置信息
  gwas_data[, c("chr_num", "pos") := tstrsplit(SNP, "_", fixed = TRUE)]
  gwas_data$chr_num <- as.numeric(gwas_data$chr_num)
  gwas_data$pos <- as.numeric(gwas_data$pos)
  
  # 转换染色体格式
  gwas_data$chr <- sapply(gwas_data$chr_num, convert_chr_format)
  
  # 整理GWAS数据格式
  gwas_final <- data.table(
    phe_id = gwas_data$SNP,
    chr = gwas_data$chr,
    pos = gwas_data$pos,
    pval = gwas_data$Pvalue,
    data_type = "GWAS",
    trait = trait_name
  )
  
  # 2. 处理TWAS数据
  cat(paste("  处理", trait_name, "的TWAS数据...\n"))
  twas_file <- paste0("gene_raw_TPM_1Mb_rn_", trait_name, "_merged.tsv")
  twas_data <- fread(twas_file)
  
  # 根据基因ID匹配位置信息
  twas_data <- merge(twas_data, gene_pos, by.x = "ID", by.y = "GeneId", all.x = TRUE)
  
  # 使用基因的start位置作为绘图位置
  twas_final <- data.table(
    phe_id = twas_data$ID,
    chr = twas_data$CHR,
    pos = twas_data$start,
    pval = twas_data$TWAS.P,
    data_type = "TWAS",
    trait = trait_name
  )
  
  # 3. 处理TCIWAS数据 - Loop
  cat(paste("  处理", trait_name, "的TCIWAS Loop数据...\n"))
  loop_file <- paste0("loop_TMMed_40kb_rn_", trait_name, "_merged.tsv")
  loop_data <- fread(loop_file)
  
  # 根据ID匹配位置信息
  loop_data <- merge(loop_data, loop_pos, by.x = "ID", by.y = "Loop_ID", all.x = TRUE)
  
  # 计算Loop位置: (start1 + end2) / 2
  loop_data$pos <- (loop_data$start1 + loop_data$end2) / 2
  
  loop_final <- data.table(
    phe_id = loop_data$ID,
    chr = loop_data$CHR,
    pos = loop_data$pos,
    pval = loop_data$TWAS.P,
    data_type = "TCIWAS_Loop",
    trait = trait_name
  )
  
  # 4. 处理TCIWAS数据 - TAD
  cat(paste("  处理", trait_name, "的TCIWAS TAD数据...\n"))
  tad_file <- paste0("tad_final_250kb_rn_", trait_name, "_merged.tsv")
  tad_data <- fread(tad_file)
  
  # 根据ID匹配位置信息
  tad_data <- merge(tad_data, tad_pos, by.x = "ID", by.y = "TAD_ID", all.x = TRUE)
  
  # 计算TAD位置: (start + end) / 2
  tad_data$pos <- (tad_data$start + tad_data$end) / 2
  
  tad_final <- data.table(
    phe_id = tad_data$ID,
    chr = tad_data$CHR,
    pos = tad_data$pos,
    pval = tad_data$TWAS.P,
    data_type = "TCIWAS_TAD",
    trait = trait_name
  )
  
  # 合并该表型的所有数据
  trait_data <- rbind(gwas_final, twas_final, loop_final, tad_final, fill = TRUE)
  
  # 移除缺失值
  trait_data <- trait_data[!is.na(chr) & !is.na(pos) & !is.na(pval) & pval > 0]
  
  all_data[[trait_name]] <- trait_data
}

# 将所有表型数据合并到一张表，以便统一处理
combined_data <- rbindlist(all_data, fill = TRUE)

# 对TWAS和TCIWAS数据进行FDR校正
cat("对TWAS和TCIWAS数据进行FDR校正...\n")
combined_data[, pval_fdr := pval] # 默认使用原始p值

# TWAS: 对每个表型单独进行FDR校正
combined_data[data_type == "TWAS", pval_fdr := p.adjust(pval, method = "BH"), by = trait]

# TCIWAS_Loop: 对每个表型单独进行FDR校正
combined_data[data_type == "TCIWAS_Loop", pval_fdr := p.adjust(pval, method = "BH"), by = trait]

# TCIWAS_TAD: 对每个表型单独进行FDR校正
combined_data[data_type == "TCIWAS_TAD", pval_fdr := p.adjust(pval, method = "BH"), by = trait]

# 准备绘制曼哈顿图
cat("准备绘制曼哈顿图...\n")

# 创建染色体顺序和累积长度
chr_order <- c(paste0("HC04_A", sprintf("%02d", 1:13)), 
               paste0("HC04_D", sprintf("%02d", 1:13)))

genome_info$chr <- factor(genome_info$chr, levels = chr_order)
genome_info <- genome_info[order(genome_info$chr)]

# 计算累积位置
genome_info$cumulative_start <- cumsum(c(0, genome_info$length[-nrow(genome_info)]))
genome_info$cumulative_mid <- genome_info$cumulative_start + genome_info$length / 2

# 绘制GWAS曼哈顿图的函数
create_gwas_manhattan <- function(data, threshold_line = NULL) {
  
  # 筛选GWAS数据
  plot_data <- data[data_type == "GWAS"]
  if (nrow(plot_data) == 0) return(NULL)
  
  # 添加累积位置
  plot_data <- merge(plot_data, genome_info[, .(chr, cumulative_start)], by = "chr")
  plot_data$cumulative_pos <- plot_data$cumulative_start + plot_data$pos
  
  # 计算-log10(p值)
  plot_data$neg_log10_p <- -log10(plot_data$pval)
  plot_data <- plot_data[is.finite(neg_log10_p)]
  
  # 设置显著性颜色
  plot_data$Color <- "gray"
  plot_data[pval < threshold_line & trait == "FS", Color := "yellow"]
  plot_data[pval < threshold_line & trait == "FL", Color := "orange"]
  
  # 绘制
  p <- ggplot(plot_data, aes(x = cumulative_pos, y = neg_log10_p, color = Color)) +
    geom_point(alpha = 0.8, size = 1.3) +
    scale_color_manual(
      values = c("gray" = "gray", "yellow" = "yellow", "orange" = "orange"),
      breaks = c("yellow", "orange", "gray"),
      labels = c("FS significant", "FL significant", "Not significant"),
      name = "GWAS"
    ) +
    scale_x_continuous(
      breaks = genome_info$cumulative_mid,
      labels = gsub("HC04_", "", genome_info$chr),
      expand = c(0.01, 0.03)
    ) +
    labs(x = "Chromosome", y = "-log10(P)") +
    geom_hline(yintercept = -log10(threshold_line), color = "black", linetype = "dashed") +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = NA, color = NA),
      plot.background = element_rect(fill = NA, color = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.title = element_text(hjust = 0.5),
      legend.justification = c(1, 1),
      legend.key.height = unit(0.3, "cm"),
      legend.background = element_rect(fill = NA),
      plot.margin = unit(c(1, 0.5, 0.2, 0.2), "cm"),
      legend.position = "right"
    )
  
  return(p)
}

# 绘制TWAS和TCIWAS曼哈顿图的函数
create_multiomics_manhattan <- function(data, threshold_line = 0.05) {
  
  # 筛选TWAS和TCIWAS数据
  plot_data <- data[data_type %in% c("TWAS", "TCIWAS_Loop", "TCIWAS_TAD")]
  if (nrow(plot_data) == 0) return(NULL)
  
  # 添加累积位置
  plot_data <- merge(plot_data, genome_info[, .(chr, cumulative_start)], by = "chr")
  plot_data$cumulative_pos <- plot_data$cumulative_start + plot_data$pos
  
  # 计算-log10(p值)
  plot_data$neg_log10_p <- -log10(plot_data$pval_fdr)
  plot_data <- plot_data[is.finite(neg_log10_p)]
  
  # 设置显著性颜色
  plot_data$Color <- "gray"
  plot_data[pval_fdr < threshold_line & data_type == "TWAS" & trait == "FS", Color := "blue"]
  plot_data[pval_fdr < threshold_line & data_type == "TWAS" & trait == "FL", Color := "red"]
  plot_data[pval_fdr < threshold_line & data_type == "TCIWAS_Loop" & trait == "FS", Color := "brown"]
  plot_data[pval_fdr < threshold_line & data_type == "TCIWAS_TAD" & trait == "FS", Color := "green"]
  plot_data[pval_fdr < threshold_line & data_type == "TCIWAS_Loop" & trait == "FL", Color := "darkmagenta"] # 深紫色
  plot_data[pval_fdr < threshold_line & data_type == "TCIWAS_TAD" & trait == "FL", Color := "cyan"] # 青色
  
  # 绘制
  p <- ggplot(plot_data, aes(x = cumulative_pos, y = neg_log10_p, color = Color)) +
    geom_point(alpha = 0.8, size = 1.3) +
    scale_color_manual(
      values = c("gray" = "gray", "blue" = "blue", "red" = "red", "brown" = "brown", "green" = "green", "darkmagenta" = "darkmagenta", "cyan" = "cyan"),
      breaks = c("blue", "red", "brown", "green", "darkmagenta", "cyan", "gray"),
      labels = c("FS-TWAS", "FL-TWAS", "FS-Loop", "FS-TAD", "FL-Loop", "FL-TAD", "Not significant"),
      name = "Multi-omics"
    ) +
    scale_x_continuous(
      breaks = genome_info$cumulative_mid,
      labels = gsub("HC04_", "", genome_info$chr),
      expand = c(0.01, 0.03)
    ) +
    labs(x = "Chromosome", y = "-log10(FDR)") +
    geom_hline(yintercept = -log10(threshold_line), color = "black", linetype = "dashed") +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = NA, color = NA),
      plot.background = element_rect(fill = NA, color = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.title = element_text(hjust = 0.5),
      legend.justification = c(1, 1),
      legend.key.height = unit(0.3, "cm"),
      legend.background = element_rect(fill = NA),
      plot.margin = unit(c(1, 0.5, 0.2, 0.2), "cm"),
      legend.position = "right"
    )
  
  return(p)
}

# 为FL和FS绘制曼哈顿图
for (trait_name in traits) {
  cat(paste("绘制表型", trait_name, "的曼哈顿图...\n"))
  
  # 这里我们不需要按表型单独筛选数据，因为绘图函数处理的是 combined_data
  # 计算GWAS阈值 (1/n)
  gwas_n <- nrow(combined_data[data_type == "GWAS"])
  gwas_threshold <- 1 / gwas_n
  
  # 创建GWAS曼哈顿图
  gwas_plot <- create_gwas_manhattan(
    combined_data, 
    gwas_threshold
  )
  
  # 创建TWAS和TCIWAS曼哈顿图
  multi_plot <- create_multiomics_manhattan(
    combined_data, 
    0.05
  )
  
  # 合并两个图
  if (!is.null(gwas_plot) && !is.null(multi_plot)) {
    combined_plot <- grid.arrange(gwas_plot, multi_plot, ncol = 1, heights = c(1, 1), top = "Manhattan Plot for FL and FS")
    
    # 保存图像
    output_file <- paste0("manhattan_combined_FL_FS.png")
    ggsave(output_file, combined_plot, width = 12, height = 8, dpi = 600)
    cat(paste("  保存图像:", output_file, "\n"))
  }
  
  # 由于我们绘制的是合并图，循环一次即可
  break
}

cat("所有曼哈顿图绘制完成！\n")

# 输出数据统计摘要
cat("\n数据处理摘要:\n")
for (trait_name in traits) {
  cat(paste("表型", trait_name, ":\n"))
  data <- combined_data[trait == trait_name]
  summary_stats <- data[, .N, by = data_type]
  print(summary_stats)
  
  # 打印显著性统计
  cat("\n显著性结果统计 (FDR < 0.05):\n")
  sig_stats <- data[pval_fdr < 0.05, .N, by = data_type]
  print(sig_stats)
  
  cat("\n")
}
