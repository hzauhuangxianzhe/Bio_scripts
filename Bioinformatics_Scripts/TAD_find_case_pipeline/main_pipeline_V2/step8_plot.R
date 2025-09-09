#!/usr/bin/env Rscript

# 脚本功能：为每个通过的三元组/四元组绘制六种图（TAD IS、基因表达、FL表型、FS表型、Loop强度、Loop有无 vs 基因型），
# 并将它们合并到一张图中 (2x3 布局)，同时标记预计算的P值。
#
# 版本 3.0 (修改版): 
# - (新) Loop状态图 (geom_tile) 现在是一个完整的 2x2 列联表热图，显示原始计数和列百分比。
# - (新) Loop强度图和Loop状态图的标题现在动态包含 Loop ID。
#
# 该脚本依赖于以下 R 包：
# - ggplot2: 用于生成高质量的统计图形。
# - tidyr: 用于数据重塑，特别是将宽格式转换为长格式。
# - dplyr: 用于数据操作和处理。
# - optparse: 用于处理命令行参数。
# - patchwork: 用于将多个 ggplot 图组合到单个图像中。
#

# 加载必要的库
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(optparse)
  library(patchwork) # 用于组合多张图
})

# --- 1. 定义和解析命令行参数 ---
option_list = list(
  make_option("--triplets", type = "character", default = NULL,
              help = "passed_triplets_results.csv file path, containing Loop-related P-values", metavar = "character"),
  make_option("--p_values", type = "character", default = NULL,
              help = "all_triplet_statistics.csv file path, containing U-test p-values for TAD, Gene, FL, and FS", metavar = "character"),
  make_option("--genotype", type = "character", default = NULL,
              help = "snp_matrix.txt file path, containing genotype information", metavar = "character"),
  make_option("--gene_exp", type = "character", default = NULL,
              help = "final_gene_TPM.tsv file path, containing gene expression levels", metavar = "character"),
  make_option("--tad_is", type = "character", default = NULL,
              help = "267sample_final_tad_is.tsv file path, containing TAD IS scores", metavar = "character"),
  make_option("--phenotypes", type = "character", default = NULL,
              help = "phenotypes_transposed.tsv file path, containing FL and FS phenotypes", metavar = "character"),
  make_option("--loop_strength", type = "character", default = NULL,
              help = "final_loop_strength.tsv file path, containing loop strength", metavar = "character"),
  make_option("--loop_status", type = "character", default = NULL,
              help = "loop_status.tsv file path, containing loop presence/absence status (0/1)", metavar = "character"),
  make_option("--output_dir", type = "character", default = NULL,
              help = "Directory to save output plots", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# 检查所有必需的文件和目录参数是否已提供
if (is.null(opt$triplets) || is.null(opt$p_values) || is.null(opt$genotype) || is.null(opt$gene_exp) ||
    is.null(opt$tad_is) || is.null(opt$phenotypes) || is.null(opt$loop_strength) || is.null(opt$loop_status) ||
    is.null(opt$output_dir)) {
  stop("所有文件和输出目录路径参数都必须提供！", call. = FALSE)
}

# 如果输出目录不存在则创建
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

# --- 2. 读取和预处理数据 ---

# 读取主三元组数据文件 (包含Loop相关的P值)
triplets_data <- read.csv(opt$triplets, check.names = FALSE)

# 读取包含U检验P值的新文件 (来自 Step 5 的 all_...csv)
p_value_data <- read.csv(opt$p_values, check.names = FALSE) %>%
  select(Gene_ID, TAD_ID, SNP_ID, Full_SNP_ID, MannWhitney_P_TAD, MannWhitney_P_Gene, MannWhitney_P_FL, MannWhitney_P_FS)

# 将TAD, Gene, FL, FS的P值数据合并到主三元组数据框中
triplets_data <- left_join(triplets_data, p_value_data, by = c("Gene_ID", "TAD_ID", "SNP_ID", "Full_SNP_ID"))

# 读取基因型数据并预处理
genotype_data <- read.table(opt$genotype, header = TRUE, sep = "\t", check.names = FALSE)
genotype_long <- genotype_data %>%
  pivot_longer(cols = -SNP_ID, names_to = "Sample_ID", values_to = "Genotype") %>%
  filter(Genotype %in% c(0, 2))

# 读取基因表达数据并转换为长格式
gene_exp_long <- read.table(opt$gene_exp, header = TRUE, sep = "\t", check.names = FALSE) %>%
  pivot_longer(cols = -Gene_ID, names_to = "Sample_ID", values_to = "Value")

# 读取TAD IS数据并转换为长格式
tad_is_long <- read.table(opt$tad_is, header = TRUE, sep = "\t", check.names = FALSE) %>%
  pivot_longer(cols = -TAD_ID, names_to = "Sample_ID", values_to = "Value")

# 读取表型数据并转换为长格式
pheno_long <- read.table(opt$phenotypes, header = TRUE, sep = "\t", check.names = FALSE) %>%
  pivot_longer(cols = -trait, names_to = "Sample_ID", values_to = "Value")

# 读取Loop强度数据并转换为长格式
loop_strength_long <- read.table(opt$loop_strength, header = TRUE, sep = "\t", check.names = FALSE) %>%
  pivot_longer(cols = -Loop_ID, names_to = "Sample_ID", values_to = "Value")

# 读取Loop有无数据并转换为长格式
loop_status_long <- read.table(opt$loop_status, header = TRUE, sep = "\t", check.names = FALSE) %>%
  pivot_longer(cols = -Loop_ID, names_to = "Sample_ID", values_to = "Status")

# --- 3. 定义绘图函数 (v3.0 MODIFIED) ---

# 定义一个返回ggplot对象的通用绘图函数（箱线图）
create_boxplot <- function(data, y_label, title_label, p_value, ref_base, alt_base) {
  
  # 调整X轴标签，将数字基因型映射
  data$Genotype_Label <- factor(data$Genotype, levels = c(2, 0), 
                                labels = c(paste0(ref_base, ref_base, "(2)"), paste0(alt_base, alt_base, "(0)")))
  
  # 绘制箱线图
  p <- ggplot(data, aes(x = Genotype_Label, y = Value, fill = Genotype_Label)) +
    geom_boxplot(outlier.shape = NA) + # 隐藏离群点
    geom_jitter(shape = 16, position = position_jitter(0.05), alpha = 0.6, size = 1) + # 增加散点
    labs(
      title = title_label,
      x = "Genotype",
      y = y_label,
      subtitle = sprintf("MWU P = %.2e", p_value) # 明确P值来源
    ) +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), # 缩小标题字体以适应更长的Loop ID
      plot.subtitle = element_text(hjust = 0.5, size = 9, face = "italic"), 
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      legend.position = "none",
      plot.margin = unit(c(5, 5, 5, 5), "mm")
    )

  return(p)
}


# *** (v3.0 MODIFIED: 重写为 2x2 列联表热图) ***
# 定义一个用于绘制Loop有无的 2x2 列联表热图函数
create_contingency_heatmap <- function(data, title_label, p_value, ref_base, alt_base) {

  # 1. 转换Status和Genotype为因子，用于正确的标签和排序
  data_to_plot <- data %>%
    mutate(
      Status_Label = factor(Status, levels = c(0, 1), labels = c("Absent (0)", "Present (1)")),
      Genotype_Label = factor(Genotype, levels = c(2, 0), 
                              labels = c(paste0(ref_base, ref_base, "(2)"), paste0(alt_base, alt_base, "(0)")))
    )

  # 2. 计算 2x2 表格中所有4个单元格的计数(N)，并使用 complete() 填充0
  count_data <- data_to_plot %>%
    group_by(Genotype_Label, Status_Label) %>%
    summarise(N = n(), .groups = 'drop') %>%
    tidyr::complete(Genotype_Label, Status_Label, fill = list(N = 0))
  
  # 3. 计算列百分比 (即在同一基因型中，0和1的比例)
  count_data <- count_data %>%
    group_by(Genotype_Label) %>%
    mutate(
      Col_Total = sum(N),
      Col_Prop = ifelse(Col_Total == 0, 0, N / Col_Total) # 避免除以零
    ) %>%
    ungroup()
  
  # 4. 创建文本标签 (显示 原始计数 和 列百分比)
  count_data$Label <- sprintf("%d\n(%.1f%%)", count_data$N, count_data$Col_Prop * 100)

  # 5. 格式化P值用于副标题
  formatted_p <- sprintf("Fisher/Chi2 P = %.4e", p_value)

  # 6. 绘制Tile图 (热图)
  p <- ggplot(count_data, aes(x = Genotype_Label, y = Status_Label, fill = N)) +
    geom_tile(color = "black", size = 0.5) + # 绘制4个色块，并添加黑色边框
    geom_text(aes(label = Label), color = "white", size = 3.5, fontface = "bold") + # 在色块内添加计数和百分比
    scale_fill_gradient(low = "#56B1F7", high = "#132B43", name = "Cell Count") + # 蓝色渐变
    labs(
      title = title_label,
      subtitle = formatted_p, # 将P值作为副标题
      x = "Genotype",
      y = "Loop Status" # 恢复Y轴标签
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), # 缩小标题字体
      plot.subtitle = element_text(hjust = 0.5, size = 9, face = "italic"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9), # 恢复Y轴文本
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = unit(c(5, 5, 5, 5), "mm")
    )
  
  return(p)
}


# --- 4. 遍历三元组并绘图 (v3.0 MODIFIED: 更新标题和函数调用) ---
for (i in 1:nrow(triplets_data)) {
  
  # 获取当前三元组信息
  row <- triplets_data[i,]
  gene_id <- as.character(row$Gene_ID)
  tad_id <- as.character(row$TAD_ID)
  snp_id <- as.character(row$SNP_ID)
  full_snp_id <- as.character(row$Full_SNP_ID)
  loop_id <- as.character(row$Loop_ID)
  
  # 获取所有相关的P值
  p_val_tad <- as.numeric(row$MannWhitney_P_TAD)
  p_val_gene <- as.numeric(row$MannWhitney_P_Gene)
  p_val_fl <- as.numeric(row$MannWhitney_P_FL)
  p_val_fs <- as.numeric(row$MannWhitney_P_FS)
  p_val_loop_strength <- as.numeric(row$Loop_Strength_P)
  p_val_loop_status <- as.numeric(row$Loop_01_P)
  
  # 确保 full_snp_id 有值且格式正确
  if (!is.na(full_snp_id) && full_snp_id != "") {
    
    # (v3.0 假设: Full_SNP_ID 格式为 ID:REF:ALT 或 CHR:POS:REF:ALT)
    # 修复了v2.0中对parts长度的错误假设，现在它取最后两个元素作为REF和ALT
    parts <- strsplit(full_snp_id, ":")[[1]]
    if (length(parts) >= 3) {
      ref_base <- parts[length(parts) - 1] # 倒数第二个元素
      alt_base <- parts[length(parts)]     # 最后一个元素
    } else {
      warning(paste("Invalid Full_SNP_ID format for triplet, skipping:", full_snp_id))
      next
    }
    
    # 根据SNP ID筛选出对应的基因型数据
    current_genotype <- genotype_long %>%
      filter(SNP_ID == full_snp_id)
    
    if (nrow(current_genotype) == 0) {
      warning(paste("No genotype data found for SNP:", full_snp_id, ", skipping triplet."))
      next
    }
    
    plot_list <- list()
    
    # --- 绘制 TAD IS 箱线图 ---
    current_tad_is <- tad_is_long %>% filter(TAD_ID == tad_id) %>% inner_join(current_genotype, by = "Sample_ID")
    if (nrow(current_tad_is) > 0) plot_list$p_tad_is <- create_boxplot(current_tad_is, "TAD IS Score", "TAD IS", p_val_tad, ref_base, alt_base)
    
    # --- 绘制 FS 表型箱线图 ---
    current_fs_pheno <- pheno_long %>% filter(trait == "FS") %>% inner_join(current_genotype, by = "Sample_ID")
    if (nrow(current_fs_pheno) > 0) plot_list$p_fs_pheno <- create_boxplot(current_fs_pheno, "FS Phenotype", "FS Phenotype", p_val_fs, ref_base, alt_base)
    
    # --- 绘制 基因表达 箱线图 ---
    current_gene_exp <- gene_exp_long %>% filter(Gene_ID == gene_id) %>% inner_join(current_genotype, by = "Sample_ID")
    if (nrow(current_gene_exp) > 0) plot_list$p_gene_exp <- create_boxplot(current_gene_exp, "Gene Expression (TPM)", "Gene Expression", p_val_gene, ref_base, alt_base)
    
    # --- 绘制 FL 表型箱线图 ---
    current_fl_pheno <- pheno_long %>% filter(trait == "FL") %>% inner_join(current_genotype, by = "Sample_ID")
    if (nrow(current_fl_pheno) > 0) plot_list$p_fl_pheno <- create_boxplot(current_fl_pheno, "FL Phenotype", "FL Phenotype", p_val_fl, ref_base, alt_base)
    
    # --- (v3.0 MODIFIED: 添加 Loop ID 到标题) ---
    if (!is.na(loop_id) && loop_id != "") {
      current_loop_strength <- loop_strength_long %>% filter(Loop_ID == loop_id) %>% inner_join(current_genotype, by = "Sample_ID")
      if (nrow(current_loop_strength) > 0) {
        # 动态创建包含Loop ID的标题
        loop_strength_title <- paste0("Loop Strength (", loop_id, ")")
        plot_list$p_loop_strength <- create_boxplot(current_loop_strength, "Loop Strength", loop_strength_title, p_val_loop_strength, ref_base, alt_base)
      } else {
        plot_list$p_loop_strength <- ggplot() + theme_void() + labs(title = "No Loop Strength Data")
      }
    } else {
      plot_list$p_loop_strength <- ggplot() + theme_void() + labs(title = "No Loop Data (B1 Pass)")
    }
    
    # --- (v3.0 MODIFIED: 调用新的 2x2 热图函数并添加 Loop ID) ---
    if (!is.na(loop_id) && loop_id != "") {
      current_loop_status <- loop_status_long %>% filter(Loop_ID == loop_id) %>% inner_join(current_genotype, by = "Sample_ID")
      if (nrow(current_loop_status) > 0) {
        # 动态创建包含Loop ID的标题
        loop_status_title <- paste0("Loop Status (", loop_id, ")")
        # *** 调用新的 2x2 热图函数 ***
        plot_list$p_loop_status_contingency <- create_contingency_heatmap(current_loop_status, loop_status_title, p_val_loop_status, ref_base, alt_base)
      } else {
        plot_list$p_loop_status_contingency <- ggplot() + theme_void() + labs(title = "No Loop Status Data")
      }
    } else {
      plot_list$p_loop_status_contingency <- ggplot() + theme_void() + labs(title = "No Loop Data (B1 Pass)")
    }
    
    # --- 组合所有子图 (v3.0 MODIFIED: 更新为新的绘图对象名称) ---
    combined_plot <- (plot_list$p_tad_is + plot_list$p_fs_pheno + plot_list$p_loop_strength) /
      (plot_list$p_gene_exp + plot_list$p_fl_pheno + plot_list$p_loop_status_contingency) +
      plot_annotation(
        title = paste0("Gene: ", gene_id, ", TAD: ", tad_id, ", SNP: ", snp_id) 
      ) & theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    
    # 保存组合图像
    output_filename_combined <- paste0(gene_id, "_", tad_id, "_", snp_id, "_", loop_id, "_Combined_Plots_v3.png")
    # 清理文件名，防止 loop_id 为 NA (B1通过时) 导致文件名中出现 "NA"
    output_filename_combined <- gsub("_NA_", "_", output_filename_combined, fixed = TRUE) 
    
    ggsave(
      filename = file.path(opt$output_dir, output_filename_combined),
      plot = combined_plot,
      width = 11,  # 宽度增加到11，以更好地容纳三列（2x2图的图例需要额外空间）
      height = 7,  
      dpi = 300 
    )
    cat(paste("Generated combined plot for:", gene_id, tad_id, snp_id, loop_id, "\n"))
  } else {
    cat(paste("Skipping triplet due to missing Full_SNP_ID:", gene_id, tad_id, snp_id, "\n"))
  }
}

cat("所有合并图已生成并保存到指定目录中。\n")