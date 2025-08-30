# ===================================================================
# CMplot思想重构版: Circular Manhattan Plot (修改版 - 陆地棉专用)
# 主要修改：
# 1. 染色体命名：使用A01-A13, D01-D13格式
# 2. 使用实际染色体长度文件：l128.genome.lst
# ===================================================================

# 加载必要的库
library(dplyr)

# ===================================================================
# 1. 数据准备
# ===================================================================
# 读取染色体长度信息
genome_info <- read.table("l128.genome.lst", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(genome_info) <- c("chr_id", "chr_length")

# 读取QTL数据文件
gene_data <- read.table("extracted_eQTL_data.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tad_data  <- read.table("extracted_TAD_QTL_data.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
loop_data <- read.table("extracted_loop_QTL_data.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 解析函数，从var_id中提取染色体和位置信息（修改版）
parse_data <- function(df) {
  df %>%
    mutate(
      chr_part = sapply(strsplit(var_id, "_"), function(x) {
        hc_part <- x[grepl("HC04", x)]
        if(length(hc_part) > 0) return(hc_part[1])
        return(x[1])
      }),
      Position = as.numeric(sapply(strsplit(var_id, "_"), function(x) {
        num_parts <- x[grepl("^[0-9]+$", x)]
        if(length(num_parts) > 0) return(num_parts[length(num_parts)])
        return(NA)
      })),
      chr_num = case_when(
        grepl("HC04_A", chr_part) ~ as.numeric(gsub("HC04_A", "", chr_part)),
        grepl("HC04_D", chr_part) ~ as.numeric(gsub("HC04_D", "", chr_part)) + 13,
        grepl("^[0-9]+$", chr_part) ~ as.numeric(chr_part),
        TRUE ~ NA_real_
      ),
      # 修改染色体命名格式：A01-A13, D01-D13
      Chromosome = case_when(
        chr_num >= 1 & chr_num <= 13 ~ paste0("A", sprintf("%02d", chr_num)),
        chr_num >= 14 & chr_num <= 26 ~ paste0("D", sprintf("%02d", chr_num - 13)),
        TRUE ~ NA_character_
      ),
      # 添加对应的HC04格式用于匹配长度信息
      HC04_format = case_when(
        chr_num >= 1 & chr_num <= 13 ~ paste0("HC04_A", sprintf("%02d", chr_num)),
        chr_num >= 14 & chr_num <= 26 ~ paste0("HC04_D", sprintf("%02d", chr_num - 13)),
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(chr_num), !is.na(Position)) %>%
    select(-chr_part)
}

# 应用解析函数
gene_data <- parse_data(gene_data)
tad_data  <- parse_data(tad_data)
loop_data <- parse_data(loop_data)

# 检查解析结果
cat("染色体解析结果：\n")
cat("基因数据染色体范围：", range(gene_data$chr_num, na.rm=TRUE), "\n")
cat("TAD数据染色体范围：", range(tad_data$chr_num, na.rm=TRUE), "\n")
cat("Loop数据染色体范围：", range(loop_data$chr_num, na.rm=TRUE), "\n")
cat("检测到的染色体：", unique(gene_data$Chromosome[!is.na(gene_data$Chromosome)]), "\n")

# ===================================================================
# 2. 核心逻辑修复区
# ===================================================================

# 步骤A: 创建全局显著SNP列表
sig_gene_snps <- gene_data %>% filter(qval < 0.05) %>% pull(var_id) %>% unique()
sig_tad_snps  <- tad_data  %>% filter(FDR_pval < 0.01 & r_squared > 0.05) %>% pull(var_id) %>% unique()
sig_loop_snps <- loop_data %>% filter(FDR_pval < 0.01 & r_squared > 0.05) %>% pull(var_id) %>% unique()

# p值转换为-log10
log10_trans <- function(p) {
  p_min_non_zero <- min(p[p > 0], na.rm = TRUE)
  p[p == 0] <- p_min_non_zero * 0.1
  -log10(p)
}

# 步骤B: 重构 prepare_plot_df 函数
prepare_plot_df <- function(df, pval_col, track_type) {

  # 定义显著性阈值
  GENE_QVAL_THRESH <- 0.05
  QTL_FDR_THRESH <- 0.01
  QTL_R2_THRESH <- 0.05

  # 内部辅助函数：仅当一个点"行显著"时，才调用此函数来确定其组合颜色
  get_combo_color <- function(snp_id) {
    gene_sig <- snp_id %in% sig_gene_snps
    tad_sig  <- snp_id %in% sig_tad_snps
    loop_sig <- snp_id %in% sig_loop_snps

    # 当前轨道显著，判断组合类型
    if (gene_sig && tad_sig && loop_sig) {
      return("Red") # 三轨道都显著
    } else if (gene_sig && tad_sig && !loop_sig) {
      if (track_type %in% c("gene", "tad")) return("Blue") else return("Grey")
    } else if (gene_sig && !tad_sig && loop_sig) {
      if (track_type %in% c("gene", "loop")) return("Green") else return("Grey")
    } else if (!gene_sig && tad_sig && loop_sig) {
      if (track_type %in% c("tad", "loop")) return("Pink") else return("Grey")
    } else {
      # 仅在当前分析的轨道类型中显著
      return("Yellow")
    }
  }

  # 使用if/else分离不同数据结构的处理逻辑
  if (track_type == "gene") {
    # 只处理 gene 数据，不涉及 r_squared
    processed_df <- df %>%
      mutate(
        is_row_significant = .data[[pval_col]] < GENE_QVAL_THRESH,
        color_group = ifelse(
          !is_row_significant,
          "Grey",
          sapply(var_id, get_combo_color)
        ),
        value = log10_trans(.data[[pval_col]])
      )
  } else {
    # 处理 tad 和 loop 数据，包含 r_squared
    processed_df <- df %>%
      mutate(
        is_row_significant = .data[[pval_col]] < QTL_FDR_THRESH & r_squared > QTL_R2_THRESH,
        color_group = ifelse(
          !is_row_significant,
          "Grey",
          sapply(var_id, get_combo_color)
        ),
        value = log10_trans(.data[[pval_col]])
      )
  }

  # 返回统一格式的结果
  processed_df %>%
    select(chr = Chromosome, x = Position, value, color_group, snp_id = var_id)
}

# 步骤C: 使用修复后的函数生成绘图数据
gene_plot_df <- prepare_plot_df(gene_data, "qval", "gene")
tad_plot_df  <- prepare_plot_df(tad_data, "FDR_pval", "tad")
loop_plot_df <- prepare_plot_df(loop_data, "FDR_pval", "loop")

# ===================================================================
# 3. 智能Y轴设置
# ===================================================================
setup_smart_yaxis <- function(df, n_ticks = 5) {
  values <- df$value[is.finite(df$value)]
  if(length(values) == 0) return(list(lim = c(0, 10), ticks = seq(0, 10, 2)))
  data_max <- max(values, na.rm = TRUE)
  y_max <- ceiling(data_max / 5) * 5
  if(y_max < data_max) y_max <- ceiling(data_max/2)*2
  ticks <- pretty(c(0, y_max), n = n_ticks)
  list(lim = c(0, max(ticks)), ticks = ticks)
}
yaxis_gene <- setup_smart_yaxis(gene_plot_df)
yaxis_tad  <- setup_smart_yaxis(tad_plot_df)
yaxis_loop <- setup_smart_yaxis(loop_plot_df)

# ===================================================================
# 4. 基因组线性化（使用实际染色体长度）
# ===================================================================
linearize_genome <- function(datasets, genome_info, inter_chr_gap = 3e7, chr_compression = 0.6, y_axis_gap_degrees = 6) {
  
  # 创建染色体顺序映射（倒序：D13-D01-A13-A01）
  chr_order_map <- data.frame(
    HC04_format = c(paste0("HC04_D", sprintf("%02d", 13:1)), 
                    paste0("HC04_A", sprintf("%02d", 13:1))),
    display_name = c(paste0("D", sprintf("%02d", 13:1)), 
                     paste0("A", sprintf("%02d", 13:1))),
    order_num = 1:26,
    stringsAsFactors = FALSE
  )
  
  # 合并基因组信息和顺序信息
  chr_lengths <- genome_info %>%
    left_join(chr_order_map, by = c("chr_id" = "HC04_format")) %>%
    filter(!is.na(order_num)) %>%
    arrange(order_num) %>%
    mutate(
      compressed_len = chr_length * chr_compression
    )
  
  cat("染色体长度信息：\n")
  print(chr_lengths)
  
  chr_info <- chr_lengths %>%
    mutate(
      chr_len_gap = compressed_len + inter_chr_gap,
      cumulative_pos = cumsum(as.numeric(chr_len_gap)),
      start_pos = lag(cumulative_pos, default = 0),
      end_pos = start_pos + compressed_len,
      mid_pos = start_pos + compressed_len / 2,
      original_len = chr_length
    ) %>%
    rename(chr = display_name) # 使用显示名称作为染色体名

  total_genome_len <- max(chr_info$cumulative_pos)
  y_axis_gap_radians <- (y_axis_gap_degrees * 2) * pi / 180
  data_angle_range <- 2 * pi - y_axis_gap_radians
  
  # 【关键修改】计算让最后一个染色体后间隔中心位于12点钟方向的起始角度
  last_chr_end <- chr_info$end_pos[nrow(chr_info)]  # 最后一个染色体的结束位置
  last_gap_center <- last_chr_end + inter_chr_gap / 2  # 最后间隔的中心位置
  
  # 目标：让 last_gap_center 对应 90度（pi/2）
  # 当前：角度 = data_start_angle + (position / total_len) * data_angle_range
  # 要求：pi/2 = data_start_angle + (last_gap_center / total_len) * data_angle_range
  # 解得：data_start_angle = pi/2 - (last_gap_center / total_len) * data_angle_range
  calculated_start_angle <- pi/2 - (last_gap_center / total_genome_len) * data_angle_range
  
  cat(sprintf("计算信息：\n"))
  cat(sprintf("  最后染色体结束位置: %.0f\n", last_chr_end))
  cat(sprintf("  最后间隔中心位置: %.0f\n", last_gap_center))
  cat(sprintf("  总基因组长度: %.0f\n", total_genome_len))
  cat(sprintf("  计算得到的起始角度: %.4f 弧度 (%.2f 度)\n", calculated_start_angle, calculated_start_angle * 180 / pi))
  
  processed_datasets <- lapply(datasets, function(df) {
    df %>%
      left_join(select(chr_info, chr, start_pos, original_len, compressed_len), by = "chr") %>%
      mutate(
        compressed_pos = (x / original_len) * compressed_len,
        abs_pos = start_pos + compressed_pos
      ) %>%
      select(-start_pos, -original_len, -compressed_len, -compressed_pos)
  })

  return(list(
    datasets = processed_datasets,
    chr_info = chr_info,
    total_len = total_genome_len,
    data_angle_range = data_angle_range,
    y_axis_gap_radians = y_axis_gap_radians,
    calculated_start_angle = calculated_start_angle  # 返回计算得到的起始角度
  ))
}

linear_data <- linearize_genome(list(gene_plot_df, tad_plot_df, loop_plot_df),
                                genome_info,
                                inter_chr_gap = 3e7,
                                chr_compression = 0.6,
                                y_axis_gap_degrees = 0)
gene_plot_df <- linear_data$datasets[[1]]
tad_plot_df  <- linear_data$datasets[[2]]
loop_plot_df <- linear_data$datasets[[3]]
chr_info     <- linear_data$chr_info
total_len    <- linear_data$total_len
data_angle_range <- linear_data$data_angle_range
y_axis_gap_radians <- linear_data$y_axis_gap_radians
calculated_start_angle <- linear_data$calculated_start_angle

# ===================================================================
# 5. 绘图函数定义
# ===================================================================

draw_sector <- function(start_angle, end_angle, r_inner, r_outer, col, border = NA) {
  theta <- seq(start_angle, end_angle, length.out = 100)
  x_outer <- r_outer * cos(theta)
  y_outer <- r_outer * sin(theta)
  x_inner <- rev(r_inner * cos(theta))
  y_inner <- rev(r_inner * sin(theta))
  polygon(c(x_outer, x_inner), c(y_outer, y_inner), col = col, border = border)
}

draw_chromosome_track <- function(chr_info, total_len, radius, angle_start_rad) {
  # 绘制所有染色体扇形
  for (i in 1:nrow(chr_info)) {
    # 使用计算得到的起始角度
    start_angle <- calculated_start_angle + (chr_info$start_pos[i] / total_len) * data_angle_range
    end_angle <- calculated_start_angle + (chr_info$end_pos[i] / total_len) * data_angle_range
    mid_angle <- calculated_start_angle + (chr_info$mid_pos[i] / total_len) * data_angle_range

    # 绘制染色体扇形（区分A和D亚基因组颜色）
    chr_color <- ifelse(grepl("^A", chr_info$chr[i]), "#4297D8", "#4297D8")  # A亚基因组深蓝，D亚基因组浅蓝
    draw_sector(start_angle, end_angle, radius[1], radius[2], chr_color, border = "black")

    # 添加染色体名称
    text_radius <- (radius[1] + radius[2]) / 2
    text_x <- text_radius * cos(mid_angle)
    text_y <- text_radius * sin(mid_angle)
    text_angle_deg <- (mid_angle * 180 / pi) - 90

    # 绘制染色体名称（使用A01-A13, D01-D13格式）
    text(text_x, text_y,
         labels = chr_info$chr[i],  # 使用新的染色体命名格式
         srt = text_angle_deg,
         adj = c(0.5, 0.5),
         cex = 0.6,
         col = "white",
         font = 2)
  }

  # 绘制外圆完整边界线
  data_end_angle <- calculated_start_angle + data_angle_range
  theta_outer <- seq(calculated_start_angle, data_end_angle, length.out = 300)
  x_outer <- radius[2] * cos(theta_outer)
  y_outer <- radius[2] * sin(theta_outer)
  lines(x_outer, y_outer, col = "black", lwd = 2)

  # 绘制内圆完整边界线
  theta_inner <- seq(calculated_start_angle, data_end_angle, length.out = 300)
  x_inner <- radius[1] * cos(theta_inner)
  y_inner <- radius[1] * sin(theta_inner)
  lines(x_inner, y_inner, col = "black", lwd = 2)
}

draw_data_points <- function(df, yaxis, radius, total_len, angle_start_rad, colors, order) {
  df_clean <- df %>% filter(is.finite(value))

  # 使用计算得到的起始角度
  angles <- calculated_start_angle + (df_clean$abs_pos / total_len) * data_angle_range
  radii <- radius[1] + (radius[2] - radius[1]) *
    (df_clean$value - yaxis$lim[1]) / (yaxis$lim[2] - yaxis$lim[1])

  valid_indices <- radii <= radius[2] & radii >= radius[1]
  df_clean <- df_clean[valid_indices, ]
  angles <- angles[valid_indices]
  radii <- radii[valid_indices]

  x_coords <- radii * cos(angles)
  y_coords <- radii * sin(angles)
  df_clean$x_plot <- x_coords
  df_clean$y_plot <- y_coords

  for (layer in order) {
    subset_data <- df_clean %>% filter(color_group == layer)
    if (nrow(subset_data) > 0) {
      points(subset_data$x_plot, subset_data$y_plot, col = colors[layer], pch = 19, cex = 0.4)
    }
  }
}

draw_grids_only <- function(tracks_info, data_angle_range, y_axis_gap_radians) {
  # 只绘制网格线，不绘制Y轴
  for (track in tracks_info) {
    for (tick in track$yaxis$ticks) {
      # 包含0值的网格线
      if (tick >= 0 && tick <= track$yaxis$lim[2]) {
        r <- track$radius[1] + (track$radius[2] - track$radius[1]) *
          (tick - track$yaxis$lim[1]) / (track$yaxis$lim[2] - track$yaxis$lim[1])

        # 使用计算得到的起始角度
        data_end_angle <- calculated_start_angle + data_angle_range
        theta <- seq(calculated_start_angle, data_end_angle, length.out = 300)

        lines(r * cos(theta), r * sin(theta), col = "#F0F0F0", lwd = 1)
      }
    }
    # 绘制阈值线
    if (!is.null(track$threshold) && track$threshold <= track$yaxis$lim[2]) {
      r <- track$radius[1] + (track$radius[2] - track$radius[1]) *
        (track$threshold - track$yaxis$lim[1]) / (track$yaxis$lim[2] - track$yaxis$lim[1])

      # 使用计算得到的起始角度
      data_end_angle <- calculated_start_angle + data_angle_range
      theta <- seq(calculated_start_angle, data_end_angle, length.out = 300)

      lines(r * cos(theta), r * sin(theta), col = "black", lty = 2, lwd = 1)
    }
  }
}

draw_yaxis_only <- function(tracks_info, data_angle_range, y_axis_gap_radians) {
  # 只绘制Y轴，不绘制网格线
  y_axis_angle <- pi/2
  for (track in tracks_info) {
    r_inner <- track$radius[1]
    r_outer <- track$radius[2]

    lines(c(r_inner * cos(y_axis_angle), r_outer * cos(y_axis_angle)),
          c(r_inner * sin(y_axis_angle), r_outer * sin(y_axis_angle)),
          col = track$axis_color, lwd = 1)

    for (tick in track$yaxis$ticks) {
      if (tick <= track$yaxis$lim[2]) {
        r_tick <- r_inner + (r_outer - r_inner) *
          (tick - track$yaxis$lim[1]) / (track$yaxis$lim[2] - track$yaxis$lim[1])

        lines(c(r_tick * cos(y_axis_angle), r_tick * cos(y_axis_angle) + 0.01),
              c(r_tick * sin(y_axis_angle), r_tick * sin(y_axis_angle)),
              col = track$axis_color, lwd = 1.5)

        text(r_tick * cos(y_axis_angle) + 0.015, r_tick * sin(y_axis_angle),
             labels = tick, adj = c(0, 0.5), cex = 0.35, col = track$axis_color)
      }
    }
  }
}

# ===================================================================
# 6. 绘图主流程
# ===================================================================

# 定义颜色和绘图顺序
colors <- c(Red = "#D43F3A", Green = "#4CAE4C", Blue = "#357EBD",
            Pink = "#FF69B4", Yellow = "#FFD700", Grey = "#CCCCCC")
point_order <- c("Grey", "Yellow", "Pink", "Blue", "Green", "Red")
gene_threshold <- -log10(0.05)
qtl_threshold <- -log10(0.01)

# 定义轨道半径
radii <- list(
  chr_track  = c(0.93, 1.0),
  loop_track = c(0.70, 0.90),
  tad_track  = c(0.45, 0.65),
  gene_track = c(0.20, 0.40)
)

# 组合轨道信息
tracks_info <- list(
  list(name = "loop", df = loop_plot_df, yaxis = yaxis_loop, radius = radii$loop_track, threshold = qtl_threshold, axis_color = "black"),
  list(name = "tad",  df = tad_plot_df,  yaxis = yaxis_tad,  radius = radii$tad_track,  threshold = qtl_threshold, axis_color = "black"),
  list(name = "gene", df = gene_plot_df, yaxis = yaxis_gene, radius = radii$gene_track, threshold = gene_threshold, axis_color = "black")
)

# 角度偏移
angle_offset_rad <- pi/2

# 初始化绘图设备
jpeg("Cotton_Circular_Manhattan_Plot.jpg", width = 6, height = 6, units = "in", res = 1800)
par(mar = c(0.2, 0.2, 0.2, 0.2), xpd = TRUE)
plot(NULL, xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1)

# 按层次绘制图形
# a. 绘制染色体轨道 (底层)
draw_chromosome_track(chr_info, total_len, radii$chr_track, angle_offset_rad)

# b. 绘制网格线 (中下层)
draw_grids_only(tracks_info, data_angle_range, y_axis_gap_radians)

# c. 绘制数据点 (中上层)
for (track in tracks_info) {
  draw_data_points(track$df, track$yaxis, track$radius, total_len, angle_offset_rad, colors, point_order)
}

# d. 绘制Y轴 (顶层)
draw_yaxis_only(tracks_info, data_angle_range, y_axis_gap_radians)

# 关闭设备
dev.off()

# 输出完成信息
cat("陆地棉圆形曼哈顿图绘制完成！\n")
cat("染色体命名格式：A01-A13, D01-D13\n")
cat("已使用实际染色体长度信息\n")