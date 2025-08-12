#!/bin/bash

# =============================================================================
#           步骤 10: 最终可视化成果汇总脚本
# =============================================================================
# 功能: 根据第9步的四元组汇总表，从第5步和第7步的输出中，
#       收集并整理所有最终相关的图片和日志文件。
# =============================================================================

# --- 脚本用法提示 ---
usage() {
    # 【核心修正】: 更新参数说明以匹配代码实现
    echo "Usage: $0 -q <quartet_file> -s <step5_dir> -d <step7_dir> -o <output_dir>"
    echo "  -q <quartet_file>     : 第9步生成的四元组汇总文件路径"
    echo "  -s <step5_dir>        : 第5步生成的 snp_pheno_boxplots 目录路径"
    echo "  -d <step7_dir>        : 第7步生成的 final_loop_analysis 目录路径"
    echo "  -o <output_dir>       : 保存最终汇总结果的总目录"
    exit 1
}

# --- 解析命令行参数 ---
QUARTET_FILE=""
STEP5_DIR=""
STEP7_DIR=""
OUTPUT_DIR=""

# 【核心修正】: getopts字符串 "q:s:d:o:" 与下面的case语句完全匹配
while getopts "q:s:d:o:" opt; do
    case ${opt} in
        q) QUARTET_FILE=$OPTARG ;;
        s) STEP5_DIR=$OPTARG ;; # -s for step5
        d) STEP7_DIR=$OPTARG ;;  # -d for step7 dir
        o) OUTPUT_DIR=$OPTARG ;;
        *) usage ;;
    esac
done

if [ -z "${QUARTET_FILE}" ] || [ -z "${STEP5_DIR}" ] || [ -z "${STEP7_DIR}" ] || [ -z "${OUTPUT_DIR}" ]; then
    echo "错误：缺少必需的参数。"
    usage
fi

# --- 脚本准备 ---
echo "--- 步骤10：最终可视化成果汇总 ---"
# 创建两个子目录
OUTPUT_BOXPLOTS="${OUTPUT_DIR}/final_snp_pheno_boxplots"
OUTPUT_LOOPS="${OUTPUT_DIR}/final_loop_analysis"
mkdir -p "${OUTPUT_BOXPLOTS}"
mkdir -p "${OUTPUT_LOOPS}"

echo "最终成果将分别保存在:"
echo "  - ${OUTPUT_BOXPLOTS}"
echo "  - ${OUTPUT_LOOPS}"
echo "---------------------------------"

# --- 主循环: 遍历四元组文件 ---
copied_s5_count=0
copied_s7_count=0
skipped_s5_count=0
skipped_s7_count=0
# 【严谨性修正】: 明确声明关联数组
declare -A seen_missing_s5
declare -A seen_missing_s7_plot

# 使用awk读取文件，跳过标题行，并用 | 作为分隔符输出
# 这样可以安全处理可能包含特殊字符的ID
awk -F',' 'NR > 1 {print $1"|"$2"|"$3"|"$4"|"$5}' "${QUARTET_FILE}" | while IFS='|' read -r gene_id tad_id loop_id snp_id full_snp_id; do

    # --- 1. 处理第5步的Boxplot图片 ---
    
    # 构建源文件路径
    # 文件名格式: GeneID_SNP-ID.png (根据draw_boxplot...R脚本)
    safe_gene_s5=$(echo "${gene_id}" | sed 's/[^A-Za-z0-9_.-]/_/g')
    safe_snp_s5=$(echo "${snp_id}" | sed 's/[^A-Za-z0-9_.-]/_/g')
    source_boxplot="${STEP5_DIR}/${tad_id}/${safe_gene_s5}_${safe_snp_s5}.png"

    if [ -f "${source_boxplot}" ]; then
        # 构建目标目录和文件路径
        target_boxplot_dir="${OUTPUT_BOXPLOTS}/${tad_id}/${gene_id}"
        mkdir -p "${target_boxplot_dir}"
        target_boxplot_path="${target_boxplot_dir}/${snp_id}.png"
        
        # 复制文件，仅当目标文件不存在时才打印日志，避免重复复制和记录
        if [ ! -f "${target_boxplot_path}" ]; then
            cp "${source_boxplot}" "${target_boxplot_path}"
            echo "  [Boxplot] 已复制: ${tad_id}/${gene_id}/${snp_id}.png"
            copied_s5_count=$((copied_s5_count + 1))
        fi
    else
        # 仅当第一次遇到此文件缺失时记录
        if [[ ! -v seen_missing_s5["${source_boxplot}"] ]]; then
              echo "  - [Boxplot] 警告: 源文件未找到，跳过: ${source_boxplot}"
              seen_missing_s5["${source_boxplot}"]=1
              skipped_s5_count=$((skipped_s5_count + 1))
        fi
    fi

    # --- 2. 处理第7步的Loop相关图片和日志 (仅当存在Loop时) ---
    if [ "${loop_id}" != "NA" ] && [ ! -z "${loop_id}" ]; then
        # 【核心修正】: 使用sed移除可能存在的浮点数后缀".0"
        cleaned_loop_id=$(echo "${loop_id}" | sed 's/\.0$//')
        
        # 文件名格式: GeneID_TADID_FullSNP-ID.png/log (根据analyze_and_plot_loops.R脚本)
        safe_snp_s7=$(echo "${full_snp_id}" | sed 's/:/_/g')
        base_filename="${gene_id}_${tad_id}_${safe_snp_s7}"
        
        # 【最终修正】: 移除路径中的tad_id和gene_id，直接访问Loop目录
        source_loop_plot="${STEP7_DIR}/Loop_${cleaned_loop_id}/${base_filename}.png"
        source_loop_log="${STEP7_DIR}/Loop_${cleaned_loop_id}/${base_filename}_chisq_test.log"

        # 构建目标目录
        target_loop_dir="${OUTPUT_LOOPS}/${tad_id}/${gene_id}/Loop_${cleaned_loop_id}"
        mkdir -p "${target_loop_dir}"
        
        # 复制图片
        if [ -f "${source_loop_plot}" ]; then
            target_loop_plot_path="${target_loop_dir}/${snp_id}.png"
            if [ ! -f "${target_loop_plot_path}" ]; then
                cp "${source_loop_plot}" "${target_loop_plot_path}"
                echo "  [Loop Plot] 已复制: ${tad_id}/${gene_id}/Loop_${cleaned_loop_id}/${snp_id}.png"
                copied_s7_count=$((copied_s7_count + 1))
            fi
        else
            if [[ ! -v seen_missing_s7_plot["${source_loop_plot}"] ]]; then
                  echo "  - [Loop Plot] 警告: 源文件未找到，跳过: ${source_loop_plot}"
                  seen_missing_s7_plot["${source_loop_plot}"]=1
                  skipped_s7_count=$((skipped_s7_count + 1))
            fi
        fi
        
        # 复制日志
        if [ -f "${source_loop_log}" ]; then
            target_loop_log_path="${target_loop_dir}/${snp_id}_chisq_test.log"
            if [ ! -f "${target_loop_log_path}" ]; then
                 cp "${source_loop_log}" "${target_loop_log_path}"
            fi
        fi
    fi
done

echo "---------------------------------"
echo "成果汇总完成。"
echo "SNP-Pheno Boxplots: 共复制 ${copied_s5_count} 个文件，跳过 ${skipped_s5_count} 个未找到的文件。"
echo "Loop 分析图表: 共复制 ${copied_s7_count} 个文件，跳过 ${skipped_s7_count} 个未找到的文件。"