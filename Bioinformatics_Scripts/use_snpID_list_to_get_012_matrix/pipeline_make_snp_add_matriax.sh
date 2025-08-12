#!/bin/bash

# --- 脚本用法提示 ---
usage() {
    echo "Usage: $0 -v <vcf_file_path> -s <snp_input_file> -c <snp_column_index> [-H] -o <final_output_matrix_name>"
    echo "       -v <vcf_file_path>         : 输入 VCF 文件的路径 (e.g., L_maf005_filtered_samples_final.vcf.gz)"
    echo "       -s <snp_input_file>        : 包含 SNP ID 的输入文件路径 (e.g., factor_15_output_sig_best_hit_formatted.txt)"
    echo "       -c <snp_column_index>      : SNP ID 在 snp_input_file 中的列索引 (从 1 开始计数)"
    echo "       -H                         : 可选。如果 snp_input_file 包含标题行，请使用此标志。"
    echo "       -o <final_output_matrix_name> : 最终输出转置矩阵的文件名 (e.g., snp_matrix_transposed.txt)"
    exit 1
}

# --- 解析命令行参数 ---
VCF_FILE=""
SNP_INPUT_FILE=""
SNP_COL_INDEX=""
HAS_HEADER_FLAG=""
FINAL_OUTPUT_MATRIX=""

while getopts "v:s:c:Ho:" opt; do
    case ${opt} in
        v) VCF_FILE=$OPTARG ;;
        s) SNP_INPUT_FILE=$OPTARG ;;
        c) SNP_COL_INDEX=$OPTARG ;;
        H) HAS_HEADER_FLAG="--has_header" ;; # 如果 -H 存在，则设置为 --has_header
        o) FINAL_OUTPUT_MATRIX=$OPTARG ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

# 检查必需参数是否提供
if [ -z "${VCF_FILE}" ] || [ -z "${SNP_INPUT_FILE}" ] || [ -z "${SNP_COL_INDEX}" ] || [ -z "${FINAL_OUTPUT_MATRIX}" ]; then
    echo "错误：缺少必需的参数。"
    usage
fi

# --- 脚本配置 ---
# 确保Python脚本在PATH中或提供完整路径
EXTRACT_SNP_SCRIPT="./extract_unique_snp_list.py" # 假设脚本在当前目录
TRANS_MATRIX_SCRIPT="./trans_matrix.py"       # 假设脚本在当前目录

# 临时文件名称定义
TEMP_UNIQUE_SNP_LIST="temp_unique_snps_to_extract.txt"
TEMP_EXTRACTED_VCF="temp_extracted_snps.vcf.gz"
TEMP_RAW_MATRIX="temp_extracted_snps_matrix.raw"

echo "==================================================="
echo "                  Coloc 数据处理流程开始           "
echo "==================================================="
echo "输入 VCF 文件: ${VCF_FILE}"
echo "SNP 输入文件: ${SNP_INPUT_FILE}"
echo "SNP 列索引 (1-based): ${SNP_COL_INDEX}"
echo "SNP 输入文件是否有标题行: ${HAS_HEADER_FLAG:---无标题行}" # 显示是否提供了 -H
echo "最终输出矩阵文件: ${FINAL_OUTPUT_MATRIX}"
echo "---------------------------------------------------"

# --- 步骤 1: 提取并去重排序 SNP 列表 ---
echo "阶段 1/4: 提取并去重排序 SNP 列表..."
python3 "${EXTRACT_SNP_SCRIPT}" \
    -i "${SNP_INPUT_FILE}" \
    -o "${TEMP_UNIQUE_SNP_LIST}" \
    -c "${SNP_COL_INDEX}" \
    ${HAS_HEADER_FLAG} # 如果 --has_header 变量为空，则不会传递此参数

if [ $? -ne 0 ]; then
    echo "错误: 提取并去重 SNP 列表失败。脚本退出。"
    exit 1
fi
echo "阶段 1/4 完成。"
echo "---------------------------------------------------"

# --- 步骤 2: 使用 PLINK2 过滤 VCF 文件 ---
echo "阶段 2/4: 使用 PLINK2 过滤 VCF 文件..."
plink2 --vcf "${VCF_FILE}" \
    --extract "${TEMP_UNIQUE_SNP_LIST}" \
    --export vcf-4.2 bgz \
    --allow-extra-chr \
    --out "${TEMP_EXTRACTED_VCF%.vcf.gz}" # PLINK2 自动添加 .vcf.gz 后缀

if [ $? -ne 0 ]; then
    echo "错误: PLINK2 VCF 过滤失败。脚本退出。"
    exit 1
fi
echo "阶段 2/4 完成。"
echo "---------------------------------------------------"

# --- 步骤 3: 使用 PLINK2 生成原始矩阵 ---
echo "阶段 3/4: 使用 PLINK2 生成原始矩阵..."
plink2 --vcf "${TEMP_EXTRACTED_VCF}" \
    --export A include-alt \
    --allow-extra-chr \
    --out "${TEMP_RAW_MATRIX%.raw}" # PLINK2 自动添加 .raw 后缀

if [ $? -ne 0 ]; then
    echo "错误: PLINK2 矩阵生成失败。脚本退出。"
    exit 1
fi
echo "阶段 3/4 完成。"
echo "---------------------------------------------------"

# --- 步骤 4: 转置并格式化矩阵 ---
echo "阶段 4/4: 转置并格式化矩阵..."
# 调用修改后的 trans_matrix.py 脚本，将最终输出文件名作为参数传入
python3 "${TRANS_MATRIX_SCRIPT}" \
    -i "${TEMP_RAW_MATRIX}" \
    -o "${FINAL_OUTPUT_MATRIX}"

if [ $? -ne 0 ]; then
    echo "错误: 矩阵转置和格式化失败。脚本退出。"
    exit 1
fi
echo "阶段 4/4 完成。"
echo "---------------------------------------------------"

# --- 清理临时文件 ---
echo "清理临时文件..."
rm -f "${TEMP_UNIQUE_SNP_LIST}" "${TEMP_EXTRACTED_VCF}" "${TEMP_EXTRACTED_VCF}.csi" "${TEMP_RAW_MATRIX}"* # PLINK2可能生成多个raw文件
echo "清理完成。"

echo "==================================================="
echo "                  所有任务完成！                   "
echo "最终结果保存在: ${FINAL_OUTPUT_MATRIX}"
echo "==================================================="
