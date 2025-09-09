#!/usr/bin/env bash

# ==================================================================================================
#
#       SNP -> Gene/TAD/Phenotype -> Loop Association Pipeline (MODIFIED V2.4)
#
# ==================================================================================================
#
# 版本: 5.4 (Production Release) - (V2脚本重构, 更新Step5逻辑, 暴露所有6个分析参数)
# 日期: 2025-09-07
# 描述: 从SNP到染色质环作用机制的完整生物信息分析流程。
#
# 特性:
#   - 8个核心分析阶段（基于V2脚本重构）。
#   - 6个可配置的关键分析参数已从子脚本暴露到主 shell (Step 5 和 Step 7)。
#   - 生产级代码结构、日志和错误处理。
#   - 流水线现在使用 Step 5 预过滤的结果。
#
# ==================================================================================================

# --- 脚本配置: 严格模式 ---
set -e
set -o pipefail

####################################################################################################
# Section 1: 核心函数库 (保持不变)
####################################################################################################

# --- 颜色与样式定义 ---
readonly C_OFF='\033[0m'
readonly C_BOLD='\033[1m'
readonly C_GREEN='\033[0;32m'
readonly C_RED='\033[0;31m'
readonly C_YELLOW='\033[0;33m'
readonly C_BLUE='\033[0;34m'
readonly C_PURPLE='\033[0;35m'
readonly C_CYAN='\033[0;36m'

# --- 脚本自身信息 ---
readonly SCRIPT_NAME=$(basename "$0")
PIPELINE_LOG_FILE=""

# @description 核心日志记录函数
log_msg() {
    local level_color="$1"
    local level_text="$2"
    local message="$3"
    echo -e "${level_color}${C_BOLD}[$(date '+%T')] [${level_text}]${C_OFF} ${message}" | tee -a "${PIPELINE_LOG_FILE}" >&2
}

log_info()    { log_msg "${C_BLUE}"   " INFO    " "$1"; }
log_success() { log_msg "${C_GREEN}"  " SUCCESS " "$1"; }
log_warn()    { log_msg "${C_YELLOW}" " WARNING " "$1"; }
log_error()   { log_msg "${C_RED}"    " ERROR   " "$1"; }

# @description 打印带边框的标题
print_header() {
    local title=" $1 "
    local len=${#title}
    local border=$(printf '─%.0s' $(seq 1 "$len"))
    echo -e "\n${C_PURPLE}╭─${border}─╮${C_OFF}" | tee -a "${PIPELINE_LOG_FILE}" >&2
    echo -e "${C_PURPLE}│ ${C_BOLD}${title}${C_OFF}${C_PURPLE} │${C_OFF}" | tee -a "${PIPELINE_LOG_FILE}" >&2
    echo -e "${C_PURPLE}╰─${border}─╯${C_OFF}" | tee -a "${PIPELINE_LOG_FILE}" >&2
}

# @description 执行一个流程阶段
run_stage() {
    local stage_num=$1
    local description=$2
    shift 2
    local command_to_run="$@"
    local log_file_name="step${stage_num}_$(echo "$description" | awk '{print $1}' | tr '[:upper:]' '[:lower:]').log"
    local stage_log_path="${LOG_DIR}/${log_file_name}"

    log_info "${C_CYAN}Stage ${stage_num}/${TOTAL_STAGES}: ${description}${C_OFF}"

    if eval "${command_to_run}" > "${stage_log_path}" 2>&1; then
        log_success "Stage ${stage_num} finished successfully."
    else
        log_error "Stage ${stage_num} failed. Check log for details: ${C_YELLOW}${stage_log_path}${C_OFF}"
        exit 1
    fi
}

####################################################################################################
# Section 2: 流程辅助函数 (用法'usage'已更新)
####################################################################################################

# @description 显示脚本用法 (V2.4 最终版 - 6个参数)
usage() {
    cat << EOF

${C_GREEN}${C_BOLD}CIWAS-Case-Analysis-Pipeline V2.4${C_OFF}

${C_BOLD}DESCRIPTION:${C_OFF}
  一个从SNP出发，建立其通过影响TAD和染色质环(Loop)来调控基因表达和表型的分子机制证据链的完整生物信息分析流程。

${C_BOLD}USAGE:${C_OFF}
  ${SCRIPT_NAME} -s <triplet_file> -o <output_prefix> -n <num_cores> [OPTIONS]

${C_BOLD}REQUIRED OPTIONS:${C_OFF}
  -s <file>     ${C_YELLOW}必需。${C_OFF}初始三元组文件路径。
  -o <prefix>   ${C_YELLOW}必需。${C_OFF}本次运行的唯一前缀，将作为主输出目录名。
  -n <cores>    ${C_YELLOW}必需。${C_OFF}并行计算时使用的CPU核心数。

${C_BOLD}PIPELINE OPTIONS:${C_OFF}
  -c            ${C_YELLOW}可选。${C_OFF}流程结束后，清理所有临时文件。

${C_BOLD}ANALYSIS PARAMETERS (Optional):${C_OFF}
  ${C_CYAN}(Step 5 Filtering):${C_OFF}
  -F <float>    ${C_YELLOW}可选。${C_OFF}(Step 5) TAD IS P值过滤 (Filter) 阈值。 (默认: 0.05)
  -G <float>    ${C_YELLOW}可选。${C_OFF}(Step 5) 基因表达 (Gene) P值过滤阈值。 (默认: 0.05)
  ${C_CYAN}(Step 7 Logic):${C_OFF}
  -T <int>      ${C_YELLOW}可选。${C_OFF}(Step 7) 分支1/2的距离阈值 (Threshold Distance) (kb)。 (默认: 80)
  -X <int>      ${C_YELLOW}可选。${C_OFF}(Step 7) 分支2的最大搜索距离 (MaX Distance) (kb)。 (默认: 250)
  -A <int>      ${C_YELLOW}可选。${C_OFF}(Step 7) 基因锚点(Anchor)搜索扩展距离 (kb)。 (默认: 40)
  -P <float>    ${C_YELLOW}可选。${C_OFF}(Step 7) 介导Loop的P值(P-value)显著性阈值。 (默认: 0.05)

EOF
    exit 1
}

# @description 检查核心外部命令依赖 (保持不变)
check_dependencies() {
    print_header "Checking Dependencies"
    log_info "Verifying required command-line tools..."
    local dependencies=("plink2" "Rscript" "python3" "awk" "tee" "mv") # 确保 mv 存在
    for cmd in "${dependencies[@]}"; do
        if ! command -v "$cmd" &>/dev/null; then
            log_error "Required command not found in PATH: ${cmd}"
            exit 1
        fi
    done
    log_success "All command-line dependencies are satisfied."
}

# @description 检查所有输入文件 (保持不变)
check_files_exist() {
    log_info "Verifying required input files and scripts..."
    for file in "$@"; do
        if [ ! -f "$file" ]; then
            log_error "Required file not found: ${file}"
            exit 1
        fi
    done
    log_success "All required files and scripts are accessible."
}

# @description 流程结束时清理临时文件 (保持不变)
cleanup() {
    if [ "$PERFORM_CLEANUP" = true ]; then
        log_info "Cleaning up temporary files based on -c flag..."
        rm -f "${TEMP_UNIQUE_SNP_LIST}" \
              "${TEMP_EXTRACTED_VCF}" \
              "${TEMP_EXTRACTED_VCF}.csi" \
              "${TEMP_RAW_MATRIX}"*
        log_success "Temporary files removed."
    else
        log_info "Skipping temporary file cleanup."
    fi
}

####################################################################################################
# Section 3: 主逻辑 (Main) - (数据库路径已更新)
####################################################################################################

main() {
    # --- 解析命令行参数 (必需参数) ---
    local TRIPLET_INPUT_FILE=""
    local OUTPUT_PREFIX=""
    local NUM_CORES=""
    local PERFORM_CLEANUP=false
    
    # --- (V2.4) 定义所有6个额外分析参数的默认值 ---
    local PARAM_TAD_IS_PVAL=0.05      # (Step 5 - F)
    local PARAM_GENE_EXP_PVAL=0.05    # (Step 5 - G)
    local PARAM_DIST_THRESHOLD=80     # (Step 7 - T)
    local PARAM_MAX_DIST=250          # (Step 7 - X)
    local PARAM_ANCHOR_DIST=40        # (Step 7 - A)
    local PARAM_LOOP_PVAL=0.05        # (Step 7 - P)

    # --- (V2.4) 更新 getopts 字符串 (添加 F: 和 G:) ---
    while getopts "s:o:n:cF:G:T:X:A:P:" opt; do # 扩展 getopts
        case ${opt} in
            s) TRIPLET_INPUT_FILE=$OPTARG ;;
            o) OUTPUT_PREFIX=$OPTARG ;;
            n) NUM_CORES=$OPTARG ;;
            c) PERFORM_CLEANUP=true ;;
            F) PARAM_TAD_IS_PVAL=$OPTARG ;;      # (新增) Step 5
            G) PARAM_GENE_EXP_PVAL=$OPTARG ;;    # (新增) Step 5
            T) PARAM_DIST_THRESHOLD=$OPTARG ;;       # (Step 7)
            X) PARAM_MAX_DIST=$OPTARG ;;             # (Step 7)
            A) PARAM_ANCHOR_DIST=$OPTARG ;;           # (Step 7)
            P) PARAM_LOOP_PVAL=$OPTARG ;;            # (Step 7)
            *) usage ;;
        esac
    done

    if [ -z "${TRIPLET_INPUT_FILE}" ] || [ -z "${OUTPUT_PREFIX}" ] || [ -z "${NUM_CORES}" ]; then
        log_error "Missing required arguments (-s, -o, -n)."
        usage
    fi

    # ========================== 1. 环境与路径设置 (路径已更新) ==========================

    # --- 核心目录与常量 ---
    readonly MAIN_OUTPUT_DIR="${OUTPUT_PREFIX}"
    readonly LOG_DIR="${MAIN_OUTPUT_DIR}/logs"
    readonly BASE_SCRIPT_PATH=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
    readonly TOTAL_STAGES=8 # V2流程共8个阶段

    # --- 初始化环境 ---
    mkdir -p "${MAIN_OUTPUT_DIR}"
    mkdir -p "${LOG_DIR}"
    PIPELINE_LOG_FILE="${LOG_DIR}/pipeline_main.log"
    > "${PIPELINE_LOG_FILE}" # 清空旧日志
    trap cleanup EXIT INT ERR

    # --- 打印欢迎横幅和配置信息 (V2.4) ---
    print_header "CIWAS-Case-Analyze-Pipeline V2.4"
    log_info "Run initiated on $(date)"
    log_info "──────────────────────────────────────────────────────────"
    log_info "Input Triplets : ${C_YELLOW}${TRIPLET_INPUT_FILE}${C_OFF}"
    log_info "Output Directory : ${C_YELLOW}${MAIN_OUTPUT_DIR}${C_OFF}"
    log_info "CPU Cores      : ${C_YELLOW}${NUM_CORES}${C_OFF}"
    log_info "Cleanup Temp   : ${C_YELLOW}${PERFORM_CLEANUP}${C_OFF}"
    log_info "──────────────────────────────────────────────────────────"
    log_info "${C_BOLD}Analysis Parameters (V2.4):${C_OFF}"
    log_info "Input Has Header: ${C_YELLOW}true (Hardcoded)${C_OFF}"
    log_info "${C_CYAN}(Step 5) Filter P-TAD (-F):  ${C_YELLOW}${PARAM_TAD_IS_PVAL}${C_OFF}"
    log_info "${C_CYAN}(Step 5) Filter P-Gene (-G): ${C_YELLOW}${PARAM_GENE_EXP_PVAL}${C_OFF}"
    log_info "${C_CYAN}(Step 7) Dist. Thresh (-T):  ${C_YELLOW}${PARAM_DIST_THRESHOLD} kb${C_OFF}"
    log_info "${C_CYAN}(Step 7) Max Distance (-X):    ${C_YELLOW}${PARAM_MAX_DIST} kb${C_OFF}"
    log_info "${C_CYAN}(Step 7) Anchor Dist (-A):   ${C_YELLOW}${PARAM_ANCHOR_DIST} kb${C_OFF}"
    log_info "${C_CYAN}(Step 7) Loop P-Value (-P):  ${C_YELLOW}${PARAM_LOOP_PVAL}${C_OFF}"
    log_info "──────────────────────────────────────────────────────────"


    # --- 检查依赖 ---
    check_dependencies

    # --- 路径注册表 (与V2流程一致) ---
    
    # 临时文件
    readonly TEMP_UNIQUE_SNP_LIST="${MAIN_OUTPUT_DIR}/temp_unique_snps.txt"
    readonly TEMP_EXTRACTED_VCF="${MAIN_OUTPUT_DIR}/temp_extracted.vcf.gz"
    readonly TEMP_RAW_MATRIX="${MAIN_OUTPUT_DIR}/temp_matrix.raw"
    
    # V2 输出文件与目录
    readonly SNP_MATRIX_FILE="${MAIN_OUTPUT_DIR}/snp_matrix.txt" # Step 4
    readonly STATS_OUTPUT_DIR_STEP5="${MAIN_OUTPUT_DIR}/step5_initial_stats"
    # (V2.4) Step 5 现在有两个输出。定义两个
    readonly STATS_FILE_STEP5_ALL="${STATS_OUTPUT_DIR_STEP5}/all_triplet_statistics.csv" # 包含所有结果 (用于Step 8绘图)
    readonly STATS_FILE_STEP5_SIG="${STATS_OUTPUT_DIR_STEP5}/significant_triplets.csv" # *** (关键) 过滤后的结果 (用于Step 6) ***
    
    readonly POS_TRIPLETS_FILE_STEP6="${MAIN_OUTPUT_DIR}/triplets_with_positions_significant.csv" # Step 6 (重命名以反映输入)
    readonly FINAL_OUTPUT_DIR_STEP7="${MAIN_OUTPUT_DIR}/step7_branch_analysis"
    readonly PASSED_TRIPLETS_FILE_STEP7="${FINAL_OUTPUT_DIR_STEP7}/passed_quadruplets_results.csv" # Step 7
    readonly PLOTS_OUTPUT_DIR_STEP8="${MAIN_OUTPUT_DIR}/final_summary_plots" # Step 8
    
    # V2 子脚本 (使用新版 step5)
    readonly EXTRACT_SNP_SCRIPT="${BASE_SCRIPT_PATH%/}/step1_extract_unique_snp_list.py"
    readonly TRANS_MATRIX_SCRIPT="${BASE_SCRIPT_PATH%/}/step4_trans_matrix.py"
    readonly STEP5_STATS_SCRIPT="${BASE_SCRIPT_PATH%/}/step5_TAD_and_Gene_analysis.py" # 假设这是新脚本的名字，请确保它匹配
    readonly STEP6_ADD_POS_SCRIPT="${BASE_SCRIPT_PATH%/}/step6_add_gene_and_tad_pos_info.py"
    readonly STEP7_BRANCH_SCRIPT="${BASE_SCRIPT_PATH%/}/step7_distance_and_loop_analysis.py"
    readonly STEP8_PLOT_SCRIPT="${BASE_SCRIPT_PATH%/}/step8_plot.R"

    # --- (V2.3) 固定依赖文件 (数据库路径已修正) ---
    readonly DATABASE_BASE_PATH="/share/home/cyzhao24/pan3D_267sample/data_for_ciwas_find_case_pipeline"

    # 1. VCF
    readonly VCF_FILE="${DATABASE_BASE_PATH}/Genotype/267filtered_13_AARR_vcf_filter_maf005.vcf.gz"
    # 2. TAD IS 矩阵
    readonly TAD_IS_FILE="${DATABASE_BASE_PATH}/TAD/267sample_final_tad_is.tsv"
    # 3. TAD 位置
    readonly TAD_POS_FILE="${DATABASE_BASE_PATH}/TAD/267sample_final_tad_pos_info.tsv"
    # 4. TAD 0-1 矩阵
    readonly TAD_01_FILE="${DATABASE_BASE_PATH}/TAD/final_TAD_01_matrix.tsv"
    # 5. Loop 0-1 矩阵
    readonly LOOP_PRESENCE_FILE="${DATABASE_BASE_PATH}/Loop/final_loop_01matrix.tsv"
    # 6. Loop 强度
    readonly LOOP_STRENGTH_FILE="${DATABASE_BASE_PATH}/Loop/final_loop_strength.tsv"
    # 7. Loop 位置
    readonly SORTED_LOOPS_FILE="${DATABASE_BASE_PATH}/Loop/final_loop_pos_info.tsv"
    # 8. Gene 位置
    readonly GENE_POS_FILE="${DATABASE_BASE_PATH}/Gene/final_gene_pos_info.tsv"
    # 9. Gene 表达 (TPM)
    readonly GENE_EXP_FILE="${DATABASE_BASE_PATH}/Gene/final_gene_TPM.tsv"
    # 10. 表型 (FL/FS)
    readonly PHENOTYPE_FILE="${DATABASE_BASE_PATH}/Phenotype/phenotypes_transposed.tsv"
    # --- (路径更新结束) ---


    # --- 检查文件存在性 (已更新) ---
    check_files_exist \
        "${TRIPLET_INPUT_FILE}" \
        "${EXTRACT_SNP_SCRIPT}" \
        "${TRANS_MATRIX_SCRIPT}" \
        "${STEP5_STATS_SCRIPT}" \
        "${STEP6_ADD_POS_SCRIPT}" \
        "${STEP7_BRANCH_SCRIPT}" \
        "${STEP8_PLOT_SCRIPT}" \
        "${VCF_FILE}" \
        "${TAD_POS_FILE}" \
        "${TAD_IS_FILE}" \
        "${TAD_01_FILE}" \
        "${GENE_POS_FILE}" \
        "${GENE_EXP_FILE}" \
        "${PHENOTYPE_FILE}" \
        "${SORTED_LOOPS_FILE}" \
        "${LOOP_STRENGTH_FILE}" \
        "${LOOP_PRESENCE_FILE}"

    # ========================== 2. 流水线阶段执行 (V2.4 逻辑) ==========================
    print_header "Starting Pipeline Execution (Rebuilt V2.4 Flow)"

    # --- 阶段 1: 提取SNP (硬编码 --has_header) ---
    run_stage 1 "Extracting & Sorting SNPs" "python3 '${EXTRACT_SNP_SCRIPT}' \
        -i '${TRIPLET_INPUT_FILE}' \
        -o '${TEMP_UNIQUE_SNP_LIST}' \
        -c 3 \
        --has_header" # <-- 根据要求硬编码

    # --- 阶段 2-4: (保持不变) ---
    run_stage 2 "Filtering VCF File" "plink2 --vcf '${VCF_FILE}' --extract '${TEMP_UNIQUE_SNP_LIST}' --export vcf-4.2 bgz --allow-extra-chr --out '${TEMP_EXTRACTED_VCF%.vcf.gz}'"
    run_stage 3 "Generating Raw Matrix" "plink2 --vcf '${TEMP_EXTRACTED_VCF}' --export A include-alt --allow-extra-chr --out '${TEMP_RAW_MATRIX%.raw}'"
    run_stage 4 "Transposing Matrix" "python3 '${TRANS_MATRIX_SCRIPT}' -i '${TEMP_RAW_MATRIX}' -o '${SNP_MATRIX_FILE}'"

    # --- 阶段 5: 初始统计与过滤 (V2.0.4 脚本) ---
    # (已更新: 传递2个新的P值阈值参数)
    run_stage 5 "Initial Stats & Filtering (Py v2.0.4)" "python3 '${STEP5_STATS_SCRIPT}' \
        '${TRIPLET_INPUT_FILE}' \
        '${SNP_MATRIX_FILE}' \
        '${GENE_EXP_FILE}' \
        '${TAD_IS_FILE}' \
        '${TAD_01_FILE}' \
        '${PHENOTYPE_FILE}' \
        '${NUM_CORES}' \
        --output_dir '${STATS_OUTPUT_DIR_STEP5}' \
        --tad_is_p_threshold ${PARAM_TAD_IS_PVAL} \
        --gene_exp_p_threshold ${PARAM_GENE_EXP_PVAL}" # <-- 两个新参数已应用

    # 检查主输出文件 (all_...csv) 是否存在
    if [ ! -f "${STATS_FILE_STEP5_ALL}" ]; then
        log_warn "Stage 5 produced no main output file (all_...csv). Pipeline terminating early."
        exit 0
    fi
    
    # (关键检查) 检查过滤后的文件 (sig_...csv) 是否存在且非空
    if [ ! -s "${STATS_FILE_STEP5_SIG}" ]; then
        log_warn "Stage 5 produced 0 significant candidates (file: ${STATS_FILE_STEP5_SIG}). Pipeline terminating."
        # 我们可以在这里正常退出，因为没有数据可以进入Step 6
        exit 0
    else
        log_info "Stage 5 found significant candidates. Proceeding with filtered list."
    fi


    # --- 阶段 6: 添加位置信息 (*** 关键修改: 输入已更改为过滤后的文件 ***) ---
    run_stage 6 "Adding Gene/TAD Positions (to Significant)" "python3 '${STEP6_ADD_POS_SCRIPT}' \
        '${STATS_FILE_STEP5_SIG}' \
        '${TAD_POS_FILE}' \
        '${GENE_POS_FILE}' && \
        mv triplets_with_positions.csv '${POS_TRIPLETS_FILE_STEP6}'" # <-- 输入已更改为 SIG 文件

    # --- 阶段 7: 核心机制分析 (输入已更改为过滤后的位置文件) ---
    run_stage 7 "Running Branch 1/2 Logic (Py)" "python3 '${STEP7_BRANCH_SCRIPT}' \
        '${POS_TRIPLETS_FILE_STEP6}' \
        '${SNP_MATRIX_FILE}' \
        '${GENE_EXP_FILE}' \
        '${TAD_IS_FILE}' \
        '${TAD_01_FILE}' \
        '${SORTED_LOOPS_FILE}' \
        '${LOOP_STRENGTH_FILE}' \
        '${LOOP_PRESENCE_FILE}' \
        '${NUM_CORES}' \
        --output_dir '${FINAL_OUTPUT_DIR_STEP7}' \
        --distance_threshold ${PARAM_DIST_THRESHOLD} \
        --max_distance ${PARAM_MAX_DIST} \
        --gene_anchor_distance ${PARAM_ANCHOR_DIST} \
        --loop_p_threshold ${PARAM_LOOP_PVAL}" # <-- 4个保留的参数已应用

    # --- 阶段 8: 生成汇总图 (*** 保持不变: 绘图仍需 'all_stats' 文件作为P值查找表 ***) ---
    if [ -s "${PASSED_TRIPLETS_FILE_STEP7}" ]; then
         log_info "Stage 7 found passed candidates. Proceeding to plotting."
         run_stage 8 "Generating Summary Plots (R)" "Rscript '${STEP8_PLOT_SCRIPT}' \
            --triplets '${PASSED_TRIPLETS_FILE_STEP7}' \
            --p_values '${STATS_FILE_STEP5_ALL}' \
            --genotype '${SNP_MATRIX_FILE}' \
            --gene_exp '${GENE_EXP_FILE}' \
            --tad_is '${TAD_IS_FILE}' \
            --phenotypes '${PHENOTYPE_FILE}' \
            --loop_strength '${LOOP_STRENGTH_FILE}' \
            --loop_status '${LOOP_PRESENCE_FILE}' \
            --output_dir '${PLOTS_OUTPUT_DIR_STEP8}'" # <-- 注意: --p_values 仍使用 ALL stats 文件
    else
        log_warn "Stage 7 produced no passed candidates. Skipping final plotting (Stage 8)."
    fi

    # ========================== 3. 流程结束报告 (已更新) ==========================
    print_header "Pipeline Finished"
    log_success "All tasks completed successfully."
    log_info "All results, logs, and temporary files are stored in: ${C_YELLOW}${MAIN_OUTPUT_DIR}${C_OFF}"

    cat << EOF | tee -a "${PIPELINE_LOG_FILE}" >&2

  Key Outputs (Rebuilt Flow V2.4):
  ──────────────────────────────────────────────────────────
  Initial Stats (All):     ${STATS_FILE_STEP5_ALL}
  Initial Stats (Filtered): ${STATS_FILE_STEP5_SIG}
  Mechanism Analysis:      ${FINAL_OUTPUT_DIR_STEP7}/
  Final Passed List:       ${PASSED_TRIPLETS_FILE_STEP7}
  Summary Plots:           ${PLOTS_OUTPUT_DIR_STEP8}/
  ──────────────────────────────────────────────────────────
  All detailed logs are in: logs

EOF
}

# --- 脚本执行入口 ---
main "$@"