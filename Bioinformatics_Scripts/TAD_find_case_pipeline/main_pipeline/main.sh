#!/usr/bin/env bash

# ==================================================================================================
#
#       SNP -> Gene/TAD/Phenotype -> Loop Association Pipeline
#
# ==================================================================================================
#
# 版本: 5.0 (Production Release)
# 日期: 2025-08-05
# 描述: 从SNP到染色质环作用机制的完整生物信息分析流程。
#
# 特性:
#   - 12个核心分析阶段。
#   - 生产级代码结构、日志和错误处理。
#   - 优雅的终端UI。
#
# ==================================================================================================

# --- 脚本配置: 严格模式 ---
set -e
set -o pipefail

####################################################################################################
# Section 1: 核心函数库
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

#
# @description 核心日志记录函数
#
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

#
# @description 打印带边框的标题
#
print_header() {
    local title=" $1 "
    local len=${#title}
    local border=$(printf '─%.0s' $(seq 1 "$len"))
    echo -e "\n${C_PURPLE}╭─${border}─╮${C_OFF}" | tee -a "${PIPELINE_LOG_FILE}" >&2
    echo -e "${C_PURPLE}│ ${C_BOLD}${title}${C_OFF}${C_PURPLE} │${C_OFF}" | tee -a "${PIPELINE_LOG_FILE}" >&2
    echo -e "${C_PURPLE}╰─${border}─╯${C_OFF}" | tee -a "${PIPELINE_LOG_FILE}" >&2
}

#
# @description 执行一个流程阶段
#
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
# Section 2: 流程辅助函数
####################################################################################################

#
# @description 显示脚本用法
#
usage() {
    cat << EOF

${C_GREEN}${C_BOLD}CIWAS-Case-Analyze-Pipeline V1.0${C_OFF}

${C_BOLD}DESCRIPTION:${C_OFF}
  一个从SNP出发，建立其通过影响TAD和染色质环(Loop)来调控基因表达和表型的分子机制证据链的完整生物信息分析流程。

${C_BOLD}USAGE:${C_OFF}
  ${SCRIPT_NAME} -s <triplet_file> -o <output_prefix> -n <num_cores> [-c]

${C_BOLD}OPTIONS:${C_OFF}
  -s <file>     ${C_YELLOW}必需。${C_OFF}初始三元组文件路径。
  -o <prefix>   ${C_YELLOW}必需。${C_OFF}本次运行的唯一前缀，将作为主输出目录名。
  -n <cores>    ${C_YELLOW}必需。${C_OFF}R脚本并行计算时使用的CPU核心数。
  -c            ${C_YELLOW}可选。${C_OFF}流程结束后，清理所有临时文件。

EOF
    exit 1
}

#
# @description 检查核心外部命令依赖
#
check_dependencies() {
    print_header "Checking Dependencies"
    log_info "Verifying required command-line tools..."
    local dependencies=("plink2" "Rscript" "python3" "awk" "tee")
    for cmd in "${dependencies[@]}"; do
        if ! command -v "$cmd" &>/dev/null; then
            log_error "Required command not found in PATH: ${cmd}"
            exit 1
        fi
    done
    log_success "All command-line dependencies are satisfied."
}

#
# @description 检查所有输入文件
#
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

#
# @description 流程结束时清理临时文件
#
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
# Section 3: 主逻辑 (Main)
####################################################################################################

main() {
    # --- 解析命令行参数 ---
    local TRIPLET_INPUT_FILE=""
    local OUTPUT_PREFIX=""
    local NUM_CORES=""
    local PERFORM_CLEANUP=false

    while getopts "s:o:n:c" opt; do
        case ${opt} in
            s) TRIPLET_INPUT_FILE=$OPTARG ;;
            o) OUTPUT_PREFIX=$OPTARG ;;
            n) NUM_CORES=$OPTARG ;;
            c) PERFORM_CLEANUP=true ;;
            *) usage ;;
        esac
    done

    if [ -z "${TRIPLET_INPUT_FILE}" ] || [ -z "${OUTPUT_PREFIX}" ] || [ -z "${NUM_CORES}" ]; then
        log_error "Missing required arguments."
        usage
    fi

    # ========================== 1. 环境与路径设置 ==========================

    # --- 核心目录与常量 ---
    readonly MAIN_OUTPUT_DIR="${OUTPUT_PREFIX}"
    readonly LOG_DIR="${MAIN_OUTPUT_DIR}/logs"
    readonly BASE_SCRIPT_PATH="/share/home/cyzhao24/workplace/01.ML/02.loop_association/my_scripts/TAD_find_case_pipeline/main_pipeline"
    readonly TOTAL_STAGES=12

    # --- 初始化环境 ---
    mkdir -p "${MAIN_OUTPUT_DIR}"
    mkdir -p "${LOG_DIR}"
    PIPELINE_LOG_FILE="${LOG_DIR}/pipeline_main.log"
    > "${PIPELINE_LOG_FILE}" # 清空旧日志
    trap cleanup EXIT INT ERR

    # --- 打印欢迎横幅和配置信息 ---
    print_header "CIWAS-Case-Analyze-Pipeline V1.0"
    log_info "Run initiated on $(date)"
    log_info "──────────────────────────────────────────────────────────"
    log_info "Input Triplets : ${C_YELLOW}${TRIPLET_INPUT_FILE}${C_OFF}"
    log_info "Output Directory : ${C_YELLOW}${MAIN_OUTPUT_DIR}${C_OFF}"
    log_info "CPU Cores      : ${C_YELLOW}${NUM_CORES}${C_OFF}"
    log_info "Cleanup Temp   : ${C_YELLOW}${PERFORM_CLEANUP}${C_OFF}"
    log_info "──────────────────────────────────────────────────────────"

    # --- 检查依赖 ---
    check_dependencies

    # --- 路径注册表 (按类别组织) ---
    # 临时文件
    readonly TEMP_UNIQUE_SNP_LIST="${MAIN_OUTPUT_DIR}/temp_unique_snps.txt"
    readonly TEMP_EXTRACTED_VCF="${MAIN_OUTPUT_DIR}/temp_extracted.vcf.gz"
    readonly TEMP_RAW_MATRIX="${MAIN_OUTPUT_DIR}/temp_matrix.raw"
    # 输出文件与目录
    readonly PLOTS_OUTPUT_DIR_STEP5="${MAIN_OUTPUT_DIR}/snp_pheno_boxplots"
    readonly FINAL_OUTPUT_DIR_STEP7="${MAIN_OUTPUT_DIR}/final_loop_analysis"
    readonly FINAL_VISUAL_DIR_STEP10="${MAIN_OUTPUT_DIR}/final_visual_summary"
    readonly SNP_MATRIX_FILE="${MAIN_OUTPUT_DIR}/snp_matrix.txt"
    readonly STATS_FILE_STEP5="${PLOTS_OUTPUT_DIR_STEP5}/filtered_plotted_statistics.csv"
    readonly LOOP_ANALYSIS_FILE_STEP6="${MAIN_OUTPUT_DIR}/gene_tad_loop_analysis.csv"
    readonly FINAL_SUMMARY_CSV_STEP7="${FINAL_OUTPUT_DIR_STEP7}/final_pvalue_summary.csv"
    readonly FILTERED_LOOP_ANALYSIS_FILE_STEP8="${MAIN_OUTPUT_DIR}/final_filtered_loop_analysis.csv"
    readonly QUARTET_SUMMARY_FILE_STEP9="${MAIN_OUTPUT_DIR}/final_quartet_summary.csv"
    readonly QUARTET_SUMMARY_FILE_STEP11="${MAIN_OUTPUT_DIR}/final_quartet_summary_with_TM1.csv"
    readonly FINAL_FILTERED_QUARTET_SUMMARY_STEP12="${MAIN_OUTPUT_DIR}/final_pvalue_filtered_quartet_summary.csv"
    # 子脚本
    readonly EXTRACT_SNP_SCRIPT="${BASE_SCRIPT_PATH%/}/step1_extract_unique_snp_list.py"
    readonly TRANS_MATRIX_SCRIPT="${BASE_SCRIPT_PATH%/}/step4_trans_matrix.py"
    readonly DRAW_BOXPLOT_SCRIPT_STEP5="${BASE_SCRIPT_PATH%/}/step5_draw_boxplot_Molecule_Phenotype_with_SNP_parallel.R"
    readonly FIND_LOOPS_SCRIPT_STEP6="${BASE_SCRIPT_PATH%/}/step6_find_gene_tad_loops.py"
    readonly ANALYZE_LOOPS_SCRIPT_STEP7="${BASE_SCRIPT_PATH%/}/step7_analyze_and_plot_loops.R"
    readonly FILTER_SCRIPT_STEP8="${BASE_SCRIPT_PATH%/}/step8_filter_loops_final.py"
    readonly QUARTET_SCRIPT_STEP9="${BASE_SCRIPT_PATH%/}/step9_create_quartet_summary.py"
    readonly COLLECT_VISUAL_SCRIPT_STEP10="${BASE_SCRIPT_PATH%/}/step10_collect_final_plots.sh"
    readonly ADD_TM1_ID_SCRIPT="${BASE_SCRIPT_PATH%/}/step11_add_tm1_gene_id.py"
    readonly FILTER_PVALUE_SCRIPT_STEP12="${BASE_SCRIPT_PATH%/}/step12_filter_by_pvalue.py"
    # 固定依赖文件
    readonly VCF_FILE="/share/home/cyzhao24/workplace/01.ML/02.loop_association/coloc/3d_with_gene_new/tad_gene/L_maf005_filtered_samples_final.vcf.gz"
    readonly TAD_POS_FILE="/share/home/cyzhao24/workplace/01.ML/02.loop_association/data/TAD/all_TAD_bin_with_pos_and_it_not_ISscore.bed"
    readonly GENE_POS_FILE="/share/home/cyzhao24/workplace/01.ML/02.loop_association/coloc/3d_with_gene_new/tad_gene_vaild/raw_data/formatted_rna_data_for_qtl.bed"
    readonly SORTED_LOOPS_FILE="/share/home/cyzhao24/workplace/01.ML/02.loop_association/coloc/3d_with_gene_new/loop_snp/sorted_loops.bed"
    readonly LOOP_STRENGTH_FILE="/share/home/cyzhao24/workplace/01.ML/02.loop_association/coloc/fusion/data_for_run_fusion/loop_strength/sorted_loops_filtered_244.bed"
    readonly LOOP_PRESENCE_FILE="/share/home/cyzhao24/workplace/01.ML/02.loop_association/coloc/3d_with_gene_new/loop_snp/loop_presence_matrix.csv"
    readonly TM1_HC04_GENE_PAIR_FILE="/share/home/cyzhao24/workplace/01.ML/02.loop_association/data/HC04_and_TM1/TM1_HC04_gene_pair.bed"

    # --- 检查文件存在性 ---
    check_files_exist "${EXTRACT_SNP_SCRIPT}" "${TRANS_MATRIX_SCRIPT}" "${DRAW_BOXPLOT_SCRIPT_STEP5}" "${FIND_LOOPS_SCRIPT_STEP6}" "${ANALYZE_LOOPS_SCRIPT_STEP7}" "${FILTER_SCRIPT_STEP8}" "${QUARTET_SCRIPT_STEP9}" "${COLLECT_VISUAL_SCRIPT_STEP10}" "${ADD_TM1_ID_SCRIPT}" "${FILTER_PVALUE_SCRIPT_STEP12}" "${TRIPLET_INPUT_FILE}" "${VCF_FILE}" "${TAD_POS_FILE}" "${GENE_POS_FILE}" "${SORTED_LOOPS_FILE}" "${LOOP_STRENGTH_FILE}" "${LOOP_PRESENCE_FILE}" "${TM1_HC04_GENE_PAIR_FILE}"

    # ========================== 2. 流水线阶段执行 ==========================
    print_header "Starting Pipeline Execution"

    run_stage 1 "Extracting & Sorting SNPs" "python3 '${EXTRACT_SNP_SCRIPT}' -i '${TRIPLET_INPUT_FILE}' -o '${TEMP_UNIQUE_SNP_LIST}' -c 3 --has_header"
    run_stage 2 "Filtering VCF File" "plink2 --vcf '${VCF_FILE}' --extract '${TEMP_UNIQUE_SNP_LIST}' --export vcf-4.2 bgz --allow-extra-chr --out '${TEMP_EXTRACTED_VCF%.vcf.gz}'"
    run_stage 3 "Generating Raw Matrix" "plink2 --vcf '${TEMP_EXTRACTED_VCF}' --export A include-alt --allow-extra-chr --out '${TEMP_RAW_MATRIX%.raw}'"
    run_stage 4 "Transposing Matrix" "python3 '${TRANS_MATRIX_SCRIPT}' -i '${TEMP_RAW_MATRIX}' -o '${SNP_MATRIX_FILE}'"
    run_stage 5 "Initial Association Analysis" "Rscript '${DRAW_BOXPLOT_SCRIPT_STEP5}' '${TRIPLET_INPUT_FILE}' '${SNP_MATRIX_FILE}' '${PLOTS_OUTPUT_DIR_STEP5}' '${NUM_CORES}'"

    if [ ! -f "${STATS_FILE_STEP5}" ]; then
        log_warn "Stage 5 produced no candidates. Pipeline terminating early."
        exit 0
    fi

    run_stage 6 "Finding Associated Loops" "python3 '${FIND_LOOPS_SCRIPT_STEP6}' --stats_file '${STATS_FILE_STEP5}' --tad_pos_file '${TAD_POS_FILE}' --gene_pos_file '${GENE_POS_FILE}' --loop_file '${SORTED_LOOPS_FILE}' --output_file '${LOOP_ANALYSIS_FILE_STEP6}'"

    loop_count=$(awk -F',' 'NR > 1 && $9 != "NA" && $9 != "" {count++} END{print count+0}' "${LOOP_ANALYSIS_FILE_STEP6}")
    if [ "${loop_count}" -gt 0 ]; then
        log_info "Found ${loop_count} records with loops. Proceeding with Stage 7."
        run_stage 7 "Loop Association Analysis" "Rscript '${ANALYZE_LOOPS_SCRIPT_STEP7}' '${LOOP_ANALYSIS_FILE_STEP6}' '${STATS_FILE_STEP5}' '${LOOP_STRENGTH_FILE}' '${SNP_MATRIX_FILE}' '${LOOP_PRESENCE_FILE}' '${FINAL_OUTPUT_DIR_STEP7}' '${FINAL_SUMMARY_CSV_STEP7}' '${NUM_CORES}'"
    else
        log_warn "No records with loops found. Skipping Stage 7 analysis."
        mkdir -p "${FINAL_OUTPUT_DIR_STEP7}"
        echo "Gene_ID,TAD_ID,Loop_ID,Full_SNP_ID,Wilcoxon_P_Value_Loop_Strength,Categorical_P_Value_Loop_Presence" > "${FINAL_SUMMARY_CSV_STEP7}"
        log_info "Created empty placeholder for p-value summary."
    fi

    run_stage 8 "Filtering by Distance & P-value" "python3 '${FILTER_SCRIPT_STEP8}' --loop_analysis_file '${LOOP_ANALYSIS_FILE_STEP6}' --pvalue_summary_file '${FINAL_SUMMARY_CSV_STEP7}' --output_file '${FILTERED_LOOP_ANALYSIS_FILE_STEP8}'"
    run_stage 9 "Generating Quartet Summary" "python3 '${QUARTET_SCRIPT_STEP9}' --filtered_loop_file '${FILTERED_LOOP_ANALYSIS_FILE_STEP8}' --stats_file '${STATS_FILE_STEP5}' --pvalue_file '${FINAL_SUMMARY_CSV_STEP7}' --output_file '${QUARTET_SUMMARY_FILE_STEP9}'"

    if [ ! -s "${QUARTET_SUMMARY_FILE_STEP9}" ]; then
        log_warn "Stage 9 did not produce a valid summary file. Skipping final steps."
    else
        run_stage 10 "Collecting Visual Results" "bash '${COLLECT_VISUAL_SCRIPT_STEP10}' -q '${QUARTET_SUMMARY_FILE_STEP9}' -s '${PLOTS_OUTPUT_DIR_STEP5}' -d '${FINAL_OUTPUT_DIR_STEP7}' -o '${FINAL_VISUAL_DIR_STEP10}'"
        run_stage 11 "Adding TM1 Gene IDs" "python3 '${ADD_TM1_ID_SCRIPT}' '${QUARTET_SUMMARY_FILE_STEP9}' '${TM1_HC04_GENE_PAIR_FILE}' '${QUARTET_SUMMARY_FILE_STEP11}'"
        run_stage 12 "Final P-value Filtering" "python3 '${FILTER_PVALUE_SCRIPT_STEP12}' --input '${QUARTET_SUMMARY_FILE_STEP11}' --output '${FINAL_FILTERED_QUARTET_SUMMARY_STEP12}'"
    fi

    # ========================== 3. 流程结束报告 ==========================
    print_header "Pipeline Finished"
    log_success "All tasks completed successfully."
    log_info "All results, logs, and temporary files are stored in: ${C_YELLOW}${MAIN_OUTPUT_DIR}${C_OFF}"

    cat << EOF | tee -a "${PIPELINE_LOG_FILE}" >&2

  Key Outputs:
  ──────────────────────────────────────────────────────────
  Initial Plots:       snp_pheno_boxplots
  Loop Analysis:       final_loop_analysis
  Visual Summary:      final_visual_summary
  Annotated Summary:   final_quartet_summary_with_TM1.csv
  Final Filtered List: final_pvalue_filtered_quartet_summary.csv
  ──────────────────────────────────────────────────────────
  All detailed logs are in: logs

EOF
}

# --- 脚本执行入口 ---
main "$@"