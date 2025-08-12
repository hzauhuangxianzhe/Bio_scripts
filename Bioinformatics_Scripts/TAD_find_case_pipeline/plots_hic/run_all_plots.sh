#!/bin/bash

# --- 参数检查 ---
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "用法: $0 <输入目录前缀> <输出基础目录> [--test]"
    echo "  --test: 可选参数，只提交前两个任务进行测试"
    exit 1
fi

input_prefix=$1
output_base_dir=$2
test_mode=false

if [ "$#" -eq 3 ] && [ "$3" == "--test" ]; then
    test_mode=true
    echo "脚本已进入测试模式，将只提交前两个任务。"
fi

job_scripts_dir="${output_base_dir}/job_scripts"
dsub_log_dir="${output_base_dir}/dsub_logs"

# 确保 dsub 日志目录存在
mkdir -p "${dsub_log_dir}"
echo "已创建 dsub 日志目录: ${dsub_log_dir}"

# --- 运行 Python 脚本生成任务脚本文件 ---
echo "正在运行 Python 脚本进行数据预处理和任务脚本生成..."
# 确保你这里的 Python 脚本文件名与实际一致
python /share/home/cyzhao24/my_scripts/TAD_find_case_pipeline/plots_hic/make_plots_scirpt.py "${input_prefix}" "${output_base_dir}"

if [ ! -d "${job_scripts_dir}" ]; then
    echo "错误: Python脚本未能生成任务脚本目录 '${job_scripts_dir}'。"
    exit 1
fi

echo "任务脚本已生成，开始提交 dsub 任务..."

# --- 遍历任务脚本文件并提交 dsub 任务 ---
file_count=0
for job_script in "${job_scripts_dir}"/*.sh; do
    if [ -f "$job_script" ]; then
        if [ "$test_mode" == true ] && [ "$file_count" -ge 2 ]; then
            echo "---"
            echo "已完成测试模式下的前两个任务，退出循环。"
            break
        fi

        filename=$(basename -- "$job_script")
        filename_no_ext="${filename%.sh}"
        job_name="${filename_no_ext}" # 使用文件名作为任务名称
        
        echo "---"
        echo "正在提交任务: ${job_name}"
        
        # 修复了 dsub 命令，将 'bash ${job_script}' 放在了最后
        dsub -n "${job_name}" \
             -q 'default' \
             -R 'cpu=2' \
             -o "${dsub_log_dir}/${job_name}_%J.out" \
             -e "${dsub_log_dir}/${job_name}_%J.err" \
             bash "${job_script}"
        
        echo "已提交任务。"
        file_count=$((file_count + 1))
    fi
done

echo "所有任务已提交。请检查 ${dsub_log_dir} 目录查看任务状态。"
