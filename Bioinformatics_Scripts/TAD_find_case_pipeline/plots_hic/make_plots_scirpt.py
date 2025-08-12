import os
import glob
import pandas as pd
import sys
from pathlib import Path
from collections import defaultdict

# --- 参数获取 ---
if len(sys.argv) < 3:
    print("用法: python 01_preprocess_and_generate_jobs.py <输入目录前缀> <输出基础目录>")
    sys.exit(1)

input_prefix = sys.argv[1]
output_base_dir = Path(sys.argv[2])
cool_base_path = Path("/share/home/xhhuang24/biodata/pan_3Dgenome/Domain/Intra_mtx/iced/4k")

# 确保输出目录和任务脚本目录存在
job_scripts_dir = output_base_dir / "job_scripts"
job_scripts_dir.mkdir(parents=True, exist_ok=True)

# --- 1. 读取染色体最大长度信息 ---
def load_chromosome_lengths(genome_file_path):
    """从文件中加载染色体长度信息"""
    chromosome_lengths = {}
    try:
        with open(genome_file_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 2:
                    chromosome_lengths[parts[0]] = int(parts[1])
    except FileNotFoundError:
        print(f"错误: 染色体长度文件不存在于 {genome_file_path}。请检查路径。")
        sys.exit(1)
    return chromosome_lengths

genome_file_path = "/share/home/cyzhao24/workplace/01.ML/02.loop_association/data/genome/l128.genome.lst"
chromosome_lengths = load_chromosome_lengths(genome_file_path)

# --- 2. 遍历目录，提取并处理 TAD-SNP 组合 ---
snp_matrix_cache = {}
target_dirs = glob.glob(f'./{input_prefix}*')
task_count = 0

print("开始遍历目录，为每个 TAD-SNP 组合生成任务脚本...")

for run_dir in target_dirs:
    run_dir_name = Path(run_dir).name
    summary_file = Path(run_dir) / "final_quartet_summary_with_TM1.csv"
    snp_matrix_file = Path(run_dir) / "snp_matrix.txt"

    if not summary_file.exists() or not snp_matrix_file.exists():
        continue
    
    if snp_matrix_file not in snp_matrix_cache:
        try:
            snp_matrix_cache[snp_matrix_file] = pd.read_csv(snp_matrix_file, sep='\t')
        except (pd.errors.EmptyDataError, FileNotFoundError):
            continue
    
    df_snp = snp_matrix_cache[snp_matrix_file]
    
    try:
        df_summary = pd.read_csv(summary_file)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        continue
    
    # 使用一个集合来跟踪当前 run_dir 中已处理的 TAD-SNP 组合
    processed_in_run_dir = set()

    for _, row in df_summary.iterrows():
        tad_id = str(row['TAD_ID'])
        full_snp_id = str(row['Full_SNP_ID'])
        tad_snp_key = (tad_id, full_snp_id)
        
        # 如果该组合在当前 run_dir 中尚未处理，则生成一个新脚本
        if tad_snp_key not in processed_in_run_dir:
            processed_in_run_dir.add(tad_snp_key)
            task_count += 1
            
            snp_data = df_snp[df_snp['SNP_ID'] == full_snp_id]
            if snp_data.empty:
                continue

            sample_cols = [col for col in df_snp.columns if col != 'SNP_ID']
            alt_alt_samples = [sample for sample in sample_cols if snp_data[sample].iloc[0] == 0]
            ref_ref_samples = [sample for sample in sample_cols if snp_data[sample].iloc[0] == 2]

            if not alt_alt_samples and not ref_ref_samples:
                continue
                
            chromosome = str(row['Chromosome'])
            tad_start = int(row['TAD_Start'])
            tad_end = int(row['TAD_End'])

            plot_start = max(0, tad_start - 480000)
            plot_end = min(chromosome_lengths.get(chromosome, tad_end + 500000), tad_end + 500000)
            
            # --- 修复 NameError：在这里定义 coordinate 变量 ---
            coordinate = f"{chromosome}:{plot_start}-{plot_end}"
            
            # 生成一个干净的文件名，替换特殊字符
            sanitized_snp_id = full_snp_id.replace(':', '_').replace('/', '_')
            job_script_path = job_scripts_dir / f"job_{run_dir_name}_{tad_id}_{sanitized_snp_id}.sh"

            with open(job_script_path, 'w') as job_file:
                job_file.write("#!/bin/bash\n\n")
                job_file.write(f"echo '正在处理 TAD-SNP 组合: {tad_id} - {full_snp_id}'\n")

                if alt_alt_samples:
                    output_dir = output_base_dir / run_dir_name / tad_id / full_snp_id / "ALT_ALT"
                    job_file.write(f"mkdir -p {output_dir}\n")
                    for sample in alt_alt_samples:
                        cool_file = cool_base_path / sample / f"{sample}_{chromosome}.cool"
                        output_file = output_dir / f"{sample}.png"
                        job_file.write(f"if [ -f {cool_file} ]; then cooler show {cool_file} {coordinate} -o {output_file}; fi\n")
                
                if ref_ref_samples:
                    output_dir = output_base_dir / run_dir_name / tad_id / full_snp_id / "REF_REF"
                    job_file.write(f"mkdir -p {output_dir}\n")
                    for sample in ref_ref_samples:
                        cool_file = cool_base_path / sample / f"{sample}_{chromosome}.cool"
                        output_file = output_dir / f"{sample}.png"
                        job_file.write(f"if [ -f {cool_file} ]; then cooler show {cool_file} {coordinate} -o {output_file}; fi\n")
            
            print(f"已生成任务脚本: {job_script_path}")

print(f"\n共生成 {task_count} 个任务脚本。请运行 Bash 脚本来提交任务。")
