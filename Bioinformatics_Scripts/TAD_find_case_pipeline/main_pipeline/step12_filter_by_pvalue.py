#!/usr/bin/env python3
import pandas as pd
import argparse
import sys
import re

def filter_by_p_values(input_file: str, output_file: str, p_threshold: float = 0.05):
    """
    根据TAD, Gene和Phenotype的P值对最终的汇总表进行筛选。

    筛选条件:
    1. (Wilcox_P_TAD < 阈值) AND (Wilcox_P_Gene < 阈值)
    2. AND
    3. (Wilcox_P_FL < 阈值) OR (Wilcox_P_FS < 阈值)
    """
    try:
        print(f"--- 开始执行最终P值筛选 (阈值 P < {p_threshold}) ---", file=sys.stderr)
        df = pd.read_csv(input_file)
        print(f"成功加载文件 '{input_file}'，共 {len(df)} 行。", file=sys.stderr)
    except FileNotFoundError:
        print(f"错误: 输入文件未找到 '{input_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"读取文件时发生错误: {e}", file=sys.stderr)
        sys.exit(1)

    initial_rows = len(df)
    condition1 = (df['Wilcox_P_TAD'] < p_threshold) & (df['Wilcox_P_Gene'] < p_threshold)
    condition2 = (df['Wilcox_P_FL'] < p_threshold) | (df['Wilcox_P_FS'] < p_threshold)
    final_mask = condition1 & condition2
    filtered_df = df[final_mask].copy() # 使用 .copy() 以避免 SettingWithCopyWarning
    final_rows = len(filtered_df)
    
    # 如果没有筛选结果，则直接退出
    if final_rows == 0:
        print("筛选后没有匹配的行，将创建一个空文件。", file=sys.stderr)
        open(output_file, 'w').close()
        sys.exit(0)

    # 提取TAD_ID和Gene_ID中的数字部分进行排序
    def extract_tad_id_number(tad_id):
        match = re.search(r'(\d+)', str(tad_id))
        return int(match.group(1)) if match else 0

    def extract_gene_id_number(gene_id):
        match = re.search(r'(\d+)$', str(gene_id))
        return int(match.group(1)) if match else 0
    
    # 应用提取函数创建用于排序的临时列
    # 注意：这里已经将 'Gene' 改为 'Gene_ID'
    filtered_df['tad_id_num'] = filtered_df['TAD_ID'].apply(extract_tad_id_number)
    filtered_df['gene_id_num'] = filtered_df['Gene_ID'].apply(extract_gene_id_number)
    
    # 进行排序：首先按TAD_ID数字，然后按Gene_ID数字
    filtered_df.sort_values(by=['tad_id_num', 'gene_id_num'], inplace=True)
    
    # 删除临时列
    filtered_df.drop(columns=['tad_id_num', 'gene_id_num'], inplace=True)

    try:
        filtered_df.to_csv(output_file, index=False)
        print(f"筛选完成。从 {initial_rows} 行中保留了 {final_rows} 行。", file=sys.stderr)
        print(f"数据已根据TAD_ID和Gene_ID排序，并保存至: '{output_file}'", file=sys.stderr)
    except Exception as e:
        print(f"保存文件时发生错误: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="根据P值对最终四元组汇总表进行筛选。")
    parser.add_argument('--input', type=str, required=True, help="第11步生成的汇总文件路径。")
    parser.add_argument('--output', type=str, required=True, help="最终筛选后输出文件的路径。")
    parser.add_argument('--p_threshold', type=float, default=0.05, help="P值显著性阈值 (默认为0.05)。")
    args = parser.parse_args()
    filter_by_p_values(args.input, args.output, args.p_threshold)