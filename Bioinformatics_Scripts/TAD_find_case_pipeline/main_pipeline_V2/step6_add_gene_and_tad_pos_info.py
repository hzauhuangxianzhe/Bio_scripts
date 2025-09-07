#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import sys
from argparse import ArgumentParser

def add_position_info(triplets_file: str, tad_pos_file: str, gene_pos_file: str):
    """
    为三元组文件添加TAD和基因的位置信息，并保存为新的CSV文件。

    Args:
        triplets_file (str): 包含三元组（Gene_ID, TAD_ID, SNP_ID）的CSV文件路径。
        tad_pos_file (str): 包含TAD位置信息的文件路径。
        gene_pos_file (str): 包含基因位置信息的文件路径。
    """
    print("开始加载位置信息文件...", file=sys.stderr)

    # 加载TAD位置信息文件
    try:
        print(f"正在读取TAD文件: {tad_pos_file}", file=sys.stderr)
        tad_df = pd.read_csv(tad_pos_file, sep='\t')
        print(f"TAD文件列名: {tad_df.columns.tolist()}", file=sys.stderr)
        
        if 'TAD_ID' not in tad_df.columns:
            print(f"错误：TAD文件中没有找到 'TAD_ID' 列", file=sys.stderr)
            return
        
        tad_df.set_index('TAD_ID', inplace=True)
        # 重命名列以避免与基因位置信息列名冲突
        tad_df.rename(columns={'chrom': 'TAD_chrom', 'start': 'TAD_start', 'end': 'TAD_end'}, inplace=True)
        print(f"TAD文件加载成功，共{len(tad_df)}行", file=sys.stderr)
        
    except Exception as e:
        print(f"加载TAD文件时出错: {e}", file=sys.stderr)
        return
    
    # 加载基因位置信息文件
    try:
        print(f"正在读取基因文件: {gene_pos_file}", file=sys.stderr)
        gene_df = pd.read_csv(gene_pos_file, sep='\t')
        print(f"基因文件列名: {gene_df.columns.tolist()}", file=sys.stderr)
        print(f"基因文件前3行:\n{gene_df.head(3)}", file=sys.stderr)
        
        # 清理列名（去除可能的空格）
        gene_df.columns = gene_df.columns.str.strip()
        
        if 'Gene_ID' not in gene_df.columns:
            print(f"错误：基因文件中没有找到 'Gene_ID' 列", file=sys.stderr)
            print(f"可用的列名: {gene_df.columns.tolist()}", file=sys.stderr)
            
            # 尝试常见的基因ID列名变体
            possible_gene_cols = ['gene_id', 'GeneID', 'Gene_Name', 'gene_name', 'ID']
            found_col = None
            for col in possible_gene_cols:
                if col in gene_df.columns:
                    found_col = col
                    break
            
            if found_col:
                print(f"找到可能的基因ID列: {found_col}，将其重命名为Gene_ID", file=sys.stderr)
                gene_df.rename(columns={found_col: 'Gene_ID'}, inplace=True)
            else:
                print("无法找到合适的基因ID列，请检查文件格式", file=sys.stderr)
                return
        
        gene_df.set_index('Gene_ID', inplace=True)
        # 重命名列以避免与TAD位置信息列名冲突
        gene_df.rename(columns={'chrom': 'Gene_chrom', 'start': 'Gene_start', 'end': 'Gene_end'}, inplace=True)
        print(f"基因文件加载成功，共{len(gene_df)}行", file=sys.stderr)
        
    except Exception as e:
        print(f"加载基因文件时出错: {e}", file=sys.stderr)
        return

    print("位置信息文件加载完成。开始加载三元组文件...", file=sys.stderr)

    # 加载三元组文件
    try:
        triplets_df = pd.read_csv(triplets_file)
        print(f"三元组文件加载成功，共{len(triplets_df)}行", file=sys.stderr)
        print(f"三元组文件列名: {triplets_df.columns.tolist()}", file=sys.stderr)
    except Exception as e:
        print(f"加载三元组文件时出错: {e}", file=sys.stderr)
        return

    print("三元组文件加载完成。开始合并数据...", file=sys.stderr)

    # 根据TAD_ID合并TAD位置信息
    try:
        merged_df = triplets_df.join(tad_df, on='TAD_ID')
        print(f"TAD信息合并完成", file=sys.stderr)
    except Exception as e:
        print(f"合并TAD信息时出错: {e}", file=sys.stderr)
        return
    
    # 根据Gene_ID合并基因位置信息
    try:
        merged_df = merged_df.join(gene_df, on='Gene_ID')
        print(f"基因信息合并完成", file=sys.stderr)
    except Exception as e:
        print(f"合并基因信息时出错: {e}", file=sys.stderr)
        return

    print("数据合并完成。开始保存结果...", file=sys.stderr)

    # 定义新的列顺序，使输出文件更具可读性
    # 只包含实际存在的列
    base_columns = [
        "Gene_ID", "TAD_ID", "SNP_ID"
    ]
    
    position_columns = [
        "Gene_chrom", "Gene_start", "Gene_end",
        "TAD_chrom", "TAD_start", "TAD_end"
    ]
    
    # 获取其他所有列
    other_columns = [col for col in merged_df.columns 
                    if col not in base_columns + position_columns]
    
    # 构建最终的列顺序
    output_columns = base_columns + position_columns + other_columns
    
    # 只保留实际存在的列
    output_columns = [col for col in output_columns if col in merged_df.columns]
    
    merged_df = merged_df.reindex(columns=output_columns)
    
    output_file = "triplets_with_positions.csv"
    # 保存为新的CSV文件，使用逗号作为分隔符
    merged_df.to_csv(output_file, index=False)
    
    print(f"结果已成功保存到文件: {output_file}", file=sys.stderr)
    print(f"输出文件包含 {len(merged_df)} 行数据，{len(merged_df.columns)} 列", file=sys.stderr)

def main():
    """
    主函数，处理命令行参数并调用核心功能。
    """
    parser = ArgumentParser(description="为三元组文件添加TAD和基因的位置信息。")
    parser.add_argument("triplets_file", help="significant_triplets.csv文件路径。")
    parser.add_argument("tad_pos_file", help="TAD位置信息文件路径。")
    parser.add_argument("gene_pos_file", help="基因位置信息文件路径。")
    
    args = parser.parse_args()
    
    add_position_info(args.triplets_file, args.tad_pos_file, args.gene_pos_file)

if __name__ == "__main__":
    main()
