#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import sys

def filter_loops_by_conditions(loop_analysis_file, pvalue_summary_file, output_file, p_threshold=0.05):
    """
    根据两阶段的复杂条件筛选Loop分析结果。
    """
    try:
        print("--- 阶段1 (Python子脚本): 正在加载输入文件... ---", file=sys.stderr)
        # Use more robust read_csv options
        loop_df = pd.read_csv(loop_analysis_file, encoding='utf-8')
        pvalue_df = pd.read_csv(pvalue_summary_file, encoding='utf-8')
        print(f"成功加载 Loop分析文件，共 {len(loop_df)} 行。", file=sys.stderr)
        print(f"成功加载 P值汇总文件，共 {len(pvalue_df)} 行。", file=sys.stderr)

    except FileNotFoundError as e:
        print(f"错误：找不到文件 {e.filename}。请检查文件路径。", file=sys.stderr)
        sys.exit(1)

    print("\n--- 阶段2 (Python子脚本): 开始执行第一阶段筛选 (基于距离)... ---", file=sys.stderr)
    
    # --- 第一阶段筛选 ---
    loop_df['distance_numeric_kb'] = np.nan
    
    kb_mask = loop_df['Distance_Relationship'].str.contains('kb', na=False)
    loop_df.loc[kb_mask, 'distance_numeric_kb'] = pd.to_numeric(
        loop_df.loc[kb_mask, 'Distance_Relationship'].str.replace('kb|\\+', '', regex=True),
        errors='coerce'
    )
    
    bp_mask = loop_df['Distance_Relationship'].str.contains('bp', na=False)
    loop_df.loc[bp_mask, 'distance_numeric_kb'] = pd.to_numeric(
        loop_df.loc[bp_mask, 'Distance_Relationship'].str.replace('bp|\\+', '', regex=True),
        errors='coerce'
    ) / 1000

    cond1_sub1 = (loop_df['Distance_Relationship'] == 'Gene in TAD')
    cond1_sub2 = (loop_df['distance_numeric_kb'] > -50) & (loop_df['distance_numeric_kb'] < 50)
    final_condition1 = cond1_sub1 | cond1_sub2
    
    df_kept_by_cond1 = loop_df[final_condition1].copy()
    df_remaining = loop_df[~final_condition1].copy()
    
    print(f"第一阶段筛选完成: {len(df_kept_by_cond1)} 行因满足距离条件而被直接保留。", file=sys.stderr)

    print("\n--- 阶段3 (Python子脚本): 开始执行第二阶段筛选 (基于Loop和P值)... ---", file=sys.stderr)
    
    # --- 第二阶段筛选 (仅对第一阶段未被保留的数据进行) ---
    
    candidates_cond2 = df_remaining[df_remaining['Matched_Loop_ID'].notna()].copy()
    
    significant_pvalues = pvalue_df[
        (pvalue_df['Wilcoxon_P_Value_Loop_Strength'] < p_threshold) |
        (pvalue_df['Categorical_P_Value_Loop_Presence'] < p_threshold)
    ].copy()
    
    # 提取显著记录中唯一的 Gene-TAD-Loop 组合
    significant_loops_to_keep = significant_pvalues[['Gene_ID', 'TAD_ID', 'Loop_ID']].drop_duplicates()
    
    # --- 最终修复: 强制清理所有匹配列 ---
    match_cols_loop = ['Gene_ID', 'TAD_ID', 'Matched_Loop_ID']
    match_cols_pvalue = ['Gene_ID', 'TAD_ID', 'Loop_ID']

    # Rename and clean pvalue_df columns
    significant_loops_to_keep.rename(columns={'Loop_ID': 'Matched_Loop_ID'}, inplace=True)
    
    # Force string conversion and strip all whitespace from relevant columns
    for col in match_cols_loop:
        if col in candidates_cond2.columns:
            candidates_cond2[col] = candidates_cond2[col].astype(str).str.strip().str.replace(' ', '')
    
    for col in match_cols_pvalue:
        if col in significant_loops_to_keep.columns:
            significant_loops_to_keep[col] = significant_loops_to_keep[col].astype(str).str.strip().str.replace(' ', '')
    
    # 修复：确保 Matched_Loop_ID 列类型一致，处理浮点数
    candidates_cond2['Matched_Loop_ID'] = pd.to_numeric(candidates_cond2['Matched_Loop_ID'], errors='coerce').astype('Int64').astype(str)
    significant_loops_to_keep['Matched_Loop_ID'] = pd.to_numeric(significant_loops_to_keep['Matched_Loop_ID'], errors='coerce').astype('Int64').astype(str)
    
    # 调试输出
    print(f"调试：candidates_cond2 条数: {len(candidates_cond2)}", file=sys.stderr)
    print(f"调试：significant_loops_to_keep 条数: {len(significant_loops_to_keep)}", file=sys.stderr)
    if len(candidates_cond2) > 0:
        print(f"调试：candidates_cond2 样本:\n{candidates_cond2[['Gene_ID', 'TAD_ID', 'Matched_Loop_ID']].head()}", file=sys.stderr)
    if len(significant_loops_to_keep) > 0:
        print(f"调试：significant_loops_to_keep 样本:\n{significant_loops_to_keep.head()}", file=sys.stderr)
            
    # 将候选行与显著的Loop组合进行内连接（inner join）
    df_kept_by_cond2 = pd.merge(
        candidates_cond2,
        significant_loops_to_keep,
        on=['Gene_ID', 'TAD_ID', 'Matched_Loop_ID'],
        how='inner'
    )

    print(f"第二阶段筛选完成: {len(df_kept_by_cond2)} 行因存在显著关联的Loop而被保留。", file=sys.stderr)
    
    print("\n--- 阶段4 (Python子脚本): 合并结果并保存... ---", file=sys.stderr)

    # --- 合并与保存 ---
    final_df = pd.concat([df_kept_by_cond1, df_kept_by_cond2], ignore_index=True)
    
    # 删除辅助列
    final_df = final_df.drop(columns=['distance_numeric_kb'])
    
    # 保存结果
    final_df.to_csv(output_file, index=False)
    
    print(f"筛选完成，最终共保留 {len(final_df)} 行。", file=sys.stderr)
    print(f"结果已保存至文件: {output_file}", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="根据距离关系和P值显著性，对Loop分析结果进行两阶段筛选。")
    parser.add_argument('--loop_analysis_file', required=True, help="第6步生成的Loop分析文件路径。")
    parser.add_argument('--pvalue_summary_file', required=True, help="第7步生成的P值汇总文件路径。")
    parser.add_argument('--output_file', required=True, help="筛选后的最终输出文件路径。")
    parser.add_argument('--p_threshold', type=float, default=0.05, help="用于判断显著性的P值阈值 (默认为0.05)。")
    
    args = parser.parse_args()
    
    filter_loops_by_conditions(
        args.loop_analysis_file,
        args.pvalue_summary_file,
        args.output_file,
        args.p_threshold
    )