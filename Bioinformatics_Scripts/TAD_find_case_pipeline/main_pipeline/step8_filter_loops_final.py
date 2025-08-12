#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import sys

def filter_loops_by_conditions(loop_analysis_file, pvalue_summary_file, output_file, p_threshold=0.05):
    """
    根据两阶段的复杂条件筛选Loop分析结果。

    Args:
        loop_analysis_file (str): 第6步生成的loop分析文件路径。
        pvalue_summary_file (str): 第7步生成的p值汇总文件路径。
        output_file (str): 最终筛选结果的输出文件路径。
        p_threshold (float): P值的显著性阈值。
    """
    try:
        print("--- 阶段1 (Python子脚本): 正在加载输入文件... ---", file=sys.stderr)
        loop_df = pd.read_csv(loop_analysis_file)
        pvalue_df = pd.read_csv(pvalue_summary_file)
        print(f"成功加载 Loop分析文件，共 {len(loop_df)} 行。", file=sys.stderr)
        print(f"成功加载 P值汇总文件，共 {len(pvalue_df)} 行。", file=sys.stderr)

    except FileNotFoundError as e:
        print(f"错误：找不到文件 {e.filename}。请检查文件路径。", file=sys.stderr)
        sys.exit(1)

    print("\n--- 阶段2 (Python子脚本): 开始执行第一阶段筛选 (基于距离)... ---", file=sys.stderr)
    
    # --- 第一阶段筛选 ---
    # 1. 创建辅助列，用于存放统一单位（kb）的数值
    loop_df['distance_numeric_kb'] = np.nan
    
    # 2. 统一单位并转换为数值
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

    # 3. 构建筛选条件①
    cond1_sub1 = (loop_df['Distance_Relationship'] == 'Gene in TAD')
    # 使用更新后的范围 (-50kb, +50kb)
    cond1_sub2 = (loop_df['distance_numeric_kb'] > -50) & (loop_df['distance_numeric_kb'] < 50)
    final_condition1 = cond1_sub1 | cond1_sub2
    
    # 4. 根据条件①，将数据拆分为两部分
    df_kept_by_cond1 = loop_df[final_condition1].copy()
    df_remaining = loop_df[~final_condition1].copy()
    
    print(f"第一阶段筛选完成: {len(df_kept_by_cond1)} 行因满足距离条件而被直接保留。", file=sys.stderr)

    print("\n--- 阶段3 (Python子脚本): 开始执行第二阶段筛选 (基于Loop和P值)... ---", file=sys.stderr)
    
    # --- 第二阶段筛选 (仅对第一阶段未被保留的数据进行) ---
    
    # 1. 从剩余数据中，筛选出存在Loop的候选行
    candidates_cond2 = df_remaining[df_remaining['Matched_Loop_ID'].notna()].copy()
    
    # 2. 从P值文件中，筛选出至少有一个P值显著的记录
    significant_pvalues = pvalue_df[
        (pvalue_df['Wilcoxon_P_Value_Loop_Strength'] < p_threshold) |
        (pvalue_df['Categorical_P_Value_Loop_Presence'] < p_threshold)
    ].copy()
    
    # 3. 提取显著记录中唯一的 Gene-TAD-Loop 组合，用于快速匹配
    # 我们只需要Gene_ID, TAD_ID, 和Loop_ID来标识一个需要保留的Loop
    significant_loops_to_keep = significant_pvalues[['Gene_ID', 'TAD_ID', 'Loop_ID']].drop_duplicates()
    # 将Matched_Loop_ID重命名以进行合并
    significant_loops_to_keep.rename(columns={'Loop_ID': 'Matched_Loop_ID'}, inplace=True)
    
    # 4. 将候选行与显著的Loop组合进行内连接（inner join）
    # 只有那些既存在Loop，又在p值文件中有显著记录的行才会被保留下来
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
