#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import sys

def create_quartet_summary(filtered_loop_file, stats_file, pvalue_file, output_file):
    """
    整合多个步骤的结果，创建一个以 Gene-TAD-Loop-SNP 四元组为核心的最终汇总表。
    """
    try:
        print("--- 阶段1 (Python子脚本): 正在加载所有输入文件... ---", file=sys.stderr)
        filtered_loop_df = pd.read_csv(filtered_loop_file)
        stats_df = pd.read_csv(stats_file)
        pvalue_df = pd.read_csv(pvalue_file)
        print(f"成功加载筛选后的Loop文件，共 {len(filtered_loop_df)} 行。", file=sys.stderr)
        print(f"成功加载初步统计文件，共 {len(stats_df)} 行。", file=sys.stderr)
        print(f"成功加载P值汇总文件，共 {len(pvalue_df)} 行。", file=sys.stderr)

    except FileNotFoundError as e:
        print(f"错误：找不到文件 {e.filename}。请检查文件路径。", file=sys.stderr)
        sys.exit(1)

    print("\n--- 阶段2 (Python子脚本): 准备关联数据... ---", file=sys.stderr)
    
    # 从初步统计文件中选取我们需要的SNP和表型P值列
    stats_subset = stats_df[[
        'Gene_ID', 'TAD_ID', 'SNP_ID', 'Full_SNP_ID', 
        'Wilcox_P_TAD', 'Wilcox_P_Gene', 'Wilcox_P_FL', 'Wilcox_P_FS',
        'Cor_R_Gene_TAD', 'Cor_P_Gene_TAD'
    ]].copy()

    # 将Loop-SNP的P值文件也准备好
    pvalue_subset = pvalue_df[[
        'Gene_ID', 'TAD_ID', 'Loop_ID', 'Full_SNP_ID',
        'Wilcoxon_P_Value_Loop_Strength', 'Categorical_P_Value_Loop_Presence'
    ]].copy()
    # 将Loop_ID重命名以匹配filtered_loop_df中的列名
    pvalue_subset.rename(columns={'Loop_ID': 'Matched_Loop_ID'}, inplace=True)

    print("--- 阶段3 (Python子脚本): 按是否存在Loop拆分并分别处理... ---", file=sys.stderr)

    # 将第8步的结果拆分为“有Loop”和“无Loop”两部分
    loops_present_df = filtered_loop_df[filtered_loop_df['Matched_Loop_ID'].notna()].copy()
    loops_absent_df = filtered_loop_df[filtered_loop_df['Matched_Loop_ID'].isna()].copy()
    
    final_present_list = []
    if not loops_present_df.empty:
        # 1. 处理有Loop的部分：首先根据Gene-TAD匹配所有关联的SNP
        merged_present_1 = pd.merge(
            loops_present_df,
            stats_subset,
            on=['Gene_ID', 'TAD_ID'],
            how='left' # 使用left join保留所有loop记录
        )
        
        # 2. 再将Loop的P值信息合并进来
        final_present = pd.merge(
            merged_present_1,
            pvalue_subset,
            on=['Gene_ID', 'TAD_ID', 'Matched_Loop_ID', 'Full_SNP_ID'],
            how='left'
        )
        final_present_list.append(final_present)
        print(f"处理完成: {len(loops_present_df)} 条有Loop的记录被扩展为 {len(final_present)} 个四元组。", file=sys.stderr)

    final_absent_list = []
    if not loops_absent_df.empty:
        # 3. 处理无Loop的部分：只需根据Gene-TAD匹配所有关联的SNP
        final_absent = pd.merge(
            loops_absent_df,
            stats_subset,
            on=['Gene_ID', 'TAD_ID'],
            how='left'
        )
        # 为这些记录手动添加空的Loop P值列
        final_absent['Wilcoxon_P_Value_Loop_Strength'] = np.nan
        final_absent['Categorical_P_Value_Loop_Presence'] = np.nan
        final_absent_list.append(final_absent)
        print(f"处理完成: {len(loops_absent_df)} 条无Loop的记录被扩展为 {len(final_absent)} 个四元组。", file=sys.stderr)
        
    print("\n--- 阶段4 (Python子脚本): 合并并格式化最终的汇总表... ---", file=sys.stderr)

    # 4. 合并所有处理好的部分
    if not final_present_list and not final_absent_list:
        print("警告：没有任何数据可以生成最终的四元组汇总表。", file=sys.stderr)
        pd.DataFrame().to_csv(output_file, index=False)
        sys.exit(0)
        
    final_df = pd.concat(final_present_list + final_absent_list, ignore_index=True)
    
    # --- 【核心修改】在此处增加对四元组的去重操作 ---
    print(f"\n--- 阶段5 (Python子脚本): 对最终四元组进行去重... ---", file=sys.stderr)
    print(f"去重前共有 {len(final_df)} 行。", file=sys.stderr)
    
    # 根据四个核心ID列进行去重
    final_df.drop_duplicates(
        subset=['Gene_ID', 'TAD_ID', 'Matched_Loop_ID', 'Full_SNP_ID'],
        inplace=True
    )
    
    print(f"去重后剩下 {len(final_df)} 行唯一的四元组。", file=sys.stderr)
    # --- 去重操作结束 ---

    print("\n--- 阶段6 (Python子脚本): 按要求排序列并保存... ---", file=sys.stderr)

    # 5. 定义最终输出的列和顺序
    final_columns_order = [
        'Gene_ID', 'TAD_ID', 'Matched_Loop_ID', 
        'SNP_ID', 
        'Full_SNP_ID',
        'Wilcox_P_TAD', 'Wilcox_P_Gene', 'Wilcox_P_FL', 'Wilcox_P_FS',
        'Cor_R_Gene_TAD', 'Cor_P_Gene_TAD',
        'Wilcoxon_P_Value_Loop_Strength', 'Categorical_P_Value_Loop_Presence',
        'Chromosome', 'Gene_Start', 'Gene_End', 'TAD_Start', 'TAD_End',
        'Distance_Relationship', 'Loop_Anchor1_Start', 'Loop_Anchor1_End',
        'Loop_Anchor2_Start', 'Loop_Anchor2_End', 'Anchor_Location',
        'Anchor1_Fine_Position', 'Anchor2_Fine_Position'
    ]
    
    # 筛选并重排序
    final_df = final_df[final_columns_order]
    
    # 保存结果
    final_df.to_csv(output_file, index=False)
    
    print(f"汇总完成，最终生成了 {len(final_df)} 条唯一的四元组记录。", file=sys.stderr)
    print(f"结果已保存至文件: {output_file}", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="整合多个流程步骤的结果，生成最终的Gene-TAD-Loop-SNP四元组汇总表。")
    parser.add_argument('--filtered_loop_file', required=True, help="第8步筛选后的Loop分析文件路径。")
    parser.add_argument('--stats_file', required=True, help="第5步生成的初步统计文件路径。")
    parser.add_argument('--pvalue_file', required=True, help="第7步生成的Loop P值汇总文件路径。")
    parser.add_argument('--output_file', required=True, help="最终四元组汇总表的输出路径。")
    
    args = parser.parse_args()
    
    create_quartet_summary(
        args.filtered_loop_file,
        args.stats_file,
        args.pvalue_file,
        args.output_file
    )