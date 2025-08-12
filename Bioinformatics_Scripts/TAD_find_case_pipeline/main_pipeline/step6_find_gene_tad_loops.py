import pandas as pd
import numpy as np
import sys
import argparse # 导入用于解析命令行参数的库

# =============================================================================
# 辅助函数定义 (这部分保持不变)
# =============================================================================

def format_distance(distance_value, prefix):
    """
    根据距离的绝对值格式化输出字符串。
    - 大于1000的转换为kb（保留两位小数）。
    - 小于等于1000的用bp表示。
    """
    if distance_value > 1000:
        return f"{prefix}{distance_value / 1000:.2f}kb"
    else:
        return f"{prefix}{int(distance_value)}bp"

def calculate_distance(gene_start, gene_end, tad_start, tad_end):
    """
    根据基因和TAD的位置关系，计算并格式化距离字符串。
    """
    if max(gene_start, tad_start) < min(gene_end, tad_end):
        return "Gene in TAD"
    elif gene_end < tad_start:
        distance = tad_start - gene_end
        return format_distance(distance, '-')
    elif gene_start > tad_end:
        distance = gene_start - tad_end
        return format_distance(distance, '+')
    else:
        return "Unknown"

def check_overlap(start1, end1, start2, end2):
    """
    检查两个区间 [start1, end1] 和 [start2, end2] 是否有重叠。
    """
    return max(start1, start2) < min(end1, end2)

def get_anchor_location(anchor_start, anchor_end, tad_loop_start, tad_loop_end, gene_loop_start, gene_loop_end):
    """
    确定单个锚点相对于TAD和Gene扩展区域的位置。
    返回 'O' (Overlap), 'T' (TAD only), 'G' (Gene only)。
    """
    in_tad = check_overlap(anchor_start, anchor_end, tad_loop_start, tad_loop_end)
    in_gene = check_overlap(anchor_start, anchor_end, gene_loop_start, gene_loop_end)
    
    if in_tad and in_gene:
        return 'O'
    elif in_tad:
        return 'T'
    elif in_gene:
        return 'G'
    else:
        return 'NA'

def fine_anchor_position(anchor_start, anchor_end, region_start, region_end, extension_start, extension_end):
    """
    计算锚点相对于实际区域的位置。
    """
    if check_overlap(anchor_start, anchor_end, region_start, region_end):
        return "True_IN"
    elif check_overlap(anchor_start, anchor_end, extension_start, extension_end):
        if anchor_end < region_start:
            distance = region_start - anchor_end
            prefix = '-'
        elif anchor_start > region_end:
            distance = anchor_start - region_end
            prefix = '+'
        else:
            return "NA" 
        return format_distance(distance, prefix)
    else:
        return "NA"

# =============================================================================
# 主逻辑
# =============================================================================
def main(stats_file, tad_pos_file, gene_pos_file, loop_file, output_file):
    """
    脚本的主执行函数。
    """
    # ------------------- 文件路径现在从函数参数传入 -------------------
    try:
        # ------------------- 2. 预加载并处理关联数据 -------------------
        # 修改点⑤: 将所有进度/日志信息输出到 stderr
        print("--- 阶段1 (Python子脚本): 正在加载所有必需的数据文件... ---", file=sys.stderr)
        
        tad_info = pd.read_csv(tad_pos_file, sep='\t').set_index('TAD_ID')
        print(f"成功加载 {len(tad_info)} 条TAD位置信息。", file=sys.stderr)
        
        gene_info_raw = pd.read_csv(gene_pos_file, sep='\t')
        gene_info_raw.rename(columns={gene_info_raw.columns[0]: gene_info_raw.columns[0].lstrip('#')}, inplace=True)
        gene_info = gene_info_raw[['chr', 'start', 'end', 'pid']].set_index('pid')
        print(f"成功加载 {len(gene_info)} 条基因位置信息。", file=sys.stderr)
        
        loop_info = pd.read_csv(loop_file, sep='\t')
        loops_by_chrom = {chrom: group for chrom, group in loop_info.groupby('chrom1')}
        print(f"成功加载 {len(loop_info)} 条Loop信息，并按染色体分组。", file=sys.stderr)
        
        main_input = pd.read_csv(stats_file)
        print(f"成功加载主输入文件，准备处理 {len(main_input)} 个基因-TAD对。", file=sys.stderr)
        
    except FileNotFoundError as e:
        # 修改点⑤: 将错误信息输出到 stderr
        print(f"错误：找不到文件 {e.filename}。请检查文件路径是否正确。", file=sys.stderr)
        sys.exit(1)

    # ------------------- 3. 循环处理每个基因-TAD对 -------------------
    print("\n--- 阶段2 (Python子脚本): 开始逐行处理基因-TAD对... ---", file=sys.stderr)
    
    results_list = []
    
    for row in main_input.itertuples(index=False):
        gene_id = row.Gene_ID
        tad_id = row.TAD_ID
        
        try:
            gene_record = gene_info.loc[gene_id]
            tad_record = tad_info.loc[tad_id]
        except KeyError:
            # 修改点⑤: 将静默跳过改为明确的警告日志
            print(f"警告: 在位置文件中找不到 Gene_ID '{gene_id}' 或 TAD_ID '{tad_id}'，已跳过此对。", file=sys.stderr)
            continue
            
        if gene_record['chr'] != tad_record['chrom']:
            # 修改点⑤: 为染色体不匹配的情况增加警告日志
            print(f"警告: Gene '{gene_id}' (CHR {gene_record['chr']}) 与 TAD '{tad_id}' (CHR {tad_record['chrom']}) 染色体不匹配，已跳过此对。", file=sys.stderr)
            continue
        
        gene_chrom = gene_record['chr']
        gene_start, gene_end = gene_record['start'], gene_record['end']
        tad_start, tad_end = tad_record['start'], tad_record['end']
        
        distance_str = calculate_distance(gene_start, gene_end, tad_start, tad_end)
        
        gene_loop_start, gene_loop_end = gene_start - 40000, gene_end + 40000
        tad_loop_start, tad_loop_end = tad_start - 40000, tad_end + 40000
        
        found_loops = []
        
        relevant_loops = loops_by_chrom.get(gene_chrom)
        
        if relevant_loops is not None:
            for loop in relevant_loops.itertuples(index=False):
                anchor1_start, anchor1_end = loop.start1, loop.end1
                anchor2_start, anchor2_end = loop.start2, loop.end2
                
                cond1 = check_overlap(anchor1_start, anchor1_end, tad_loop_start, tad_loop_end) and \
                        check_overlap(anchor2_start, anchor2_end, gene_loop_start, gene_loop_end)
                
                cond2 = check_overlap(anchor1_start, anchor1_end, gene_loop_start, gene_loop_end) and \
                        check_overlap(anchor2_start, anchor2_end, tad_loop_start, tad_loop_end)

                if cond1 or cond2:
                    loc_A1 = get_anchor_location(anchor1_start, anchor1_end, tad_loop_start, tad_loop_end, gene_loop_start, gene_loop_end)
                    loc_A2 = get_anchor_location(anchor2_start, anchor2_end, tad_loop_start, tad_loop_end, gene_loop_start, gene_loop_end)
                    location_desc = f"A1:{loc_A1}_A2:{loc_A2}"
                    found_loops.append((loop, location_desc))
        
        base_info = [
            gene_chrom,
            gene_id, gene_start, gene_end,
            tad_id, tad_start, tad_end,
            distance_str
        ]
        
        if not found_loops:
            results_list.append(base_info + ['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'])
        else:
            for loop, location_desc in found_loops:
                loop_info_cols = [
                    loop.Loop_ID,
                    loop.start1, loop.end1,
                    loop.start2, loop.end2,
                    location_desc
                ]
                
                anchor1_pos = "NA"
                anchor2_pos = "NA"
                
                if location_desc == "A1:O_A2:O":
                    anchor1_pos = "Attention"; anchor2_pos = "Attention"
                elif location_desc == "A1:O_A2:T":
                    anchor1_pos = fine_anchor_position(loop.start1, loop.end1, gene_start, gene_end, gene_loop_start, gene_loop_end)
                    anchor2_pos = fine_anchor_position(loop.start2, loop.end2, tad_start, tad_end, tad_loop_start, tad_loop_end)
                elif location_desc == "A1:O_A2:G":
                    anchor1_pos = fine_anchor_position(loop.start1, loop.end1, tad_start, tad_end, tad_loop_start, tad_loop_end)
                    anchor2_pos = fine_anchor_position(loop.start2, loop.end2, gene_start, gene_end, gene_loop_start, gene_loop_end)
                elif location_desc == "A1:T_A2:O":
                    anchor1_pos = fine_anchor_position(loop.start1, loop.end1, tad_start, tad_end, tad_loop_start, tad_loop_end)
                    anchor2_pos = fine_anchor_position(loop.start2, loop.end2, gene_start, gene_end, gene_loop_start, gene_loop_end)
                elif location_desc == "A1:T_A2:G":
                    anchor1_pos = fine_anchor_position(loop.start1, loop.end1, tad_start, tad_end, tad_loop_start, tad_loop_end)
                    anchor2_pos = fine_anchor_position(loop.start2, loop.end2, gene_start, gene_end, gene_loop_start, gene_loop_end)
                elif location_desc == "A1:G_A2:T":
                    anchor1_pos = fine_anchor_position(loop.start1, loop.end1, gene_start, gene_end, gene_loop_start, gene_loop_end)
                    anchor2_pos = fine_anchor_position(loop.start2, loop.end2, tad_start, tad_end, tad_loop_start, tad_loop_end)
                elif location_desc == "A1:G_A2:O":
                    anchor1_pos = fine_anchor_position(loop.start1, loop.end1, gene_start, gene_end, gene_loop_start, gene_loop_end)
                    anchor2_pos = fine_anchor_position(loop.start2, loop.end2, tad_start, tad_end, tad_loop_start, tad_loop_end)
                
                results_list.append(base_info + loop_info_cols + [anchor1_pos, anchor2_pos])

    # ------------------- 4. 生成并保存最终文件 -------------------
    print(f"\n--- 阶段3 (Python子脚本): 处理完成，共生成 {len(results_list)} 条记录。正在生成输出文件... ---", file=sys.stderr)
    
    output_columns = [
        'Chromosome', 'Gene_ID', 'Gene_Start', 'Gene_End',
        'TAD_ID', 'TAD_Start', 'TAD_End', 'Distance_Relationship',
        'Matched_Loop_ID', 'Loop_Anchor1_Start', 'Loop_Anchor1_End',
        'Loop_Anchor2_Start', 'Loop_Anchor2_End', 'Anchor_Location',
        'Anchor1_Fine_Position', 'Anchor2_Fine_Position'
    ]
    
    output_df = pd.DataFrame(results_list, columns=output_columns)
    output_df.to_csv(output_file, index=False)
    print(f"分析完成！结果已保存到: {output_file}", file=sys.stderr)

# =============================================================================
# 脚本执行入口
# =============================================================================
if __name__ == "__main__":
    # --- **新功能**：使用 argparse 解析命令行参数 ---
    # 创建一个解析器对象，并提供脚本的描述信息
    parser = argparse.ArgumentParser(description="根据SNP-TAD-Gene统计结果，寻找关联的染色质环(Loop)并进行详细的位置分析。")
    
    # 添加所有必需的文件路径参数
    # 主输入文件，由步骤5生成
    parser.add_argument('--stats_file', required=True, help="输入的统计文件路径 (e.g., filtered_plotted_statistics.csv)")
    # 关联文件
    parser.add_argument('--tad_pos_file', required=True, help="TAD位置信息文件路径 (e.g., all_TAD_bin_pos_with_id.bed)")
    parser.add_argument('--gene_pos_file', required=True, help="基因位置信息文件路径 (e.g., formatted_rna_data_for_qtl.bed)")
    parser.add_argument('--loop_file', required=True, help="Loop位置信息文件路径 (e.g., sorted_loops.bed)")
    # 输出文件
    parser.add_argument('--output_file', required=True, help="输出分析结果的文件路径 (e.g., gene_tad_loop_analysis_results.csv)")
    
    # 解析命令行中实际传入的参数
    args = parser.parse_args()
    
    # 使用从命令行解析到的参数来调用主函数
    main(
        stats_file=args.stats_file,
        tad_pos_file=args.tad_pos_file,
        gene_pos_file=args.gene_pos_file,
        loop_file=args.loop_file,
        output_file=args.output_file
    )