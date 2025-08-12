import pandas as pd
import sys

# 该脚本通过命令行参数接收输入和输出文件路径
if len(sys.argv) != 4:
    print("用法: python3 add_tm1_gene_id_step11.py <quartet_file> <gene_pair_file> <output_file>", file=sys.stderr)
    sys.exit(1)

quartet_input = sys.argv[1]
gene_pair_input = sys.argv[2]
output = sys.argv[3]

def add_tm1_gene_id_to_file():
    """
    对quartet_file文件添加TM1_Gene_ID列。
    """
    
    # 1. 读取基因对文件，以备查询
    try:
        # 使用sep='\s+'来稳健地处理空格或制表符分隔
        gene_map_df = pd.read_csv(gene_pair_input, sep=r'\s+', header=0)
        if 'TM1' not in gene_map_df.columns or 'HC04' not in gene_map_df.columns:
            print(f"错误: 基因对文件 '{gene_pair_input}' 的列名不符合预期 (需包含 'TM1' 和 'HC04')。", file=sys.stderr)
            sys.exit(1)
        print("成功加载基因对文件。")
    except FileNotFoundError:
        print(f"错误: 基因对文件未找到: {gene_pair_input}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"读取基因对文件时发生错误: {e}", file=sys.stderr)
        sys.exit(1)
    
    # 2. 读取第9步生成的四元组汇总文件
    try:
        quartet_df = pd.read_csv(quartet_input)
        if 'Gene_ID' not in quartet_df.columns:
            print(f"错误: 四元组文件 '{quartet_input}' 缺少必需的列 'Gene_ID'。", file=sys.stderr)
            sys.exit(1)
        print(f"成功加载四元组汇总文件，共 {len(quartet_df)} 行。")
    except FileNotFoundError:
        print(f"错误: 四元组文件未找到: {quartet_input}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"读取四元组文件时发生错误: {e}", file=sys.stderr)
        sys.exit(1)

    # 3. 将两个DataFrame进行左连接（Left Merge），根据HC04基因ID进行匹配
    # 'how=left' 确保保留 quartet_df 的所有行。
    # 如果找不到匹配项，新列 'TM1' 将自动标记为 NaN（在 CSV 中显示为 NA）。
    merged_df = pd.merge(quartet_df, gene_map_df, left_on='Gene_ID', right_on='HC04', how='left')

    # 4. 重命名新添加的'TM1'列为'TM1_Gene_ID'
    merged_df = merged_df.rename(columns={'TM1': 'TM1_Gene_ID'})

    # 5. 删除合并过程中产生的重复列，即原来的'HC04'列
    merged_df = merged_df.drop(columns=['HC04'])

    # 6. 将'TM1_Gene_ID'列移动到第一列
    tm1_gene_id_col = merged_df.pop('TM1_Gene_ID')
    merged_df.insert(0, 'TM1_Gene_ID', tm1_gene_id_col)
    
    # 7. 保存最终的CSV文件，不包含索引列
    merged_df.to_csv(output, index=False)

    print(f"处理完成，已生成包含TM1基因ID的最终汇总表: {output}")

if __name__ == "__main__":
    add_tm1_gene_id_to_file()