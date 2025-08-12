import argparse
import pandas as pd
import re # 导入正则表达式模块，用于解析 chr_pos 格式

def extract_deduplicate_and_sort_snps(input_filepath: str, output_filepath: str, snp_column_one_indexed: int, has_header: bool):
    """
    从输入文件中提取指定列的 SNP ID，进行去重，并按照 chr_pos 顺序（染色体从小到大，位置从小到大）
    保存到新文件中，不包含标题行。
    列索引从 1 开始计数，即第一列的索引是 1。

    Args:
        input_filepath (str): 输入文件的路径。
        output_filepath (str): 提取、去重并排序后的 SNP ID 的输出文件路径。
        snp_column_one_indexed (int): 包含 SNP ID 的列的索引（从 1 开始计数）。
        has_header (bool): 输入文件是否包含标题行 (True 表示有，False 表示没有)。
    """
    print(f"开始从文件 '{input_filepath}' 提取、去重并排序 SNP ID...")
    print(f"SNP ID 位于列索引: {snp_column_one_indexed} (从1开始计数)")
    print(f"输入文件 {'有' if has_header else '没有'} 标题行。")

    # 将用户提供的 1-based 索引转换为 Python 中使用的 0-based 索引
    snp_column_zero_indexed = snp_column_one_indexed - 1

    # 检查转换后的索引是否有效
    if snp_column_zero_indexed < 0:
        print(f"错误：列索引必须大于等于 1。您输入了 {snp_column_one_indexed}。")
        return

    try:
        # 根据 has_header 参数设置 pandas 的 header 参数
        header_param = 0 if has_header else None

        # 使用 pandas 读取 TSV 文件，强制所有列为字符串
        df = pd.read_csv(
            input_filepath,
            sep='\t',
            header=header_param,
            dtype=str
        )

        # 确保 DataFrame 有足够的列数
        if snp_column_zero_indexed >= df.shape[1]:
            print(f"错误：指定的列索引 {snp_column_one_indexed} (0-based: {snp_column_zero_indexed}) 超出了文件 '{input_filepath}' 的实际列数 ({df.shape[1]})。")
            return

        # 提取 SNP ID 列的数据
        snp_ids = df.iloc[:, snp_column_zero_indexed]
        print(f"原始提取了 {len(snp_ids)} 个 SNP ID。")

        # 对提取的 SNP ID 进行去重操作
        unique_snp_ids = snp_ids.drop_duplicates().reset_index(drop=True)
        print(f"去重后剩下 {len(unique_snp_ids)} 个唯一 SNP ID。")

        # --- 排序核心逻辑 ---
        # 创建一个列表来存储解析后的 SNP ID (chr, pos, original_id)
        parsed_snps = []
        # 正则表达式用于匹配 'chr_pos' 格式
        # group(1) 匹配染色体部分 (数字)
        # group(2) 匹配位置部分 (数字)
        snp_pattern = re.compile(r'(\d+)_(\d+)')

        for snp_id in unique_snp_ids:
            match = snp_pattern.match(snp_id)
            if match:
                # 提取染色体和位置，并转换为整数进行排序
                chrom = int(match.group(1))
                pos = int(match.group(2))
                parsed_snps.append((chrom, pos, snp_id))
            else:
                # 如果 SNP ID 格式不符合 chr_pos，则打印警告并跳过排序
                # 实际上，根据提供的样本，所有SNP ID都应符合格式。
                # 但为了严谨性，可以加入此处理。这里简化为跳过。
                print(f"警告：SNP ID 格式 '{snp_id}' 不符合 'chr_pos'，将跳过排序该ID。")
                parsed_snps.append((float('inf'), float('inf'), snp_id)) # 将不符合格式的放到最后

        # 对解析后的 SNP 列表进行排序
        # 默认会先按第一个元素 (chrom) 排序，然后按第二个元素 (pos) 排序
        parsed_snps.sort()

        # 提取排序后的原始 SNP ID
        sorted_snp_ids = [snp_id for chrom, pos, snp_id in parsed_snps]

        # 将排序后的 SNP ID 保存到新的文件，不包含标题行
        # index=False 表示不写入行索引
        # header=False 表示不写入列标题
        pd.Series(sorted_snp_ids).to_csv(output_filepath, sep='\t', index=False, header=False)

        print(f"去重并排序后的 SNP ID 已成功保存到: '{output_filepath}'")
        print(f"最终文件中包含 {len(sorted_snp_ids)} 个唯一且排序的 SNP ID。")

    except FileNotFoundError:
        print(f"错误：文件未找到，请检查路径: {input_filepath}")
    except Exception as e:
        print(f"处理过程中发生错误: {e}")

if __name__ == "__main__":
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser(
        description="从指定文件中提取 SNP ID，进行去重，并按照 chr_pos 顺序（染色体从小到大，位置从小到大）保存到新文件。"
    )

    # 添加命令行参数
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help="输入文件的路径。"
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help="提取、去重并排序后的 SNP ID 的输出文件路径。"
    )
    parser.add_argument(
        '-c', '--column',
        type=int,
        required=True,
        help="包含 SNP ID 的列的索引（从 1 开始计数）。"
    )
    parser.add_argument(
        '-H', '--has_header',
        action='store_true', # 这是一个布尔标志，如果存在则为 True，否则为 False
        help="如果输入文件包含标题行，请使用此标志。此标志仅用于防止读取时跳过第一行数据，提取仍按列索引进行。"
    )

    # 解析命令行参数
    args = parser.parse_args()

    # 调用函数执行操作，传入解析后的参数
    extract_deduplicate_and_sort_snps(args.input, args.output, args.column, args.has_header)
