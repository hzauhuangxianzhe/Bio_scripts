import argparse
import pandas as pd
import re # 导入正则表达式模块，用于解析 chr_pos 格式
import sys # 导入sys模块，用于在出错时退出

def extract_deduplicate_and_sort_snps(input_filepath: str, output_filepath: str, snp_column_one_indexed: int, has_header: bool):
    """
    从输入文件中提取指定列的 SNP ID，进行去重，并按照 chr_pos 顺序
    （染色体从小到大，位置从小到大）保存到新文件中。
    """
    # 修改点⑤: 将所有进度/日志信息输出到 stderr
    print(f"--- 开始从文件 '{input_filepath}' 提取 SNP ID ---", file=sys.stderr)
    print(f"目标列: 第 {snp_column_one_indexed} 列", file=sys.stderr)
    print(f"文件包含标题行: {'是' if has_header else '否'}", file=sys.stderr)

    # 将用户提供的 1-based 索引转换为 Python 中使用的 0-based 索引
    snp_column_zero_indexed = snp_column_one_indexed - 1

    if snp_column_zero_indexed < 0:
        print(f"错误：列索引必须大于等于 1。您提供了 {snp_column_one_indexed}。", file=sys.stderr)
        sys.exit(1)

    try:
        # 根据 has_header 参数设置 pandas 读取时的 header 参数
        # 如果有标题行，header=0 会使用第一行作为列名；否则 header=None 将不使用标题行
        header_param = 0 if has_header else None

        # --- 核心读取逻辑 ---
        # 使用 pandas 读取制表符分隔的文件。
        # **关键点**: `dtype=str` 强制将所有列都读取为字符串类型，
        # 这可以从根本上防止纯数字或科学记数法格式的SNP ID被错误地解析为数字（float或int），
        # 从而解决了 'expected string or bytes-like object, got float' 的错误。
        df = pd.read_csv(
            input_filepath,
            sep='\t',
            header=header_param,
            dtype=str,
            comment=None # 确保即使有#号开头的行也能被正常读取（如果header=None）
        )

        # 确保指定的列存在
        if snp_column_zero_indexed >= df.shape[1]:
            print(f"错误：指定的列索引 {snp_column_one_indexed} 超出了文件 '{input_filepath}' 的实际列数 ({df.shape[1]})。", file=sys.stderr)
            sys.exit(1)

        # 提取 SNP ID 列的数据，并移除所有可能的前后空格
        snp_ids = df.iloc[:, snp_column_zero_indexed].str.strip()
        # 修改点⑤: 将所有进度/日志信息输出到 stderr
        print(f"原始提取了 {len(snp_ids)} 个 SNP ID。", file=sys.stderr)

        # 去除NA/空值并进行去重
        unique_snp_ids = snp_ids.dropna().unique()
        # 修改点⑤: 将所有进度/日志信息输出到 stderr
        print(f"去重后剩下 {len(unique_snp_ids)} 个唯一 SNP ID。", file=sys.stderr)

        # --- 排序核心逻辑 ---
        # 正则表达式用于匹配 '数字_数字' 格式
        # 第一个括号 (\d+) 捕获染色体号，第二个括号 (\d+) 捕获位置号
        snp_pattern = re.compile(r'^(\d+)_(\d+)$')
        
        parsed_snps = []
        malformed_snps = []

        for snp_id in unique_snp_ids:
            match = snp_pattern.match(str(snp_id)) # 再次确保是字符串
            if match:
                # 如果匹配成功，提取染色体和位置，并转换为整数以便进行数值排序
                chrom = int(match.group(1))
                pos = int(match.group(2))
                parsed_snps.append((chrom, pos, snp_id))
            else:
                # 如果格式不匹配，记录下来并在最后处理
                malformed_snps.append(snp_id)
        
        if malformed_snps:
            # 修改点⑤: 将所有进度/日志信息输出到 stderr
            print(f"警告：发现 {len(malformed_snps)} 个格式不符合 'chr_pos' 的 SNP ID，它们将被放置在输出文件的末尾。", file=sys.stderr)
            # 打印前5个不符合格式的ID作为例子
            print(f"  例如: {malformed_snps[:5]}", file=sys.stderr)

        # 对解析成功的 SNP 列表进行排序
        # Python的sort会先按元组的第一个元素(chrom)排序，如果相同，再按第二个元素(pos)排序
        parsed_snps.sort()

        # 从排序后的元组列表中提取出原始的 SNP ID 字符串
        sorted_snp_ids = [snp_id for chrom, pos, snp_id in parsed_snps]
        
        # 将格式不正确的SNP ID附加到列表末尾
        final_sorted_list = sorted_snp_ids + malformed_snps

        # 将最终排序好的 SNP ID 列表写入输出文件
        # 不写入pandas的索引和标题行
        pd.Series(final_sorted_list).to_csv(output_filepath, index=False, header=False)

        # 修改点⑤: 将所有进度/日志信息输出到 stderr
        print(f"去重并排序后的 SNP ID 已成功保存到: '{output_filepath}'", file=sys.stderr)
        print(f"最终文件中包含 {len(final_sorted_list)} 个 SNP ID。", file=sys.stderr)

    except FileNotFoundError:
        print(f"错误：文件未找到 '{input_filepath}'。", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"处理过程中发生未知错误: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    # 创建一个命令行参数解析器对象
    parser = argparse.ArgumentParser(
        description="从指定文件中提取 SNP ID，去重，并按照染色体和位置进行排序。"
    )

    # 添加输入文件参数
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help="输入文件的路径 (制表符分隔)。"
    )
    # 添加输出文件参数
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help="排序去重后的 SNP ID 输出文件路径。"
    )
    # 添加列索引参数
    parser.add_argument(
        '-c', '--column',
        type=int,
        required=True,
        help="包含 SNP ID 的列索引（从 1 开始计数）。"
    )
    # 添加是否有标题行的标志参数
    parser.add_argument(
        '--has_header',
        action='store_true', # 如果使用此参数，其值为True，否则为False
        help="指定输入文件包含标题行。"
    )

    # 解析命令行传入的实际参数
    args = parser.parse_args()

    # 使用解析到的参数调用主函数
    extract_deduplicate_and_sort_snps(args.input, args.output, args.column, args.has_header)
