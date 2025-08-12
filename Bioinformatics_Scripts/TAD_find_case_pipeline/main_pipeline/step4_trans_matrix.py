#!/usr/bin/env python3
import pandas as pd
import argparse
import re
import sys # 修改点⑤: 引入sys模块

def transform_snp_matrix(input_raw_matrix_path: str, output_matrix_path: str):
    """
    读取 PLINK2 生成的原始 SNP 矩阵，删除不需要的列，转置，
    并重新格式化 SNP ID，最后保存结果。

    Args:
        input_raw_matrix_path (str): PLINK2 生成的原始矩阵文件路径 (e.g., extracted_snps_matrix.raw)。
        output_matrix_path (str): 转置并格式化后的矩阵输出文件路径 (e.g., snp_matrix_transposed.txt)。
    """
    # 修改点⑤: 将所有进度/日志信息输出到 stderr
    print(f"开始转置和格式化矩阵文件: '{input_raw_matrix_path}'...", file=sys.stderr)

    try:
        # 读取文件
        # sep='\s+' 表示以一个或多个空白字符作为分隔符
        df = pd.read_csv(input_raw_matrix_path, sep='\s+')
        # 修改点⑤: 将所有进度/日志信息输出到 stderr
        print(f"已读取原始矩阵文件，原始维度: {df.shape[0]} 行, {df.shape[1]} 列。", file=sys.stderr)

        # 删除不需要的列，只保留IID和SNP列
        # 保留 'IID' 列以及那些包含 '_' 但不是 PLINK2 默认元数据列的列
        cols_to_keep = ['IID'] + [col for col in df.columns if '_' in col and col not in ['FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE']]
        df_filtered = df[cols_to_keep]
        # 修改点⑤: 将所有进度/日志信息输出到 stderr
        print(f"过滤掉不必要的列后，维度: {df_filtered.shape[0]} 行, {df_filtered.shape[1]} 列。", file=sys.stderr)


        # 转置 DataFrame，并将 'IID' 列设置为索引
        df_transposed = df_filtered.set_index('IID').T
        # 修改点⑤: 将所有进度/日志信息输出到 stderr
        print(f"矩阵转置完成，新维度: {df_transposed.shape[0]} 行 (SNP), {df_transposed.shape[1]} 列 (样本)。", file=sys.stderr)

        # 重新格式化SNP ID (根据您的要求，此部分逻辑不变)
        new_index = []
        snp_id_pattern = re.compile(r'(\d+_\d+)_([A-Za-z]+)(?:\(/([A-Za-z]+)\)|_([A-Za-z]+))')

        for snp in df_transposed.index:
            match = snp_id_pattern.match(snp)
            if match:
                chr_pos = match.group(1)
                ref_allele = match.group(2)
                alt_allele = match.group(3) if match.group(3) else match.group(4)
                
                if alt_allele:
                    new_snp = f"{chr_pos}:{ref_allele}:{alt_allele}"
                else:
                    new_snp = f"{chr_pos}:{ref_allele}"
            else:
                # 修改点⑤: 将所有进度/日志信息输出到 stderr
                print(f"警告：SNP ID 格式 '{snp}' 不符合预期，将保留原样。", file=sys.stderr)
                new_snp = snp
            new_index.append(new_snp)

        df_transposed.index = new_index
        df_transposed.index.name = 'SNP_ID'

        # 保存结果
        df_transposed.to_csv(output_matrix_path, sep='\t')
        # 修改点⑤: 将所有进度/日志信息输出到 stderr
        print(f"转置和格式化完成，最终结果已保存到: '{output_matrix_path}'", file=sys.stderr)
        print(f"最终矩阵样本数: {df_transposed.shape[1]}, SNP数: {df_transposed.shape[0]}", file=sys.stderr)

    except FileNotFoundError:
        # 修改点⑤: 将所有进度/日志信息输出到 stderr
        print(f"错误：输入矩阵文件未找到，请检查路径: {input_raw_matrix_path}", file=sys.stderr)
        # 增加退出码，更规范
        sys.exit(1)
    except Exception as e:
        # 修改点⑤: 将所有进度/日志信息输出到 stderr
        print(f"处理过程中发生未知错误: {e}", file=sys.stderr)
        # 增加退出码，更规范
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="转置 PLINK2 原始矩阵并格式化 SNP ID。"
    )
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help="PLINK2 生成的原始矩阵文件路径 (e.g., extracted_snps_matrix.raw)"
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        required=True,
        help="转置并格式化后的矩阵输出文件路径 (e.g., snp_matrix_transposed.txt)"
    )
    args = parser.parse_args()
    transform_snp_matrix(args.input, args.output)