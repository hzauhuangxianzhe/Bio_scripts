#!/usr/bin/env python3
import pandas as pd
import argparse
import re

def transform_snp_matrix(input_raw_matrix_path: str, output_matrix_path: str):
    """
    读取 PLINK2 生成的原始 SNP 矩阵，删除不需要的列，转置，
    并重新格式化 SNP ID，最后保存结果。

    Args:
        input_raw_matrix_path (str): PLINK2 生成的原始矩阵文件路径 (e.g., extracted_snps_matrix.raw)。
        output_matrix_path (str): 转置并格式化后的矩阵输出文件路径 (e.g., snp_matrix_transposed.txt)。
    """
    print(f"开始转置和格式化矩阵文件: '{input_raw_matrix_path}'...")

    try:
        # 读取文件
        # sep='\s+' 表示以一个或多个空白字符作为分隔符
        df = pd.read_csv(input_raw_matrix_path, sep='\s+')
        print(f"已读取原始矩阵文件，原始维度: {df.shape[0]} 行, {df.shape[1]} 列。")

        # 删除不需要的列，只保留IID和SNP列
        # 保留 'IID' 列以及那些包含 '_' 但不是 PLINK2 默认元数据列的列
        cols_to_keep = ['IID'] + [col for col in df.columns if '_' in col and col not in ['FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE']]
        df_filtered = df[cols_to_keep]
        print(f"过滤掉不必要的列后，维度: {df_filtered.shape[0]} 行, {df_filtered.shape[1]} 列。")


        # 转置 DataFrame，并将 'IID' 列设置为索引
        df_transposed = df_filtered.set_index('IID').T
        print(f"矩阵转置完成，新维度: {df_transposed.shape[0]} 行 (SNP), {df_transposed.shape[1]} 列 (样本)。")

        # 重新格式化SNP ID
        new_index = []
        # 正则表达式用于匹配 'chr_pos_ref_alt' 或 'chr_pos_allele(ref/alt)' 等复杂格式
        # 目标是将其统一为 chr:pos:ref:alt 或 chr_pos:ref:alt
        # 例如 1_228282_A_G 或 1_228282_A(/G) 转换为 1:228282:A:G
        # 考虑到PLINK2 --export A include-alt 的输出格式通常是 CHROM:POS_A1/A2
        # 我们需要更通用的匹配或直接拆分
        
        # 尝试匹配 'chr_pos_allele1_allele2' 或 'chr_pos_allele1(/allele2)' 形式
        # 这里假设SNP ID通常为 CHROM_POS_REF_ALT 或者 CHROM_POS_REF(ALT)
        # PLINK2 --export A include-alt 默认输出格式通常是 CHROM:POS_A1_A2 (例如 1:123_A_G)
        # 但如果VCF的ID列有值，它会用ID列，这需要注意
        
        # 定义一个更健壮的正则表达式，以适应 PLINK2 --export A include-alt 的输出格式
        # 常见的格式可能是 CHROM:POS_ALLELE1/ALLELE2 或者 CHROM_POS_ALLELE1_ALLELE2 (如果原始VCF ID列是CHROM_POS)
        # 但根据例子，SNP ID的原始格式是 chr_pos (例如 1_228282)
        # 而转置后 PLINK2 可能会生成如 1_228282_A_G 或者 1_228282_A(/G) 这样的列名
        snp_id_pattern = re.compile(r'(\d+_\d+)_([A-Za-z]+)(?:\(/([A-Za-z]+)\)|_([A-Za-z]+))')

        for snp in df_transposed.index:
            match = snp_id_pattern.match(snp)
            if match:
                chr_pos = match.group(1) # 例如 '1_228282'
                ref_allele = match.group(2) # 例如 'A'
                # 尝试获取 / 括号内的ALT或下划线后的ALT
                alt_allele = match.group(3) if match.group(3) else match.group(4)
                
                if alt_allele:
                    new_snp = f"{chr_pos}:{ref_allele}:{alt_allele}"
                else: # 如果没有alt等位基因（比如原始ID就只有chr_pos_ref）
                    new_snp = f"{chr_pos}:{ref_allele}"
            else:
                # 如果不匹配预期的模式，则保持原样，并打印警告
                # 这可能是因为原始VCF的ID列就有复杂格式，或者PLINK2输出了其他形式
                print(f"警告：SNP ID 格式 '{snp}' 不符合预期，将保留原样。")
                new_snp = snp
            new_index.append(new_snp)

        df_transposed.index = new_index
        df_transposed.index.name = 'SNP_ID'

        # 保存结果
        df_transposed.to_csv(output_matrix_path, sep='\t')
        print(f"转置和格式化完成，最终结果已保存到: '{output_matrix_path}'")
        print(f"最终矩阵样本数: {df_transposed.shape[1]}, SNP数: {df_transposed.shape[0]}")

    except FileNotFoundError:
        print(f"错误：输入矩阵文件未找到，请检查路径: {input_raw_matrix_path}")
    except Exception as e:
        print(f"处理过程中发生错误: {e}")

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